#!/usr/bin/env python

r'''
a code generator for serializing objects

Copyright (c) 2010 Cheng Zhang

It uses a template, locates "typedef struct" code
and try to generate io functions, 
 * cfgopen(): read parameters from configuration file
 * close():   recursively free call resources
 * clear():   clear data
 * binread(), binwrite(): read and write data in binary format

Commands syntax
  $cmd = args;
or
  $cmd;  (lazy command mainly for switches)
  + cmd contains characters, numbers, and _ or -
  + the operator = can be replaced by :, it can also be := for
    persistent commands
  + currently args cannot contain ; even within a string
    it should be replaced by \\;
  + $# ... ; means a comment on commands
  + if the terminal ; is missing, args extend to the end of line

Commands of an object:
  * $skip, $skipme: skip this object, do not generate anything

  * $fprefix:       function prefix for the object
  * $ptrname:       a common variable used as pointers to the object
  
  * $cfgopen:       opens configuration file in its initializer,
                    instead of accepting an existing handle

Commands of an item:
  * $def:     default value for a variable, 
              if the variable is a dynamic/static array, this is
              value for array members

  * $cnt:     declare a pointer variable as a dynamic array,
              ineffective for non-pointer variables
              for 2d array:
                $cnt: m, n;
              allocates m*n elements, but should be understood,
              when needed, as m x n array
              if $cnt = 0; no allocation is performed during
              initialization, a free() is however performed
              during close() if the pointer is not NULL.

  * $io:      declare I/O type, can be a combination of cbt,
              for configuration, binary and text respectively
              it can also be `none' or `all'

  * $key:     key to get from configuration file
  * $kprefix: prefix to be added when it is deduced from 
              variable name, usually as a persistent command
  * $kdepfx:  if variable name starts with $key_unprefix, that part 
              is first removed before applying $kprefix
  * $kargs:   printf arguments in constructing key
  
  * $must:       a critial key that must present in configuration file

  * $prereq:      prerequisite to be tested *before* reading a variable
                  from configuration file, but after assigning the 
                  given by $def; this is the default order because it
                  ensures a variable is assigned at a reasonable value
  * $tfirst:      if this is set, we do not assign the default value $def 
                  unless $prereq is true, useful for dummy variables
  * $valid:       a condition to be tested for validity *after* reading
                  a variable from configuration file

  * $bin_prereq:  a prerequisite that applies to reading/writing binary files
                  it can contain a variable `ver', a version number of
                  binary data.
                  binread()/binwrite() rely on the this condition.

  * $com_prereq:  a general prerequisite that applies to accessing data
                  usually related to MPI rank.
                  clear()/close() rely on this condition

  * $obj:     a member object, a function needs to called to properly 
              initialize it; 
              for an object array, also set $cnt;

  * $flag:    two fields, name of flag and its value, something like 
              $flag: ABC_FLAG  0x000010;

  * $usr:     the variable should be passed as a parameter, 
              not read from the configuration file
              $usr:parent; declare the variable represents 
              a pointer to the parent object
              $usr:bintmp; temporary variable used in bin

A fold is an embed object that has its own i/o routines
  * $fold:          follow by a fold-identifier
  * $fold_bin_add:  add some variables (separated by ,) during writing binary
  * $fold_prereq:   skip io if this condition is not true
                    applies to:
                      binread(), binwrite(), clear()
  * $fold_valid:    raise an error if condition is not true

In a stand-alone comment

  * $assert:  asserts a condition and aborts the program if 
              the condition fails
  * $call:    means a direct translation of the argument to code

In arguments of the command, `@' has special meanings
  * @@:       this item
  * @var:     ptrname->var
  * @-func:   means: fprefix ## func
  * @:        ptrname
  * @^:       parent (parent is the variable name that has $usr:parent;)
  * @~:       the corresponding member under parent, e.g.,
                Dad_t  *d; /* $usr:parent; */
                int    a;  /* $def: @~;  */
              the first line declares `d' as a pointer to parent variable,
              which is a parameter of cfgopen(), @^ means d;
              on the second line, the default value of `a' is set to 
              the corresponding variable in `d', that is, d->a.
               
TODO:
  * multiple-line support, e.g., typedef\nstruct
  * preprocessor
  * skip struct in comment or string
'''


import os, sys, re, filecmp
from copy    import copy, deepcopy 
from objpos  import P 
from objcmt  import CComment
from objcow  import CCodeWriter
from objcmd  import Commands
from objitm  import Item
from objext  import *

USE_MPI = "USE_MPI"  # macro indicating we have MPI

class Object:
  '''
  a struct defined as
  typedef struct {
    ...
  } abc_t; /* cmds */
  '''
  def __init__(self, src, pos = None):
    ''' initialize myself from a bunch of src '''
    if len(src) == 0: return # return an empty object

    if type(pos) == int:  # integer input
      pos = P(pos, 0)
    self.items = []
    self.empty = 1
    self.parse(src, pos)

  def parse(self, src, p):
    ''' parse src starting from 'p' '''
   
    # print "start searching object from %s" % p
    if not self.find_beginning(src, p):
      return -1

    while 1: # parse the content within { ... },
      # aggressively search for an item, and update p
      item = Item(src, p)
      if item.isempty(): break

      #print "line %3d has %d item(s): %s" % (p.row + 1, len(item.itlist), item)
      self.items += [item]  # do not expand multiple declarators, 
                            # since we need to merge comments first
    # parse }
    self.find_ending(src, p, 0)
    if self.cmds and ("skipme" in self.cmds
        or "skip" in self.cmds):
      self.empty = -1
      return -1

    self.merge_comments() # merging multiple-line C++ comments
    self.expand_multidecl_items() # must be done after merging comments
    self.fill_cmds()   # get cmds from comment, must be after merging comments
    self.add_dummies() # add dummy declarations for flags and alts
    self.get_gtype()   # generic type, needs cmds, must be after get_cmds
    self.get_prefix()  # determine fprefix, ptrname and parent
    self.init_folds()  # initialize different folds, get prefices
    self.enrich_cmds() # assign default values, etc.
    return 0

  def find_beginning(self, src, p, aggr = 1):
    '''
    search the block-beginning marker {
    aggressive search: repeat until found
    position 'p' is move to the end of file
      if nothing is found
    limited to single-line case
    '''
    while 1:
      line = p.getline(src)
      if line == None: break

      pattern = r"\s*(typedef)\s+struct\s*\w*\s*\{"
      m = re.search(pattern, line)
      if m: # a match is found
        self.begin = P(p.row, p.col + m.start(1))
        p.col += m.end(0)
        #print "object-beginning is found in %s, %s" % (p, src[p.row])
        return 1
      if aggr:
        p.nextline()
      else:
        break
    return 0  # not found
  
  def search_ending_comment(self, src, p):
    cmt = CComment(src, p, 2) # search a nearby comment
    if cmt.isempty(): 
      self.cmds = Commands("")
      return
    # print "found a comment %s after object %s" % (cmt.raw, self.name)
    self.cmt = cmt
    if cmt: self.end = copy(p)
    self.cmds = Commands(cmt.raw)

  def find_ending(self, src, p, aggr = 0):
    '''
    search the block-ending mark
    passive search: since we assume aggressive search has been
    performed for items, we only search the current line
    position 'p' is unchanged, if nothing is found
    '''
    while 1:
      line = p.getline(src)
      if line.strip() == "": # blank line 
        p.nextline()
        continue
      pattern = r'\s*(\}\s*)(\w+)(\s*;).*'
      m = re.match(pattern, line)
      if m:
        p.col += m.end(3)
        self.name = m.group(2)
        self.end = copy(p)
        self.empty = 0
        #print "object-ending is found in %s, %s" % (p, src[p.row])
        self.search_ending_comment(src, p)
        return 1
      if aggr:
        p.nextline()
      else:
        break
    print "no object-ending is found, p at", p
    return 0  # not found
  
  def __str__(self): return self.name
  def __repr__(self): return self.name
  def isempty(self): return self.empty

  def merge_comments(self):
    ''' 
    merge multiple comments 
    since we merge comments before expand multiple declarators,
    we should test it.dl instead of it.decl
    '''
    i = 1
    items = self.items
    while i < len(items):
      it = items[i]
      itp = items[i-1]
      if (it.cmt and not it.dl  # note, test dl instead of decl
          and i > 0 and itp.cmt  # stand-alone comment allowing another
          and itp.cmt.type == "line"
          and it.cmt.type == "line"
          and it.cmt.begin.col == itp.cmt.begin.col): # starting at the same column
        #print "merging commands from item %d (%s) and item %d (%s)" % (
        #    i-1, itp.cmt.begin, i, it.cmt.begin)
        #print "'%s'\n+\n'%s'" % (itp.cmt.raw, it.cmt.raw)
        #raw_input()
        itp.cmt.raw += ' ' + it.cmt.raw
        items = items[:i] + items[i+1: ] # create a new list
        continue
      i += 1
    self.items = items
    #print "%d %d" % (len(items), len(self.items))
    #raw_input()

  def expand_multidecl_items(self):
    ''' expand lines that declares multiple items, 
    e.g.,   int a, b, c;  '''
    its = []
    for it in self.items: 
      its += it.itlist
    self.items = its

  def fill_cmds(self):
    ''' 
    obtain cmds from comments and apply persistent commands
    persistent commands apply to this item and all items afterwards
    after this function it.cmds will not be None (at least 
    it is an empty set {})
    '''
    p_cmds = Commands("") # persistent commands
    p_cmds["fold"] = ""   # apply an empty fold
    for it in self.items:
      if it.cmt:
        cmds = Commands(it.cmt.raw)
      else:
        cmds = Commands("")
        if it.pre:
          cmds.addpre(it.pre.pif, it.pre.pelse, it.pre.pendif)
      #if "cnt" in p_cmds:
      #  print "cnt is persistent when analysing %s in %s" % (it, self)
      #  raw_input()

      # update global persistent commands
      for key in cmds.persist:
        persist = cmds.persist[key]
        if persist > 0:
          p_cmds[key] = cmds[key]
          # print "add global command [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
        elif key in p_cmds:
          pval0 = p_cmds[key]
          if persist < 0:
            #print "remove global command [%s:%s] for %s" % (key, cmds[key], it)
            #raw_input()
            del p_cmds[key]
          if persist == -2: # preprocessor #else
            assert (key == "#if")
            newval = "!(%s)" % pval0
            if cmds[key]: newval += " && (%s)" % cmds[key] 
            p_cmds[key] = newval
            cmds[key] = newval
            #print "new[%s] is %s" % (key, newval)
          
      # merge with global persistent commands
      for key in p_cmds:
        if key not in cmds:
          cmds.cmds[key] = p_cmds[key]
          # print "inherit [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
      it.cmds = cmds
 
  def var2item(self, nm):  # for debug use
    if nm.startswith("@"):
      nm = nm[1:]
    elif hasattr(self, 'ptrname'):
      ptrpfx = self.ptrname + "->"
      if nm.startswith(ptrpfx): # remove the prefix, if any
        nm = nm[len(ptrpfx):]

    for it in self.items:
      if it.decl and it.decl.name == nm: 
        return it
    else: return None

  def add_dummies(self):
    ''' 
    add a declaration field for dummy variables, 
    such as flags, and alternative input
    '''
    # first flag all variable with declaration as non-dummy
    for it in self.items: 
      it.isdummy = (not it.decl)

    # now address dummy variables
    for it in self.items:
      if not it.isdummy: continue
      
      # determine the reference variable
      if "flag" in it.cmds: 
        it.prepare_flag()
        varref = it.cmds["flag_var"]
      elif "altvar" in it.cmds:
        varref = it.cmds["altvar"]
      else: continue

      # now copy the declaration of the reference
      it.decl = deepcopy(self.var2item(varref).decl)
      it.isdummy = 1
      # print ("copied decl from %s for %s\n" % (varref, it))
      # raw_input()

  def get_gtype(self):
    ''' get generic type, it uses information from commands '''
    for it in self.items:
      if it.decl: 
        it.gtype = it.get_gtype()
  
  def init_folds(self, verbose = 0):
    ''' 
    initialize folds 
    '''
    # add a dictionary of Fold's, each of which is an object
    # the key the fold name
    folds = {}
    folds[''] = Fold(self.fprefix, '')
    for it in self.items:
      f = it.cmds["fold"]
      if f not in folds:
        folds[f] = Fold(self.fprefix, f)
      folds[f].additem(it)
    if verbose:
      print "%s has %d folds" % (self.name, len(folds))
    
    for f in folds:
      if verbose: print "fold %-8s has %3d items" % (f, len(folds[f]))
      folds[f].prereq = 1  # prereq is always met
      folds[f].valid  = 1  # always valid
      #for it in folds[f].items: print it.getraw()
      #raw_input()
    self.folds = folds

  def init_folds2(self):
    ''' stuff need to be done after @'s are substituted '''
    for it in self.items: # fill fold-specific prerequisite
      f = it.cmds["fold"]
      if len(f) == 0: continue
      pr = it.cmds["fold_prereq"]
      if pr: self.folds[f].prereq = pr
      vd = it.cmds["fold_valid"]
      if vd: self.folds[f].valid = vd
    
  def get_prefix(self):
    '''
    determine ptrname, variable name for a pointer to the object
    and fprefix, the name attached to functions of the object
    also find variable's name is parent
    '''
    cmds = self.cmds
    name = self.name
    if name.endswith("_t"): name = name[:-2]
    self.ptrname = cmds["ptrname"] if (cmds and "ptrname" in cmds) else name
    self.fprefix = cmds["fprefix"] if (cmds and "fprefix" in cmds) else name + "_"
    self.parent = "parent"
    for it in self.items:
      if it.decl and it.cmds["usr"] == "parent":
        self.parent = it.decl.name

  def enrich_cmds(self):
    ''' refine and enrich commands '''
    state = 0
    for it in self.items:
      it.fill_def()  # make sure we have default values
      it.fill_key()  # fill in default keys
      it.fill_dim()  # test array `cnt' and fill `dim'
      it.fill_io()
      it.fill_test()
      it.sub_prefix(self.ptrname, self.fprefix, 
          self.parent) # @ to $ptrname
    self.init_folds2()
 
  def sort_items(self, items, tag):
    ''' 
    sort items according to ordering tags 
    `tag' can be bin or cfg
    '''
    sprev = "%s_prev" % tag
    nitems = len(items)
    ideps = [None] * nitems
    # for each item build a dependency list
    for i in range(nitems):
      it = items[i]

      # parse dependency string
      vars = it.cmds[sprev]
      deps = [s.strip() for s in vars.split(",")] if vars else []

      idep = [] # build dependency list
      for dep in deps:
        itdep = self.var2item(dep)
        if not itdep:
          print "dependencies variable %s of %s is not found!" % (dep, it.decl.name)
          raise Exception
        if itdep not in items: continue
        idep += [items.index(itdep)]
      ideps[i] = idep
    idx = sortidx(ideps)
    #print '\n\n'.join(["%3d: %s" % (i, items[idx[i]].getraw()) for i in range(nitems)] )
    #print idx; raw_input()
    return [items[idx[i]] for i in range(nitems)]

  def gen_code(self):
    ''' 
    generate code for an object to decl, header, source 
    '''
    funclist = []
    funclist += self.gen_func_cfgopen2()
    funclist += [self.gen_func_close()]
    # fold-specific functions
    for f in self.folds:
      funclist += [self.gen_func_clear(f)]
      funclist += self.gen_func_binrw2(f, "r")
      funclist += self.gen_func_binrw2(f, "w")
    funclist += [self.gen_func_manifest()]
    if not disabled(self.cmds["mpi"]):
      funclist += [self.gen_func_initmpi()]

    decl = self.gen_decl().rstrip()
    desc = self.cmds["desc"]
    if desc and len(desc):
      decl += " /* " + self.cmds["desc"] + " */"
    decl += "\n\n"
    decl += self.gen_flags_def()
    self.decl = decl;

    source = ""
    header = ""
    for f in funclist: # add functions
      header += f[0]
      source += "\n" + f[1] 
    self.header = header
    self.source = source

  def gen_decl(self):
    '''
    write code for object declaration
    '''
    tab    = 4
    offset = 2
    cow = CCodeWriter()
    cow.addln("typedef struct {")
    block = [] # block list
    for it in self.items:
      block = it.declare_var(cow, block, tab, offset)
    cow.dump_block(block, tab, offset)
    block = []
    if not disabled(self.cmds["mpi"]):
      cow.addln("int mpi_rank;")
      cow.addln("#ifdef %s\nMPI_Comm mpi_comm;\n#endif\n", USE_MPI)
    cow.addln("} " + self.name + ";")
    
    return cow.gets()

  def gen_flags_def(self):
    tab = 4
    offset = 2
    cow = CCodeWriter()
    block = []
    for it in self.items:
      if "flag" in it.cmds:
        flag = it.cmds["flag"]
        if "flag_val" not in it.cmds:
          print "missing value for flag [%s]" % it.cmds["flag"]
          raise Exception
        bval = it.cmds["flag_val"]
        cmt = it.cmds["desc"]
        cmt = "/* %s */" % cmt if cmt else ""
        # print flag, bval, cmt; raw_input()
        block += [("#define", flag, bval, cmt)]
      else:
        cow.dump_block(block, tab, offset)
        block = []
    cow.dump_block(block, tab, offset)
    return cow.gets()

  def get_usr_vars(self, tag, f = None):
    ''' tag can be "cfg", "bin", "txt", "close" '''
    if tag not in ("cfg", "bin", "txt", "close"): raise Exception
    param_list = ""
    var_list = ""
    for it in self.items:
      if it.decl and "usr" in it.cmds:
        decl0 = it.decl.datatype
        decl1 = it.decl.raw
        name = it.decl.name
        usrval = it.cmds["usr"]
        usr_name = name
        # replace the name by the proper usr_name in the declarator
        if usrval == "parent" and tag == "cfg":
          pass
        elif usrval == tag or (tag == "cfg" and usrval == 1):
          if f and it not in self.folds[f].items:
            continue  # wrong fold
          if tag == "cfg": usr_name = "usr_%s" % name
          decl1 = re.sub(name, usr_name, decl1) 
        else:
          continue # irrelavent variables
          #print "strange: tag:%s, usr:%s," % (tag, usrval); raw_input()
        param_list += ", %s %s" % (decl0, decl1)
        var_list += ", %s" % usr_name
    return param_list, var_list

  def gen_func_cfgopen_low(obj):
    ''' 
    write a function that initialize an object from configuration file
    Note: it loads configuration parameters only, 
          it does not load previous data 
    '''
    objtp   = obj.name
    objptr  = obj.ptrname
    fread   = "%scfgopen_low"  % obj.fprefix
    fdesc   = '''initialize members of %s from configuration 
    file `cfg', or if unavailable, from default values ''' % objtp
    
    cfgdecl = "cfgdata_t *cfg"
    usrvars = obj.get_usr_vars("cfg")
    fdecl = "int *%s(%s *%s, %s%s)" % (fread, objtp, objptr, cfgdecl, usrvars[0])
    cow = CCodeWriter()
    cow.begin_function(fread, fdecl, fdesc)

    # read configuration file
    for it in obj.items:
      it.cfgget_var(cow, obj.ptrname)

    cow.end_function("return 0;", silence = ["cfg"])
    return cow.prototype, cow.function

  def gen_func_cfgopen(obj):
    ''' 
    write a function that returns an new object 
    initialized from configuration file
    Note: it loads configuration parameters only, 
          it does not load previous data 
    '''
    cfgopen = (obj.cmds["cfgopen"] != None)
    objtp   = obj.name
    objptr  = obj.ptrname
    fopen   = "%scfgopen"  % obj.fprefix
    fread   = fopen + "_low"
    fdesc   = '''return a pointer of an initialized %s
              if possible, initial values are taken from configuration
              file `cfg', otherwise default values are assumed ''' % objtp
    
    cfgdecl = "cfgdata_t *cfg"
    usrvars = obj.get_usr_vars("cfg")
    fdecl = "%s *%s(%s%s)" % (obj.name, fopen, 
            "const char *cfgname" if cfgopen else cfgdecl, 
            usrvars[0])
    cow = CCodeWriter()
    cow.begin_function(fopen, fdecl, fdesc)
    cow.declare_var("%s *%s" % (objtp, objptr), objptr)

    # open a configuration file
    if cfgopen:
      cow.add_comment("open configuration file")
      cow.declare_var("cfgdata_t *cfg", "cfg")
      scfg = "(cfg = cfgopen(cfgname)) == NULL"
      cow.die_if(scfg,
          "%s: cannot open config. file %%s." % objtp,
          "cfgname", onerr = "return NULL;")
      cow.addln()

    # allocate memory
    cow.add_comment("allocate memory for %s" % obj.name)
    sobj = "(%s = calloc(1, sizeof(*%s))) == NULL" % (objptr, objptr)
    cow.die_if(sobj, "no memory for %s (size %%u)" % obj.name,
        "(unsigned) sizeof(*%s)" % objptr);
    cow.addln()
   
    cow.add_comment("call low level function")
    s = "0 != %s(%s, cfg%s)" % (fread, objptr, usrvars[1])
    ermsg = objtp + ": error while reading configuration file"
    onerr = "free(%s);\nreturn NULL;" % objptr
    if cfgopen:
      cow.die_if(s, ermsg + " %s", "cfgname", onerr = onerr)
    else:
      cow.die_if(s, ermsg, onerr = onerr) 
    cow.addln()
    
    if cfgopen:
      cow.add_comment("close handle to configuration file")
      cow.addln("cfgclose(cfg);")
    cow.end_function("return %s;" % objptr)
    return cow.prototype, cow.function

  def gen_func_cfgopen2(self):
    list = [self.gen_func_cfgopen_low()]
    if not self.cmds["private"]:
      list += [self.gen_func_cfgopen()]
    return list

  def gen_func_binrw_low(self, f, rw):
    ''' 
    write a function for reading/writing data in binary form 
    low level functions have no overhead from file-opening, checking-bytes
    and thus can be recursively applied to member objects
    '''
    readwrite = "read" if rw == "r" else "write"
    fprefix = self.folds[f].fprefix
    cow = CCodeWriter()
    funcnm = "%sbin%s_low" % (fprefix, readwrite)
    usrs = self.get_usr_vars("bin", f)[0]
    fdecl = "int %s(%s *%s, FILE *fp, int ver%s%s)" % (
        funcnm, self.name, self.ptrname, 
        ", unsigned flags, int endn" if rw == "r" else "", usrs)
    title = self.name + (("/%s" % f) if len(f) else "")
    cow.begin_function(funcnm, fdecl, 
        "%s %s data as binary" % (readwrite, title))
    callclear = "%sclear(%s);" % (fprefix, self.ptrname)

    if rw == "r":
      #cow.declare_var("int err = 0")
      #cow.declare_var("int verify")
      #cow.addln("verify = !(flags & IO_NOVERIFY);")
      cow.add_comment("clear data before reading")
      cow.addln(callclear)
      cow.addln()
    
    bin_items = self.sort_items(self.folds[f].items, "bin")
    for it in bin_items:
      if it.pre: continue
      if it.cmds["fold"] != f: continue
      #print "%s" % it.cmds; raw_input()
     
      xvars = it.cmds["fold_bin_add"]
      if xvars: # additional vars to the fold variables
        xvl = [s.strip() for s in xvars.split(",")]
        for xv in xvl:
          xit = self.var2item(xv)
          if not xit: raise Exception
          xit.binrw_var(cow, xv, rw)
        continue
      if not it.decl or it.isdummy:
        call = it.cmds[rw + "b_call"]
        if call: cow.addln(call)
        continue
      if not it.cmds["io_bin"]: continue

      varname = (self.ptrname+"->" if (it.cmds["usr"] != "bintmp")
          else "") + it.decl.name
      it.binrw_var(cow, varname, rw)

    unused_vars = ["ver"] + (["flags"] if rw == "r" else [])
    cow.end_function("return 0;\nERR:\n%sreturn -1;" % (
      callclear+"\n" if rw == "r" else ""), # clear data 
      silence = unused_vars)
    return cow.prototype, cow.function

  def gen_func_binrw(self, f, rw):
    ''' 
    write a wrapper function for reading/writing data in binary form 
    '''
    readwrite = "read" if rw == "r" else "write"
    cow = CCodeWriter()
    funcnm = "%sbin%s" % (self.folds[f].fprefix, readwrite)
    funcnm_ = funcnm + "_low"
    usrs = self.get_usr_vars("bin", f)
    fdecl = "int %s(%s *%s, const char *fname, int %s%s)" % (
        funcnm, self.name, self.ptrname, 
        "ver" if rw == "w" else "*pver, unsigned flags", 
        usrs[0])
    cow.begin_function(funcnm, fdecl, 
        "%s %s data as binary" % (readwrite, 
          self.name + ("/"+f if len(f) else "") ) )

    self.folds[f].add_fold_tests(cow, funcnm, "0")

    cow.declare_var("FILE *fp")
    condf = '(fp = fopen(fname, "%sb")) == NULL' % rw
    cow.die_if(condf, r"cannot %s binary file [%%s]." % readwrite, 
      "fname", onerr = "goto ERR;")
    cow.addln()

    # add checking bytes
    if rw == "r": 
      cow.rb_checkbytes();
      cow.declare_var("int ver")
      cow.rb_var("ver", "int")
      cow.addln("*pver = ver;")
    else:
      cow.wb_checkbytes();
    cow.addln()
    #cow.declare_var("int verify")
    #cow.addln("verify = !(flags & IO_NOVERIFY)")

    # call the actual low-level function
    cow.add_comment("call low level %s function for members" % readwrite)
    cow.declare_var("int i")
    cow.addln("i = %s(%s, fp, ver%s%s);", funcnm_, self.ptrname, 
        ", flags, endn" if rw == "r" else "", usrs[1]);
    cow.end_function("fclose(fp);\nreturn i;\nERR:\nfclose(fp);\nreturn -1;");
    return cow.prototype, cow.function

  def gen_func_binrw2(self, f, rw):
    list = [self.gen_func_binrw_low(f, rw)]
    private =  self.cmds["private"]
    if not private:
      list += [self.gen_func_binrw(f, rw)]
    return list

  def gen_func_clear(self, f):
    ''' write a function for clearing data '''
    cow = CCodeWriter()
    funcnm = "%sclear" % (self.folds[f].fprefix)
    fdecl = "void %s(%s *%s)" % (funcnm, self.name, self.ptrname)
    cow.begin_function(funcnm, fdecl, 
        "clear %s data" % self.name + ("/" + f if len(f) else ""))

    # add fold prereq and validation
    self.folds[f].add_fold_tests(cow, funcnm, "")
    
    ## for f = "", we clear everything including contents in folds
    #items = self.folds[f].items if len(f) else self.items
    items = self.folds[f].items
    items = [it for it in items if it.cmds["clr"]]
    #print "fold: %s\nitems: %s" % (f, items)
    for it in items:
      it.clear_var(cow, self.ptrname)
    cow.end_function("")
    return cow.prototype, cow.function

  def gen_func_close(self):
    cow = CCodeWriter()
    usrs = self.get_usr_vars("close")

    macroname = "%sclose" % self.fprefix
    funcnm = macroname + "_low"
    macro = "#define %s(%s%s) { %s(%s%s); free(%s); %s = NULL; }" % (
      macroname, self.ptrname, usrs[1],
      funcnm, self.ptrname, usrs[1], 
      self.ptrname, self.ptrname)
    fdecl = "void %s(%s *%s%s)" % (funcnm, self.name, self.ptrname, usrs[0])
    cow.begin_function(funcnm, fdecl, "close a pointer to %s" % self.name, 
        None if self.cmds["private"] else macro)

    # compute the longest variable
    calc_len = lambda it: len(it.decl.name)+2+len(self.ptrname) if (
        it.decl and it.gtype in ("char *", "dynamic array") ) else 0
    maxwid = max(calc_len(it) for it in self.items) if len(self.items) else 0
    for it in self.items:
      it.close_var(cow, self.ptrname, maxwid)

    cow.addln("memset(%s, 0, sizeof(*%s));", self.ptrname, self.ptrname)
    cow.end_function("")
    return cow.prototype, cow.function

  def gen_func_manifest(self):
    ''' write a function to present current data '''
    cow = CCodeWriter()
    funcnm = "%smanifest" % (self.fprefix)
    fdecl = "void %s(%s *%s, FILE *fp)" % (funcnm, self.name, self.ptrname)
    cow.begin_function(funcnm, fdecl, "clear %s data" % self.name)

    for it in self.items:
      it.manifest_var(cow, self.ptrname)
    cow.end_function("")
    return cow.prototype, cow.function

  def gen_func_initmpi(self):
    '''
    write a function to initialize MPI 
    every node calls this function
    '''
    cow = CCodeWriter()
    funcnm = "%sinitmpi" % (self.fprefix)
    obj = self.name
    ptr = self.ptrname
    fdecl = "int %s(%s *%s, MPI_Comm comm)" % (funcnm, obj, ptr)
    cow.begin_function(funcnm, fdecl, 
        "initialize MPI for %s" % obj,
        pp = "#ifdef %s" % USE_MPI)

    cow.die_if(ptr+" == NULL", "null pointer %s to %s" % (ptr, obj), 
        onerr = "return -1;")
    cow.addln("MPI_Bcast(%s, sizeof(*%s), MPI_BYTE, %s, comm);\n",
        ptr, ptr, MASTERID);
    # assign rank and communicator
    # must be done after Bcast to avoid being overwritten 
    cow.addln("%s->mpi_comm = comm;", ptr)
    cow.addln("MPI_Comm_rank(%s->mpi_comm, &%s->mpi_rank);", ptr, ptr)

    for it in self.items:
      it.initmpi_var(cow, self.ptrname)
    cow.addln("return 0;")
    cow.end_function("")
    return cow.prototype, cow.function

class Parser:
  '''
  handle objects in a file
  '''
  def __init__(self, src):
    '''
    `src' is a list of lines 
    '''
    self.src = src
    self.change_use_mpi()
    self.parse_lines()

  def parse_lines(self):
    '''
    process the input template 
    '''
    objs = self.objs = []
    p = P()

    # search for objects
    while 1:
      obj = Object(self.src, p) # try to get an object
      if obj.empty == 1: break
      if obj.empty == -1: continue  # skip

      print "found a new object '%s' with %d items, from %s to %s" % (
          obj.name, len(obj.items), obj.begin, obj.end)
      objs += [obj]
    #print "I got %d objects" % (len(objs))

  def output(self):
    for obj in self.objs:
      obj.gen_code()

    # add header
    s = self.change_objs_decl()
    
    # append functions
    for obj in self.objs:
      s += obj.source
    return s

  def change_objs_decl(self):
    '''
    imbed header to the corresponding places in the template
    '''
    objs = self.objs
    if objs == []: return ""
    # combine all headers
    allheaders = ""
    for obj in objs:
      allheaders += obj.header
    
    src = self.src
    p0 = P(0,0)
    s = p0.gettext(src, objs[0].begin) # to the beginning of the first object
    nobjs = len(objs)
    for i in range(nobjs):
      obj  = objs[i]
      s += obj.decl
      if i == nobjs - 1: # add function headers after the last the object declaration
        s += allheaders
      pnext = objs[i+1].begin if i < nobjs - 1 else P(len(src), 0)
      s += obj.end.gettext(src, pnext) # this object's end to the next object's beginning
    return s

  def change_use_mpi(self):
    # find `#define USE_MPI xxx'
    for i in range(len(self.src)):
      line = self.src[i]
      m = re.match(r"\#define\s+USE_MPI\s+(.*)$", line)
      if not m: continue
      # found a match
      global USE_MPI
      USE_MPI = m.group(1).strip()
      print "USE_MPI is now %s" % (USE_MPI)
      self.src[i] = "/* %s */\n" % line.strip()
      idef = i
      break
    else:
      return

    # substitute every line with USE_MPI 
    for i in range(len(self.src)):
      if i == idef: continue
      line = self.src[i]
      if line.find("USE_MPI") < 0: continue
      self.src[i] = re.sub("USE_MPI", USE_MPI, line)
