#!/usr/bin/env python
r'''
a code generator for serializing objects

It uses a template, locates "typedef struct" code
and try to generate io functions, 
e.g., read parameters from configuration file

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

  * $key:           key to get from configuration file
  * $key_prefix:    prefix to be added when it is deduced from 
                    variable name, usually as a persistent command
  * $key_unprefix:  if variable name starts with $key_unprefix, that part 
                    is first removed before applying $key_prefix
  * $key_args:      printf arguments in constructing key
  
  * $must:       a critial key that must present in configuration file

  * $prereq:      prerequisite to be tested *before* reading a variable
                  from configuration file, but after assigning the 
                  given by $def; this is the default order because it
                  ensures a variable is assigned at a reasonable value
                  
  * $test_first:  if this is set, we do not assign the default value $def 
                  unless $prereq is true, useful for dummy variables
  * $valid:       a condition to be tested for validity *after* reading
                  a variable from configuration file

  * $obj:     a member object, a function needs to called to properly 
              initialize it; 
              for an object array, also set $cnt;

  * $flag:    two fields, name of flag and its value, something like 
              $flag: ABC_FLAG  0x000010;

  * $usr:     the variable should be passed as a parameter, 
              not read from the configuration file
              $usr:parent; declare the variable represents 
              a pointer to the parent object

  * $fold:            follow by a fold-identifier
  * $fold_bin_add:    add some variables (separated by ,) during writing binary
  * $fold_bin_prereq: quit writing binary if this condition is true
  * $fold_bin_valid:  test if condition is true

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
from objccw  import CCodeWriter
from objcmd  import Commands
from objitm  import Item
from objext  import *

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
            print "new[%s] is %s" % (key, newval)
          
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
  
  def init_folds(self):
    ''' 
    initialize folds 
    '''
    # add a dictionary of Fold's, each of which is an object
    # the key the fold name
    folds = {}
    for it in self.items:
      f = it.cmds["fold"]
      if f not in folds:
        folds[f] = Fold(self.fprefix, f)
      folds[f].additem(it)
    print "%s has %d folds" % (self.name, len(folds))
    
    for f in folds:
      print "fold %-8s has %3d items" % (f, len(folds[f]))
      folds[f].prereq = 1
      folds[f].valid  = 1
      #for it in folds[f].items: print it.getraw()
      #raw_input()
    self.folds = folds

  def init_folds2(self):
    ''' stuff need to be done after @'s are substituted '''
    for it in self.items: # fill fold-specific prerequisite
      f = it.cmds["fold"]
      if len(f) == 0: continue
      pr = it.cmds["fold_bin_prereq"]
      if pr: self.folds[f].prereq = pr
      vd = it.cmds["fold_bin_valid"]
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
    funclist += [self.gen_func_cfgopen_()]
    funclist += [self.gen_func_cfgopen()]
    funclist += [self.gen_func_close()]
    funclist += [self.gen_func_binwrite_("")]
    funclist += [self.gen_func_binwrite("")]
    # fold-specific functions
    for f in self.folds:
      if f == "": continue
      funclist += [self.gen_func_binwrite_(f)]
      funclist += [self.gen_func_binwrite(f)]

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
    cw = CCodeWriter()
    cw.addln("typedef struct {")
    block = [] # block list
    for i in range(len(self.items)):
      item = self.items[i]
      decl0 = decl1 = scmt = ""
      # 1. handle pre
      if item.pre:
        dump_block(block, cw)
        block = [] # empty the block
        cw.addln("#" + item.pre.raw)
        continue
      
      # 2. handle declaration
      if item.cmds["usr"] not in (None, 1, "cfg"):
        continue
      if item.decl:
        decl0 = item.decl.datatype
        if item.isdummy: continue
        decl1 = item.decl.raw
      else:
        # skip a comment if it contains instructions
        if item.cmds["assert"]: continue
        if item.cmds["call"]:   continue
        if item.cmds["flag"]:   continue # skip flags
        #raw_input ("%s %s" % (decl0, decl1))

      # 3. handle comment
      desc = item.cmds["desc"]
      if len(desc) > 0:
        scmt = "/* " + desc + " */"
      
      #print "decl: %s, cmds: %s" % (item.decl, item.cmds)
      #raw_input()

      # 4. output the content
      if len(decl0) > 0:  # put it to a code block
        block += [(decl0, decl1+';', scmt)] # just buffer it
      else:
        dump_block(block, cw, tab, offset)
        block = [] # empty the block
        # also print this line
        if len(scmt) > 0: 
          cw.addln(scmt)
    dump_block(block, cw, tab, offset)
    block = []
    cw.addln("} " + self.name + ";")
    
    return cw.gets()

  def gen_flags_def(self):
    tab = 4
    offset = 2
    cw = CCodeWriter()
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
        dump_block(block, cw, tab, offset)
        block = []
    dump_block(block, cw, tab, offset)
    return cw.gets()

  def gen_func_close(self):
    ow = CCodeWriter()
    macroname = "%sclose" % self.fprefix
    funcnm = macroname + "_"
    macro = "#define %s(%s) { %s(%s); free(%s); %s = NULL; }" % (
      macroname, self.ptrname, funcnm, self.ptrname, self.ptrname, self.ptrname)
    fdecl = "void %s(%s *%s)" % (funcnm, self.name, self.ptrname)
    ow.begin_function(funcnm, fdecl, "close a pointer to %s" % self.name, macro)

    # compute the longest variable
    calc_len = lambda it: len(it.decl.name)+2+len(self.ptrname) if (
        it.decl and it.gtype in ("char *", "dynamic array") ) else 0
    maxwid = max(calc_len(it) for it in self.items)
    for it in self.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if not it.decl or it.isdummy: continue
      if it.cmds["usr"] == "parent": continue

      varnm = it.decl.name;
      varname = "%s->%s" % (self.ptrname, varnm)
      
      #if varnm == "cache":
      #  print ("destroying %s, gtype: %s" % (it.decl.name, it.gtype)); 
      #  raw_input()
      if it.gtype == "char *": 
        funcfree = "ssdelete"
      elif it.gtype == "dynamic array":
        funcfree = "free"
      elif it.gtype == "object pointer":
        fpfx = it.get_obj_fprefix()
        funcfree = "%sclose" % fpfx
      elif it.gtype == "object array": 
        fpfx = it.get_obj_fprefix()
        funcfree = "%sclose_" % fpfx
        cnt = it.cmds["cnt"]
        ow.declare_var("int i")

        cond = "%s != NULL" % varname
        ow.begin_if(cond)
        ow.addln("for (i = 0; i < %s; i++) {" % cnt)
        ow.addln("%s(%s + i);", funcfree, varname) 
        ow.addln("}");
        ow.addln("free(%s);", varname);
        ow.end_if(cond)
        continue
      else:
        continue
  
      ow.addln("if (%-*s != NULL) %s(%s);",
               maxwid, varname, funcfree, varname)

    ow.addln("memset(%s, 0, sizeof(*%s));", self.ptrname, self.ptrname)
    ow.end_function("")
    return ow.prototype, ow.function

  def get_usr_vars(self, tag, f = None):
    ''' tag can be "cfg", "bin", "txt" '''
    if tag not in ("cfg", "bin", "txt"): raise Exception
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

  def gen_func_cfgopen_(obj):
    ''' 
    write a function that initialize an object from configuration file
    Note: it loads configuration parameters only, 
          it does not load previous data 
    '''
    objtp   = obj.name
    objptr  = obj.ptrname
    fread   = "%scfgopen_"  % obj.fprefix
    fdesc   = '''initialize members of %s from configuration 
    file `cfg', or if unavailable, from default values ''' % objtp
    
    cfgdecl = "cfgdata_t *cfg"
    usrvars = obj.get_usr_vars("cfg")
    fdecl = "int *%s(%s *%s, %s%s)" % (fread, objtp, objptr, cfgdecl, usrvars[0])
    ow = CCodeWriter()
    ow.begin_function(fread, fdecl, fdesc)

    # read configuration file
    for it in obj.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if it.cmds["usr"] == "parent":
        continue

      defl    = it.cmds["def"]
      key     = it.cmds["key"]
      prereq  = it.cmds["prereq"]
      tfirst  = it.cmds["test_first"]
      valid   = it.cmds["valid"]
      must    = it.cmds["must"]
      desc    = it.cmds["desc"]
      cnt     = it.cmds["cnt"]

      if not it.decl: # stand-alone comment
        ow.add_comment(desc)
        ow.insist(it.cmds["assert"], desc)
        call = it.cmds["call"]
        if call: ow.addln(call + ";")
        continue
      
      varnm = it.decl.name
      varname = obj.ptrname + "->" + varnm

      # add a remark before we start
      ow.add_comment(desc)

      if "flag" in it.cmds and "key" in it.cmds: # input flag
        flag = it.get_flag()
        fmt = it.type2fmt(it.gtype)
        ow.cfgget_flag(varname, key, flag, it.gtype, fmt,
            defl, prereq, desc)
      
      elif it.gtype == "object pointer":
        fpfx = it.get_obj_fprefix()
        s = "(%s = %scfgopen(cfg%s)) == NULL" % (varname, 
          fpfx, it.get_init_args())
        ow.die_if(s, r"failed to initialize %s\n" % varname)

      elif it.gtype == "object array":
        fpfx = it.get_obj_fprefix()
        etp = it.get_gtype(offset = 1)
        ow.alloc_darr(varname, etp, cnt)
        ow.declare_var("int i;")
        ow.addln("for (i = 0; i < %s; i++) {" % cnt)
        s = "0 != %scfgopen_(%s+i, cfg%s)" % (
          fpfx, varname, it.get_init_args())
        ow.die_if(s, r"failed to initialize %s[%%d]\n" % varname, "i")
        ow.addln("}\n")

      elif it.gtype == "dynamic array":
        ow.init_darr(varname, it.get_gtype(offset = 1),
            defl, cnt, desc)
      
      elif it.gtype == "static array":
        etp = it.get_gtype(offset = 1)
        # print "%s: static array of element = [%s] io = %s, defl = [%s]" % (varnm, etp, it.cmds["io_cfg"], defl); raw_input()
        if not it.cmds["io_cfg"]: 
          ow.init_sarr(varname, etp, defl, cnt, desc)
        else:
          fmt  = it.type2fmt(etp)
          cmpl = it.cmds["complete"]
          ow.cfgget_sarr(varname, key, etp, fmt, defl, cnt, 
              must, valid, cmpl, desc)

      else: # regular variable
        usrval = it.cmds["usr"]
        if usrval:
          if usrval in ("cfg", 1): # usr variables
            ow.assign(varname, "usr_%s" % varnm, it.gtype)
          else: pass # ignore other usr variables
        elif not it.cmds["io_cfg"]: # assign default value
          if len(defl): 
            ow.assign(varname, defl, it.gtype);
        else:
          fmt = it.type2fmt(it.gtype)
          ow.cfgget_var(varname, key, it.gtype, fmt, 
            defl, must, prereq, tfirst, valid, desc)
    
    ow.end_function("return 0;")
    return ow.prototype, ow.function

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
    fread   = fopen + "_"
    fdesc   = '''return a pointer of an initialized %s
              if possible, initial values are taken from configuration
              file `cfg', otherwise default values are assumed ''' % objtp
    
    cfgdecl = "cfgdata_t *cfg"
    usrvars = obj.get_usr_vars("cfg")
    fdecl = "%s *%s(%s%s)" % (obj.name, fopen, 
            "const char *fname" if cfgopen else cfgdecl, 
            usrvars[0])
    ow = CCodeWriter()
    ow.begin_function(fopen, fdecl, fdesc)
    ow.vars.addln("%s *%s;" % (objtp, objptr))

    # open a configuration file
    if cfgopen:
      ow.vars.addln("cfgdata_t *cfg;")
      scfg = "(cfg = cfgopen(fname)) == NULL"
      ow.begin_if(scfg)
      ow.errmsg(r"%s: cannot open config. file %%s.\n" % objtp,
          "fname");
      ow.addln("return NULL;")
      ow.end_if(scfg)

    # allocate memory
    sobj = "(%s = calloc(1, sizeof(*%s))) == NULL" % (objptr, objptr)
    ow.die_if(sobj, "no memory for %s" % obj.name);
   
    # ow.declare_var("int ret")
    s = "0 != %s(%s, cfg%s)" % (fread, objptr, usrvars[1])
    ow.begin_if(s)
    ow.errmsg("%s: failed to read config." % objtp)
    ow.addln("free(%s);" % objptr)
    ow.addln("return NULL;")
    ow.end_if(s)
    
    if cfgopen:
      ow.addln("cfgclose(cfg);")
    ow.end_function("return %s;" % objptr)
    return ow.prototype, ow.function

  def gen_func_binwrite_(self, f):
    ''' write a function of writing data in binary form '''
    ow = CCodeWriter()
    funcnm = "%sbinwrite_" % self.folds[f].fprefix
    fdecl = "int %s(%s *%s, FILE *fp, int ver)" % (
        funcnm, self.name, self.ptrname)
    title = self.name + (("/%s" % f) if len(f) else "")
    ow.begin_function(funcnm, fdecl, 
        "write %s data as binary" % title)
    
    usrs = self.get_usr_vars("bin", f)[0]
    if len(usrs):
      usrs = [s.strip() for s in usrs.lstrip(", ").split(",")]
      for usr in usrs:
        ow.declare_var(usr)
        usrtp, usrvar = usr.split()
        uit = self.var2item(usrvar)
        defl = uit.cmds["def"]
        if defl:
          ow.assign("usrvar", "defl", usrtp)

    bin_items = self.sort_items(self.folds[f].items, "bin")
    for it in bin_items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if it.cmds["fold"] != f: continue
      #print "%s" % it.cmds; raw_input()
     
      xvars = it.cmds["fold_bin_add"]
      if xvars: # additional vars to the fold variables
        xvl = [s.strip() for s in xvars.split(",")]
        for xv in xvl:
          xit = self.var2item(xv)
          if not xit: raise Exception
          xit.binwrite_var(ow, xv)
        continue
      if not it.decl or it.isdummy: continue
      if not it.cmds["io_bin"]: continue

      varname = "%s->%s" % (self.ptrname, it.decl.name)
      it.binwrite_var(ow, varname)

    ow.addln("return 0;")
    ow.addln("ERR:")
    ow.addln("return -1;")
    ow.end_function("")
    return ow.prototype, ow.function

  def gen_func_binwrite(self, f):
    ''' write a function of writing data in binary form '''
    ow = CCodeWriter()
    funcnm = "%sbinwrite" % self.folds[f].fprefix
    funcnm_ = funcnm + "_"
    fdecl = "int %s(%s *%s, const char *fname, int ver)" % (
        funcnm, self.name, self.ptrname)
    ow.begin_function(funcnm, fdecl, 
        "write %s/%s data as binary" % (self.name, f))

    # add prerequisite/validation tests
    if self.folds[f].prereq != 1:
      ow.addlnraw("if (!(%s)) return 0;" % self.folds[f].prereq)
    if self.folds[f].valid != 1:
      cond = "%s" % self.folds[f].valid
      ow.insist(cond, "validation in write binary")

    ow.declare_var("FILE *fp")
    condf = '(fp = fopen(fname, "wb")) == NULL'
    ow.begin_if(condf)
    ow.addlnraw(r'printf("cannot write binary history file [%s].\n", fname);')
    ow.addln("return -1;")
    ow.end_if(condf)

    # add checking bytes
    ow.wb_checkbytes();
    ow.wb_var("ver", "int")
    ow.declare_var("int i")
    ow.addln("i = %s(%s, fp, ver);", funcnm_, self.ptrname);
    ow.addln("fclose(fp);")
    ow.addln("return i;")
    ow.addln("ERR:")
    ow.addln("fclose(fp);")
    ow.addln("return -1;")
    ow.end_function("")
    return ow.prototype, ow.function



class Parser:
  '''
  handle objects in a file
  '''
  def __init__(self, src):
    self.src = src
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

      print "found a new object '%s' with %d items, from %s to %s\n" % (
          obj.name, len(obj.items), obj.begin, obj.end)
      objs += [obj]
    print "I got %d objects" % (len(objs))

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

def handle(file, template = ""):
  ''' generate a C file according to a template '''
  # guess the template name if it's empty
  if template == "":
    pt = os.path.splitext(file)
    template = pt[0] + '.0' + pt[1]
    ref = pt[0] + '1' + pt[1]

  # read in the template
  src = open(template, 'r').readlines()
  # construct a parser
  psr = Parser(src)
  # write the processed file
  open(file, 'w').write(psr.output())

  if os.path.exists(ref):
    if not filecmp.cmp(file, ref):
      print "File %s is different from %s" % (file, ref)
      raw_input()

def main():
  files = ["spb.c", "at.c", "mb.c"]
  #files = ["mb.c"]
  if len(sys.argv) > 1: files = [sys.argv[1]]
  for file in files:
    handle(file)

if __name__ == "__main__":
  try:
    import psyco
    psyco.full()
  except ImportError: pass
  main()

