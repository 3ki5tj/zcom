#!/usr/bin/env python
r'''
a code generator for serializing objects

It uses a template, locates "typedef struct" code
and try to generate io functions, 
e.g., read parameters from configuration file

Commands syntax
  $cmd = args;
or
  $cmd[;]  (lazy command for switches)
  + cmd contains characters, numbers, and _ or -
  + the operator = can be replaced by :, it can also be := for
    persistent commands
  + currently args cannot contain ; even within a string
    it should be replaced by \\;
  + $# ... ; means a comment commands
  + if the terminal ; is missing, args extend to the end of line

Commands of an object:
  * $skip, $skipme: skip this object, do not generate anything
  * $fprefix:       function prefix for the object
  * $ptrname:       a common variable used as pointers to the object

Commands of an item:
  * $def:     default value for a variable, 
              if the variable is a dynamic/static array, this is
              value for array members

  * $cnt:     if not 0, declare a pointer variable as a dynamic array,
              ineffective for non-pointer variables

  * $io:      declare I/O type, can be a combination of cbt,
              for configuration, binary and text respectively
              it can also be `none' or `all'

  * $key:           key to get from configuration file
  * $key_prefix:    prefix to be added when it is deduced from 
                    variable name, usually as a persistent command
  * $key_unprefix:  if variable name starts with $key_unprefix, that part 
                    is first removed before applying $key_prefix
  * key_must:       a critial key that must present in configuration file

  * $prereq:      a condition to be tested *before* reading a varible
                  from configuration file, but after assigning the 
                  given by $def; this is the default order because it
                  ensures a variable is assigned at a reasonable value
                  
  * $test_first:  if this is set, we do not assign the default value $def 
                  unless $prereq is true, useful for dummy variables
  * $valid:       a condition to be tested for validity *after* reading
                  a variable from configuration file

  * $obj:     a member object, a function needs to called to properly 
              initialize it
  * $objarr:  object array

  * $flag:    two fields, name of flag and its value, something like 
              $flag: ABC_FLAG  0x000010;

  * $usr:     the variable should be passed as a parameter, 
              not read from the configuration file

In a stand-alone comment

  * $assert:  asserts a condition and aborts the program if 
              the condition fails
  * $call:    means a direct translation of the argument to code

In arguments of the command, `@' has special meanings
  * @@:       this item
  * @var:     ptrname->var
  * @_func:   fprefix_func
  * @:        ptrname

TODO:
  * multiple-line support, e.g., typedef\nstruct
  * preprocessor
  * skip struct in comment or string
'''

import os, sys, re, filecmp
from copy    import copy, deepcopy 
from objpos  import P 
from objdcl  import CDeclaratorList
from objcmt  import CComment
from objpre  import CPreprocessor
from objccw  import CCodeWriter
from objcmd  import Commands
from objitm  import Item

class Fold:
  ''' a sub-object embedded inside an object '''
  def __init__(self, obj, name):
    # generate a name
    self.name = obj.name
    if len(name): self.name += "-%s" % name

    # generate a function prefix
    if len(name) and not name.endswith("_"):
      name += "_"
    self.fprefix = obj.fprefix + name
    self.items = []
  def additem(self, item):
    self.items += [item]
  def __len__(self):
    return len(self.items)


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
    self.get_prefix()  # determine fprefix and ptrname
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
    nitems = []
    for it in self.items:
      nitems += it.itlist
    self.items = nitems

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
      if not it.cmt:
        it.cmds = copy(p_cmds)
        it.cmds["desc"] = ""  # add an empty description
        continue
      #if "cnt" in p_cmds:
      #  print "cnt is persistent when analysing %s in %s" % (it, self)
      #  raw_input()

      cmds = Commands(it.cmt.raw)
      # update global persistent commands
      for key in cmds.persist:
        persist = cmds.persist[key]
        if persist > 0:
          p_cmds[key] = cmds[key]
          # print "add global command [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
        elif persist < 0 and key in p_cmds:
          # print "remove global command [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
          del p_cmds[key]
          
      # merge with global persistent commands
      for key in p_cmds:
        if key not in cmds:
          cmds.cmds[key] = p_cmds[key]
          # print "inherit [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
      it.cmds = cmds
 
  def getdeclfrom(obj, var):
    ''' 
    for dummy variable, e.g., flag or alternative input,
    copy declaration of var to me 
    '''
    if var == None:
      print "variable name is None"
      raise Exception
    # print "calling substitution for [%s]" % var; raw_input
    
    # we first strip away the prefix, if any
    if var.startswith("@"): # in case substitution
      var = var[1:]

    # search over items
    for it in obj.items:
      if (not it.isdummy and it.decl.name == var):
        return copy(it.decl)
    else:
      print "cannot determine the type of [%s]" % var
      raise Exception
    
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
      it.decl = self.getdeclfrom(varref)
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
    self.folds = {}
    for it in self.items:
      f = it.cmds["fold"]
      if f not in self.folds:
        self.folds[f] = Fold(self, f)
      self.folds[f].additem(it)
    print "%s has %d folds" % (self.name, len(self.folds))
    for f in self.folds:
      print "fold %-8s has %3d items" % (f, len(self.folds[f]))
    

  def get_prefix(self):
    '''
    determine ptrname, variable name for a pointer to the object
    and fprefix, the name attached to functions of the object
    '''
    cmds = self.cmds
    name = self.name
    if name.endswith("_t"): name = name[:-2]
    self.ptrname = cmds["ptrname"] if (cmds and "ptrname" in cmds) else name
    self.fprefix = cmds["fprefix"] if (cmds and "fprefix" in cmds) else name + "_"
    
  def enrich_cmds(self):
    ''' refine and enrich commands '''
    for it in self.items:
      it.fill_def()  # make sure we have default values
      it.fill_key()  # fill in default keys
      it.test_cnt()
      it.fill_io()
      it.fill_test()
      it.sub_prefix(self.ptrname, self.fprefix) # @ to $ptrname
  
  def gen_code(self):
    funclist = []
    funclist += [self.gen_func_cfgread()]
    funclist += [self.gen_func_cfgopen()]
    funclist += [self.gen_func_close()]
    # fold-specific functions
    for f in self.folds:
      if f == "": continue
      pass

    # assemble everything to header and source
    header = self.gen_decl().rstrip()
    desc = self.cmds["desc"]
    if desc and len(desc):
      header += " /* " + self.cmds["desc"] + " */"
    header += "\n\n"
    source = ""

    header += self.gen_flags_def()
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
        self.dump_block(block, cw)
        block = [] # empty the block
        cw.addln("#" + item.pre.raw)
        continue

      # 2. handle declaration
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
        self.dump_block(block, cw, tab, offset)
        block = [] # empty the block
        # also print this line
        if len(scmt) > 0: 
          cw.addln(scmt)
    self.dump_block(block, cw, tab, offset)
    block = []
    cw.addln("} " + self.name + ";")
    
    return cw.gets()

  def dump_block(self, block, cw, tab = 4, offset = 2):
    if len(block) == 0: return
    nfields = len(block[0])
    wid = [8] * nfields
    #print "block of %s\n%s" % (len(block), block)
    #raw_input()

    # align and format all items in the bufferred block
    # and append the result to string s
    for j in range(nfields):
      w = max(len(item[j]) for item in block) + 1 # +1 for the following blank
      wid[j] = ((w + offset + tab - 1) // tab)*tab - offset

    for item in block:
      s = ""
      for j in range(nfields-1):
        s += "%-*s" % (wid[j], item[j])
      s += item[nfields - 1] + "\n"
      cw.addln(s)

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
        self.dump_block(block, cw, tab, offset)
        block = []
    self.dump_block(block, cw, tab, offset)
    return cw.gets()

  def gen_func_close(self):
    ow = CCodeWriter()
    macroname = "%sclose" % self.fprefix
    fname = macroname + "_"
    macro = "#define %s(%s) { %s(%s); free(%s); %s = NULL; }" % (
      macroname, self.ptrname, fname, self.ptrname, self.ptrname, self.ptrname)
    fdecl = "void %s(%s *%s)" % (fname, self.name, self.ptrname)
    ow.begin_function(fname, fdecl, "close a pointer to %s" % self.name, macro)

    # compute the longest variable
    calc_len = lambda it: len(it.decl.name) if (it.decl and 
        it.gtype in ("char *", "dynamic array") ) else 0
    maxwid = max(calc_len(it) for it in self.items)
    for it in self.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if not it.decl or it.isdummy: continue

      # print ("destroying %s, gtype: %s" % (it.decl.name, it.gtype)); 
      # raw_input()
      if it.gtype == "char *": 
        funcfree = "ssdelete"
      elif it.gtype == "dynamic array":
        funcfree = "free"
      else: continue
      
      ow.addln("if (%s->%-*s != NULL) %s(%s->%s);",
               self.ptrname, maxwid, it.decl.name,
               funcfree, self.ptrname, it.decl.name)

    ow.addln("memset(%s, 0, sizeof(*%s));", self.ptrname, self.ptrname)
    ow.end_function("")
    return ow.prototype, ow.function

  def get_usr_vars(self):
    param_list = ""
    var_list = ""
    for it in self.items:
      if it.decl and "usr" in it.cmds:
        decl0 = it.decl.datatype
        decl1 = it.decl.raw
        name = it.decl.name
        usr_name = "usr_%s" % name
        decl1 = re.sub(name, usr_name, decl1)
        param_list += ", %s %s" % (decl0, decl1)
        var_list += ", %s" % decl1
    return param_list, var_list

  def gen_func_cfgopen(obj):
    ''' 
    write a function that returns an new object 
    initialized from configuration file
    '''
    objtp  = obj.name
    objptr = obj.ptrname
    fopen  = "%scfgopen" % obj.fprefix
    fread  = "%scfgread" % obj.fprefix
    fdesc  = "return an initialized %s, wrapper of %s" % (objtp, fread)
    usrvars = obj.get_usr_vars()
    fdecl = "%s *%s(cfgdata_t *cfg%s)" % (obj.name, fopen, usrvars[0])
    ow = CCodeWriter()
    ow.begin_function(fopen, fdecl, fdesc)
    ow.vars.addln("%s %s;" % (objtp, objptr))

    # allocate memory
    s = "(%s = calloc(1, sizeof(*%s))) == NULL" % (objptr, objptr)
    ow.die_if (s, "no memory for %s" % obj.name);
    ow.addln("");
    
    # read configuration 
    s = "0 != %s(%s, cfg%s)" % (fread, objptr, usrvars[1])
    ow.begin_if (s)
    ow.addln(r'fprintf(stderr, "failed to open %s\n");', objtp);
    ow.addln("return NULL;")
    ow.end_if (s)
    ow.end_function("return %s;" % objptr)
    return ow.prototype, ow.function
 
  def gen_func_cfgread(obj):
    ''' write a function that reads configuration file '''
    fname = "%scfgread" % obj.fprefix
    fdesc = '''initialize %s by reading members from configuration file `cfg'
               if `cfg' is NULL, default values are assigned ''' % obj.name
    fdecl = "int %s(%s *%s, cfgdata_t *cfg%s)" % (
        fname, obj.name, obj.ptrname, obj.get_usr_vars()[0])
    ow = CCodeWriter()
    ow.begin_function(fname, fdecl, fdesc)

    for it in obj.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue

      defl    = it.cmds["def"]
      key     = it.cmds["key"]
      prereq  = it.cmds["prereq"]
      tfirst  = it.cmds["test_first"]
      valid   = it.cmds["valid"]
      must    = it.cmds["key_must"]
      desc    = it.cmds["desc"]

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

      if it.gtype == "pointer to object":
        # TODO: call the appropriate initializer
        '''
        obj = search_the_object_name (it.decl.datatype)
        '''
        obj_ptrname = it.decl.datatype[:-2]
        ow.addln("%s = %sinit(%s);", varname, 
          obj_ptrname, obj.ptrname)
      
      elif "flag" in it.cmds and "key" in it.cmds: # input flag
        flag = it.get_flag()
        fmt = it.type2fmt(it.gtype)
        ow.cfgget_flag(varname, key, flag, it.gtype, fmt,
            defl, prereq, desc)
      
      elif it.gtype == "dynamic array":
        gtype = it.get_gtype(offset = 1) 
        ow.init_dynamic_array(varname, gtype,
            defl, it.cmds["cnt"], desc)
      
      elif it.gtype == "static array":
        gtype = it.get_gtype(offset = 1)
        # print "%s: static array of %s" % (varnm, it.gtype)
        ow.init_static_array(varname, gtype, 
            defl, it.cmds["cnt"], desc)
      
      else: # regular variable
        if "usr" in it.cmds: # usr variables
          ow.assign(varname, "usr_%s" % varnm, it.gtype)
        elif not it.cmds["io_cfg"]: # assign default value
          if len(defl): 
            ow.assign(varname, defl, it.gtype);
        else:
          fmt = it.type2fmt(it.gtype)
          ow.cfgget_var(varname, key, it.gtype, fmt, 
            defl, must, prereq, tfirst, valid, desc)

    ow.end_function("return 0;")
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
    src = self.src
    p0 = P(0,0)
    s = p0.gettext(src, objs[0].begin) # to the beginning of the first object
    nobjs = len(objs)
    for i in range(nobjs):
      obj  = objs[i]
      s += obj.header
      pnext = objs[i+1].begin if i < nobjs - 1 else P(len(src), 0)
      s += obj.end.gettext(src, pnext) # this object's end to the next object's beginning
    return s
   

def handle(file, template = ""):
  ''' generate a C file according to a template '''
  # guess the template name if it's empty
  if template == "":
    pt = os.path.splitext(file)
    template = pt[0] + '.0' + pt[1]

  # read in the template
  src = open(template, 'r').readlines()
  # construct a parser
  psr = Parser(src)
  # write the processed file
  open(file, 'w').write(psr.output())

def main():
  files = ("at.c", "mb.c")
  if len(sys.argv) > 1: files = (sys.argv[1])
  for file in files:
    handle(file)

if __name__ == "__main__":
  main()

