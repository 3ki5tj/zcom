#!/usr/bin/env python
'''
a code generator for serializing objects

It uses a template, locates "typedef struct" code
and try to generate io functions, 
e.g., read parameters from configuration file

Commands syntax
  $cmd = args;
or
  $cmd[;]  (lazy command)
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

  * $key:               key to get from configuration file
  * $key_prefix:        prefix to be added when it is deduced from 
                        variable name, usually as a persistent command
  * $key_prefix_minus:  prefix of the variable name to be removed 
                        before applying $key_prefix

  * $test:      a condition to be tested for validity *after* reading
                a variable from configuration file
  * $pre_test:  a condition to be tested *before* reading a varible
                from configuration file

  * $obj:     a member object, a function needs to called to properly initialize it
  * $objarr:  object array


TODO:
  * multiple-line support, e.g., typedef\nstruct
  * preprocessor
  * skip struct in comment or string
'''

import os, sys, re, filecmp
from copy    import * 
from objpos  import P 
from objdcl  import CDeclaratorList
from objcmt  import CComment
from objpre  import CPreprocessor
from objccw  import CCodeWriter
from objcmd  import Commands
from objitm  import Item

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
    self.get_gtype()   # generic type, needs cmds, must be after get_cmds
    self.get_prefix()  # default ptrname and fprefix
    self.sub_ptrname() # substitute @ by ptrname
    self.enrich_cmds() # assign default values
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
      self.cmds = None
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
    obtain cmds from comments
    persist commands apply to this item and all items afterwards
    after this function it.cmds will not be None (at least 
    it is an empty set {})
    '''
    p_cmds = {} # persistent commands
    for it in self.items:
      if not it.cmt:
        it.cmds = copy(p_cmds)
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
          del p_cmds[key]
          # print "remove global command [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
          
      # merge with global persistent commands
      for key in p_cmds:
        if key not in cmds:
          cmds.cmds[key] = p_cmds[key]
          # print "inherit [%s:%s] for %s" % (key, cmds[key], it)
          # raw_input()
      it.cmds = cmds

  def get_gtype(self):
    ''' get generic type, it uses information from commands '''
    for it in self.items:
      if it.decl: 
        it.gtype = it.get_generic_type()
  
  def get_prefix(self):
    '''
    determine ptrname, variable name for a pointer to the object
    and fprefix, the name attached to functions of the object
    '''
    cmds = self.cmds
    name = self.name
    if name.endswith("_t"): name = name[:-2]
    self.fprefix = cmds["fprefix"] if (cmds and "fprefix" in cmds) else name + "_"
    self.ptrname = cmds["ptrname"] if (cmds and "ptrname" in cmds) else name

  def sub_ptrname(self):
    '''
    substitute @ by $ptrname
    '''
    for it in self.items:
      if not it.cmds: continue
      it.sub_ptrname(self.ptrname)
    
  def enrich_cmds(self):
    ''' refine and enrich commands '''
    for it in self.items:
      it.fill_def()  # make sure we have default values
      it.fill_key()  # fill in default keys
      it.test_cnt()
      it.fill_io()
  
  def gen_code(self):
    funclist = []
    funclist += [self.gen_func_cfgread()]
    funclist += [self.gen_func_close()]

    # assemble everything to header and source
    header = self.gen_decl().rstrip()
    if self.cmds and "desc" in self.cmds:
      header += " /* " + self.cmds["desc"] + " */"
    header += "\n\n"
    source = ""
    for f in funclist: # add functions
      header += f[0]
      source += "\n" + f[1] 
    self.header = header
    self.source = source

  def gen_decl(self):
    '''
    write code for object declaration
    '''
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
        decl1 = item.decl.raw
        #raw_input ("%s %s" % (decl0, decl1))
      # 3. handle comment
      cmds = item.cmds
      if (cmds and "desc" in cmds
          and len(cmds["desc"]) > 0):
        scmt = "/* " + cmds["desc"] + " */"
      # 4. output the content
      if len(decl0) > 0:  # put it to a code block
        block += [(decl0, decl1, scmt)] # just buffer it
      else:
        self.dump_block(block, cw)
        block = [] # empty the block
        # also print this line
        if len(scmt) > 0: 
          cw.addln(scmt)
    self.dump_block(block, cw)
    block = []
    cw.addln("} " + self.name + ";")
    return cw.gets()
  def dump_block(self, block, cw):
    if len(block) == 0: return
    # align and format all items in the bufferred block
    # and append the result to string s
    wid0 = max(len(item[0]) for item in block)
    wid1 = max(len(item[1]) for item in block) + 1 # for ;
    for item in block:
      line = "%-*s %-*s %s\n" % (wid0, item[0], wid1, item[1]+';', item[2])
      cw.add(line)

  def gen_func_close(self):
    ow = CCodeWriter()
    owh = CCodeWriter()
    macroname = "%sclose" % self.fprefix
    funcname = macroname + "_"
    owh.addln("#define %s(%s) { %s(%s); free(%s); %s = NULL; }",
      macroname, self.ptrname, funcname, self.ptrname, self.ptrname, self.ptrname)
    funcdesc = "/* %s: close a pointer to %s */" % (funcname, self.name)
    #owh.addln(funcdesc)
    ow.addln(funcdesc)
    fdecl = "void %s(%s *%s)" % (funcname, self.name, self.ptrname)
    owh.addln(fdecl + ";") # prototype
    ow.addln(fdecl)
    ow.addln("{")

    # compute the longest variable
    calc_len = lambda it: len(it.decl.name) if (it.decl and 
        it.gtype in ("char *", "dynamic array") ) else 0
    maxwid = max(calc_len(it) for it in self.items)
    for it in self.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if not it.decl: continue

      funcfree = None
      if it.gtype == "char *": 
        funcfree = "ssdelete"
      elif it.gtype == "dynamic array":
        funcfree = "free"
      if not funcfree: continue
      
      ow.addln("if (%s->%-*s != NULL) %s(%s->%s);",
               self.ptrname, maxwid, it.decl.name,
               funcfree, self.ptrname, it.decl.name)

    ow.addln("memset(%s, 0, sizeof(*%s));", self.ptrname, self.ptrname)
    ow.addln("}")
    ow.addln("")
    return owh.gets(), ow.gets()

  def gen_func_cfgread(self):
    ow = CCodeWriter()
    owh = CCodeWriter()
    funcname = "%scfgread" % self.fprefix
    funcdesc = "/* %s: initialize %s by reading members from a configuration file */" % (
        funcname, self.name)
    #owh.addln(funcdesc)
    ow.addln(funcdesc)
    fdecl = "void %s(%s *%s, cfgdata_t *cfg)" % (
        funcname, self.name, self.ptrname)
    owh.addln(fdecl + ";") # prototype
    ow.addln(fdecl)
    ow.addln("{")

    for it in self.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if not it.decl or not it.cmds["io_cfg"]: continue

      if it.gtype == "pointer to object":
        # TODO: call the appropriate initialize
        '''
        obj = search_the_object_name (it.decl.datatype)
        '''
        obj_ptrname = it.decl.datatype[:-2]
        ow.addln("%s->%s = %sinit(%s);", self.ptrname, it.decl.name, 
          obj_ptrname, self.ptrname)
        continue

      if "def" not in it.cmds:
        print "missing default value for %s, %s" % (it, it.gtype)
        raw_input()
        continue
      if "key" not in it.cmds:
        print "missing key for %s, %s" % (it, it.gtype)
        raw_input()
        continue
      desc = it.cmds["desc"] if "desc" in it.cmds else ""
      if not it.gtype in ("static array", "dynamic array"):
        type = it.gtype if ("pointer" not in it.gtype) else "void *"
        ow.addln('XM_CFGGETC_(cfg, %s->%s, "%s", %s, %s, 1, 1, , "%s", "missing", );', 
          self.ptrname, it.decl.name, 
          it.cmds["key"], type, it.cmds["def"], 
          desc)
      elif it.gtype == "dynamic array":
        try:
          type = it.get_generic_type(offset = 1) + " *"
        except:
          print "name: %s" % it.name(); raw_input()

        ow.addln('XM_DARR_(%s->%s, %s, %s, %s, , "%s", "error", exit(1));', 
          self.ptrname, it.decl.name, 
          type, it.cmds["def"], it.cmds["cnt"], desc)
      elif it.gtype == "static array":
        type = it.get_generic_type(offset = 1)
        ow.addln('XM_SARR_(%s->%s, %s, %s, %s, , "%s");', 
          self.ptrname, it.decl.name, 
          type, it.cmds["def"], it.cmds["cnt"], desc)

    ow.addln("}")
    ow.addln("")
    return owh.gets(), ow.gets()


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
  handle("test2.c")

if __name__ == "__main__":
  main()

