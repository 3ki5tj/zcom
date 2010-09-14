#!/usr/bin/env python
import os, sys, re, filecmp
from copy    import copy
from objpos  import P 
from objdcl  import CDeclaratorList
from objcmt  import CComment
from objpre  import CPreprocessor
from objccw  import CCodeWriter
from objcmd  import Commands
'''
intend to be a code generator
TODO:
  * multiple-line support, e.g., typedef\nstruct
  * preprocessor
  * skip struct in comment or string
'''



class Item:
  '''
  item:
    preprocessor
    declaration
    declaration comment
    comment
  '''

  def __init__(self, src, p):
    '''
    initialize the structure from a piece of code
    '''
    self.empty = 1
    self.begin = copy(p)
    if src == None: return
    self.src = src
    self.parse(src, p)

  def parse(self, src, p):
    '''
    parse an item aggressively
    '''
    s = p.skipspace(src)
    if s[0] == "}": return -1
 
    # try to see it's a preprocessor
    pre = CPreprocessor(src, p)
    if not pre.isempty():
      #print "pos %s, preprocessor %s" % (p, pre.raw)
      #raw_input()
      self.empty = 0
      self.pre = pre
      self.dl = self.cmt = None
    else:
      self.pre = None

      # try to get a declaration
      self.decl = None
      dl = CDeclaratorList(src, p)
      self.dl = None if dl.isempty() else dl
      
      # try to get a comment
      cmt = CComment(src, p)
      self.cmt = cmt if not cmt.isempty() else None
    
    if self.pre or self.cmt or self.dl:
      #print "successfully got one item, at %s" % p
      self.end = copy(p)
      self.empty = 0
    else:
      print "bad line: at %s, s = [%s]" % (p, s.strip())
      raise Exception
   
    self.expand_multiple_declarators()

  def __str__(self):
    return "pre: %5s, dl: %5s, cmt: %5s, raw = [%s]" % (
        self.pre != None, self.dl != None, self.cmt != None, 
        self.begin.gettext(self.src, self.end).strip() )
    
  def get_generic_type(self):
    '''
    get a generic type, and type name
    also need commands
    '''
    types = self.decl.types
    cmds = self.cmds if self.cmds else []
    nlevels = len(types)
    if nlevels == 1:
      return types[0]
    elif types[0].startswith("array"):
      self.arrcnt = types[0][5:]
      return "static array"
    elif types[0] == "pointer":
      #towhat = ' '.join(types[1:])
      if nlevels == 2 and types[1] == "function":
        return "function pointer"
      elif nlevels == 2 and types[1] == "char":
        return "string"
      elif 'cnt' in cmds:
        return "dynamic array"
      elif 'obj' in cmds:
        return "object"
      else:
        return "pointer"

  def expand_multiple_declarators(self):
    '''
    create a list of items, each with only one declarator
    since an item may contain a list of declarators (variable names), e.g.
      int a, b, c;
    it is inconvenient to us
    here we create a bunch of items, each of which is devoted for a single
    declarator, such as a, b, and c, the list is saved to item_list,
    Object should use items in the list instead of this item itself
    '''
    if self.dl == None:
      self.decl = None
      self.itlist = [self]
      return
    dclist = self.dl.dclist
    n = len(dclist)
    if n <= 0:
      self.decl = None
      self.itlist = [self]
      return

    self.itlist = []
    for i in range(n):
      decl = dclist[i]
      decl.types += [self.dl.datatype]
      it = copy(self)  # make a shallow copy of myself
      it.decl = decl
      #it.gtype = it.get_generic_type()
      self.itlist += [it]
    
  def isempty(self):
    return self.empty

class Object:
  '''
  a struct defined as
  typedef struct {
    ...
  } abc_t;
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
      return

    while 1: # parse the content within { ... },
      # aggressively search for an item, and update p
      item = Item(src, p)
      if item.isempty(): break

      #print "line %3d has %d item(s): %s" % (p.row + 1, len(item.itlist), item)
      self.items += item.itlist
   
    # parse }
    self.find_ending(src, p, 0)

    self.merge_comments() # merging multiple-line C++ comments
    self.get_cmds_gtype() # should be done after merging comments
    self.determine_prefix()  #

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
        self.name = line[m.start(2) : m.end(2)]
        self.end = copy(p)
        self.empty = 0
        #print "object-ending is found in %s, %s" % (p, src[p.row])
        return 1
      if aggr:
        p.nextline()
      else:
        break
    print "no object-ending is found, p at", p
    return 0  # not found

  def isempty(self): return self.empty

  def determine_prefix(self):
    '''
    determine the nickname, variable name for a pointer to the object
    and the prefix, the name attached to functions of the object
    '''
    name = self.name
    if name.endswith("_t"): name = name[:-2]
    if name.endswith("data") and len(name) > 4: # abcdata --> abc
      name = name[:-4]
    name = name.lower()
    self.prefix   = name + "_"
    self.nickname = name

  def merge_comments(self):
    ''' merge multiple comments '''
    i = 1
    items = self.items
    while i < len(items):
      it = items[i]
      itp = items[i-1]
      if (it.cmt and not it.decl 
          and i > 0 and itp.cmt  # stand-alone comment allowing another
          and it.cmt.begin.col == itp.cmt.begin.col): # starting at the same column
        #print "merging commands from item %d and item %d" % (i-1, i)
        #print "'%s'\n+\n'%s'" % (itp.cmt.raw, it.cmt.raw)
        #raw_input()
        itp.cmt.raw += it.cmt.raw
        items = items[:i] + items[i+1: ] # create a new list
        continue
      i += 1
    self.items = items
    #print "%d %d" % (len(items), len(self.items))
    #raw_input()

  def get_cmds_gtype(self):
    ''' get cmds and determine generic type ''' 
    for it in self.items:
      if it.cmt:
        it.cmds = Commands(it.cmt.raw)
      else:
        it.cmds = None
    
    # gtype uses information from commands
    for it in self.items:
      if it.decl:
        it.gtype = it.get_generic_type()
      
  def gen_code(self):
    funclist = []
    funclist += [self.gen_func_close()]

    # assemble everything to header and source
    header = self.gen_decl()
    source = ""
    for f in funclist:
      header += f[0]
      source += "\n" + f[1] 
    header = header[:-1] # strip the final '\n'
    self.header = header
    self.source = source

  def gen_decl(self):
    '''
    return a string for formatted declaration
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
          and len(cmds["desc"][0]) > 0):
        scmt = "/* " + cmds["desc"][0] + " */"
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
    macroname = "%sclose" % self.prefix
    funcname = macroname + "_"
    owh.addln("#define %s(%s) { %s(%s); free(%s); %s = NULL; }",
      macroname, self.nickname, funcname, self.nickname, self.nickname, self.nickname)
    funcdesc = "/* %s: close a pointer to %s */" % (funcname, self.name)
    #owh.addln(funcdesc)
    ow.addln(funcdesc)
    fdecl = "void %s(%s *%s)" % (funcname, self.name, self.nickname)
    owh.addln(fdecl + ";") # prototype
    ow.addln(fdecl)
    ow.addln("{")

    # compute the longest variable
    calc_len = lambda it: len(it.decl.name) if (it.decl and 
        it.gtype in ("string", "dynamic array") ) else 0
    maxwid = max(calc_len(it) for it in self.items)
    for it in self.items:
      if it.pre:
        ow.addln("#" + it.pre.raw)
        continue
      if not it.decl: continue

      funcfree = None
      if it.gtype == "string": 
        funcfree = "ssdelete"
      elif it.gtype == "dynamic array":
        funcfree = "free"
      if not funcfree: continue
      
      ow.addln("if (%s->%-*s != NULL) %s(%s->%s);",
               self.nickname, maxwid, it.decl.name,
               funcfree, self.nickname, it.decl.name)

    ow.addln("memset(%s, 0, sizeof(*%s));", self.nickname, self.nickname)
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
      if obj.isempty(): break

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

