#!/usr/bin/env python
import os, sys, re 
from copy    import * 
from objpos  import P 
from objdcl  import CDeclaratorList
from objcmt  import CComment
from objpre  import CPreprocessor
from objccw  import CCodeWriter
from objcmd  import Commands

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
      cmt = CComment(src, p, 2)
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
    return "pre: %5s, dl: %5s, decl: %5s, cmt: %5s, raw = [%s]" % (
        self.pre != None, self.dl != None, 
        self.decl != None, self.cmt != None, 
        self.begin.gettext(self.src, self.end).strip() )
    
  def get_generic_type(self, offset = 0):
    '''
    get a generic type, and type name
    it need commands
    `offset' is the starting offset from the bottom
    '''
    types = self.decl.types
    try:
      cmds = self.cmds if self.cmds else {}
    except AttributeError:
      print "missing commmands! %s" % self
      raw_input()

    n = len(types)
    if offset >= n: return None

    nlevels = n - offset
    if nlevels == 1:
      return types[offset]
    elif types[offset].startswith("array"):
      ''' add a command for array count '''
      arrcnt = types[offset][5:]
      if offset == 0:
        self.cmds["cnt"] = arrcnt # override cnt
        #print "static array cnt has been set to %s for %s" % (arrcnt, self)
        #raw_input()
      return "static array"
    elif types[offset] == "pointer":
      #towhat = ' '.join(types[1:])
      if nlevels == 2 and types[offset+1] == "function":
        return "function pointer"
      elif nlevels == 2 and types[offset+1] == "char":
        return "char *"
      elif 'cnt' in cmds:
        return "dynamic array"
      elif 'obj' in cmds:
        return "pointer to object"
      else:
        return "pointer"

  def fill_default(self):
    ''' set 'def' value '''
    if not self.decl or "def" in self.cmds: return
    defval = None
    gtype = self.gtype
    # for arrays, 'def' means the default value of its elements
    if len(self.decl.types) == 2 and gtype in (
        "static array", "dynamic array"):
      gtype = self.get_generic_type(offset = 1)
    
    if gtype in ("int", "unsigned", "unsigned int", 
        "long", "unsigned long"):
      defval = "0"
    elif gtype in ("float", "real"):
      defval = "0.0f"
    elif gtype in ("double"):
      defval = "0.0"
    elif gtype in ("pointer", "function pointer", "char *", 
        "pointer to object"):
      defval = "NULL"
    else: 
      defval = "(%s) 0" % self.decl.datatype
    if defval != None: 
      #print "set default %s for gtype %s, var %s" % (defval, gtype, self.decl.name)
      #raw_input()
      self.cmds["def"] = defval

  def add_item_to_list(self, dcl): 
    ''' add an item with decl being dcl '''
    it = copy(self) # make a shallow copy of myself
    it.decl = dcl
    it.dl = None   # destory the list
    self.itlist += [it]

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
    self.itlist = []
    if self.dl == None:
      self.add_item_to_list(None)
      return
    dclist = self.dl.dclist
    n = len(dclist)
    for i in range(n):
      decl = dclist[i]
      decl.types += [self.dl.datatype]
      self.add_item_to_list(decl)
    
  def isempty(self):
    return self.empty

