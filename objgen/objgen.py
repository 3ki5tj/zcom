#!/usr/bin/env python
import os, sys, re, filecmp
from copy import copy


cmt0 = "/*"
cmt1 = "*/"

class Decl:
  '''
  C declaration, adapted from K & R
  '''
  def __init__(self, s):
    if type(s) != type(""): raise Exception
    self.raw = s.strip()
    self.nraw = len(self.raw)
    self.types = []
    self.dcl_wrapper()

  def __str__(self): 
    return ' '.join(self.types)

  def dcl_wrapper(self):
    '''
    get the data type, then call dcl
    '''
    self.pos = 0
    self.token = self.ttype = None
    self.gettoken()
    if self.ttype != "abc":
      print "complex type is yet supported, [%s]" % self.raw
      print self.ttype, self.token, self.pos
      raise Exception
    self.datatype = self.token
    self.dcl()
    self.types += [self.datatype]

  def gettoken(self):
    '''
    get a token starting from the current pos
    '''
    s = self.raw[self.pos:]
    # skip leading space
    m = re.match(r"(\s*).*", s)
    if not m: raise Exception
    self.pos += m.end(1)
    s = self.raw[self.pos:]

    if len(s) == 0: 
      #print "gettok: string exhausted"
      self.ttype = self.token = None
      return 
    m = re.match(r"([a-zA-Z_]\w*).*", s)
    if m: 
      self.ttype, self.token = "abc", s[:m.end(1)]
      self.pos += m.end(1)
      return
    m = re.match(r"(\[)(.*?)(\]).*", s)
    if m: 
      self.ttype, self.token = "[]", s[m.start(2):m.end(2)]
      self.pos += m.end(3)
      return
    if s.startswith("()"):
      self.ttype = self.token = "()"
      self.pos += 2;
      return
    self.ttype = self.token = s[0]
    #print self.pos, self.raw[self.pos]
    self.pos += 1
    #print self.pos, self.raw[self.pos]

  def dcl(self):
    '''
    declaration
    '''
    ns = 0
    while self.pos < self.nraw:
      #print "while:", self.raw, "pos:", self.pos
      self.gettoken()
      #print "while:", self.raw, "pos:", self.pos
      #print "while:  current:", self.raw[self.pos], "ttype:", self.ttype, "pos:", self.pos
      if self.ttype != "*": break
      ns += 1
    #print "dcl:    current:", self.raw[self.pos], "ttype:", self.ttype, "pos:", self.pos
    self.dirdcl()
    self.types += ["pointer"] * ns

  def dirdcl(self):
    '''
    direct declaration
    '''
    if self.ttype == '(': # ( dcl )
      #print "dirdcl: current:", self.raw[self.pos], "ttype:", self.ttype, "pos:", self.pos
      self.dcl()
      if self.ttype != ')':
        print "missing ), s: %s, pos %d" % (self.raw, self.pos)
        raise Exception
    elif self.ttype == "abc":
      self.name = self.token
    else:
      print "expected name or ( at pos: %d, s: %s" % (self.pos, self.raw)
    #print "dirdc2: current:", self.raw[self.pos], "ttype:", self.ttype, "pos:", self.pos
    while 1:
      self.gettoken()
      #print "dirtok: ttype:", self.ttype, "pos:", self.pos
      if self.ttype == "()":
        self.types += ["function"]
      elif self.ttype == "[]":
        self.types += ["array" + self.token]
      else: break


class P:
  ''' 
  simple position in the code (line, column)
  unfortunately, the most vulnerable class
  '''
  def __init__(self, row = -1, col = -1):
    self.row = row
    self.col = col
    self.check()
  def __str__(self):
    self.check()
    return "line: %d, col %d" % (self.row + 1, self.col + 1)
  def __repr__(self):
    self.check();
    return "%d, %d" % (self.row, self.col)
  def __cmp__(a, b):
    a.check()
    if type(a) != type(b):
      print "trying to compare %s with %s" % (type(a), type(b))
      raise Exception 
    b.check()
    return (a.row - b.row) * 10000 + a.col - b.col
  def check(self):
    if type(self.row) != int:
      print "row %s is not a number" % (self.row)
      raise Exception  
    if type(self.col) != int:
      print "col %s is not a number" % (self.col)
      raise Exception  
  def nextline(self):
    self.check()
    self.row += 1
    self.col = 0



class Comment:
  '''
  a C comment block 
  '''
  def __init__(self, src, pos):
    if type(pos) == int: pos = P(pos, 0)
    self.span(src, pos)     # .begin, .end
    self.extract_text(src)  # .text
    self.get_commands()     # .cmds

  def span(self, src, pos):
    ''' 
    find the beginning and end of a comment,
    starting from pos
    '''
    self.begin = copy(pos)
    self.end = copy(self.begin)
    line = src[self.begin.row]
    j = line.find(cmt0, self.begin.col)
    if j < 0:
      print "missing comment from %s" % (self.begin)
      self.checkbe()
      return -1
    else:
      self.begin.col = j
    now = P(self.begin.row, j + len(cmt0)) # create a new position

    # search the end of comment
    n = len(src)
    while now.row < n:
      j = src[now.row].find(cmt1, now.col)
      if j >= 0:
        self.end = P(now.row, j + 2)
        break
      now.nextline() # not the next line
    else:
      print "Warning: cannot find the end of comment from %s, now = %s" % (
          self.begin, now)
      raise Exception
      #return -1
    self.checkbe()
    return 0

  def extract_text(self, src):
    ''' 
    extract text within the comment
    '''
    self.checkbe()
    p = P(self.begin.row, self.begin.col + len(cmt0))
    #print "creating P at %s, self.begin: %s, self.end: %s" % ( p, self.begin, self.end)
    #raw_input()
    s = ""
    while p.row < self.end.row:
      s += src[p.row][p.col:]
      p.nextline()
    if p.col < self.end.col:
      s += src[p.row][p.col : self.end.col - len(cmt1)]
    self.text = s

  def checkbe(self):
    self.begin.check()
    self.end.check()

  def get_commands(self):
    '''
    parse text to commands, save results to .cmds
    '''
    self.cmds = {}
    s = self.text
    pos = 0
    while 1:
      # note: $$ or \$ means a literal $
      # : or = means set the current
      # := or :: means set the parser's property
      pattern = r"[^\$\\](\$)(\w+)\s*([\:\=]{1,2})\s*(.*?)(\;)"
      m = re.search(pattern, s, re.MULTILINE | re.DOTALL)
      if m == None: break
      # a command is found
      cmd = s[m.start(2) : m.end(2)]
      persist = m.end(3) - m.start(3) - 1 # the operator
      param = s[m.start(4) : m.end(4)]
      self.cmds[cmd] = [param, persist] # add to dictionary
      s = s[:m.start(1)] + s[m.end(5):]; # remove the command
      #print "a pattern is found [%s]: [%s], rest: %s" % (
      #    cmd, param, s)
      #raw_input()

    # merge the remaining string to description
    if not self.cmds.has_key("desc"):
      s = re.sub(r"[\$\\]\$", "$", s)
      sa = [a.strip() for a in s.splitlines()] # split to lines
      s = ' '.join(sa).strip().rstrip(";,.")
      self.cmds["desc"] = [s, 0] # description is not persistent
      #print "the remain string is [%s]" % self.cmds["desc"]

  def __str__(self): return self.text
  def __repr__(self): return __str__(self)


class Item:
  ''' a variable defined in a struct '''

  def __init__(self, src, pos):
    '''
    initialize the structure from a piece of code
    '''
    if type(pos) == int: pos = P(pos, 0)
    self.parse(src, pos)

  def parse(self, src, pos):
    self.begin = copy(pos)
    self.end = P()
    line = src[pos.row]
    pivot = line.find(";", pos.col) # search end of a statement
    assert pivot >= 0
    self.decl_end = P(pos.row, pivot)
    self.cmt_begin = P(pos.row, pivot+1)
    
    # everything before ';' is a declaration
    self.decl = line[pos.col:pivot].strip()
    self.find_type()
    #print "decl:", decl, "\nrest:", line[pivot+1:]
    #raw_input()

    # search the associated comment
    #print "start searching comment with from: %s cmt_begin: %s" % (
    #    src[self.cmt_begin.row][self.cmt_begin.col:], self.cmt_begin)
    cmt = Comment(src, self.cmt_begin)
    self.end = copy(cmt.end)
    self.cmds = cmt.cmds  # copy ?
    #raw_input("associated comment over: " + str(self.cmds))

    self.get_generic_type()

  def find_type(self):
    '''
    parse the declaration string, to obtain the name and type
    '''
    s = self.decl
    
    chmap = {} # use dictionary

    # check the validity of the declaration
    unwanted = "~`!@#$%^&-+={}|\\:;'\"<>,.?/"
    for c in s:
      if c in unwanted:
        print "bad declaration [%s], has character [%s]" % (s, c)
        raise Exception
      else: 
        chmap[c] = (1 if not chmap.has_key(c) else (chmap[c]+1))

    d = Decl(s)
    self.types = d.types
    self.name = d.name
    #print "type stack:", self.types

  def get_generic_type(self):
    '''
    get a generic type, and type name
    '''
    types = self.types
    nlevels = len(types)
    if nlevels == 1:
      self.typeg = ("simple", types[0])
    elif types[0].startswith("array"):
      self.typeg = ("static array", ' '.join(types[1:]))
      self.arrcnt = types[0][5:]
    elif types[0] == "pointer":
      towhat = ' '.join(types[1:])
      if nlevels == 2 and types[1] == "function":
        self.typeg = ("function", ' '.join(types[2:]))
      elif 'arr' in self.cmds:
        if 'obj' in self.cmds and nlevels == 2:
          self.typeg = ("object array", towhat)
        else:
          self.typeg = ("dynamic array", towhat)
      elif 'obj' in self.cmds:
        self.typeg = ("object", towhat) 
      else:
        self.typeg = ("pointer", towhat)


class Object:
  '''
  a struct defined as
  typedef struct {
  } abc_t
  '''
  def __init__(self, src, pos = None):
    ''' initialize myself from a bunch of src '''
    if len(src) == 0: return # return an empty object

    if type(pos) == int:  # integer input
      pos = P(pos, 0)
    self.begin = copy(pos)
    self.items = []
    self.parse(src, pos)

  def has_definition(self, line):
    ''' check if line i has a definition of Object '''
    return line.strip().startswith("typedef struct")

  def parse(self, src, pos):
    ''' parse src starting from 'pos' '''
    n = len(src)
    #print 'src:', src, '\nn:', n, '\npos:', pos
    #raw_input()

    # search the beginning {
    for i in range(pos.row, n):
      line = src[i]
      j = line.find("{") 
      if j >= 0:
        p = P(i, j+1)
        break
    else:
      print ("cannot find the structure beginner from line", 
          self.begin.row + 1)
      exit(1)
   
    #print 'starting block from', p
    #raw_input()

    # parse the content with { ... } 
    # until } is encountered
    while p.row < n:
      line_orig = src[p.row]
      line = line_orig[p.col:]

      # search the block-ending mark
      if line.strip().startswith("}"):
        pattern = r'(}\s*)(\w+)\s*;'
        m = re.search(pattern, line)
        if m == None:
          print "no structure end:", line
          exit(1)
        self.name = line[m.start(2) : m.end(2)]
        self.end = P(p.row, p.col + m.end(0))
        return
      elif line.strip() == "": # current line is exhausted or empty 
        p.nextline()
      elif line.strip().startswith(cmt0): # a stand-alone comment
        cmt = Comment(src, p)
        if len(cmt.cmds) > 1: # we will always have `desc'
          print "comment has %d commands: %s" % (len(cmt.cmds), cmt.cmds)
          self.items += [cmt]
        p = copy(cmt.end)
      elif line_orig.find(";", p.col): # search for item definition
        item = Item(src, p)
        print "new item: %s, type: %s, cmds: %s" % (
            item.name, item.typeg, item.cmds)
        #raw_input()
        self.items += [item]
        p = copy(item.end) # move to the end of item
      else:
        raw_input("suspetious line %s, continue?" % line)
        p.nextline()


def complete_lines(src):
  ''' process the input src '''
  n    = len(src)
  dest = []
  objs = []

  # first pass: get objects
  p = P(0, 0)
  obj = Object("") # get a handle 
  while p.row < n:
    if obj.has_definition(src[p.row][p.col:]):
      obj  = Object(src, p.row)
      raw_input ("found a new object '%s' with %d items, from %s to %s\n" 
          % (obj.name, len(obj.items), obj.begin, obj.end))
      objs += [obj]
      p = copy(obj.end) # go to the end of object definition
    else:
      p.nextline()

  # second pass: do stuff
  for obj in objs: pass

  return dest

def complete(file, template = ""):
  ''' generate a C file according to a template '''
  # guess the template name if it's empty
  if template == "":
    pt = os.path.splitext(file)
    template = pt[0] + '.0' + pt[1]

  # read in the template
  src = open(template, 'r').readlines()
  # write the processed file
  open(file, 'w').writelines(complete_lines(src))

def main():
  complete("test.c")
  #print Decl(sys.argv[1])
  

if __name__ == "__main__":
  main()

