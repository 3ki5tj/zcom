#!/usr/bin/env python
import os, sys, re, filecmp

'''
python implementation of K & R's  dcl program
see K & R appendix A8 and $5.12
'''

class CDecl:
  def __init__(self, s):
    if type(s) != type(""): raise Exception
    self.s = s.strip()
    self.nraw = len(self.s)
    self.types = ""
    self.dcl_wrapper()

  def __str__(self):
    return self.types

  def dcl_wrapper(self):
    '''
    get the data type, then call dcl
    '''
    self.pos = 0
    self.token = self.ttype = None
    self.param_level = 0
    self.datatype = self.dclspec()
    self.dcl()
    self.types += self.datatype

  def dclspec(self, allow_empty = 0):
    # storage classifiers: auto, register, static, extern, typedef
    # type qualifier: const, volatile
    # we do not support structure yet
    alltypes = ("auto", "register", "static", "extern", "typedef",
      "const", "volatile", 
      "void", "char", "short", "int", "long", "float", "double", "signed", "unsigned")
    type = ""
    while 1:
      self.gettoken()
      if self.ttype != "word":
        if allow_empty or len(type) > 0:
          self.ungettok()
          break
        print "expected data type, %s" % (self.dbg())
        raise Exception
      if self.token not in alltypes: 
        self.ungettok()
        break
      type += self.token + " "
    return type[:-1]

  def dcl(self):
    '''
    declarator
    wrapper: counter the leading # of stars, and lead to dirdcl()
    '''
    ns = 0
    while self.pos < self.nraw:
      if self.gettoken() != "*": 
        self.ungettok()
        break
      ns += 1
    self.dirdcl()
    self.types += "pointer to " * ns

  def dirdcl(self):
    '''
    direct declarator
    '''
    self.gettoken()
    if self.ttype == '(': # ( dcl )
      self.dcl() 
      if self.gettoken() != ')':
        print "missing ) %s" % (self.dbg())
        raise Exception
    elif self.ttype == "word":
      self.name = self.token
    elif self.param_level > 0: # allow missing name
      self.ungettok()
    else:
      print "expected name or (, %s" % (self.dbg())
      raise Exception
    while 1:
      ttype = self.gettoken()
      if ttype == "(":
        self.types += "function of ("
        self.ptlist()
        if self.gettoken() != ')':
          print "missing ) after parameter type list, %s" % (self.dbg())
        self.types += ") returning "
      elif ttype == "[]":
        self.types += "array " + self.token + " of "
      else:
        self.ungettok()
        break 

  def ptlist(self):
    ''' parameter type list '''
    while 1:
      self.pdecl()
      if self.gettoken() != ',':
        self.ungettok()
        break
      self.types += ", "

  def pdecl(self):
    ''' parameter declarator '''
    self.param_level += 1
    type = self.dclspec(allow_empty = 1)
    self.dcl()
    self.types += type
    self.param_level -= 1

  def gettoken(self):
    '''
    get a token starting from the current pos
    '''
    s = self.s[self.pos:]
    m = re.match(r"(\s*).*", s)
    if not m: raise Exception
    self.pos += m.end(1)
    s = self.s[self.pos:]

    if len(s) == 0:
      self.ttype = self.token = None
      return self.ttype
    m = re.match(r"([a-zA-Z_]\w*).*", s) # a variable
    if m:
      self.ttype, self.token = "word", m.group(1)
      self.pos += m.end(1)
      return self.ttype
    m = re.match(r"(\[)(.*?)(\]).*", s)  # [...]
    if m:
      self.ttype, self.token = "[]", m.group(2)
      self.pos += m.end(3)
      return self.ttype
    self.ttype = self.token = s[0]
    self.pos += 1
    return self.ttype

  def ungettok(self):
    if self.ttype == None: return
    n = len(self.token)
    if n > self.pos:
      print "unable to unget %s" % (self.dbg())
    self.pos -= n

  def dbg(self):
    return "s: [%s], pos: %s, token: [%s], ttype: [%s]" % (
        self.s, self.pos, self.token, self.ttype)


def test(s):
  d = CDecl(s)
  print "s = %-36s desc = %s, name = %s" % (s+',', d, d.name)

def main():
  #test(sys.argv[1])
  if len(sys.argv) > 1:
    test(sys.argv[1])
    return
  test("char **argv")
  test("int (*daytab)[13]")
  test("int *daytab[13]")
  test("void *comp()")
  test("void (*comp)()")
  test("char (*(*x())[])()")
  test("char (*(*x[3])())[5]")
  test("char (*(*x(double, int *, char *s, char (*)(float)))[5])()")

if __name__ == "__main__":
  main()

