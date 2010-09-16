#!/usr/bin/env python
import os, sys, re 
from copy import copy

'''
C declarator list
'''

class CDeclarator:
  '''
  a single C declarator, adapted from K & R
  '''
  def __init__(self, src, p):
    self.cnt = 0
    self.parse(src, p)

  def parse(self, src, p):
    s = p.skipspace(src)
    self.begin = copy(p)

    self.param_level = 0 # not inside a function parameter list
    self.types = []
    
    # find the declarator, e.g. *a, (*a)(), ...
    if 0 != self.dcl(src, p): 
      return -1
    self.cnt = 1
    
    self.end = copy(p)
    self.raw = self.begin.gettext(src, self.end)
    #print "complete a decl {%s} from %s to %s\n{%s}" % (self, self.begin, self.end, self.raw)
    return 0

  def dclspec(self, src, p, allow_empty = 0):
    # storage classifiers: auto, register, static, extern, typedef
    # type qualifier: const, volatile
    p0 = copy(p)
    alltypes = ("auto", "register", "static", "extern",
          "const", "volatile", 
          "void", "char", "short", "int", "long", 
          "float", "double", "signed", "unsigned", 
          "FILE", "size_t", "ssize_t", "fpos_t",
          "div_t", "ldiv_t", "clock_t", "time_t", "tm",
          "va_list", "jmp_buf", "__jmp_buf_tag",
          "lconv", "fenv_t", "fexcept_t"
          "wchar_t", "wint_t", "wctrans_t", "wctype_t",
          "errno_t", "sig_atomic_t", 
          "int8_t", "uint8_t", "int16_t", "uint16_t",
          "int32_t", "uint32_t", "int64_t", "uint64_t",
          "bool", "_Bool", "complex",
          "real")
    # MPI function types are yet added
    alltypes += ("MPI_Comm", 
          "MPI_Datatype", "MPI_Group", "MPI_Win",
          "MPI_File", "MPI_Op", "MPI_Topo_type",
          "MPI_Errhandler", "MPI_Request", "MPI_Info",
          "MPI_Aint", "MPI_Fint", "MPI_Offset", 
          "MPI_Status")
    type = ""
    while 1:
      token, ttype = p.gettoken(src)
      if ttype != "word":
        if allow_empty or len(type) > 0:
          p.ungettok(src)
          break
        print "expected data type from %s, type=%s, %s" % (p0, type, self.dbg(src, p))
        raise Exception
        return None
      # guess if the token is a type specifier
      if (token not in alltypes and
          token[:2]  not in ("t_", "T_") and
          token[:3]  not in ("_t_", "_T_") and
          token[:4]  not in ("__t_", "__T_") and
          token[-2:] not in ("_t", "_T") and 
          token[-3:] not in ("_t_", "_T_") and 
          token[-4:] not in ("_t__", "_T__") ):
        p.ungettok(src)
        break
      type += token + " "
      #print "datatype is [%s], pos %s" % (self.type, p)
    return type[:-1]

  def dcl(self, src, p):
    '''
    declarator
    wrapper: counter the leading # of stars, and lead to dirdcl()
    '''
    ns = 0
    while 1:
      token, type = p.gettoken(src)
      if token != "*":
        p.ungettok(src)
        break
      ns += 1
    if self.dirdcl(src, p) != 0: return -1
    if self.param_level == 0: self.types += ["pointer"] * ns
    return 0

  def dirdcl(self, src, p):
    '''
    direct declarator
    '''
    token, ttype = p.gettoken(src)
    if ttype == '(': # ( dcl )
      self.dcl(src, p) 
      if p.gettoken(src) != (')', ')'):
        print "missing ) %s" % (self.dbg(src, p))
        raise Exception
    elif ttype == "word":
      self.name = token
    elif self.param_level > 0:
      p.ungettok(src)
    else:
      print "expected name or (, %s" % (self.dbg(src, p))
      raise Exception
    while 1:
      token, ttype = p.gettoken(src)
      if ttype == "(":
        self.ptlist(src, p)
        token, type = p.gettoken(src)
        if token != ')':
          print "missing ) after parameter list, token: [%s], %s" % (token, self.dbg(src, p))
          raise Exception
        if self.param_level == 0: self.types += ["function"]
      elif ttype == "[": 
        # array dimension can only span a one line
        pattern = r"(.*)\]"
        s = p.getline(src)
        m = re.match(pattern, s)
        if m:
          token = m.group(1)
          if self.param_level == 0: self.types += ["array" + token]
          p.col += m.end(0)
        else:
          print "missing ], %s" % (self.dbg(src, p))
          raise Exception
      else:
        p.ungettok(src)
        break
    return 0

  def ptlist(self, src, p):
    ''' parameter type list '''
    while 1:
      self.pdecl(src, p)
      token, ttype = p.gettoken(src)
      if token != ',':
        p.ungettok(src)
        break

  def pdecl(self, src, p):
    ''' parameter declarator '''
    self.param_level += 1
    self.dclspec(src, p, allow_empty = 1)
    self.dcl(src, p)
    self.param_level -= 1

  def __str__(self):
    return ' '.join(self.types) + ' ' + self.name
  def __repr__(self):
    return self.types
  def dbg(self, src, p):
    ''' return the debug message '''
    return "line: [%s], pos: %s" % (p.getline(src).rstrip(), p)
  def isempty(self): return 0 if self.cnt else 1

class CDeclaratorList(CDeclarator): 
  '''
  C declarations of multiple variables
  we inherit `dclspec()' and `datatype' from CDeclarator
  '''
  def __init__(self, src, p):
    self.dclist = []
    if 0 != self.parse(src, p): self.cnt = 0

  def parse(self, src, p):
    p.skipspace(src)
    self.begin = copy(p)

    # pre-filter out the case where the first token is not a word
    token, type = p.peektok(src)
    if type != "word": return -1

    # find the type, e.g. int, void, ...
    self.param_level = 0 # dclspec needs it
    self.datatype = self.dclspec(src, p)
    if self.datatype == None: return -1

    #print "start to find variables from %s" % (p)
    while 1: # repeatedly get variables
      # we find declarators (variables) through CDeclarator's
      # can s
      d = CDeclarator(src, p)
      if d.isempty(): return -1
      d.datatype = self.datatype # add datatype for individual CDeclarator
      self.dclist += [d]
      
      # handle the last semicolon
      token, ttype = p.gettoken(src)
      if token == ";": 
        self.end = copy(p)
        break
      if token != ",":
        print "missing ; or , token [%s] %s" % (token, self.dbg(src, p))
        raise Exception
        return -1
    
    self.raw = self.begin.gettext(src, self.end)

  def isempty(self): return 0 if len(self.dclist) > 0 else 1

