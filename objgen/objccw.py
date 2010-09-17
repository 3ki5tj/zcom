#!/usr/bin/env python
import os, sys, re
from copy import copy

class CCodeWriter:
  '''
  output code in an indented fashion
  '''
  sindent = "  "
  def __init__(self, nindents = 0):
    self.nindents = nindents
    self.s = ""
    self.funcname = None

  def inc(self): self.nindents += 1
  def dec(self): 
    self.nindents = self.nindents-1 if self.nindents > 0 else 0

  def addraw(self, t):
    who = self.body if self.funcname else self
    if t.lstrip().startswith("}"): who.dec()
    if (who.s == "" or who.s.endswith("\n")) and not t.lstrip().startswith("#"):
      who.s += who.sindent * who.nindents
    who.s += t
    if t.rstrip().endswith("{"): who.inc()

  def add(self, t, *args):
    who = self.body if self.funcname else self
    who.addraw( t % args)

  def addln(self, t, *args):
    ''' we write to self.body, if we are in a function '''
    who = self.body if self.funcname else self
    who.addraw( (t.rstrip() % args) + '\n' )

  def remove_empty_pp(self):
    ''' remove empty preprocessor blocks 
    #if
    #else
    #endif
    '''
    lines = self.s.splitlines()
    i = 1
    while i < len(lines):
      if lines[i].startswith("#endif") and i > 0:
        lnprev = lines[i-1]
        if (lnprev.startswith("#else") or 
            lnprev.startswith("#elif")):
          # remove an empty #else/#elif branch
          lines = lines[:i-1] + lines[i:]
          continue
        elif lnprev.startswith("#if"):
          lines = lines[:i-1] + lines[i+1:]
          continue
      i += 1
    self.s = '\n'.join(lines) + '\n'

  def gets(self):
    self.remove_empty_pp()
    return self.s
 
  def print_file_line(self):
    self.addln(r'fprintf(stderr, "FILE: %%s, LINE: %%s\n", __FILE__, __LINE__);')

  def start_function(self, name, funcdef, desc):
    self.function = ""

    decl = self.decl = CCodeWriter(self.nindents) # start a declaration writer
    decl.addln("/* %s: %s */", name, desc)
    decl.addln(funcdef)
    decl.addln("{")
    
    self.vars = CCodeWriter(self.nindents + 1)

    hdr = self.hdr = CCodeWriter(0)  # prototype writer
    hdr.addln(funcdef + ";")
    self.prototype = hdr.gets()

    self.body = CCodeWriter(self.nindents + 1)  # body writer
    
    tail = self.tail = CCodeWriter(self.nindents)
    tail.addln("}")
    tail.addln("")
    
    self.funcname = name
    self.used_array_index = 0

  def init_static_array(self, var, type, default, cnt, desc):
    ''' write code for initializing a static array '''
    if self.funcname != "":
      self.used_array_index = 1
    self.addln("for (i = 0; i < %s; i++)" % cnt)
    self.addln(self.sindent + "%s[i] = %s;" % (var, default))
    self.addln("")

  def init_dynamic_array(self, var, type, default, cnt, desc):
    ''' write code for initializing a dynamic array '''
    cntptn = r"(.*)\s*\+\s*([0-9]+)$"
    m = re.match(cntptn, cnt)
    if m:
      #print "before: cnt = %s" % cnt
      cnta = "(%s + %s)" % (m.group(1).strip(), str(int(m.group(2))+1))
      #print "after:  cnt = %s" % cnt
      #raw_input()
    else:
      cnta = "(%s + 1)" % cnt.strip()
      
    self.addln("if ((%s = calloc(%s, sizeof(%s))) == NULL) {", var, cnta, type)
    self.addln(r'fprintf(stderr, "no memory! var: %s, type: %s\n");', var, type)
    self.print_file_line();
    self.addln("exit(1)")
    self.addln("}")
    self.init_static_array(var, type, default, cnt, desc)

  def cfgget_var(self, var, key, type, fmt, default, must, test, valid, desc):
    has_test = 0 if test in ("TRUE", "YES", "1", 1) else 1
    has_valid = 0 if valid in ("TRUE", "YES", "1", 1) else 1
    
    self.addln("%s = %s;", var, default)
    if has_test:
      self.addln("if (%s) {", test)

    self.addln(r'if (0 != cfgget(cfg, &%s, "%s", "%s")) {', var, key, fmt)
    if must:
      self.addln(r'fprintf(stderr, "missing var: %s, key: "%s", fmt: %%%s\n");', 
          var, key, fmt)
      self.print_file_line()
      self.addln(r'exit(1)')
    else:
      self.addln(r'fprintf(stderr, "assuming default value  %s = %s\n");', 
          var, default)
    self.addln("}")
        
    if has_valid:
      self.addln("if (!(%s)) {", valid)
      self.addln(r'fprintf(stderr, "%s failed validation:\n"', var); 
      self.addln(r'                "\t%s\n");', valid); 
      self.print_file_line();
      self.addln("exit(1);")
      self.addln("}")
    if has_test:
      self.addln("}")


  def finish_function(self):
    if self.used_array_index:
      self.vars.addln("int i;")
    if self.vars.s != "":
      self.vars.addln("");

    s  = self.decl.gets()
    s += self.vars.gets()
    s += self.body.gets()
    s += self.tail.gets()
    self.function = s
    self.funcname = None  # no longer inside the function
    self.used_array_index = 0
