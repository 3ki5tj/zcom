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
    self.addln(r'fprintf(stderr, "FILE: %%s, LINE: %%d\n", __FILE__, __LINE__);')

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
    
    self.funcname = name
    self.used_array_index = 0
  
  def simple_if(self, tag, cond, msgs):
    self.addln("%s_if (%s,", tag, cond)
    n = len(msgs)
    lead = self.sindent
    for i in range(n-1):
      self.addln(lead + '"%s"', msgs[i])
    self.addln(lead + '"%s");', msgs[n-1])

  def die_if(self, cond, *msgs):
    self.simple_if("die", cond, msgs)

  def msg_if(self, cond, *msgs):
    self.simple_if("msg", cond, msgs)
    ''' without using msg_if
    self.addln("if (%s) {", cond)
    n = len(msgs)
    for i in range(n):
      self.addln(r'fprintf(stderr, "%s");', msgs[i])
    self.addln("}")
    '''

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
      cnta = "(%s + %s)" % (m.group(1).strip(), str(int(m.group(2))+1))
    else:
      cnta = "(%s + 1)" % cnt.strip()
    
    cond = "(%s = calloc(%s, sizeof(%s))) == NULL" % (var, cnta, type)
    msg = r"no memory! var: %s, type: %s\n" % (var, type)
    self.die_if(cond, msg)
    self.init_static_array(var, type, default, cnt, desc)

  def assign(self, var, value, type):
    if type == "char *":
      self.addln('%s = ssdup("%s");', var, value)
    else:
      self.addln("%s = %s;", var, value)

  def cfgget_var(self, var, key, type, fmt, default, 
      must, test, test_first, valid, desc, addnl = 1):
    has_test  = 0 if test  in ("TRUE", "YES", "1", 1) else 1
    has_valid = 0 if valid in ("TRUE", "YES", "1", 1) else 1
    
    if not test_first:
      self.assign(var, default, type)

    if has_test:
      self.addln("if (%s) {", test)

    if test_first:
      self.assign(var, default, type)

    cond = r'0 != cfgget(cfg, &%s, "%s", "%s")' % (var, key, fmt)
    if must:
      # do not die if cfg is NULL, because this is the 'lazy' mode
      cond = 'cfg != NULL && ' + cond 
      self.die_if(cond, 
        r"missing var: %s, key: [%s], fmt: %s\n" % (var, key, fmt))
    else:
      # for optional variables, we print the message immediately
      cond = "cfg == NULL || " + cond
      self.msg_if(cond, 
        r'assuming default value\n',
        r'var: %s, key: [%s], def: %s\n' % (var, key, default) )
        
    if has_valid:
      self.die_if("!(%s)" % valid, r"failed validation: %s\n" % valid)
    if has_test:
      self.addln("}")
    if addnl:  self.addln("")

  def cfgget_flag(self, var, key, flag, type, fmt, default, 
      test, test_first, desc):
    # get flag to a temporary variable i
    if default not in ("1", "0"):
      print "flag default can only be 0 or 1, var = %s" % var
      raise Exception
    ival = "i"
    self.cfgget_var(ival, key, type, fmt, default, 
        "FALSE", # flag is usually optional
        test, test_first, "i == 0 || i == 1", desc, 0)
    self.addln("if (i) {")
    self.addln("%s |= %s;", var, flag)
    self.addln("} else {")
    self.addln("%s &= ~%s;", var, flag)
    self.addln("}")
    self.addln("")

  def finish_function(self, sret = ""):
    if self.used_array_index:
      self.vars.addln("unsigned i;")
    if self.vars.s != "":
      self.vars.addln("");

    s  = self.decl.gets()
    s += self.vars.gets()
    s += self.body.gets()
    
    tail = self.tail = CCodeWriter(self.nindents + 1)
    if len(sret): tail.addln(sret);
    tail.addln("}")
    tail.addln("")
    
    s += self.tail.gets()
    self.function = s
    self.funcname = None  # no longer inside the function
    self.used_array_index = 0
