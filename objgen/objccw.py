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

  def err_msg(self, msg, *args):
    msg = r'fprintf(stderr, "error: %s\n");' % (msg % args)
    self.addln(msg)
    self.print_file_line()

  def err_ret(self, msg, ret):
    self.err_msg(msg)
    self.addln("return %s;", ret)

  def add_comment(self, text):
    ''' format text and put them in comment '''
    if not text or len(text) == 0: return
    arr = text.splitlines()
    n = len(arr)
    for i in range(n):
      s = "/* " if i == 0 else " * "
      s += arr[i].strip()
      if i == n - 1: s += " */"
      self.addln(s)

  def notalways(self, cond):
    ''' test if a condition is missing or always true '''
    return (cond and cond not in ("TRUE", "1", 1))

  def begin_if(self, cond):
    if self.notalways(cond):
      self.addln("if (%s) {" % cond)

  def begin_else(self):
    self.addln("} else {")

  def end_if(self, cond):
    if self.notalways(cond):
      self.addln("}")

  def simple_if(self, tag, cond, msgs):
    self.addln("%s_if (%s,", tag, cond)
    n = len(msgs)
    lead = self.sindent
    for i in range(n-1):
      self.addln(lead + '"%s"', msgs[i])
    self.addln(lead + '"%s");', msgs[n-1])

  def msg_if(self, cond, *msgs):
    self.simple_if("msg", cond, msgs)

  def die_if(self, cond, *msgs):
    self.simple_if("die", cond, msgs)

  def insist(self, cond, *msgs):
    ''' die if cond is *not* true '''
    if self.notalways(cond):
      self.simple_if("die", " !(%s)" % cond, msgs)

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

  def cfgget_var_low(self, var, key, type, fmt, default, 
      must, valid, desc):
    ''' get a variable and make a validation '''
    cond = r'0 != cfgget(cfg, &%s, "%s", "%s")' % (var, key, fmt)
    if must:
      # do not die if cfg is NULL, because this is the 'lazy' mode
      cond = 'cfg != NULL && ' + cond 
      self.die_if(cond, 
        r"missing var: %s, key: %s, fmt: %s\n" % (var, key, fmt))
    else:
      # for optional variables, we print the message immediately
      cond = "cfg == NULL || " + cond
      self.msg_if(cond, 
        r'assuming default value\n',
        r'var: %s, key: %s, def: %s\n' % (var, key, default) )
    self.insist(valid, r"failed validation: %s\n" % valid)

  def cfgget_var(self, var, key, type, fmt, default,
      must, prereq, tfirst, valid, desc):
    # if prereq is missing, it defaults to 1 
    if not tfirst: self.assign(var, default, type)
    self.begin_if(prereq)    
    if tfirst: self.assign(var, default, type)
    
    self.cfgget_var_low(var, key, type, fmt, default,
        must, valid, desc)
    
    self.end_if(prereq)
    self.addln("")

  def cfgget_flag(self, var, key, flag, type, fmt, default, 
      prereq, desc):
    ''' 
    assign a flag to var if prereq is true
    '''
    # get flag to a temporary variable i
    if (default not in ("1", "0") or 
        type not in ("unsigned", "unsigned int", "unsigned long")):
      print "flag %s (key %s) with bad default %s or type %s" % (var, key, default, type)
      raise Exception
    
    self.begin_if(prereq)
    ivar = "i"
    self.assign(ivar, default, type)
    self.cfgget_var_low(ivar, key, type, fmt, default, 
        "FALSE", # flag is usually optional
        "%s == 0 || %s == 1" % (ivar, ivar), desc)
    self.begin_if(ivar)
    self.addln("%s |= %s;", var, flag)
    self.begin_else()
    self.addln("%s &= ~%s;", var, flag)
    self.end_if(ivar)
    self.end_if(prereq)
    self.addln("")

  def begin_function(self, name, fdecl, desc, macro = ""):
    ''' start four code-writers:  decl, vars, hdr, body '''
    self.function = ""

    decl = self.decl = CCodeWriter(self.nindents) # start a declaration writer
    decl.add_comment("%s: %s" % (name, desc))
    decl.addln(fdecl)
    decl.addln("{")
    
    self.vars = CCodeWriter(self.nindents + 1)

    hdr = self.hdr = CCodeWriter(0)  # prototype writer
    if len(macro): hdr.addln(macro)
    hdr.addln(fdecl + ";")
    self.prototype = hdr.gets()

    self.body = CCodeWriter(self.nindents + 1)  # body writer
    
    self.funcname = name
    self.used_array_index = 0
 
  def end_function(self, sret = ""):
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
    self.decl = self.vars = self.body = self.hdr = self.tail = None
