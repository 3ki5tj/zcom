#!/usr/bin/env python
import os, sys, re
from copy import copy
  
def dm(s):
  ''' %d --> %%d in a string '''
  return re.sub(r"\%", "%%", s)

class CCodeWriter:
  '''
  output code in an indented fashion
  '''
  enabletpcheck = 1
  sindent = "  "

  def __init__(self, nindents = 0):
    self.nindents = nindents
    self.s = ""
    self.funcname = None

  def inc(self): self.nindents += 1
  def dec(self): 
    self.nindents = self.nindents-1 if self.nindents > 0 else 0

  def addraw(self, t):
    ''' write a plain string `t' without format conversion '''
    # we write to self.body, if we are in a function
    who = self.body if self.funcname else self

    if t.lstrip().startswith("}"): who.dec()
    if (who.s == "" or who.s.endswith("\n")) and (
        not t.lstrip().startswith("#") and 
        not t.strip().endswith(":")):  # label
      who.s += who.sindent * who.nindents
    who.s += t
    if t.rstrip().endswith("{"): who.inc()

  def addlnraw(self, t):
    self.addraw(t + '\n')

  def add(self, t, *args):
    self.addraw( t % args)

  def addln(self, t, *args):
    self.addraw( (t % args) + '\n' )

  def remove_empty_pp(self):
    ''' remove empty preprocessor blocks 
    #if
    #else
    #endif
    '''
    endl = '\n' if self.s.endswith('\n') else ''
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
    self.s = '\n'.join(lines) + endl

  def gets(self):
    self.remove_empty_pp()
    return self.s
 
  def print_file_line(self):
    self.addln(r'fprintf(stderr, "FILE: %%s, LINE: %%d\n", __FILE__, __LINE__);')

  def errmsg(self, msg, args = ""):
    ''' print an error message with FILE/LINE 
    NOTE: args are arguments for fprintf, not formatting args for msg '''
    if len(args):
      args = ", " + args.lstrip(" ,")
    msg = r'fprintf(stderr, "error: %s\n"%s);' % (msg, args)
    self.addraw(msg + '\n')
    self.print_file_line()

  def err_ret(self, ret, msg, args = ""):
    self.errmsg(msg, args)
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
    '''
    start an if block
    note: it starts only if cond is not always true
    '''
    if self.notalways(cond):
      self.addln("if (%s) {", cond)

  def begin_else_if(self, cond):
    self.addln("} else if (%s) {", cond)

  def begin_else(self):
    self.addln("} else {")

  def end_if(self, cond):
    if self.notalways(cond):
      self.addln("}")

  def simple_if(self, tag, cond, msg, args):
    ''' 
    `args' are parameters to be passed to 
    printf-like functions: msg_if and die_if
    '''
    self.addln("%s_if (%s,", tag, cond)

    # break msg into multiple line messages
    msgs = msg.splitlines()
    n = len(msgs)
    lead = self.sindent
    for i in range(n-1):
      self.addln('%s"%s"', lead, msgs[i])
    s = '%s"%s"' % (lead, msgs[n-1])

    # build an argument list
    if args:
      s += ",\n"
      self.addraw(s)  # stop the last message line
      al = [s.strip() for s in args.split(",")]
      s = lead + ', '.join(al)
    s += ");\n"
    self.addraw(s)

  def msg_if(self, cond, msg, args = None):
    self.simple_if("msg", cond, msg, args)

  def die_if(self, cond, msg, args = None):
    self.simple_if("die", cond, msg, args)

  def insist(self, cond, msg, args = None):
    ''' die if cond is *not* true '''
    if self.notalways(cond):
      self.simple_if("die", " !(%s)" % cond, msg, args)

  def init_sarr(self, var, type, default, cnt, desc):
    ''' write code for initializing a static array '''
    self.declare_var("int i") # declare index i
    self.addln("for (i = 0; i < %s; i++)", cnt)
    self.addln(self.sindent + "%s[i] = %s;", var, default)
    self.addln("")

  def alloc_darr(self, var, type, cnt):
    self.checktp('*'+var, type)

    # modify the size to allow some margin
    cntptn = r"(.*)\s*\+\s*([0-9]+)$"
    m = re.match(cntptn, cnt)
    if m:
      cnta = "(%s + %s)" % (m.group(1).strip(), str(int(m.group(2))+1))
    else:
      cnta = "(%s + 1)" % cnt.strip()

    # allocate array
    cond = "(%s = calloc(%s, sizeof(%s))) == NULL" % (var, cnta, type)
    msg = r"no memory! var: %s, type: %s\n" % (var, type)
    self.die_if(cond, msg)

  def init_darr(self, var, type, default, cnt, desc):
    ''' write code for initializing a dynamic array '''
    if cnt.strip() != "0":
      self.alloc_darr(var, type, cnt)
      self.init_sarr(var, type, default, cnt, desc)

  def assign(self, var, value, type):
    ''' assignment (for simple types) '''
    if type == "char *" and value != "NULL":
      self.addln('%s = ssdup("%s");', var, value)
    else:
      self.addln("%s = %s;", var, value)

  def checktp(self, var, type):
    # add type checking
    if self.enabletpcheck:
      csiz = "sizeof(%s) == sizeof(%s)" % (var, type)
      msgz = r"wrong size: %s is not %s" % (var, type)
      self.insist(csiz, msgz)


  def getkeystr(self, key):
    ''' breaking key '''
    keyslen = 128
    # support dynamic key
    pivot = key.find("#")
    if pivot >= 0: # special key
      keystr = "keystr"
      fmt  = key[:pivot]
      args = key[pivot+1:]
      # attentatively check the length of the resulting keystr
      if len(fmt) + len(args.split(",")) * 20 > keyslen:
        print "key [%s] seems too long" % key
        raise Exception
      self.declare_var("char keystr[%d];" % keyslen, "keystr")
      self.addraw('sprintf(keystr, "%s", %s);\n' % (fmt, args));
    else:
      keystr = '"%s"' % key
      fmt = key
      args = ""
    return keystr, fmt, args


  def cfgget_var_low(self, var, key, type, fmt, default, 
      must, valid, desc):
    ''' get a variable and make a validation '''

    self.checktp(var, type)

    keystr, keyfmt, keyargs = self.getkeystr(key)

    cond = r'0 != cfgget(cfg, &%s, %s, "%s")' % (var, keystr, fmt)
    if must:
      # do not die if cfg is NULL, because this is the 'lazy' mode
      cond = 'cfg != NULL && ' + cond 
      self.die_if(cond, 
        r"missing var: %s, key: %s, fmt: %s\n" % (var, dm(key), fmt))
    else:
      # for optional variables, we print the message immediately
      cond = "cfg == NULL || " + cond
      self.msg_if(cond, 
        (r'assuming default value\n' + '\n' +
         r'var: %s, key: %s, def: %s\n') % (var, dm(key), default) )
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

  def str_to_arr(self, sbuf, var, type, fmt, cnt, 
      valid, cmpl, desc):
    ''' parse a string to array item '''
    ps = "ps"

    sbufgood = "%s != NULL" % sbuf
    self.begin_if(sbufgood)
    # declare local variables
    self.addln("char *%s = %s;", ps, sbuf)
    self.addln("int nps = 0;\n") 
    # self.declare_var("char *%s" % ps, ps)
    # self.declare_var("int nps")
    
    self.addln("for (i = 0; i < %s; i++) {" % cnt)
    s = r'1 != sscanf(ps, "%s%%n", %s+i, &nps)' % (fmt, var)
    if cmpl:  # allow incomplete
      self.die_if(s, "incomplete array %s at i = %%s" % var, "i")
    else:
      self.begin_if(s)
      self.addln("break;")
      self.end_if(s)
    if valid:
      self.insist(valid, 
        r"%s: validation %s failed at i = %%s\n" % (var, valid),
        "i")
    self.addln("ps += nps;")
    self.addln("}")
    self.end_if(sbufgood)

  def cfgget_sarr(self, var, key, type, fmt, default, cnt, 
      must, valid, cmpl, desc):
    # make initial assignment first 
    self.init_sarr(var, type, default, cnt, desc)

    # read a string from configuration file
    sbuf = "sbuf" # declare a string buffer
    self.declare_var("char *%s" % sbuf, sbuf)

    self.addln("%s = NULL;", sbuf)
    self.cfgget_var_low(sbuf, key, type, "%s", default,
        must, 1, desc)

    # parse the string to get individual items
    self.str_to_arr(sbuf, var, type, fmt, cnt, valid, cmpl, desc)
    self.addln("")

  def cfgget_flag(self, var, key, flag, type, fmt, default, 
      prereq, desc):
    ''' 
    assign a flag to var if prereq is true
    '''
    # get flag to a temporary variable i
    if (default not in ("1", "0") or 
        type not in ("unsigned", "unsigned int", "unsigned long")):
      print "flag %s (key %s) with bad default %s or type %s" % (var, dm(key), default, type)
      raise Exception
    
    self.begin_if(prereq)
    ivar = "i"
    self.assign(ivar, default, type)
    self.cfgget_var_low(ivar, key, "int", fmt, default, 
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
    ''' 
    start a function
    construct four code-writers:  decl, vars, hdr, body 
    set `funcname' to the function name, indicating
      we are in a function
    start an empty variable declaration list, `vls'
    '''
    if not name:
      print "need a function name"
      raise Exception

    # declare the function
    decl = self.decl = CCodeWriter(self.nindents) # start a declaration writer
    decl.add_comment("%s: %s" % (name, desc))
    decl.addln(fdecl)
    decl.addln("{")
    
    # create a variable list `vls'
    self.vars = CCodeWriter(self.nindents + 1)
    self.vls = {}

    # prototype
    hdr = self.hdr = CCodeWriter(0)  # prototype writer
    if len(macro): hdr.addln(macro)  # extra macro
    hdr.addln(fdecl + ";")
    self.prototype = hdr.gets()

    # body writer
    self.body = CCodeWriter(self.nindents + 1)  # body writer
    
    self.funcname = name
    self.function = ""  # the content of the function
 
  def end_function(self, sret = ""):
    # write variable list
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
    self.vls = None  # kill variable list
    self.decl = self.vars = self.body = self.hdr = self.tail = None

  def declare_var(self, dcl, var = None):
    ''' 
    declare a variable within a function
    `var': variable name
    `dcl': declarator
    e.g. in 
      double *arr;
    `var' is arr; `double *arr' is the declarator
    '''
    if self.funcname == "":
      print "need to be in a function to declare var. %s of type %s" % (var, tp)
      raise Exception

    dcl = dcl.rstrip(";")
    if var == None: # make a cheap guess of variable
      arr = dcl.split()
      var = arr[len(arr) - 1].strip("*()")
      if var == "" or var[0].isdigit() or not re.match("\w+$", var):
        print "cannot determine the var [%s]" % var
        raise Exception

    if var in self.vls:
      newdcl = self.vls[var]
      if dcl != newdcl:  # reject if the new type is different from the existing one
        print "var. %s is already declared as %s, diff. from %s" % (var, dcl, newdcl)
        raise Exception
    else:
      self.vls[var] = dcl
      self.vars.addln(dcl + ";")
  
  def wb_var(self, var, type):
    if type == "int":
      s = "BIO_WI(%s);" % var
    elif type == "double":
      s = "BIO_WD(%s);" % var
    else:
      print "don't know how to write type %s for %s" % (type, var)
      raise Exception;
    self.addraw(s + "\n")

  def wb_arr1d(self, arr, cnt, type):
    if type == "int":
      s = "BIO_WIARR(%s, %s);" % (arr, cnt)
    elif type == "double":
      s = "BIO_WDARR(%s, %s);" % (arr, cnt)
    else:
      print "don't know how to write type %s for %s" % (type, arr)
      raise Exception
    self.addraw(s + "\n")

  def wb_checkbytes(self):
    ''' write some checking bytes '''
    self.declare_var("int size");
    self.addln("size = sizeof(int);")
    self.wb_var("size", "int")
    self.addln("size = sizeof(int);")
    self.wb_var("size", "int")

  def wb_arr2d(self, arr, dim, tp, trim):
    self.declare_var("int j")
    self.declare_var("int i")
    self.declare_var("int size")
    
    pb = "pb%s" % tp[:1]
    self.declare_var("%s *%s" % (tp, pb))
    if trim:
      self.declare_var("int imax")
      self.declare_var("int imin")

    self.addln("for (j = 0; j < %s; j++) {" % dim[0])
    self.addln("%s = %s + j * %s;" % (pb, arr, dim[1]))
    self.addln("for (imin = 0, i = %s-1; i >= 0 && %s[i] > 0.0; i--) ;",
      dim[1], pb)
    self.addln("imax = i + 1;")
    self.addln("for (i = 0; i < imax && %s[i] > 0.0; i++) ;",
      pb)
    self.addln("imin = i;")
    self.addln("if ((size = imax - imin) <= 0) continue;")
    self.wb_var("j",    "int")  # current row index
    self.wb_var("imin", "int")  # lowest
    self.wb_var("size", "int")
    self.wb_arr1d("%s+imin" % pb, "size", tp)
    self.addln("}")

  def wb_arr(self, arr, dim, tp, trim):
    ndim = len(dim)
    if ndim == 2:
      self.wb_arr2d(arr, dim, tp, trim)
    elif ndim == 1:
      self.wb_arr1d(arr, dim[0], tp)
    else: raise Exception

