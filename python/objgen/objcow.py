#!/usr/bin/env python
import os, sys, re
from copy   import copy
from objext import *

def dm(s):
  ''' %d --> %%d in a string '''
  return re.sub(r"\%", "%%", s)

class CCodeWriter:
  '''
  output code in an indented fashion
  '''
  enabletpcheck = 0  # disable it by default
  sindent = "  "

  def __init__(self, nindents = 0):
    self.nindents = nindents
    self.s = ""
    self.funcname = None

  def inc(self):
    self.nindents += 1
  def dec(self):
    self.nindents = self.nindents-1 if self.nindents > 0 else 0

  def addraw(self, t):
    r'''
    write a plain string `t' without format conversion
    indents are added if the current text ends with `\n'
    no newline is appended after `t'
    '''
    # we write to self.body, if we are in a function
    who = self.body if self.funcname else self

    if t.lstrip().startswith("}"): who.dec()
    if (who.s == "" or who.s.endswith("\n")) and (
        t.strip() != "" and   # no need to indent for a blank line
        not t.lstrip().startswith("#") and   # preprocessor
        not t.strip().endswith(":")):        # label
      who.s += who.sindent * who.nindents
      while t.startswith("\\t"):
        t = t[2:]
        who.s += who.sindent
      while t.startswith("\t"):
        t = t[1:]
        who.s += who.sindent
    who.s += t
    if t.rstrip().endswith("{"): who.inc()

  def addlnraw(self, t):
    if not '\n' in t:
      self.addraw(t + '\n') # single line
    else:
      for line in t.splitlines():
        self.addraw(line + '\n')

  def addln(self, t = None, *args):
    r'''
    add text with an appending '\n'
    now support multiple line
    '''
    if not t: return self.addraw('\n')
    self.addlnraw(t % args)

  def gets(self):
    return trimcode(self.s)

  def dump_block(self, block, tab = 4, offset = 2):
    '''
    dump a block of lines that are broken to fields
    '''
    if len(block) == 0: return
    nfields = len(block[0])
    wid = [8] * nfields
    #print "block of %s\n%s" % (len(block), block)
    #raw_input()

    # align and format all items in the bufferred block
    # and append the result to string s
    for j in range(nfields):
      w = max(len(item[j]) for item in block) + 1 # +1 for the following blank
      wid[j] = ((w + offset + tab - 1) // tab)*tab - offset

    for item in block:
      s = ""
      for j in range(nfields-1):
        s += "%-*s" % (wid[j], item[j])
      s += item[nfields - 1]
      self.addln(s.rstrip())

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

  def begin_if(self, cond):
    '''
    start an if block
    note: it starts only if cond is not always true
    '''
    if notalways(cond):
      self.addln("if (%s) {", cond)

  def begin_else_if(self, cond):
    self.addln("} else if (%s) {", cond)

  def begin_else(self):
    self.addln("} else {")

  def end_if(self, cond):
    if notalways(cond):
      self.addln("}")

  def print_file_line(self):
    self.addln(r'fprintf(stderr, "FILE: %%s, LINE: %%d\n", __FILE__, __LINE__);')

  def errmsg(self, msg, args = "", addfl = 0, pfx = None):
    '''
    print an error message with FILE/LINE
    NOTE: args are arguments for fprintf, not formatting args for msg
    arguments are separated by `,'
    '''
    if msg in (None, ""): return
    lead = self.sindent * 2

    # break msg into multiple line messages
    msgs = msg.splitlines()
    n = len(msgs)
    sarr = [""] * n
    sarr[0] = r'fprintf(stderr, "%s%s\n"' % (
        (pfx + ": " if pfx else ""), msgs[0])
    for i in range(1, n): # additional lines
      sarr[i] = r'%s"%s\n"' % (lead, msgs[i])

    # build an argument list
    if args:
      al = [a.strip() for a in args.split(",")] if args else []
      sarg = ', '.join(al)
      if len(sarg) + len(sarr[n-1]) >= 78: # last line is too long
        sarr[n-1] += ','
        # start a new line for arguments
        sarr += [lead + sarg + ');']
      else:
        sarr[n-1] += ', ' +sarg + ');'
    else:
      sarr[n-1] += ');'

    # print each line
    for l in sarr: self.addraw(l + '\n')

    if addfl: self.print_file_line()

  def foo_if(self, cond, msg, args, addfl, onerr = None, pfx = None):
    '''
    `args' are parameters to be passed to
    printf-like functions: msg_if and die_if
    '''
    self.begin_if(cond)
    self.errmsg(msg, args, addfl, pfx)
    if onerr: self.addlnraw(onerr)
    self.end_if(cond)

  def msg_if(self, cond, msg, args = None, onerr = None):
    self.foo_if(cond, msg, args, addfl = 0, onerr = onerr)

  def die_if(self, cond, msg, args = None, onerr = "exit(1);"):
    self.foo_if(cond, msg, args, addfl = 1, onerr = onerr)

  def insist(self, cond, msg, args = None, onerr = "exit(1);"):
    ''' die if cond is *not* true '''
    if notalways(cond):
      self.foo_if(" !(%s) " % cond, msg, args, addfl = 1, onerr = onerr)

  def init_sarr(self, var, default, cnt, pp = None):
    ''' write code for initializing a static array '''
    if default == None:
      print "Warning: don't know how to initialize array %s" % var
      return
    self.declare_var("int i", pp = pp) # declare index i
    self.addln("for (i = 0; i < %s; i++)", cnt)
    self.addln(self.sindent + "%s[i] = %s;", var, default)

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
    msg = r"no memory! var: %s, type: %s" % (var, type)
    self.die_if(cond, msg)

  def init_darr(self, var, type, default, cnt, valid, desc, pp = None, onerr = "goto ERR;"):
    ''' write code for initializing a dynamic array '''
    # since we allocate array space by default (as long as $cnt exists
    # assign cnt to 0 is the trick to avoid dynamic allocation
    if cnt.strip() != "0":
      self.alloc_darr(var, type, cnt)
      self.init_sarr(var, default, cnt, pp)
      self.validate(valid, var, onerr = onerr)
    else:
      self.addln(var+" = NULL;")

  def assign(self, var, val, type, strspec = 0):
    ''' assignment (for simple types) '''
    callcmd = "$>"
    if val.startswith(callcmd):
      self.addln(val[len(callcmd):].lstrip() + ";")
    else:
      if strspec and type == "char *" and re.match(r'".*"$', val):
        val = "ssdup(%s)" % val # apply ssdup for a bare string
      self.addln("%s = %s;", var, val)

  def checktp(self, var, type):
    # add type checking
    if self.enabletpcheck:
      sz1 = "sizeof(%s)" % var
      sz2 = "sizeof(%s)" % type
      msgz = r"wrong size: %s=%%d, %s=%%d" % (sz1, sz2)
      self.insist("%s == %s" % (sz1, sz2), msgz,
          "(int) %s, (int) %s" % (sz1, sz2))

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

  def validate(self, valid, var = None,
      fmt = "", args = None, onerr = "exit(1);"):
    if notalways(valid):
      msg = (var + ": ") if var else ""
      msg += "failed validation: %s" % escape(valid)
      if fmt: msg += ", " + fmt
      self.insist(valid, msg, args, onerr)

  def cfgget_var_low(self, var, key, type, fmt, default,
      must, valid, desc, var_ = None, onerr = "exit(1);"):
    '''
    get a variable and make a validation
    var_ is the actual variable, used for printing
    '''

    if not var_: var_ = var;
    self.checktp(var, type)
    keystr, keyfmt, keyargs = self.getkeystr(key)
    cond = '0 != cfg_get(cfg, &%s, %s, "%s")' % (var, keystr, fmt)
    if must in (1, "1"):
      # do not die if cfg is NULL, because this is the 'lazy' mode
      cond = 'cfg != NULL && ' + cond
      self.die_if(cond,
        "missing var: %s, key: %s, fmt: %%%s" % (var_, dm(key), fmt),
        onerr = onerr)
    else:
      # for optional variables, we print the message immediately
      cond = "cfg == NULL || " + cond
      self.msg_if(cond,
        'assuming default: %s = %s, key: %s'
        % (var_, escape(default), dm(key)))

    self.validate(valid, var_, onerr = onerr)

  def cfgget_var(self, var, key, type, fmt, default,
      must, prereq, tfirst, valid, desc, onerr = "goto ERR;"):

    # char * is to be allocated dynamically, so be careful
    if (not tfirst and notalways(prereq)
        and type == "char *" and default != "NULL"):
      self.assign(var, "NULL", type)
      tfirst = 1  # force test first

    if not tfirst: self.assign(var, default, type, 1)
    self.begin_if(prereq) # not effective if prereq == None
    if tfirst: self.assign(var, default, type, 1)

    self.cfgget_var_low(var, key, type, fmt, default,
        must, valid, desc, onerr = onerr)

    self.end_if(prereq)
    self.addln()

  def str_to_arr(self, sbuf, var, type, fmt, cnt,
      valid, cmpl, desc, onerr = "exit(1);"):
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
      self.die_if(s, "incomplete array %s at i = %%d" % var, "i", onerr = onerr)
    else:
      self.begin_if(s)
      self.addln("break;")
      self.end_if(s)
    self.validate(valid, var, fmt = "i = %d", args = "i", onerr = onerr)
    self.addln("ps += nps;")
    self.addln("}")
    self.end_if(sbufgood)

  def cfgget_sarr(self, var, key, type, fmt, default, cnt,
      must, valid, cmpl, desc, onerr = "goto ERR;"):
    # make initial assignment first
    self.init_sarr(var, default, cnt)

    # read a string from configuration file
    sbuf = "sbuf" # declare a string buffer
    self.declare_var("char *%s" % sbuf, sbuf)

    self.addln("%s = NULL;", sbuf)
    self.cfgget_var_low(sbuf, key, type, "%s", default,
        must, 1, desc, onerr = onerr)

    # parse the string to get individual items
    self.str_to_arr(sbuf, var, type, fmt, cnt, valid, cmpl,
        desc, onerr = onerr)
    self.addln("ssdelete(%s);", sbuf)

  def cfgget_flag(self, var, key, flag, type, fmt, default,
      prereq, desc, onerr = "goto ERR;"):
    '''
    assign a flag to var if prereq is true
    '''
    # get flag to a temporary variable i
    if (default not in ("1", "0") or
        type not in ("unsigned", "unsigned int", "unsigned long")):
      print "flag %s (key %s) with bad default [%s] or type [%s]" % (var, dm(key), default, type)
      raise Exception

    self.begin_if(prereq)
    ivar = "i"
    self.assign(ivar, default, type)
    self.cfgget_var_low(ivar, key, "int", fmt, default,
        "FALSE", # flag is usually optional
        "%s == 0 || %s == 1" % (ivar, ivar), desc, var_ = flag, onerr = onerr)
    self.begin_if(ivar)
    self.addln("%s |= %s;", var, flag)
    self.begin_else()
    self.addln("%s &= ~%s;", var, flag)
    self.end_if(ivar)
    self.end_if(prereq)
    self.addln()

  def begin_function(self, name, fdecl, desc, macro = None, pp = None):
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

    # save the preprocessing condition, like #if XXX
    self.pp = pp

    # declare the function
    decl = self.decl = CCodeWriter(self.nindents) # start a declaration writer
    if pp: decl.addln(pp)
    decl.add_comment("%s: %s" % (name, desc))
    decl.addln(fdecl)
    decl.addln("{")

    # create a variable list `vls'
    self.vls = {}

    # prototype, placed in the header
    hdr = self.hdr = CCodeWriter(0)  # prototype writer
    if macro: hdr.addln(macro)  # extra macro
    if pp: hdr.addln(pp)
    hdr.addln(fdecl + ";")
    if pp: hdr.addln("#endif")
    self.prototype = hdr.gets()

    # body writer
    self.body = CCodeWriter(self.nindents + 1)  # body writer

    self.funcname = name
    self.function = ""  # the content of the function

  def print_vars(self, dopp):
    for v in self.vls:
      dcl,pp = self.vls[v]
      if dopp:
        if not pp: continue
        self.decl.addln("#if %s" % pp)
        self.decl.addln(dcl+";")
        self.decl.addln("#endif")
        #print dcl, pp
      else:
        if pp: continue
        self.decl.addln(dcl+";")
        #print dcl

  def end_function(self, sret = "", silence = None):
    '''
    end a function, return header and prototype
    `silence' is a list variables to be silenced if unused, e.g.,
      ["cfg", "err!"]
    `!' at the variable end is used to force a silencing statement
    '''
    # write variable list
    if len(self.vls):
      # write normal variables first, then those with pp's
      self.print_vars(0)
      self.print_vars(1)
      self.decl.addln()
      #raw_input()

    s  = self.decl.gets()
    sb = self.body.gets()

    tail = self.tail = CCodeWriter(self.nindents + 1)
    # add void statements to avoid warnings for unused variables
    if silence: # silence unused variables
      if type(silence) == str: silence = [silence] # in case a single string
      for var in silence:
        force = var.endswith("!")
        var.rstrip("!")
        #print "var %s, in decl: %s; in body: %s;" % (var, findvar(var, s), findvar(var, sb)); raw_input()
        if force or (not findvar(var, sb) # not found in the body
            and findvar(var, s)): # but can be found in the declaration
          tail.addln("(void) %s;", var);

    s += sb

    if len(sret): tail.addln(sret);
    tail.addln("}")
    if self.pp: tail.addln("#endif")
    tail.addln()

    s += self.tail.gets()
    self.function = s
    self.funcname = None  # no longer inside the function
    self.vls = None  # kill variable list
    self.pp = None
    self.decl = self.body = self.hdr = self.tail = None

  def declare_var(self, dcl, var = None, pp = None):
    '''
    declare a variable within a function
    `dcl': declarator
    `var': variable name
    e.g. in
      double *arr;
    `var' is arr; `double *arr' is the declarator
    '''
    if self.funcname == "":
      print "need to be in a function to declare var. %s of type %s" % (var, tp)
      raise Exception

    dcl = dcl.rstrip(";")
    if var == None: # make a cheap guess of variable
      dcl0 = dcl
      pivot = dcl.find("=")
      if pivot >= 0: dcl0 = dcl[:pivot]
      arr = dcl0.split()
      var = arr[len(arr) - 1].strip("*()")
      if var == "" or var[0].isdigit() or not re.match("\w+$", var):
        print "cannot determine the var [%s], decl0 = [%s]" % (var, dcl0)
        raise Exception

    if var in self.vls: # existing variable
      olddcl, oldpp = self.vls[var]
      if dcl != olddcl:  # reject if the new type is different from the existing one
        print "var. %s is already declared as [%s], diff. from [%s]" % (var, olddcl, dcl)
        raise Exception

      # update preprocessing condition
      if oldpp: # an existing condition
        if not pp: # relieve the pp
          self.vls[var] = (olddcl, None)
        elif oldpp != pp:
          self.vls[var] = (olddcl, oldpp + " || " +pp)
        else: pass # same old condition
    else: # new variable
      self.vls[var] = dcl, pp

  def wb_arr1d(self, arr, cnt, type):
    if type not in ("int", "double"):
      print "don't know how to write type %s for %s" % (type, arr)
      raise Exception
    self.rwb_atom("w", arr, cnt = cnt)

  def wb_checkbytes(self):
    ''' write some checking bytes '''
    self.declare_var("int size");
    self.addln("size = (int) sizeof(int);")
    self.wb_var("size", "int")
    self.addln("size = (int) sizeof(double);")
    self.wb_var("size", "int")

  def wb_arr2d(self, arr, dim, tp, trim):
    major = "i"
    self.declare_var("int "+major)
    minor = "j"
    if trim:
      self.declare_var("int size")

    pb = "p%s" % tp[:1]
    self.declare_var("%s *%s" % (tp, pb))
    if trim:
      jmin = minor+"min"
      jmax = minor+"max"
      self.declare_var("int "+jmin)
      self.declare_var("int "+jmax)

    self.addln("for (%s = 0; %s < %s; %s++) {", major, major, dim[0], major)
    self.addln("%s = %s + %s * %s;", pb, arr, major, dim[1])
    if trim:
      self.addln("for (%s = %s; %s > 0 && %s[%s - 1] <= 0.0; %s--) ;",
        jmax, dim[1], jmax, pb, jmax, jmax)
      self.addln("for (%s = 0; %s < %s && %s[%s] <= 0.0; %s++) ;",
        jmin, jmin, jmax, pb, jmin, jmin)
      self.addln("if ((size = %s - %s) <= 0) continue;", jmax, jmin)
    self.wb_var(major, "int")  # current row index
    if trim:
      self.wb_var(jmin, "int")  # lowest
      self.wb_var("size", "int")
      self.wb_arr1d(pb+"+"+jmin, "size", tp)
    else:
      self.wb_arr1d(pb, dim[1], tp)
    self.addln("}")

  def rwb_obj(self, tag, var, funcall):
    cond = "0 != " + funcall % var
    self.die_if(cond, "error %s object %s" % ("reading" if tag == "r" else "writing", var),
        onerr = "goto ERR;")

  def rwb_objarr1d(self, tag, arr, dim, funcall, pp, widx, imin, imax, idx="i"):
    self.declare_var("int "+idx, pp = pp)
    if imin or imax:
      pass
    else:
      imin = 0
      imax = dim[0]
    self.addln("for (%s = %s; %s < %s; %s++) {", idx, imin, idx, imax, idx)
    if widx[0]: self.rwb_var(tag, idx, "int", match = 1) # major index
    self.rwb_obj(tag, arr+"+"+idx, funcall)
    self.addln("}")

  def rwb_objarr2d(self, tag, arr, dim, funcall, pp, widx, j0, j1):
    major = "i"
    self.declare_var("int "+major, pp = pp)
    minor = "j"
    self.addln("for (%s = 0; %s < %s; %s++) {", major, major, dim[0], major)
    if widx[0]: # write major index
      self.rwb_var(tag, major, "int", match = 1)
    arr1d = "%s+%s*%s" % (arr, major, dim[1])
    if j0 or j1:
      jmin = minor+"min"
      jmax = minor+"max"
      self.declare_var("int "+jmin, pp = pp)
      self.declare_var("int "+jmax, pp = pp)
      #self.declare_var("int size", pp = pp)
      self.addln(jmin+" = " + (j0 if j0 else "0") + ";")
      self.addln(jmax+" = " + (j1 if j1 else dim[1]) + ";")
      #self.addln("size = %s - %s;", jmax, jmin)
      size = "size"
      arr1d += "-"+jmin
    else:
      jmin = 0
      jmax = dim[1]
      size = dim[1]
    self.rwb_objarr1d(tag, arr1d, [size], funcall, pp,
        widx[1:], jmin, jmax, idx = minor)
    self.addln("}")

  def rwb_objarr(self, tag, arr, dim, funcall, pp, widx, jmin, jmax):
    ''' wrapper for reading/writing an object array '''
    self.die_if("%s == NULL" % arr, "%s is null" % arr)
    ndim = len(dim)
    if ndim == 2:
      self.rwb_objarr2d(tag, arr, dim, funcall, pp, widx, jmin, jmax)
    elif ndim == 1:
      self.rwb_objarr1d(tag, arr, dim, funcall, pp, widx, jmin, jmax)


  def rwb_atom(self, tag, var, cnt = None, onerr = None):
    ''' read a variable or array '''
    default_endian = 1 # big endian

    if type(cnt) == int and cnt <= 0: raise Exception

    # endian checked reading
    cond1 = None
    if cnt: # array
      ptr = var
      val = "*(%s)" % var
      if isint(cnt):
        icnt = int(cnt)
        if icnt <= 0: raise Exception
      else:
        cond1 = "%s > 0" % cnt
    else: # single variable
      cnt = 1
      ptr = "&"+var
      val = var
      cond1 = None
    cond2 = "endn_f%s(%s, sizeof(%s), %s, fp, %s) != %s" % (
        ("read" if tag == "r" else "write"), ptr, val, cnt,
        ("endn" if tag == "r" else default_endian),
        cnt if (type(cnt) != str or cnt.isdigit()) else "(size_t) (%s)" % cnt)
    if cond1:
      cond = "(%s)\n  && (%s)" % (cond1, cond2)
    else:
      cond = cond2
    self.begin_if(cond)
    msg = "error in "+("reading" if tag == "r" else "writing")+" "+var
    if cnt == 1:
      self.errmsg(msg)
    else:
      self.errmsg(msg+", n = %s(%%d)" % cnt, cnt)
    self.addln(onerr if onerr else "goto ERR;")
    self.end_if(cond)

  def match_var(self, x, var, prec = None):
    ''' raise an error if x is not equal to var '''
    cint = "%s != %s" % (x, var)
    cdbl = "fabs(%s - %s) > %s" % (x, var, prec)
    cond = cdbl if prec else cint
    fmt = "%g" if prec else "%d"
    self.die_if (cond,
       var+" mismatch, expect: "+fmt+", read: "+fmt+", pos: %#lx",
       args = var+", "+x+", (unsigned long) ftell(fp)",
       onerr = "goto ERR;")


  def wb_var(self, var, type): self.rwb_var("w", var, type)
  def rb_var(self, var, type, match=0, prec="1e-5", onerr=0, valid=None):
    self.rwb_var("r", var, type, match, prec, onerr, valid)

  def rwb_var(self, tag, var, tp, match=0, prec="1e-5", onerr=0, valid=None):
    ''' match: read to temporary space '''
    tpint = ("int", "unsigned", "unsigned int")
    tpdbl = ("double",)
    #if match: print "var = %s, match = %s, tag = %s, type = %s" % (var, match, tag, type(match))
    if tag == "w": match = 0
    if match: onerr = 0 # default: onerr = goto ERR;
    if tp not in (tpint + tpdbl):
      print "don't know how to read type %s for %s" % (tp, var)
      raise Exception
    if tag == "r":
      isdbl = tp in tpdbl
      tmp = "dtmp" if isdbl else "itmp"
      rvar = tmp if match else var # read to tmp var. if we need to match
      if match:
        self.declare_var("%s %s" % (tp, tmp)) # declare tmp var
        if type(match) == str: # nasty, start an branch
          #print "start if for %s" % match
          self.begin_if(match)
    else: rvar = var

    # read the variable
    self.rwb_atom(tag, rvar, onerr = onerr)

    if tag == "r":
      if match:
        self.match_var(tmp, var, prec if isdbl else None)
        if type(match) == str:  # finishing the nonmatching part
          self.begin_else()
          self.rwb_atom(tag, var, onerr = 0)
          self.end_if(match)
      self.validate(valid, var, onerr = "goto ERR;")

  def rb_arr1d(self, arr, cnt, tp, match = None, valid = None, onerr = "goto ERR;"):
    if tp not in ("int", "double"):
      print "don't know how to write type %s for %s" % (tp, arr)
      raise Exception
    if match:
      if match != 1:
        self.begin_if(match)
      self.declare_var("int i")
      self.addln("for (i = 0; i < %s; i++) {", cnt)
      self.rb_var(arr+"[i]", tp, match=1)
      self.addln("}")
      if match != 1: # conditional match
        self.begin_else()
        self.rwb_atom("r", arr, cnt = cnt)
        self.validate(valid, onerr = onerr) # only validate if we don't want to match
        self.errmsg("rb: update array %s, %%d of %s"%(arr,tp),
            "%s"%cnt )
        self.end_if(match)
    else:
      self.rwb_atom("r", arr, cnt = cnt)
      self.validate(valid, onerr = onerr)

  def rb_arr2d(self, arr, dim, tp, trim):
    '''
    read two dimensional array
    allow a compact format for array of many zeroes, `trim':
    but need an offset and size for each 1D array
    '''
    major = "i"
    self.declare_var("int "+major)
    if trim:
      self.declare_var("int size")

    pb = "p%s" % tp[0]
    self.declare_var("%s *%s" % (tp, pb))
    if trim:
      jmin = "jmin"
      self.declare_var("int %s" % jmin)

    self.addln("for (%s = 0; %s < %s; %s++) {", major, major, dim[0], major)

    # read major index
    self.rb_var("itmp", "int", match = 0, onerr = "if (feof(fp)) break;")

    # handle major index mismatch
    cidx = "itmp > %s && itmp < %s" % (major, dim[0])
    self.begin_if(cidx)
    self.addln("%s = itmp;", major)
    self.begin_else_if("%s != itmp" % major)
    self.errmsg("%s bad major index, %s: %%d, read %%d" % (arr, major),
        major+", itmp")
    self.addln("goto ERR;")
    self.end_if(cidx)

    self.addln("%s = %s + %s * %s;", pb, arr, major, dim[1])
    if trim:
      self.rb_var(jmin, "int")  # lowest
      self.insist("%s >= 0 && %s < %s" % (jmin, jmin, dim[1]),
          r"%s: base index %%d out of boudary [0, %s=%%d)" % (arr, dim[1]),
          "%s, %s" % (jmin, dim[1]), onerr = "goto ERR;")
      self.rb_var("size", "int")
      self.insist("size > 0 && %s + size <= %s" % (jmin, dim[1]),
          r"%s: invalid size %%d, %s=%%d, [0, %s=%%d)" % (arr, jmin, dim[1]),
          "size, %s, %s" % (jmin, dim[1]), onerr = "goto ERR;")
      self.rb_arr1d(pb+"+"+jmin, "size", tp)
    else:
      self.rb_arr1d(pb, dim[1], tp)

    self.addln("}") # end of the loop

  def rwb_arr(self, rw, arr, dim, tp, trim = 1,
      match = None, valid = None, onerr = "goto ERR;"):
    ndim = len(dim)
    if ndim == 2:
      if rw == "w":
        self.wb_arr2d(arr, dim, tp, trim)
      else:
        self.rb_arr2d(arr, dim, tp, trim)
        self.validate(valid, onerr = onerr)
    elif ndim == 1:
      if rw == "w":
        self.wb_arr1d(arr, dim[0], tp)
      else:
        self.rb_arr1d(arr, dim[0], tp, match, valid)
    else: raise Exception

  def rb_checkbytes(self):
    tmp = "itmp"
    self.declare_var("int %s" % tmp)
    self.declare_var("int endn")

    self.add_comment("determine file endian")
    ref = "sizeof(int)"
    cond = "(endn = endn_rmatchi(&%s, %s, fp)) < 0" % (tmp, ref)
    self.die_if(cond,
        "%s 0x%%X cannot match %s 0x%%X" % (tmp, ref),
        "(unsigned) %s, (unsigned) %s" % (tmp, ref),
        onerr = "goto ERR;")
    self.rb_var(tmp, "int")
    self.match_var(tmp, "(int) sizeof(double)")

  def test_arrempty(self, var, cnt, tp, isobj = 0):
    ''' test if an array is empty '''

    if tp in ("double", "float"):
      test = "fabs(%s[i]) > 1e-30" % var
    elif tp in ("int", "unsigned", "long"):
      test = var + "[i]"
    elif isobj:
      test = "*((char *) %s + i)" % var
      cnt = "%s * sizeof(%s)" % (cnt, tp)
    else: return None
    self.declare_var("int i")
    self.addln("for (i = %s - 1; i >= 0; i--) if (%s) break;", cnt, test)
    return "i >= 0"

  def mpibcast(self, mpirank, mpisize, var, cnt, tp, master, comm, onerr = "exit(1);"):
    ''' bcast from master to others '''
    cond = "%s > 1" % mpisize
    if cnt in (1, "1"):
      bytecnt = "sizeof(%s)" % tp
    else:
      bytecnt = "(%s)*sizeof(%s)" % (cnt, tp)
    self.begin_if(cond)
    self.die_if("MPI_SUCCESS != MPI_Bcast(%s, %s, MPI_BYTE, %s, %s)"
        % (var, bytecnt, master, comm),
        "%%3d/%%3d: failed to bcast %s (%%p), type = %s, size = %s (%%d), comm = 0x%%lX"
        % (var, tp, cnt),  # msg
        "%s, %s, %s, %s, (unsigned long) %s" % (mpirank, mpisize, var, cnt, comm), # args
        onerr = onerr)
    self.end_if(cond)

  def mpireduce(self, mpirank, mpisize, var, tmp, cnt, tp, master, comm, onerr = "exit(1);"):
    ''' sum local array to temporary array '''
    self.declare_var("int i")
    cond = "%s > 1" % mpisize
    self.begin_if(cond)
    self.die_if("MPI_SUCCESS != MPI_Reduce(%s, %s, %s, %s, MPI_SUM, %s, %s)"
        % (var, tmp, cnt, mpitype(tp), master, comm), # cond
        "%%3d/%%3d: failed to reduce %s to %s (%%p), type = %s, size = %s (%%d), comm = 0x%%lX"
        % (var, tmp, tp, cnt),  # msg
        "%s, %s, %s, %s, (unsigned long) %s" % (mpirank, mpisize, var, cnt, comm), # args
        onerr = onerr)
    self.begin_else()
    self.addln("for (i = 0; i < %s; i++) %s[i] = %s[i];", cnt, tmp, var)
    self.end_if(cond)
    # clear the local variable
    self.addln("for (i = 0; i < %s; i++) %s[i] = 0.0;", cnt, var)

