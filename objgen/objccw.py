#!/usr/bin/env python
import os, sys, re
from copy   import copy
from objext import notalways

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

  def inc(self): self.nindents += 1
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
    who.s += t
    if t.rstrip().endswith("{"): who.inc()

  def addlnraw(self, t):
    if not t.find('\n'):
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

  def remove_idle_pp(self):
    ''' remove empty preprocessor blocks 
    #if
    #else
    #endif
    '''
    endl = '\n' if self.s.endswith('\n') else ''
    lines = self.s.splitlines()
    i = 1
    # remove empty pp
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

    # merge neighoring preprocessing block
    # if their conditions are the same
    cond = None
    i = 1
    while i < len(lines):
      line = lines[i].strip()
      if line.startswith("#if "):
        if not cond: cond = line[4:].strip()
      elif line.startswith("#endif"):
        if (not cond or i >= len(lines) - 2
            or not lines[i+1].startswith("#if ")):
          cond = None
          i += 1
          continue
        ncond = lines[i+1].strip()[4:].strip()
        if cond == ncond:
          # print "removing\n%s\n%s" % (lines[i], lines[i+1]); raw_input()
          lines = lines[:i] + lines[i+2:]
        else:
          cond = None
      elif line.startswith("#el"):
        cond = None
      i += 1
    self.s = '\n'.join(lines) + endl


  def merge_if_blocks(self):
    '''
    merge neighboring if-blocks with the same condition
        if (abc) {
          ...
    X   }
    X   if (abc) {
          ...
        }
    '''
    endl = '\n' if self.s.endswith('\n') else ''
    lines = self.s.splitlines()
    cond = None
    i = 1
    # remove empty pp
    while i < len(lines):
      line = lines[i]
      pattern = r"if\s*\(.*\)\s*\{$"
      m = re.match(pattern, line.strip())
      if m:
        #print "i:%4d, %s" % (i, line); raw_input()
        if not cond: cond = line
      elif (line.strip() == "}" and cond
          and line.index("}") == cond.index("if")):
        #print "match! i=%d\n%s\n%s\n" % (i, cond, line); raw_input()
        if i < len(lines)-1 and lines[i+1] == cond: 
          #print "removing i: %d\n%s\n%s\n...%s\n" % (i, lines[i], lines[i+1], lines[i+2]); raw_input()
          lines = lines[:i] + lines[i+2:]
          continue
        cond = None  # terminate condition
      i += 1
    self.s = '\n'.join(lines) + endl

  def gets(self):
    self.remove_idle_pp()
    self.merge_if_blocks()
    return self.s
 
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
    self.declare_var("int i", pp = pp) # declare index i
    self.addln("for (i = 0; i < %s; i++)", cnt)
    self.addln(self.sindent + "%s[i] = %s;", var, default)
    self.addln()

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

  def init_darr(self, var, type, default, cnt, desc, pp = None):
    ''' write code for initializing a dynamic array '''
    if cnt.strip() != "0":
      self.alloc_darr(var, type, cnt)
      self.init_sarr(var, default, cnt, pp)

  def assign(self, var, value, type):
    ''' assignment (for simple types) '''
    if type == "char *" and value != "NULL":
      self.addln('%s = ssdup("%s");', var, value)
    else:
      self.addln("%s = %s;", var, value)

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
        r"missing var: %s, key: %s, fmt: %%%s\n" % (var, dm(key), fmt))
    else:
      # for optional variables, we print the message immediately
      cond = "cfg == NULL || " + cond
      self.msg_if(cond, 
        ('assuming default value\n' +
         'var: %s, key: %s, def: %s') % (var, dm(key), default) )
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
    self.addln()

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
      self.die_if(s, "incomplete array %s at i = %%d" % var, "i")
    else:
      self.begin_if(s)
      self.addln("break;")
      self.end_if(s)
    if valid:
      self.insist(valid, 
        r"%s: validation %s failed at i = %%d" % (var, valid),
        "i")
    self.addln("ps += nps;")
    self.addln("}")
    self.end_if(sbufgood)

  def cfgget_sarr(self, var, key, type, fmt, default, cnt, 
      must, valid, cmpl, desc):
    # make initial assignment first 
    self.init_sarr(var, default, cnt)

    # read a string from configuration file
    sbuf = "sbuf" # declare a string buffer
    self.declare_var("char *%s" % sbuf, sbuf)

    self.addln("%s = NULL;", sbuf)
    self.cfgget_var_low(sbuf, key, type, "%s", default,
        must, 1, desc)

    # parse the string to get individual items
    self.str_to_arr(sbuf, var, type, fmt, cnt, valid, cmpl, desc)
    self.addln()

  def cfgget_flag(self, var, key, flag, type, fmt, default, 
      prereq, desc):
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
        "%s == 0 || %s == 1" % (ivar, ivar), desc)
    self.begin_if(ivar)
    self.addln("%s |= %s;", var, flag)
    self.begin_else()
    self.addln("%s &= ~%s;", var, flag)
    self.end_if(ivar)
    self.end_if(prereq)
    self.addln()

  def begin_function(self, name, fdecl, desc, macro = None):
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
    self.vls = {}

    # prototype
    hdr = self.hdr = CCodeWriter(0)  # prototype writer
    if macro: hdr.addln(macro)  # extra macro
    hdr.addln(fdecl + ";")
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

  def end_function(self, sret = ""):
    # write variable list
    if len(self.vls):
      # write normal variables first, then those with pp's
      self.print_vars(0)
      self.print_vars(1)
      self.decl.addln()
      #raw_input()
    
    s  = self.decl.gets()
    s += self.body.gets()
    
    tail = self.tail = CCodeWriter(self.nindents + 1)
    if len(sret): tail.addln(sret);
    tail.addln("}\n")
    
    s += self.tail.gets()
    self.function = s
    self.funcname = None  # no longer inside the function
    self.vls = None  # kill variable list
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
      arr = dcl.split()
      var = arr[len(arr) - 1].strip("*()")
      if var == "" or var[0].isdigit() or not re.match("\w+$", var):
        print "cannot determine the var [%s]" % var
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
  
  def wb_var(self, var, type):
    if type in ("int", "unsigned", "unsigned int"):
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
    if trim:
      self.declare_var("int size")
    
    pb = "pb%s" % tp[:1]
    self.declare_var("%s *%s" % (tp, pb))
    if trim:
      self.declare_var("int imax")
      self.declare_var("int imin")

    self.addln("for (j = 0; j < %s; j++) {" % dim[0])
    self.addln("%s = %s + j * %s;" % (pb, arr, dim[1]))
    if trim:
      self.addln("for (imin = 0, i = %s-1; i >= 0 && %s[i] > 0.0; i--) ;",
        dim[1], pb)
      self.addln("imax = i + 1;")
      self.addln("for (i = 0; i < imax && %s[i] > 0.0; i++) ;",
        pb)
      self.addln("imin = i;")
      self.addln("if ((size = imax - imin) <= 0) continue;")
    self.wb_var("j",    "int")  # current row index
    if trim:
      self.wb_var("imin", "int")  # lowest
      self.wb_var("size", "int")
      self.wb_arr1d("%s+imin" % pb, "size", tp)
    else:
      self.wb_arr1d(pb, dim[1], tp)
    self.addln("}")
  
  def wb_obj(self, var, funcall):
    cond = "0 != " + funcall % var
    self.die_if(cond, "error writing object %s" % var,
        onerr = "goto ERR;")

  def wb_objarr1d(self, arr, dim, funcall, pp, widx, imin, imax):
    self.declare_var("int i", pp = pp)
    self.addln("for (i = 0; i < %s; i++) {", dim[0])
    if widx[0]: self.wb_var("i", "int") # write index
    self.wb_obj(arr + "+i", funcall)
    self.addln("}")

  def wb_objarr2d(self, arr, dim, funcall, pp, widx, imin, imax):
    self.declare_var("int j", pp = pp)
    self.addln("for (j = 0; j < %s; j++) {", dim[0])
    if widx[0]: self.wb_var("j", "int") # write index j
    self.declare_var("int imin", pp = pp)
    self.declare_var("int imax", pp = pp)
    self.declare_var("int size", pp = pp)
    arr1d = "%s+j*%s" % (arr, dim[1])
    if imin or imax:
      self.addln("imin = " + (imin if imin else "0") + ";")
      self.addln("imax = " + (imax if imax else dim[1]) + ";")
      self.addln("size = imax - imin;")
      size = "size"
      arr1d += "+imin"
    else:
      size = dim[1]
    self.wb_objarr1d(arr1d, [size], funcall, pp, 
        widx[1:], None, None)
    self.addln("}")

  def wb_objarr(self, arr, dim, funcall, pp, widx, imin, imax):
    ''' wrapper for writing an object array '''
    self.die_if("%s == NULL" % arr, "%s is null" % arr)
    ndim = len(dim)
    if ndim == 2:
      self.wb_objarr2d(arr, dim, funcall, pp, widx, imin, imax)
    elif ndim == 1:
      self.wb_objarr1d(arr, dim, funcall, pp, widx, imin, imax)

  def rb_var(self, var, type, match = 0, prec = "1e-5", tolerr = 0):
    if type in ("int", "unsigned", "unsigned int"):
      if match:
        self.declare_var("int itmp")
        self.addln("BIO_RMI(itmp, %s);", var)
      else:
        if tolerr:
          self.addln("BIO_RI_ERR(%s);", var)
        else:
          self.addln("BIO_RI(%s);", var)
    elif type == "double":
      if match:
        self.declare_var("double dtmp")
        self.addln("BIO_RMD(dtmp, %s, %s);", var, prec)
      else:
        if tolerr:
          self.addln("BIO_RD_ERR(%s);", var)
        else:
          self.addln("BIO_RD(%s);", var)
    else:
      print "don't know how to read type %s for %s" % (type, var)
      raise Exception;
  
  def rb_arr1d(self, arr, cnt, type):
    if type == "int":
      self.addraw( "BIO_RIARR(%s, %s);\n" % (arr, cnt) )
    elif type == "double":
      self.addraw( "BIO_RDARR(%s, %s);\n" % (arr, cnt) )
    else:
      print "don't know how to write type %s for %s" % (type, arr)
      raise Exception

  def rb_arr2d(self, arr, dim, tp, trim):
    self.declare_var("int j")
    self.declare_var("int i")
    if trim:
      self.declare_var("int size")
    
    pb = "pb%s" % tp[0]
    self.declare_var("%s *%s" % (tp, pb))
    if trim:
      self.declare_var("int imin")

    self.addln("for (j = 0; j < %s; j++) {" % dim[0])

    # read major index
    self.rb_var("itmp", "int", match = 0, tolerr = 1)
    ceof = "feof(fp)"
    self.begin_if(ceof) # EOF encountered, normal quit */
    self.addln("break;")
    self.begin_else_if("err")
    self.errmsg(r"%s error: j = %%d, itmp = %%d\n" % arr,
        "j, itmp")
    self.addraw("goto ERR;\n")
    self.end_if(ceof)

    # handle major index mismatch
    cidx = "itmp > i && itmp < %s" % dim[0]
    self.begin_if(cidx)
    self.addln("i = itmp;")
    self.begin_else()
    self.die_if("i != itmp", 
        "%s major index mismatch at i = %%d, (read) %%d" % arr,
        "i, itmp", onerr = "goto ERR;")
    self.end_if(cidx)

    self.addln("%s = %s + j * %s;" % (pb, arr, dim[1]))
    if trim:
      self.rb_var("imin", "int")  # lowest
      self.insist("imin >= 0 && imin < %s" % dim[1],
          r"%s: base index %%d out of boudary [0, %s=%%d)" % (arr, dim[1]),
          "imin, %s" % dim[1], onerr = "goto ERR;")
      self.rb_var("size", "int")
      self.insist("size > 0 && imin + size <= %s" % dim[1],
          r"%s: invalid size %%d, imin=%%d, [0, %s=%%d)" % (arr, dim[1]),
          "size, imin, %s" % dim[1], onerr = "goto ERR;")
      self.rb_arr1d("%s+imin" % pb, "size", tp)
    else:
      self.rb_arr1d(pb, dim[1], tp)

    self.addln("}") # end of the loop

  def rwb_arr(self, arr, dim, tp, trim, rw):
    ndim = len(dim)
    if ndim == 2:
      if rw == "w":
        self.wb_arr2d(arr, dim, tp, trim)
      else:
        self.rb_arr2d(arr, dim, tp, trim)
    elif ndim == 1:
      if rw == "w":
        self.wb_arr1d(arr, dim[0], tp)
      else:
        self.rb_arr1d(arr, dim[0], tp)
    else: raise Exception

  def rb_checkbytes(self):
    self.declare_var("int size")
    self.declare_var("int endn")
    self.declare_var("int err")

    self.add_comment("determine file endian")
    val = "size"
    ref = "sizeof(int)"
    cond = "(endn = endn_rmatchi(&%s, %s, fp)) < 0" % (val, ref)
    self.die_if(cond,
        "%s 0x%%X cannot match %s 0x%%X" % (val, ref),
        "(unsigned) %s, (unsigned) %s" % (val, ref), 
        onerr = "goto ERR;")
    self.addraw("BIO_RMI(size, sizeof(double));\n")

  def rb_obj(self, var, funcall):
    cond = "0 != " + funcall % var
    self.die_if(cond, "error reading object %s" % var,
        onerr = "goto ERR;")

  def rb_objarr1d(self, arr, dim, funcall, pp, widx, imin, imax):
    self.declare_var("int i", pp = pp)
    self.addln("for (i = 0; i < %s; i++) {", dim[0])
    if widx[0]: self.rb_var("i", "int", match = 1)
    self.rb_obj(arr + "+i", funcall)
    self.addln("}")

  def rb_objarr2d(self, arr, dim, funcall, pp, widx, imin, imax):
    self.declare_var("int j", pp = pp)
    self.addln("for (j = 0; j < %s; j++) {", dim[0])
    if widx[0]: # write index j
      self.rb_var("j", "int", match = 1) 
    arr1d = "%s+j*%s" % (arr, dim[1])
    size = dim[1]
    self.rb_objarr1d(arr1d, [size], funcall, pp, 
        widx[1:], None, None)
    self.addln("}")

  def rb_objarr(self, arr, dim, funcall, pp, widx, imin, imax):
    ''' wrapper for reading an object array '''
    self.die_if("%s == NULL" % arr, "%s is null" % arr)
    ndim = len(dim)
    if ndim == 2:
      self.rb_objarr2d(arr, dim, funcall, pp, widx, imin, imax)
    elif ndim == 1:
      self.rb_objarr1d(arr, dim, funcall, pp, widx, imin, imax)

  '''
  def rb_atom(self, var, cnt, onerr, onsuc):
    # endian checked reading
    cond = "(%s > 0) && end_fread(%s, sizeof(*(%s)), %s, fp, endn) != %s" % (
        cnt, var, var, cnt, cnt)
    self.begin_if(cond)
    self.errmsg("error while reading %s, n = %s" % (var, cnt))
    self.end_if(cond)
  '''

