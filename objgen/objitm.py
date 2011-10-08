#!/usr/bin/env python
import os, sys, re 
from copy    import * 
from objdcl  import CDeclaratorList
from objcmt  import CComment
from objpre  import CPreprocessor
from objcow  import CCodeWriter
from objext  import *

simple_types = ("int", "unsigned int", "unsigned", 
  "long", "unsigned long", "real", "double", "float", 
  "char *")

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
    if src == None: return
    self.src = src
    self.begin = copy(p)
    self.parse(src, p)

  def duckcopy(s, dcl = 1):
    ''' 
    cheaper version of deepcopy 
    if dcl = 0, .decl and .dl are not copied
    '''
    c = Item(None, None)
    c.empty = s.empty
    c.src = s.src
    c.begin = copy(s.begin)
    c.end = copy(s.end)
    c.pre = deepcopy(s.pre) # this is cheap, optimized
    c.cmt = deepcopy(s.cmt) # also cheap
    if hasattr(s, "gtype"): c.gtype = s.gtype
    if hasattr(s, "cmds"): c.cmds = deepcopy(s.cmds)
    c.itlist = s.itlist  # only a reference
    if dcl:
      if hasattr(s, "decl"): c.decl = deepcopy(s.decl)
      if hasattr(s, "dl"):   c.dl   = deepcopy(s.dl)
    return c

  def __deepcopy__(s, memo):
    c = Item(None, None)
    memo[id(s)] = c
    for n, v in s.__dict__.iteritems():
      setattr(c, n, deepcopy(v, memo))
    return c

  def parse(self, src, p):
    '''
    parse an item aggressively
    '''
    s = p.skipspace(src)
    if s[0] == "}": return -1
    dodcl = docmt = 1
 
    # try to see it's a preprocessor
    pre = CPreprocessor(src, p)
    if not pre.isempty():
      #print "pos %s, preprocessor %s" % (p, pre.raw)
      #raw_input()
      self.empty = 0
      self.pre = pre
      self.dl = self.cmt = None
      dodcl = docmt = 0
      if pre.cmt: docmt = 1
    else:
      self.pre = None

    if dodcl:
      # print "try to get a declaration from %s" % p; raw_input()
      self.decl = None
      dl = CDeclaratorList(src, p)
      self.dl = None if dl.isempty() else dl
      
    if docmt:
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

  def getraw(self):
    return self.begin.gettext(self.src, self.end).strip()

  def __repr__(self):
    return self.getraw()

  def __str__(self):
    return "pre: %5s, dl: %5s, decl: %5s, cmt: %5s, raw = [%s]" % (
        self.pre != None, self.dl != None, 
        self.decl != None, self.cmt != None, 
        self.getraw() )

  def name(self):
    if self.decl: 
      return "item: %s" % self.decl.name
    elif self.pre:
      return "prep: %s" % self.pre.raw
    else:
      return "(aux)"


  def get_gtype(self, offset = 0):
    '''
    get a generic type, and type name
    it need commands
    `offset' is the starting offset from the bottom
    '''
    types = self.decl.types
    cmds = self.cmds

    n = len(types)
    if offset >= n: return None

    nlevels = n - offset
    if nlevels == 1:
      ret = types[offset]
    elif types[offset].startswith("array"):
      ''' add a command for array count '''
      dim = []
      for tp in types[offset:]:
        if tp.startswith("array"):
          dim += [ tp[5:] ]
        else: break
      dim.reverse()  # reverse the order 
      if offset == 0: # 
        self.cmds["dim"] = dim
        self.cmds["cnt"] = arrcnt = '*'.join(dim) # override cnt
        #print "static array cnt has been set to %s (%s) for %s" % (arrcnt, dim, self)
        #raw_input()
      ret = "static array"
    elif types[offset] == "pointer":
      #towhat = ' '.join(types[1:])
      if types[offset+1] == "function":
        ret = "function pointer"
      elif nlevels == 2 and types[offset+1] == "char":
        #print "a string is encountered %s" % self.decl.name; raw_input
        ret = "char *"
      elif "obj" in cmds:
        ret = "object array" if cmds["cnt"] else "object pointer"
      elif "cnt" in cmds:
        ret = "dynamic array"
      else:
        if (offset == 0 and "cnt" not in cmds 
          and nlevels == 2 
          and types[offset+1] in simple_types):
          print "Warning: `%s' is considered as a pointer" % self.decl.name,
          print "  set $cnt if it is actually an array"
        ret = "pointer"

    # debug code
    if self.decl.name == "XXXmbYYY":  
      print "gtype [%s] with %s, %s" % (ret, offset, types); 
      print "commands: %s" % self.cmds
      raw_input()        
    return ret

  def get_elegtype(self):
    ''' return gtype for simple, or gtype of elements for arrays '''
    gtype = self.gtype
    if gtype in ("static array", "dynamic array"):
      return self.get_gtype(offset = 1)
    else:
      return gtype

  def fill_def(self):
    ''' set 'def' value '''
    if not self.decl or self.cmds["def"]: return
    self.cmds["def"] = type2zero(self.get_elegtype())
    #print "set default %s for gtype %s, var %s" % (defval, gtype, self.decl.name)
    #raw_input()

  def add_item_to_list(self, dcl): 
    ''' add an item with decl being dcl '''
    it = self.duckcopy(dcl = 0) # make a copy cheaper than deepcopy
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

  def fill_key(it):
    ''' 
    add default key for configuration file reading
    unless "key" is specified, it is calculated as
    1. it.decl.name is the first guess
    2. if it starts with kdepfx, remove it
    3. add kprefix
    4. append kargs, after a `#'
    '''

    # skip if it is not a declaration
    #if not it.decl or "key" in it.cmds: return
    if not it.decl: return

    key0 = it.cmds["key"]
    cfgkey = key0 if key0 else it.decl.name
    
    if not key0 or it.cmds["flag"]:
      pm = it.cmds["kdepfx"]
      if pm and cfgkey.startswith(pm):
        cfgkey = cfgkey[len(pm):]
      pfx  = it.cmds["kprefix"]
      if pfx:  cfgkey = pfx + cfgkey
    # append key argument
    args = it.cmds["kargs"]
    if args: cfgkey += "#" + args
    it.cmds["key"] = cfgkey

  def fill_dim(it):
    ''' 
    test if cnt is empty 
    for two dimensional array parse `cnt', which looks like "n1, n2"
    into dim[2] = ["n1", "n2"]
    '''
    if not it.decl or not it.gtype.endswith("array"): return
    if "cnt" not in it.cmds:
      print "item %s: array requires a count" % it.decl.name
      raise Exception
    cnt = it.cmds["cnt"].strip()
    if len(cnt) == 0:
      print "item %s: $cnt command cannot be empty! comment: [%s]" % (it.decl.name, it.cmt)
      raise Exception
    if it.gtype == "static array": # static type
      return
    dim = [s.strip() for s in cnt.split(",")] # split into dimensions
    
    if len(dim) > 2: #
      print "item %s: do not support > 2 dimensional array!"
      raise Exception
    elif len(dim) > 1:
      for i in range(len(dim)):
        if (re.search(r"[\+\-]", dim[i])
            and not (dim[i][0] == "(" and dim[i][-1:] == ")") ):
          print "for %s add parentheses around %s, %d" % (
              it.decl.name, dim[i], i); # raw_input()
          dim[i] = "(" + dim[i] + ")"
    cnt = '*'.join(dim)
    it.cmds["cnt"] = cnt
    it.cmds["dim"] = dim

  def fill_io(it):
    ''' expand the io string '''
    if not it.decl: return
    # default I/O mode, bt for array, cbt for others
    if it.cmds["usr"]:
      iodef = ""
    elif it.gtype in ("dynamic array", "object array", "object pointer"):
      iodef = "bt"
    elif it.gtype in simple_types:
      iodef = "c" 
    else:
      iodef = ""

    # assume iodef if necessary
    io0 = io = it.cmds["io"]
    io = io if io != None else iodef
    io = io.strip()

    # a human shortcut
    if io.startswith("no"): io = ""
    elif io in ("all", "everything"): io = "cbt"

    usr = it.cmds["usr"]

    # create itemized io
    it.cmds["io_cfg"] = ("c" in io)
    if usr == None or io0 != None: # normal var or user var with explicit $io command
      it.cmds["io_bin"] = ("b" in io)
    else:
      it.cmds["io_bin"] = (usr in ("bintmp", "rbtmp", "wbtmp"))
    it.cmds["io_txt"] = ("t" in io)

  def fill_test(it):
    ''' 
    fill test conditions 
    note: we assign commands even if there's no declaration
    becomes a stand-alone comment can contain assignment
    or other commands
    '''
    # turn on $tfirst for dummy variables
    if ("tfirst" not in it.cmds and it.isdummy):
      it.cmds["tfirst"] = 1

  def prepare_flag(it):
    ''' 
    parse the original field into FLAG and flagval 
    attach prefix
    '''
    flag = it.cmds["flag"]
    if not flag: return
    # attach flag prefix
    if it.cmds["flagprefix"]:
      prefix = it.cmds["flagprefix"]
      if not flag.startswith(prefix):
        flag = prefix + flag
    it.cmds["flag"] = flag

    # assign the default flagvar
    if "flagvar" not in it.cmds:
      it.cmds["flagvar"] = "@flags"
    #print "%s %s" % (it.cmds["flag"], it.cmds["flagval"])
    #raw_input()
    
  def getflag(it):
    if "flag" not in it.cmds: return None
    return it.cmds["flag"]

  def get_obj_fprefix(it):
    ''' get prefix of an object '''
    if not it.cmds["obj"]:
      print "item is not an object %s" % it
      raise Exception
    pfx = it.cmds["objfprefix"]
    if pfx == None: # make a guess
      tp = it.decl.datatype
      if tp.endswith("_t"): tp = tp[:-2]
      pfx = tp + "_"
    return pfx

  def getargs(it, tag):
    ''' get additional arguments to be passed to an object initializer '''
    if not it.cmds["obj"]:
      print "item is not an object %s" % it
      raise Exception
    args = it.cmds[("%sargs" % tag)]
    if not args: return ""
    arr = args.strip().split(",")
    arr = [""] + [s.strip() for s in arr]
    args = ", ".join(arr)
    # print args; raw_input()
    return args

  def declare_var(it, cow, block, tab, offset):
    ''' declare a variable '''
    decl0 = decl1 = scmt = ""
    # 1. handle pre
    if it.pre and not it.cmds["key"]:
      cow.dump_block(block)
      block = [] # empty the block
      cow.addln("#" + it.pre.raw)
      return block
    
    # 2. handle declaration
    usr = it.cmds["usr"]
    # we need to declare $usr:cfg (and possibly $usr:xxxdup) 
    if usr != None and usr != 1 and usr != "cfg":
      return block
    if it.decl:
      decl0 = it.decl.datatype
      if it.isdummy: return block
      decl1 = it.decl.raw
    else:
      # skip a comment if it contains instructions
      if it.cmds["assert"] or it.cmds["call"] or it.cmds["flag"]:
        return block
      #raw_input ("%s %s" % (decl0, decl1))

    # 3. handle comment
    desc = it.cmds["desc"]
    if len(desc) > 0:
      scmt = "/* " + desc + " */"
    
    #print "decl: %s, cmds: %s" % (it.decl, it.cmds)
    #raw_input()

    # 4. output the content
    if len(decl0) > 0:  # put it to a code block
      block +=  [(decl0, decl1+';', scmt)] # just buffer it
    else:
      cow.dump_block(block, tab, offset)
      block = [] # empty the block
      # also print this line
      if len(scmt) > 0: 
        cow.addln(scmt)
    return block

  def cfgget_var(it, cow, ptrname):
    ''' read a variable from config. '''
    if it.pre and not it.cmds["key"]: return
  
    usr = it.cmds["usr"]
    if usr == "parent" or (usr not in (None, 0, 1) and 
        usr.endswith("tmp") and usr != "cfgtmp"):
      return
    
    key     = it.cmds["key"]
    defl    = it.cmds["def"]
    prereq  = it.cmds["cfgprereq"]
    if not prereq:  prereq = it.cmds["prereq"]
    tfirst  = it.cmds["tfirst"]
    valid   = it.cmds["cfgvalid"]
    if valid == None: valid = it.cmds["valid"]
    must    = it.cmds["must"]
    cnt = it.cmds["cfgcnt"]
    if cnt == None: cnt = it.cmds["cnt"]
    desc    = it.mkdesc(1)
    flag    = it.cmds["flag"]
    usrval  = it.cmds["usr"]
    pp = it.cmds["#if"]
    if pp: cow.addln("#if %s", pp)
    ###raw_input("if pp=[%s], item=%s"%(pp, it))

    passme = 0
    if not it.decl: # stand-alone comment
      insist = it.cmds["assert"]
      call = it.cmds["call"]
      if insist or call:
        cow.add_comment(desc)
        cow.insist(insist, it.cmds["desc"])
        if call: cow.addln(call + ";")
      passme = 1
    if flag and (not key or key == "flags"):
      passme = 1
    if usrval and usrval not in ("cfg", "cfgdup", 1): # usr variables 
      passme = 1

    if not passme: 
      varnm = it.decl.name
      varname = ptrname + "->" + varnm
      # add a remark before we start
      cow.add_comment(desc)

    if passme:
      pass
    elif flag: # input flag
      #raw_input("flag with key %s" % key)
      flag = it.getflag()
      fmt = type2fmt_s(it.gtype)
      cow.cfgget_flag(varname, key, flag, it.gtype, fmt,
          defl, prereq, desc)
  
    elif it.gtype == "object pointer":
      if notalways(prereq): cow.addln(varname+" = NULL;")
      cow.begin_if(prereq)
      fpfx = it.get_obj_fprefix()
      s = "(%s = %scfgopen(cfg%s)) == NULL" % (varname, 
        fpfx, it.getargs("cfg"))
      cow.die_if(s, r"failed to initialize %s\n" % varname, 
          onerr = "goto ERR;")
      cow.validate(valid, varname)
      cow.end_if(prereq)

    elif it.gtype == "object array":
      if notalways(prereq): cow.addln(varname+" = NULL;")
      cow.begin_if(prereq)
      fpfx = it.get_obj_fprefix()
      etp = it.get_gtype(offset = 1)
      cow.alloc_darr(varname, etp, cnt)
      cow.declare_var("int i;", pp = pp)
      cow.addln("for (i = 0; i < %s; i++) {" % cnt)
      s = "0 != %scfgopen_low(%s+i, cfg%s)" % (
        fpfx, varname, it.getargs("cfg"))
      cow.die_if(s, r"failed to initialize %s[%%d]\n" % varname, "i", 
          onerr = "goto ERR;")
      cow.addln("}\n")
      cow.validate(valid, varname)
      cow.end_if(prereq)

    elif it.gtype == "dynamic array":
      if notalways(prereq): cow.addln(varname+" = NULL;")
      cow.begin_if(prereq)
      cow.init_darr(varname, it.get_gtype(offset = 1),
          defl, cnt, valid, desc, pp)
      cow.end_if(prereq)
    
    elif it.gtype == "static array":
      etp = it.get_gtype(offset = 1)
      #print "%s: static array of element = [%s] io = %s, defl = [%s]" % (varnm, etp, it.cmds["io_cfg"], defl); raw_input()
      if not it.cmds["io_cfg"]: 
        cow.init_sarr(varname, defl, cnt, pp)
      else:
        fmt  = type2fmt_s(etp)
        cmpl = it.cmds["complete"]
        cow.cfgget_sarr(varname, key, etp, fmt, defl, cnt, 
            must, valid, cmpl, desc)

    else: # regular variable
      if usrval:
        cow.assign(varname, varnm, it.gtype)
      elif not it.cmds["io_cfg"]: # assign default value
        if defl and len(defl): 
          cow.assign(varname, defl, it.gtype);
      else:
        fmt = type2fmt_s(it.gtype)
        cow.cfgget_var(varname, key, it.gtype, fmt, 
          defl, must, prereq, tfirst, valid, desc)

    if pp: cow.addln("#endif")
    

  def rwb_var(it, cow, rw, varname):
    desc = it.mkdesc(1)
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    bincnt = it.cmds[rw+"bcnt"]
    if bincnt == None: bincnt = it.cmds["bincnt"]
    prereq = it.cmds[rw+"bprereq"]
    if prereq == None: prereq = it.cmds["binprereq"]
    defl = it.cmds[rw+"bdef"] # try rbdef and rwbdef first
    if defl == None: defl = it.cmds["def"]
    usrval = it.cmds["usr"]
    if rw == "r":
      readwrite = "read"
      verify = it.cmds["rbverify"]
      valid = it.cmds["rbvalid"]
      if not valid:
        valid = it.cmds["valid"]
    else:
      readwrite = "write"
      verify = 0
      valid = None
    prep = it.cmds[rw+"bprep"]
    iobin = it.cmds["io_bin"]

    # skip normal non-io variables
    if not iobin and not usrval: return

    pp = it.cmds["#if"]
    if pp: cow.addln("#if %s", pp)
    if notalways(prereq): cow.begin_if(prereq)

    if not it.decl:
      print "no declaration gtype = %s" % it.gtype
      raise Exception

    # init. temporary variables before writing
    passme = 0
    usrinit = 0
    if usrval in ("bintmp", rw+"btmp"):
      cow.declare_var(it.decl.datatype+" "+it.decl.raw, it.decl.name)
      usrinit = 1
    elif usrval != None and usrval.endswith("tmp"):
      passme = 1
    elif usrval == "cfg":
      #print "usrval %s, %s" % (it.decl.name, usrval); raw_input()
      passme = 1

    if not iobin: passme = 1

    if not passme:
      if desc: cow.add_comment(desc);
      if prep: cow.addln(prep + ";")

    if usrinit and defl: 
      # we always assign default values before reading / writing 
      if it.gtype == "static array":
        cow.init_sarr(varname, defl, cnt)
      else:
        cow.addln("%s = %s;", varname, defl)
      #print "varname = %s" % varname; raw_input()

    if passme:
      pass
    elif it.gtype == "dynamic array":
      # support nasty flag $bincnt
      if len(dim) == 1 and bincnt:
        dim[0] = bincnt
      cow.rwb_arr(rw, varname, dim, it.decl.datatype, trim = 1, 
          match = verify, valid = valid)
    elif it.gtype in ("object array", "object pointer"):
      fpfx = it.get_obj_fprefix();
      extra = ", endn" if rw == "r" else ""
      extra += it.getargs(rw+"b")
      #print "extra args = %s" % extra; raw_input()
      funcall = "%s%sbin_low(%%s, fp, ver%s)" % (fpfx, readwrite, extra);
      if it.gtype == "object array":
        imin = it.cmds["bin_jmin"]
        imax = it.cmds["bin_jmax"]
        cow.rwb_objarr(rw, varname, dim, funcall, pp, [1, 1], imin, imax)
      else:
        cow.rwb_obj(rw, varname, funcall)
    elif it.decl.datatype == "char" and it.gtype in ("char *", "static array"):
      if bincnt: 
        cnt = bincnt
      cow.die_if ("(size_t) (%s) != f%s(%s, 1, %s, fp)" % (cnt, readwrite, varname, cnt),
        "cannot "+readwrite+" string of %d for "+varname, cnt, 
        onerr = "goto ERR;");
      cow.validate(valid);
    else:
      cow.rwb_var(rw, varname, it.gtype, verify, valid = valid)

    if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")

   
  def clear_var(it, cow, ptr):
    if not it.decl or it.isdummy: return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    defl = it.cmds["def"]
    #defl = it.get_zero()  zero is unreliable, e.g., sth needs to be 1.0
    prereq = it.cmds["prereq"]

    pp = it.cmds["#if"]
    if pp: cow.addln("#if %s", pp)
    if notalways(prereq): cow.begin_if(prereq)

    if it.gtype in ("dynamic array", "static array"):
      cow.init_sarr(varname, defl, cnt, pp = pp)
    elif it.gtype in ("object array", "object pointer"):
      fpfx = it.get_obj_fprefix();
      cond = "%s != NULL" % varname
      cow.begin_if(cond)      
      funcall = "%sclear(%%s)" % fpfx      
      if it.gtype == "object array":
        cow.declare_var("int i", pp = pp) # declare index i
        cow.addln("for (i = 0; i < %s; i++)", cnt)
        cow.addln(cow.sindent + funcall % (varname+"+i") + ";")
      else:
        cow.addln(funcall % varname + ";")
      cow.end_if(cond)
    else:
      cow.assign(varname, defl, it.gtype)

    if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")

  def close_var(it, cow, ptrname, maxwid):
    ''' free a pointer / close an object '''
    if not it.decl or it.isdummy: return
    usr = it.cmds["usr"]
    if usr == "parent": return
    if type(usr) == str and usr.endswith("ref"): return
    if it.gtype not in ("char *", "dynamic array",
        "object array", "object pointer"): return

    call = it.cmds["closecall"]
    varnm = it.decl.name;
    if varnm == "iopref":
      print "usr : %s" % usr; raw_input
    varname = "%s->%s" % (ptrname, varnm)

    pp = it.cmds["#if"]
    prereq = it.cmds["prereq"]
    if pp: cow.addln("#if %s", pp)
    if notalways(prereq): cow.begin_if(prereq)    
 
    if call: cow.addln(call+";")

    handled = 0
    #if varnm == "cache":
    #  print ("destroying %s, gtype: %s" % (it.decl.name, it.gtype)); 
    #  raw_input()
    if it.gtype == "char *": 
      funcfree = "ssdelete"
    elif it.gtype == "dynamic array":
      funcfree = "free"
    elif it.gtype == "object pointer":
      fpfx = it.get_obj_fprefix()
      funcfree = "%sclose" % fpfx
    elif it.gtype == "object array": 
      fpfx = it.get_obj_fprefix()
      funcfree = "%sclose_low" % fpfx
      cnt = it.cmds["cnt"]
      cow.declare_var("int i")

      cond = "%s != NULL" % varname
      cow.begin_if(cond)
      cow.addln("for (i = 0; i < %s; i++) {" % cnt)
      cow.addln("%s(%s + i);", funcfree, varname) 
      cow.addln("}");
      cow.addln("free(%s);", varname);
      cow.end_if(cond)
      handled = 1

    if not handled:
      cow.addln("if (%-*s != NULL) %s(%s);",
          maxwid, varname, funcfree, varname)
    
    if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")
    

  def manifest_var(it, cow, ptr):
    arrnum = 5 # elements for printing arrays

    if not it.decl: return
    if it.cmds["usr"] not in (None, 1): return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    desc = it.mkdesc(1)
    prereq = it.cmds["prereq"]
    flag = it.cmds["flag"]
    key  = it.cmds["key"]
    if not flag and it.isdummy: return

    pp = it.cmds["#if"]
    if pp: cow.addln("#if %s", pp)
    if notalways(prereq): cow.begin_if(prereq)

    cow.add_comment(desc)
    if it.gtype in ("dynamic array", "static array"):
      isarr = 1
      etype = it.get_elegtype()
      try:
        fmt = type2fmt_p(etype, varname)
        cow.addln(r'fprintf(fp, "%s: %s of %s: ");', 
          varname, it.gtype, cnt)
      except TypeError:
        cow.addln(r'fprintf(fp, "%s: %s of %s %%p (unknown type: %s)", %s);',
          varname, it.gtype, cnt, etype, varname)
        isarr = 0
      if cnt == "0": isarr = 0
      if isarr:
        emptest = cow.test_arrempty(varname, cnt, etype)
        if emptest: 
          cow.begin_if(emptest)
          cow.addln('if ((arrmax < 0 || arrmax > 3) && %s > 6)\n\tfprintf(fp, "\\n");',
              cnt)
        cow.declare_var("int i", pp = pp) # declare index i
        cow.declare_var("int pacnt", pp = pp)
        
        cow.addln("for (pacnt = 0, i = 0; i < %s; i++) {", cnt)
        cond = 'i == arrmax && i < %s-arrmax' % cnt
        cow.begin_if(cond)
        cow.addln(r'if (arrmax > 3 && pacnt %% 10 != 0) fprintf(fp, "\n");')
        cow.addln(r'fprintf(fp, "..., ");')
        cow.addln(r'if (arrmax > 3) fprintf(fp, "\n");');
        cow.end_if(cond)
        cow.addln(r'if (arrmax >= 0 && i >= arrmax && i < (%s-arrmax)) continue;', cnt)
        cow.addln('fprintf(fp, "%s, ", %s);', fmt, varname + "[i]")
        cow.addln(r'if (++pacnt %% 10 == 0) fprintf(fp, "\n");')
        cow.addln("}")
        cow.addln(r'if (pacnt %% 10 != 0) fprintf(fp, "\n");')
        if emptest:
          cow.begin_else()
          cow.addln(r'fprintf(fp, " {0}\n");')
          cow.end_if(emptest)
      else:
        cow.addln(r'fprintf(fp, "\n");')

    elif it.gtype in ("object array", "object pointer"):
      fpfx = it.get_obj_fprefix();
      cond = "%s != NULL" % varname
      cow.begin_if(cond)      
      funcall = "%smanifest(%%s, fp, arrmax)" % fpfx
      if it.gtype == "object array":
        cow.addln(r'fprintf(fp, "%s: %s array of %s:");', 
          varname, it.decl.datatype, cnt)
        #raw_input("%s: %s x %s" % (varname, it.decl.datatype, cnt))
        emptest = cow.test_arrempty(varname, cnt, it.decl.datatype, isobj = 1)
        if emptest: # test if the object array is empty
          cow.begin_if(emptest)
          cow.addln(r'fprintf(fp, "\n");')
        cow.declare_var("int i", pp = pp) # declare index i

        cow.addln("for (i = 0; i < %s; i++) {", cnt)
        cond = 'i == arrmax && i < %s-arrmax' % cnt
        cow.begin_if(cond)
        cow.addln(r'fprintf(fp, "\n...\n");')
        cow.end_if(cond)
        cow.addln(r'if (arrmax >= 0 && i >= arrmax && i < (%s-arrmax)) continue;', cnt)
        cow.addln(r'fprintf(fp, "%s[%%d]:\n", i);', varname)
        cow.addln(funcall % (varname+"+i") + ";")
        cow.addln("}")
        if emptest:
          cow.begin_else()
          cow.addln(r'fprintf(fp, " {0}\n");')
          cow.end_if(emptest)
      else:
        cow.addln(r'fprintf(fp, "%s: %s to %s\n");',
          varname, it.gtype, it.decl.datatype)
        cow.addln(funcall % varname + ";")
      cow.end_if(cond)

    elif flag:
      flag = it.getflag()
      cow.addln('fprintf(fp, "%s & %s (0x%%X, %s): %%s\\n",\n'
          + '\t%s, (%s & %s) ? "on" : "off");',
        varname, flag, escape(key), flag, varname, flag)

    elif it.gtype in ("pointer",):
      print "skip var. [%s] of type [%s]" % (varname, it.gtype); # raw_input()
    
    else:
      try:
        fmt = type2fmt_p(it.gtype, varname)
        cow.addln(r'fprintf(fp, "%s: %s, %s\n", %s);',
          varname, it.gtype, fmt, varname)
      except TypeError:
        cow.addln(r'fprintf(fp, "%s: %s, 0x%%X\n", (unsigned) %s);',
          varname, it.gtype, varname)
    cow.addln()

    if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")


  def mpipp(it):
    ''' basically cmds["#if"], but avoids defined(USE_MPI) '''
    pp = it.cmds["#if"]
    if not pp: return None
    from objgen import USE_MPI
    return pp if pp != "defined(%s)"%USE_MPI else None

  def initmpi_var(it, cow, ptr):
    if not it.decl or it.isdummy: return
    if it.cmds["usr"] not in (None, 1): return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    rank = ptr + "->mpi_rank"
    size = ptr + "->mpi_size"
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    desc = it.mkdesc()
    mpicnt = it.cmds["mpicnt"]
    if mpicnt: cnt = mpicnt;
    prereq = it.cmds["prereq"]
    mpi = it.cmds["mpi"]
    etype = it.get_elegtype()

    pp = it.mpipp()
    if pp: cow.addln("#if %s", pp)
    #if notalways(prereq): cow.begin_if(prereq)

    notmaster = "%s->mpi_rank != %s" % (ptr, MASTERID)
    if mpi in ("alloc", "2"): # allocate memory only, temporary variables
      if it.gtype not in ("dynamic array",):
        print "cannot just allocate space for %s (%s)" % (it.decl.name, it.gtype)
        raise Exception
      cow.begin_if(notmaster)
      cow.add_comment(desc)
      cow.alloc_darr(varname, etype, cnt)
      cow.end_if(notmaster)
    
    # an array to be bcasted, allocate space on nonmasters
    # receive from master
    elif mpi in (1, "1", "bcast"):
      if it.gtype == "static array":
        pass  # no need to do anything
      elif it.gtype == "dynamic array":
        cow.add_comment(desc)
        cow.begin_if(notmaster)
        cow.alloc_darr(varname, etype, cnt)
        cow.end_if(notmaster)
        cow.mpibcast(rank, size, varname, cnt, etype, MASTERID, "comm")
        cow.addln()
      elif it.gtype in ("object array", "object pointer"):
        isarr = (it.gtype == "object array")
        if not isarr: cnt = "1"
        cow.begin_if(notmaster)
        cow.add_comment(desc)
        cow.alloc_darr(varname, it.decl.datatype, cnt)
        cow.end_if(notmaster)
        if isarr: 
          cow.declare_var("int i", pp=pp)
          cow.addln('for (i = 0; i < %s; i++) {' % cnt)
        fpfx = it.get_obj_fprefix()
        cow.addln("%sinitmpi(%s, comm);", fpfx, varname+("+i" if isarr else ""))
        if isarr: cow.addln('}\n')
      elif it.gtype == "char *":
        cow.add_comment(desc)
        cow.declare_var("int i", pp=pp)
        cow.addln("i = strlen(%s);", varname)
        cow.mpibcast(rank, size, "&i", 1, "int", MASTERID, "comm")
        cow.begin_if(notmaster)
        cow.addln("%s = ssnew(i);" % varname)
        cow.end_if(notmaster)
        cow.mpibcast(rank, size, varname, "i", "char", MASTERID, "comm")
        cow.addln()
      else:
        print "cannot bcast %s of type %s" % (it.decl.name, it.gtype)
        raise Exception

    elif mpi in (0, "0", None): # mpi == None
      if it.gtype in ("dynamic array", "char *", "object array", "object pointer"):
        cow.begin_if(notmaster)
        cow.addln('%s = NULL;', varname)
        cow.end_if(notmaster)
      elif it.gtype == "static array" and etype == "char *":
        cow.begin_if(notmaster)
        cow.declare_var("int i", pp=pp)
        cow.addln('for (i = 0; i < %s; i++) {\n%s[i] = NULL;\n}', cnt, varname)
        cow.end_if(notmaster)
      # do nothing otherwise
    
    else:
      print "unknown mpi tag %s" % mpi
      raise Exception

    #if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")

  def mpitask_var(it, cow, tag, ptr):
    call = it.cmds[tag+"call"]
    valid = it.cmds[tag+"valid"]
    cond="%s->mpi_rank == %s"%(ptr, MASTERID)
    if valid:
      cow.validate(valid)
      return
    if call:
      cow.begin_if(cond)
      cow.addln(call+";"); 
      cow.end_if(cond)
      return
    if not it.decl or it.isdummy: return
    if it.cmds["usr"] not in (None, 1): return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    rank = ptr + "->mpi_rank"
    size = ptr + "->mpi_size"
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    desc = it.mkdesc(1)
    prereq = it.cmds["prereq"]
    etype = it.get_elegtype()
    task = it.cmds[tag]
    comm = ptr+"->mpi_comm"
    if tag == "reduce":
      redtmp = it.cmds["redtmp"]
    if not task: return

    pp = it.mpipp()
    #if pp: print pp; raw_input()
    if pp: cow.addln("#if %s", pp)

    if it.gtype in ("object array", "object pointer"):
      isarr = (it.gtype == "object array")
      if not isarr: cnt = "1"
      cow.add_comment(desc)
      if isarr: 
        cow.declare_var("int i", pp=pp)
        cow.addln('for (i = 0; i < %s; i++) {' % cnt)
      fpfx = it.get_obj_fprefix()
      varnm1 = varname+("+i" if isarr else "")
      cond = "%s%s(%s) != 0" % (fpfx, tag, varnm1)
      cow.die_if(cond, "cannot %s %s" % (tag, varname + (" %d" if isarr else "")), 
         args = "i" if isarr else None,
         onerr = "goto ERR;")
      if isarr: cow.addln('}\n') 
    
    elif it.gtype == "dynamic array":
      cow.add_comment(desc)
      if tag == "reduce":
        if not redtmp:
          print "need temporary array for reduce, use $redtmp:"
          raise Exception
        cow.mpireduce(rank, size, varname, redtmp, cnt, etype, MASTERID, comm)
        cow.begin_if(cond)
        # add synchonized variable on the master
        cow.declare_var("int i")
        cow.addln("for (i = 0; i < %s; i++) %s[i] += %s[i];",
          cnt, task, redtmp)
        cow.end_if(cond)
      elif tag == "bcast":
        cow.mpibcast(rank, size, varname, cnt, etype, MASTERID, comm)
      
      cow.addln()
    else:
      print "cannot %s %s of type %s" % (tag, it.decl.name, it.gtype)
      raise Exception

    if pp: cow.addln("#endif")

  def mkdesc(it, force = 0):
    ''' friendly description of a variable '''
    if not it.decl:
      name = None
    elif it.cmds["flag"]:
      name = it.getflag()
    else:
      name = it.decl.name

    desc = it.cmds["desc"]
    if desc:
      if name:
        desc = "%s: %s" % (name, desc)
    elif force:
      desc = name
    return desc

