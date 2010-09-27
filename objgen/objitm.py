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
  "char *",
  "dummy_t")

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
 
    # try to see it's a preprocessor
    pre = CPreprocessor(src, p)
    if not pre.isempty():
      #print "pos %s, preprocessor %s" % (p, pre.raw)
      #raw_input()
      self.empty = 0
      self.pre = pre
      self.dl = self.cmt = None
    else:
      self.pre = None

      # print "try to get a declaration from %s" % p; raw_input()
      self.decl = None
      dl = CDeclaratorList(src, p)
      self.dl = None if dl.isempty() else dl
      
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

  def get_zero(self):
    if not self.decl:
      raise Exception
    gtype = self.get_elegtype()
    zval = None
    if gtype in ("int", "unsigned", "unsigned int", 
        "long", "unsigned long"):
      zval = "0"
    elif gtype in ("float", "real"):
      zval = "0.0f"
    elif gtype in ("double"):
      zval = "0.0"
    elif gtype in ("pointer", "function pointer", "char *", 
        "pointer to object"):
      zval = "NULL"
    elif gtype == "MPI_Comm":
      zval = "MPI_COMM_NULL"    
    # if not zval: zval = "(%s) 0" % self.decl.datatype
    return zval

  def fill_def(self):
    ''' set 'def' value '''
    if not self.decl or "def" in self.cmds: return
    self.cmds["def"] = self.get_zero()
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

  def subpfx(it, val, ptrname, fprefix, parent):
    ''' change `@' by ptrname-> in command arguments '''
    # @- means function prefix
    pattern = r"(?<![\@\\])\@\-(?=\w)"
    val = re.sub(pattern, fprefix, val)
    
    # @^ means function prefix
    pattern = r"(?<![\@\\])\@\^(?=\w)"
    val = re.sub(pattern, "parent", val)

    # @@ means this variable 
    pattern = r"(?<![\@\\])\@\@(?!\w)"
    if re.search(pattern, val):
      if not it.decl:
        print "use @@ for %s without decl." % it
        raise_exception
      val = re.sub(pattern, ptrname + "->" + it.decl.name, val)

    # @~ means parent's corresponding name 
    pattern = r"(?<![\@\\])\@\~(?!\w)"
    if re.search(pattern, val):
      if not it.decl:
        print "use @~ for %s without decl." % it
        raise_exception
      val = re.sub(pattern, "%s->" % parent + it.decl.name, val)

    # @< means keyname
    kprefix = it.cmds["key_prefix"]
    if kprefix:
      pattern = r"(?<![\@\\])\@\<(?=\w)"
      val = re.sub(pattern, kprefix, val)

    # @var
    pattern = r"(?<![\@\\])\@(?=\w)" # exclude @@, \@
    val = re.sub(pattern, ptrname + "->", val)
    
    # for a single hanging @, means ptrname
    pattern = r"(?<![\@\\])\@(?!\w)"
    if re.search(pattern, val):
      val = re.sub(pattern, ptrname, val)

    return val

  def sub_prefix(it, ptrname, fprefix, parent):
    ''' wrapper of subpfx '''
    for key in it.cmds:
      val = it.cmds[key]
      if type(val) == str:
        it.cmds[key] = it.subpfx(val, ptrname, fprefix, parent)
      elif type(val) == list:
        for i in range(len(val)):
          val[i] = it.subpfx(val[i], ptrname, fprefix, parent)
      else: continue


  def fill_key(it):
    ''' 
    add default key for configuration file reading
    unless "key" is specified, it is calculated as
    1. it.decl.name is the first guess
    2. if it starts with key_unprefix, remove it
    3. add key_prefix
    4. append key_args, after a `#'
    '''

    # skip if it is not a declaration
    #if not it.decl or "key" in it.cmds: return
    if not it.decl: return

    key0 = it.cmds["key"]
    cfgkey = key0 if key0 else it.decl.name
    
    if not key0 or it.cmds["flag"]:
      pm = it.cmds["key_unprefix"]
      if pm and cfgkey.startswith(pm):
        cfgkey = cfgkey[len(pm):]
      pfx  = it.cmds["key_prefix"]
      if pfx:  cfgkey = pfx + cfgkey
    # append key argument
    args = it.cmds["key_args"]
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
    io = it.cmds["io"] if it.cmds["io"] else iodef
    io = io.strip()

    # a human shortcut
    if io in ("no", "none"): io = ""
    elif io in ("all", "everything"): io = "cbt"

    # create itemized io
    it.cmds["io_cfg"] = 1 if "c" in io else 0
    it.cmds["io_bin"] = 1 if ("b" in io or it.cmds["usr"]=="bintmp") else 0
    it.cmds["io_txt"] = 1 if "t" in io else 0

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
    if "flag" not in it.cmds:
      return

    # split into two parts
    flag = it.cmds["flag"].strip()
    pt = flag.split()
    if len(pt) < 2:
      print "flag is invalid! flag %s" % flag
      raise Exception
    flag = pt[0].strip()
    it.cmds["flag_val"] = pt[1].strip()

    # attach flag_prefix
    if "flag_prefix" in it.cmds:
      prefix = it.cmds["flag_prefix"]
      if not flag.startswith(prefix):
        flag = prefix + flag
    it.cmds["flag"] = flag

    # assign the default flag_var
    if "flag_var" not in it.cmds:
      it.cmds["flag_var"] = "@flags"
    #print "%s %s" % (it.cmds["flag"], it.cmds["flag_val"])
    #raw_input()
    
  def get_flag(it):
    if "flag" not in it.cmds: return None
    return it.cmds["flag"]

  def get_obj_fprefix(it):
    ''' get prefix of an object '''
    if not it.cmds["obj"]:
      print "item is not an object %s" % it
      raise Exception
    pfx = it.cmds["obj_fprefix"]
    if pfx == None: # make a guess
      tp = it.decl.datatype
      if tp.endswith("_t"): tp = tp[:-2]
      pfx = tp + "_"
    return pfx

  def get_init_args(it):
    return it.get_args("init")

  def get_args(it, tag):
    ''' get additional arguments to be passed to an object initializer '''
    if not it.cmds["obj"]:
      print "item is not an object %s" % it
      raise Exception
    args = it.cmds[("%s_args" % tag)]
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
    if it.pre:
      cow.dump_block(block)
      block = [] # empty the block
      cow.addln("#" + it.pre.raw)
      return block
    
    # 2. handle declaration
    if it.cmds["usr"] not in (None, 1, "cfg"):
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
    if it.pre:
      cow.addln("#" + it.pre.raw)
      return
  
    usr = it.cmds["usr"]
    if usr == "parent" or usr in ("tmp", "bintmp", "txttmp"):
      return
    
    key     = it.cmds["key"]
    defl    = it.cmds["def"]
    prereq  = it.cmds["prereq"]
    tfirst  = it.cmds["tfirst"]
    valid   = it.cmds["valid"]
    must    = it.cmds["must"]
    cnt     = it.cmds["cnt"]
    desc    = it.cmds["desc"]
    pp      = it.cmds["#if"]
    flag    = it.cmds["flag"]

    if not it.decl: # stand-alone comment
      cow.add_comment(desc)
      cow.insist(it.cmds["assert"], desc)
      call = it.cmds["call"]
      if call: cow.addln(call + ";")
      return
    if flag and (not key or key == "flags"): return
    
    varnm = it.decl.name
    varname = ptrname + "->" + varnm

    # add a remark before we start
    cow.add_comment(desc)

    if flag: # input flag
      #raw_input("flag with key %s" % key)
      flag = it.get_flag()
      fmt = type2fmt_s(it.gtype)
      cow.cfgget_flag(varname, key, flag, it.gtype, fmt,
          defl, prereq, desc)
  
    elif it.gtype == "object pointer":
      fpfx = it.get_obj_fprefix()
      s = "(%s = %scfgopen(cfg%s)) == NULL" % (varname, 
        fpfx, it.get_init_args())
      cow.die_if(s, r"failed to initialize %s\n" % varname)

    elif it.gtype == "object array":
      fpfx = it.get_obj_fprefix()
      etp = it.get_gtype(offset = 1)
      cow.alloc_darr(varname, etp, cnt)
      cow.declare_var("int i;", pp = pp)
      cow.addln("for (i = 0; i < %s; i++) {" % cnt)
      s = "0 != %scfgopen_low(%s+i, cfg%s)" % (
        fpfx, varname, it.get_init_args())
      cow.die_if(s, r"failed to initialize %s[%%d]\n" % varname, "i")
      cow.addln("}\n")

    elif it.gtype == "dynamic array":
      cow.init_darr(varname, it.get_gtype(offset = 1),
          defl, cnt, desc, pp)
    
    elif it.gtype == "static array":
      etp = it.get_gtype(offset = 1)
      # print "%s: static array of element = [%s] io = %s, defl = [%s]" % (varnm, etp, it.cmds["io_cfg"], defl); raw_input()
      if not it.cmds["io_cfg"]: 
        cow.init_sarr(varname, defl, cnt, pp)
      else:
        fmt  = type2fmt_s(etp)
        cmpl = it.cmds["complete"]
        cow.cfgget_sarr(varname, key, etp, fmt, defl, cnt, 
            must, valid, cmpl, desc)

    else: # regular variable
      usrval = it.cmds["usr"]
      if usrval:
        if usrval in ("cfg", 1): # usr variables
          cow.assign(varname, "usr_%s" % varnm, it.gtype)
        else: pass # ignore other usr variables
      elif not it.cmds["io_cfg"]: # assign default value
        if defl and len(defl): 
          cow.assign(varname, defl, it.gtype);
      else:
        fmt = type2fmt_s(it.gtype)
        cow.cfgget_var(varname, key, it.gtype, fmt, 
          defl, must, prereq, tfirst, valid, desc)


  def binrw_var(it, cow, varname, rw):
    if rw == "w": 
      it.binwrite_var(cow, varname)
    elif rw == "r":
      it.binread_var(cow, varname)
    else:
      raise Exception
   
  def binread_var(it, cow, varname):
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    bincnt = it.cmds["bin_cnt"]
    cond = it.cmds["bin_prereq"]
    defl = it.cmds["def"]
    usrval = it.cmds["usr"]
    verify = it.cmds["verify"]
    valid = it.cmds["valid"]
    pp = it.cmds["#if"]
    if pp:
      cow.addln("#if %s", pp)
    if notalways(cond):
      cow.begin_if(cond)

    if not it.decl:
      print "no declare gtype = %s" % it.gtype
      raw_input()

    if usrval == "bintmp": 
      cow.declare_var(it.decl.datatype + " " + it.decl.raw, it.decl.name)
    
    if it.gtype == "dynamic array":
      # support nasty flag $bin_cnt
      if len(dim) == 1 and bincnt:
        dim[0] = bincnt
      cow.rwb_arr(varname, dim, it.decl.datatype, 1, "r")
    elif it.gtype in ("object array", "object pointer"):
      fpfx = it.get_obj_fprefix();
      funcall = "%sbinread_low(%%s, fp, ver, flags, endn%s)" % (
          fpfx, it.get_args("binwrite"));
      
      if it.gtype == "object array":
        imin = it.cmds["bin_imin"]
        imax = it.cmds["bin_imax"]
        cow.rb_objarr(varname, dim, funcall, it.cmds["#if"], 
            [1, 1], imin, imax)
      else:
        cow.rb_obj(varname, funcall)

    elif it.decl.datatype == "char" and (
        it.gtype in ("char *", "static array")):
      if it.gtype == "static array" and bincnt: 
        cnt = bincnt
      cond1 = "%s != fread(%s, 1, %s, fp)" % (cnt, varname, cnt)
      cow.begin_if(cond1)
      cow.addln(r'fprintf(stderr, "cannot read string of %%d for %s\n", %s);', 
          varname, cnt)
      cow.addln("goto ERR;")
      cow.end_if(cond1)
    else:
      cow.rb_var(varname, it.gtype, verify, valid = valid)

    if notalways(cond):
      cow.end_if(cond)
    if pp:
      cow.addln("#endif")
   
  def binwrite_var(it, cow, varname):
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    bincnt = it.cmds["bin_cnt"]
    cond = it.cmds["bin_prereq"]
    defl = it.cmds["def"]
    usrval = it.cmds["usr"]
    pp = it.cmds["#if"]
    if pp:
      cow.addln("#if %s", pp)
    if notalways(cond):
      cow.begin_if(cond)

    if usrval == "bintmp": 
      cow.declare_var(it.decl.datatype + " " + it.decl.raw, it.decl.name)
      if defl:
        if it.gtype == "static array":
          cow.init_sarr(varname, defl, cnt)
        else:
          cow.addln("%s = %s;", varname, defl)
    
    if it.gtype == "dynamic array":
      # support nasty flag $bin_cnt
      if len(dim) == 1 and bincnt:
        dim[0] = bincnt
      cow.rwb_arr(varname, dim, it.decl.datatype, 1, "w")
    
    elif it.gtype in ("object array", "object pointer"):
      fpfx = it.get_obj_fprefix();
      funcall = "%sbinwrite_low(%%s, fp, ver%s)" % (
          fpfx, it.get_args("binwrite"))
      if it.gtype == "object array":
        imin = it.cmds["bin_imin"]
        imax = it.cmds["bin_imax"]
        cow.wb_objarr(varname, dim, funcall, pp,
            [1, 1], imin, imax)
      else:
        cow.wb_obj(varname, funcall)
    
    elif it.decl.datatype == "char" and (
        it.gtype in ("char *", "static array")):
      if it.gtype == "static array" and bincnt: 
        cnt = bincnt
      cond1 = "%s != fwrite(%s, 1, %s, fp)" % (cnt, varname, cnt)
      cow.begin_if(cond1)
      cow.addln(r'fprintf(stderr, "cannot write string of %%d for %s\n", %s);', 
          varname, cnt)
      cow.addln("goto ERR;")
      cow.end_if(cond1)
    
    else:
      cow.wb_var(varname, it.gtype)

    if notalways(cond):
      cow.end_if(cond)
    if pp:
      cow.addln("#endif")

  def clear_var(it, cow, ptr):
    if not it.decl or it.isdummy: return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    prereq = it.cmds["com_prereq"]
    defl = it.get_zero()

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
    if it.cmds["usr"] == "parent": return
    if it.gtype not in ("char *", "dynamic array",
        "object array", "object pointer"): return

    varnm = it.decl.name;
    varname = "%s->%s" % (ptrname, varnm)

    pp = it.cmds["#if"]
    prereq = it.cmds["com_prereq"]
    if pp: cow.addln("#if %s", pp)
    if notalways(prereq): cow.begin_if(prereq)    
 
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
    if not it.decl or it.isdummy: return
    if it.cmds["usr"] not in (None, 1): return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    prereq = it.cmds["com_prereq"]

    pp = it.cmds["#if"]
    if pp: cow.addln("#if %s", pp)
    if notalways(prereq): cow.begin_if(prereq)

    if it.gtype in ("dynamic array", "static array"):
      parr = 1
      try:
        fmt = type2fmt_p(it.get_elegtype(), varname)
        cow.addln(r'fprintf(fp, "%s: %s of %s\n");',
            varname, it.gtype, cnt)
      except TypeError:
        cow.addln(r'fprintf(fp, "%s: %s of %s %%p\n", %s);',
          varname, it.gtype, cnt, varname)
        parr = 0

      if cnt == "0": parr = 0
      if parr:
        cow.declare_var("int i", pp = pp) # declare index i
        cow.addln("for (i = 0; i < %s; i++) {", cnt)
        cow.addln('fprintf(fp, "%s, ", %s);', fmt, varname + "[i]")
        cow.addln(r'if ((i+1) %% 10 == 0) printf("\n");')
        cow.addln("}")
        
        cow.addln(r'if ((%s) %% 10 != 0) fprintf(fp, "\n");', cnt)
    elif it.gtype in ("object array", "object pointer"):
      fpfx = it.get_obj_fprefix();
      cond = "%s != NULL" % varname
      cow.begin_if(cond)      
      funcall = "%smanifest(%%s, fp)" % fpfx      
      if it.gtype == "object array":
        cow.addln(r'fprintf(fp, "%s: %s of %s %s\n");',
          varname, it.gtype, cnt, it.decl.datatype)
        cow.declare_var("int i", pp = pp) # declare index i
        cow.addln("for (i = 0; i < %s; i++)", cnt)
        cow.addln(cow.sindent + funcall % (varname+"+i") + ";")
      else:
        cow.addln(r'fprintf(fp, "%s: %s to %s\n");',
          varname, it.gtype, it.decl.datatype)
        cow.addln(funcall % varname + ";")
      cow.end_if(cond)
    elif it.gtype in ("pointer"):
      pass
    else:
      try:
        fmt = type2fmt_p(it.gtype, varname)
        cow.addln(r'fprintf(fp, "%s: %s, %s\n", %s);',
          varname, it.gtype, fmt, varname)
      except TypeError:
        cow.addln(r'fprintf(fp, "%s: %s, 0x%%X\n", (unsigned) %s);',
          varname, it.gtype, varname)

    if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")


  def initmpi_var(it, cow, ptr):
    if not it.decl or it.isdummy: return
    if it.cmds["usr"] not in (None, 1): return

    varnm = it.decl.name
    varname = ptr + "->" + varnm
    dim = it.cmds["dim"]
    cnt = it.cmds["cnt"]
    desc = it.cmds["desc"]
    mpicnt = it.cmds["mpi_cnt"]
    if mpicnt: cnt = mpicnt;
    prereq = it.cmds["com_prereq"]
    mpi = it.cmds["mpi"]
    etype = it.get_elegtype()

    pp = it.cmds["#if"]
    if pp: cow.addln("#if %s", pp)
    #if notalways(prereq): cow.begin_if(prereq)

    notmaster = "%s->rank != %s" % (ptr, MASTERID)
    if mpi in ("alloc", "2"): # allocate memory only, temporary variables
      if it.gtype not in ("dynamic array"):
        print "cannot just allocate space for %s (%s)" % (it.decl.name, it.gtype)
        raise Exception
      cow.add_comment(desc)
      cow.begin_if(notmaster)
      cow.alloc_darr(varname, etype, cnt)
      cow.end_if(notmaster)
    
    elif mpi in (1, "1", "bcast"):  # allocate and bcast
      if it.gtype == "dynamic array":
        cow.add_comment(desc)
        cow.begin_if(notmaster)
        cow.alloc_darr(varname, etype, cnt)
        cow.end_if(notmaster)
        cow.addln("MPI_Bcast(%s, (%s) * sizeof(%s[0]), MPI_BTYE, %s, comm);",
            varname, cnt, varname, MASTERID)
        cow.addln()
      elif it.gtype in ("object array", "object pointer"):
        isarr = (it.gtype == "object array")
        if not isarr: cnt = "1"
        cow.add_comment(desc)
        cow.begin_if(notmaster)
        cow.alloc_darr(varname, it.decl.datatype, cnt)
        cow.end_if(notmaster)
        if isarr: 
          cow.declare_var("int i")
          cow.addln('for (i = 0; i < %s; i++) {' % cnt)
        fpfx = it.get_obj_fprefix()
        cow.addln("%sinitmpi(%s, comm);", fpfx, varname+("+i" if isarr else ""))
        if isarr: cow.addln('}\n')
      elif it.gtype == "char *":
        cow.declare_var("int i")
        cow.add_comment(desc)
        cow.addln("i = strlen(%s);", varname)
        cow.addln("MPI_Bcast(&i, sizeof(int), MPI_BTYE, %s, comm);", MASTERID)
        cow.begin_if(notmaster)
        cow.addln("%s = ssnew(i);" % varname)
        cow.end_if(notmaster)
        cow.addln("MPI_Bcast(%s, i*sizeof(char), MPI_BTYE, %s, comm);", 
            varname, MASTERID)
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
        cow.addln('for (i = 0; i < %s; i++) {\n%s[i] = NULL;\n}', cnt, varname)
        cow.end_if(notmaster)
      # do nothing otherwise
    
    else:
      print "unknown mpi tag %s" % mpi
      raise Exception

    #if notalways(prereq): cow.end_if(prereq)
    if pp: cow.addln("#endif")

