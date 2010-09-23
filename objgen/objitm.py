#!/usr/bin/env python
import os, sys, re 
from copy    import * 
from objpos  import P 
from objdcl  import CDeclaratorList
from objcmt  import CComment
from objpre  import CPreprocessor
from objccw  import CCodeWriter
from objcmd  import Commands

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

  def __str__(self):
    return "pre: %5s, dl: %5s, decl: %5s, cmt: %5s, raw = [%s]" % (
        self.pre != None, self.dl != None, 
        self.decl != None, self.cmt != None, 
        self.begin.gettext(self.src, self.end).strip() )

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

    if self.decl.name == "000":  # debug code
      print "gtype [%s] with %s, %s" % (ret, offset, types); 
      print "commands: %s" % self.cmds
      raw_input()        
    return ret

  def type2fmt(self, tp):
    fmt = ""
    if   tp == "int":      
      fmt = "%d"
    elif tp == "long":
      fmt = "%ld"
    elif tp == "unsigned": 
      fmt = "%u"
    elif tp == "unsigned long": 
      fmt = "%ul"
    elif tp == "float":
      fmt = "%f"
    elif tp == "double":
      fmt = "%lf"
    elif tp == "char *":
      fmt = "%s"
    else:
      print "no format string for type [%s] item: %s" % (tp, self)
      raise Exception
    return fmt

  def fill_def(self):
    ''' set 'def' value '''
    if not self.decl or "def" in self.cmds: return
    defval = None
    gtype = self.gtype
    # for arrays, 'def' means the default value of its elements
    if gtype in ("static array", "dynamic array"):
      gtype = self.get_gtype(offset = 1)
      # print "gtype of %s is %s" % (self.decl.name, gtype); raw_input()
    
    if gtype in ("int", "unsigned", "unsigned int", 
        "long", "unsigned long"):
      defval = "0"
    elif gtype in ("float", "real"):
      defval = "0.0f"
    elif gtype in ("double"):
      defval = "0.0"
    elif gtype in ("pointer", "function pointer", "char *", 
        "pointer to object"):
      defval = "NULL"
    elif gtype == "MPI_Comm":
      defval = "MPI_COMM_NULL"
    else: 
      defval = "(%s) 0" % self.decl.datatype
    if defval != None: 
      #print "set default %s for gtype %s, var %s" % (defval, gtype, self.decl.name)
      #raw_input()
      self.cmds["def"] = defval

  def add_item_to_list(self, dcl): 
    ''' add an item with decl being dcl '''
    it = deepcopy(self) # make a copy
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

    if "key" in it.cmds:
      cfgkey = it.cmds["key"]
    else:
      cfgkey = it.decl.name
      pm = it.cmds["key_unprefix"]
      if pm and cfgkey.startswith(pm):
        cfgkey = cfgkey[len(pm):]
      pfx = it.cmds["key_prefix"]
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
    cnt = '*'.join(dim)
    it.cmds["cnt"] = cnt
    it.cmds["dim"] = dim

  def fill_io(it):
    ''' expand the io string '''
    if not it.decl: return
    # default I/O mode, bt for array, cbt for others
    if it.gtype in ("static array", "dynamic array"):
      iodef = "bt"
    elif it.gtype == "char *":
      iodef = "c"
    elif it.gtype in simple_types:
      iodef = "cbt"
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
    it.cmds["io_bin"] = 1 if "b" in io else 0
    it.cmds["io_txt"] = 1 if "t" in io else 0

  def fill_test(it):
    ''' 
    fill test conditions 
    note: we assign commands even if there's no declaration
    becomes a stand-alone comment can contain assignment
    or other commands
    '''
    # turn on $test_first for dummy variables
    if ("test_first" not in it.cmds and it.isdummy):
      it.cmds["test_first"] = 1

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
    ''' get additional arguments to be passed to an object initializer '''
    if not it.cmds["obj"]:
      print "item is not an object %s" % it
      raise Exception
    args = it.cmds["init_args"]
    if not args: return ""
    arr = args.strip().split(",")
    arr = [""] + [s.strip() for s in arr]
    args = ", ".join(arr)
    # print args; raw_input()
    return args



