#!/usr/bin/env python
import os, sys, re, filecmp
from copy import copy

'''
intend to be a code generator
TODO:
  * multiple-line support with regular expression
  * preprocessor
  * interpretor
  * skip struct in comment or string
'''


cmt0 = "/*"
cmt1 = "*/"
lncmt = "//"

class P:
  '''
  simple position in the code (line, column)
  unfortunately, the most vulnerable class
  '''
  def __init__(self, row = 0, col = 0):
    self.row = row
    self.col = col
    self.check() # see if input are integers
  def __str__(self):
    self.check()
    return "line: %d, col %d" % (self.row + 1, self.col + 1)
  def __repr__(self):
    self.check();
    return "%d, %d" % (self.row, self.col)
  def __cmp__(a, b):
    a.check()
    if type(a) != type(b):
      print "trying to compare %s with %s" % (type(a), type(b))
      raise Exception
    b.check()
    return (a.row - b.row) * 10000 + a.col - b.col
  def check(self):
    if type(self.row) != int:
      print "row %s is not a number" % (self.row)
      raise Exception
    if type(self.col) != int:
      print "col %s is not a number" % (self.col)
      raise Exception
  def getline(p, src):
    ''' get a line starting from the current position '''
    return src[p.row][p.col:] if p.row < len(src) else None
  def nextline(self, src = None):
    '''
    go to the next line
    we don't really need src, just make it look like other methods
    '''
    self.check()
    self.row += 1
    self.col = 0

  def skipspace(p, src):
    '''
    skip spaces (including multiple blank lines)
    return the line
    '''
    p.check()
    while 1:
      s = p.getline(src)
      if s == None:
        self.ttype = self.token = None
        break
      # skip leading space
      m = re.match(r"(\s*).*", s)
      if not m: raise Exception
      p.col += m.end(1)
      s = p.getline(src)
      if len(s) != 0: break # found a token
      # this line is exhausted
      p.nextline()
      continue
    return s

  def gettoken(p, src):
    '''
    aggressively get a token starting from the current `p'
    `p' is updated to the end of the current token
    return the token only
    '''
    s = p.skipspace(src)
    
    # try to get a word
    if s != None:
      m = re.match(r"([a-zA-Z_]\w*).*", s)
      if m:
        p.col += m.end(1)
        p.ttype, p.token = "word", s[:m.end(1)]
      else:
        p.col += 1
        p.ttype = p.token = s[0]
    return p.token, p.ttype

  def ungettok(p, src):
    '''
    try to unget a token, it will not unget space
    '''
    p.check()
    tok = p.token
    if not tok:
      print "cannot unget None"
      raise Exception
    p.col -= len(tok)
    if (p.col < 0 or
        not p.getline(src).startswith(tok)):
      raise Exception
    p.token = p.ttype  = None

  def peektok(p, src):
    '''
    peek at the next token
    '''
    p.check()
    q = copy(p) # create a backup
    return q.gettoken(src)

  def gettext(self, src, q):
    '''
    get entire string from p to q
    '''
    p = copy(self) # so self is intact
    if q <= P(0,0):
      q = P(len(src), 0) # end of string
    s = ""
    while p.row < q.row:
      s += p.getline(src)
      p.nextline()
    if p.row == q.row and p.col < q.col:
      s += src[p.row][p.col: q.col]
    return s

class Preprocessor:
  '''
  preprocessor, taken literally, not understood
  '''
  def __init__(self, src, p):
    self.empty = 1
    self.parse(src, p)

  def parse(self, src, p):
    s = p.skipspace(src)
    self.begin = copy(p)

    if not s.startswith("#"): return -1
    self.raw = s[1:].strip()
    self.empty = 0
    p.col += len(s)

  def isempty(self): return self.empty

class Declarator:
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
          token[-2:] not in ("_t", "_T")):
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
          token = s[m.start(1) : m.end(1)]
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

class DeclaratorList(Declarator): 
  '''
  C declarations of multiple variables
  we inherit `dclspec()' and `datatype' from Declarator
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
      # we find declarators (variables) through Declarator's
      # can s
      d = Declarator(src, p)
      if d.isempty(): return -1
      d.datatype = self.datatype # add datatype for individual Declarator
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

class Comment:
  '''
  C block-comment or one-line comment
  '''
  def __init__(self, src, p):
    self.empty = 1
    if 0 != self.span(src, p):       # .begin, .end, and .raw
      return
    self.get_commands()     # .cmds with .cmds['desc'] being the description
    self.empty = 0

  def span(self, src, p):
    '''
    find the beginning and end of a comment,
    starting from pos
    '''
    s = p.skipspace(src) # note, its skips multiple blank lines
    if s == None: return -1

    if s.startswith(cmt0) or s.startswith(lncmt):
      self.begin = copy(p)
    else:
      #print "missing comment beginning mark %s" % p
      return -1

    p.col += 2 # skip `/*' or `//'

    if s.startswith(lncmt): # line comment
      self.raw = s[2:].strip()
      p.nextline()
      self.end = copy(p)
    else: # block comment
      # aggressively search the end of comment
      while 1:
        s = p.skipspace(src)
        if s.startswith(cmt1):
          p.col += len(cmt1)
          self.end = copy(p)
          #print "comment found from %s to %s" % (self.begin, self.end)
          break
        p.gettoken(src) # get token by token
      else:
        print "Warning: cannot find the end of comment from %s, now = %s" % (
            self.begin, now)
        raise Exception
      self.raw = self.extract_text(src)
    return 0

  def extract_text(self, src):
    '''
    extract text within block comment
    '''
    self.checkbe()
    s = self.begin.gettext(src, self.end).strip()
    if s.startswith(cmt0) and s.endswith(cmt1):
      s = s[len(cmt0) : len(s) - len(cmt1)]
    #print "comment text is [%s] from %s to %s" % (s, self.begin, self.end)
    return s

  def checkbe(self):
    self.begin.check()
    self.end.check()

  def get_commands(self):
    '''
    parse raw text in comment to commands, save results to .cmds
    '''
    self.cmds = {}
    s = self.raw
    pos = 0
    
    ''' 
    for multiple-line comments
    remove the leading`* ' for subsequent lines
    '''
    sa = [a.strip() for a in s.splitlines()]
    if len(sa) > 1: # sa[:1] is the first line
      sa = sa[:1] + [(a[2:].rstrip() if a.startswith("* ") else a) for a in sa[1:] ] 
    s = ' '.join(sa)
    while 1:
      # note: $$ or \$ means a literal $
      # command line:
      #   cmd op args;
      # op can be one of ":", "=", ":=", "::" or "" (nothing)
      #   : or = means set the current variable only
      #   := or :: means also set the parser's current state
      pattern = r"[^\$\\]*(\$)(\w+)\s*(\:|\=|\:\=|\:\:|)\s*(.*?)\;"
      m = re.search(pattern, s, re.MULTILINE | re.DOTALL)
      if m == None: break
      # a command is found
      # group 2 is the cmd
      # group 3 is the operator
      # group 4 is the argument
      cmd = s[m.start(2) : m.end(2)]
      persist = 1 if m.end(3) - m.start(3) == 2 else 0 # the operator
      param = s[m.start(4) : m.end(4)]
      self.cmds[cmd] = [param, persist] # add to dictionary
      s = s[:m.start(1)] + s[m.end(0):]; # remove the command from comment
      #print "a pattern is found [%s]: [%s], rest: %s" % (
      #    cmd, param, s)
      #raw_input()

    # merge the remaining string to description
    if not self.cmds.has_key("desc"):
      s = re.sub(r"[\$\\]\$", "$", s) # literal $
      sa = [a.strip() for a in s.splitlines()] # split to lines
      s = ' '.join(sa).strip()
      self.cmds["desc"] = [s, 0] # description is not persistent
      #print "the remain string is [%s]" % self.cmds["desc"]

  def isempty(self): return self.empty
  def __str__(self): return self.raw
  def __repr__(self): return __str__(self)


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
    self.begin = copy(p)
    if src == None: return
    self.src = src
    self.parse(src, p)

  def parse(self, src, p):
    '''
    parse an item aggressively
    '''
    s = p.skipspace(src)
    if s[0] == "}": return -1
 
    # try to see it's a preprocessor
    pp = Preprocessor(src, p)
    if not pp.isempty():
      #print "pos %s, preprocessor %s" % (p, pp.raw)
      #raw_input()
      self.empty = 0
      self.pp = pp
      self.dl = self.cmt = None
    else:
      self.pp = None

      # try to get a declaration
      self.decl = None
      dl = DeclaratorList(src, p)
      self.dl = None if dl.isempty() else dl
      
      # try to get a comment
      cmt = Comment(src, p)
      self.cmt = cmt if not cmt.isempty() else None
    
    if self.pp or self.cmt or self.dl:
      #print "successfully got one item, at %s" % p
      self.end = copy(p)
      self.empty = 0
    else:
      print "bad line: at %s, s = [%s]" % (p, s.strip())
      raise Exception
    
    self.expand_multiple_declarators()

  def __str__(self):
    return "pp: %5s, dl: %5s, cmt: %5s, raw = [%s]" % (
        self.pp != None, self.dl != None, self.cmt != None, 
        self.begin.gettext(self.src, self.end).strip() )
    
  def get_generic_type(self):
    '''
    get a generic type, and type name
    '''
    types = self.decl.types
    cmds = self.cmt.cmds if self.cmt else []
    nlevels = len(types)
    if nlevels == 1:
      return types[0]
    elif types[0].startswith("array"):
      self.arrcnt = types[0][5:]
      return "static array"
    elif types[0] == "pointer":
      #towhat = ' '.join(types[1:])
      if nlevels == 2 and types[1] == "function":
        return "function pointer"
      elif nlevels == 2 and types[1] == "char":
        return "string"
      elif 'arr' in cmds:
        return "dynamic array"
      elif 'obj' in cmds:
        return "object"
      else:
        return "pointer"

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
    if self.dl == None:
      self.decl = None
      self.itlist = [self]
      return
    dclist = self.dl.dclist
    n = len(dclist)
    if n <= 0:
      self.decl = None
      self.itlist = [self]
      return

    self.itlist = []
    for i in range(n):
      decl = dclist[i]
      decl.types += [self.dl.datatype]
      it = copy(self)  # make a shallow copy of myself
      it.decl = decl
      it.gtype = it.get_generic_type()
      self.itlist += [it]
    
  def isempty(self):
    return self.empty

class Object:
  '''
  a struct defined as
  typedef struct {
  } abc_t;
  '''
  def __init__(self, src, pos = None):
    ''' initialize myself from a bunch of src '''
    if len(src) == 0: return # return an empty object

    if type(pos) == int:  # integer input
      pos = P(pos, 0)
    self.items = []
    self.empty = 1
    self.parse(src, pos)

  def parse(self, src, p):
    ''' parse src starting from 'p' '''
   
    # print "start searching object from %s" % p
    if not self.find_beginning(src, p):
      return

    while 1: # parse the content within { ... },
      # aggressively search for an item, and update p
      item = Item(src, p)
      if item.isempty(): break

      #print "line %3d has %d item(s): %s" % (p.row + 1, len(item.itlist), item)
      self.items += item.itlist
   
    # parse }
    self.find_ending(src, p, 0)

    self.determine_prefix();

  def find_beginning(self, src, p, aggr = 1):
    '''
    search the block-beginning marker {
    aggressive search: repeat until found
    position 'p' is move to the end of file
      if nothing is found
    limited to single-line case
    '''
    while 1:
      line = p.getline(src)
      if line == None: break

      pattern = r"\s*(typedef)\s+struct\s*\w*\s*\{"
      m = re.search(pattern, line)
      if m: # a match is found
        self.begin = P(p.row, p.col + m.start(1))
        p.col += m.end(0)
        #print "object-beginning is found in %s, %s" % (p, src[p.row])
        return 1
      if aggr:
        p.nextline()
      else:
        break
    return 0  # not found
  
  def find_ending(self, src, p, aggr = 0):
    '''
    search the block-ending mark
    passive search: since we assume aggressive search has been
    performed for items, we only search the current line
    position 'p' is unchanged, if nothing is found
    '''
    while 1:
      line = p.getline(src)
      if line.strip() == "": # blank line 
        p.nextline()
        continue
      pattern = r'\s*(\}\s*)(\w+)(\s*;).*'
      m = re.match(pattern, line)
      if m:
        p.col += m.end(3)
        self.name = line[m.start(2) : m.end(2)]
        self.end = copy(p)
        self.empty = 0
        #print "object-ending is found in %s, %s" % (p, src[p.row])
        return 1
      if aggr:
        p.nextline()
      else:
        break
    print "no object-ending is found, p at", p
    return 0  # not found

  def isempty(self): return self.empty

  def write_decl(self):
    '''
    return a string for formatted declaration
    '''
    cw = CodeWriter()
    cw.addln("typedef struct {")
    block = [] # block list
    for i in range(len(self.items)):
      item = self.items[i]
      decl0 = decl1 = scmt = ""
      # 1. handle pp
      if item.pp:
        self.dump_block(block, cw)
        block = [] # empty the block
        cw.addln("#" + item.pp.raw)
        continue
      # 2. handle declaration
      if item.decl:
        decl0 = item.decl.datatype
        decl1 = item.decl.raw
        #raw_input ("%s %s" % (decl0, decl1))
      # 3. handle comment
      cmt = item.cmt
      if (cmt and cmt.cmds.has_key("desc")
          and len(cmt.cmds["desc"][0]) > 0):
        scmt = "/* " + cmt.cmds["desc"][0] + " */"
      # 4. output the content
      if len(decl0) > 0:  # put it to a code block
        block += [(decl0, decl1, scmt)] # just buffer it
      else:
        self.dump_block(block, cw)
        block = [] # empty the block
        # also print this line
        if len(scmt) > 0: 
          cw.addln(scmt)
    self.dump_block(block, cw)
    block = []
    cw.add("} " + self.name + ";")
    return cw.gets()
  def dump_block(self, block, cw):
    if len(block) == 0: return
    # align and format all items in the bufferred block
    # and append the result to string s
    wid0 = max(len(item[0]) for item in block)
    wid1 = max(len(item[1]) for item in block) + 1 # for ;
    for item in block:
      line = "%-*s %-*s %s\n" % (wid0, item[0], wid1, item[1]+';', item[2])
      cw.add(line)

  def determine_prefix(self):
    name = self.name
    if name.endswith("_t"): name = name[:-2]
    # guess a short nickname
    if len(name) > 3: name = name[:2].lower()
    self.nickname = name
    self.prefix = self.nickname + "_"

  def write_programs(self):
    s = ""
    s += self.write_prog_close()
    return s

  def write_prog_close(self):
    ow = CodeWriter()
    funcname = "%sclose" % self.prefix
    ow.addln("/* %s: close a pointer to %s */", funcname, self.name)
    ow.addln("#define %s(%s) { %s_(%s); free(%s); %s = NULL; }",
      funcname, self.nickname, funcname, self.nickname, self.nickname, self.nickname)
    ow.addln("void %s_(%s *%s)", funcname, self.name, self.nickname)
    ow.addln("{")
    
    # free strings 
    s_len = lambda it: len(it.decl.name) if (it.decl and it.gtype == "string") else 0
    s_width = max(s_len(it) for it in self.items)
    for it in self.items:
      if it.decl and it.gtype == "string":
        ow.addln("if (%s->%-*s != NULL) ssdelete(%s->%s);",
                 self.nickname, s_width, it.decl.name,
                 self.nickname, it.decl.name)

    # free pointers
    p_len = lambda it: len(it.decl.name) if (it.decl and it.gtype == "dynamic array") else 0
    p_width = max(p_len(it) for it in self.items)
    for it in self.items:
      if it.decl and it.gtype == "dynamic array":
        ow.addln("if (%s->%-*s != NULL) free(%s->%s);", 
                 self.nickname, p_width, it.decl.name,
                 self.nickname, it.decl.name)

    ow.addln("memset(%s, 0, sizeof(*%s));", self.nickname, self.nickname)
    ow.addln("}")
    ow.addln("")
    return ow.gets()

class CodeWriter:
  '''
  output code in an indented fashion
  '''
  sindent = "  "
  def __init__(self, nindents = 0):
    self.nindents = nindents
    self.s = ""

  def inc(self): self.nindents += 1
  def dec(self): 
    self.nindents = self.nindents-1 if self.nindents > 0 else 0

  def add(self, t, *args):
    t = t % args
    if t.lstrip().startswith("}"): self.dec()
    if self.s.endswith("\n") and not t.startswith("#"):
      self.s += self.sindent * self.nindents
    self.s += t
    if t.rstrip().endswith("{"): self.inc()

  def addln(self, t, *args):
    self.add( (t.strip() % args) + '\n' )

  def gets(self): return self.s
    

class Parser:
  '''
  handle objects in a file
  '''
  def __init__(self, src):
    self.src = src
    self.parse_lines()

  def parse_lines(self):
    '''
    process the input template 
    '''
    objs = self.objs = []
    p = P()

    # search for objects
    while 1:
      obj = Object(self.src, p) # try to get an object
      if obj.isempty(): break

      print "found a new object '%s' with %d items, from %s to %s\n" % (
          obj.name, len(obj.items), obj.begin, obj.end)
      objs += [obj]
    print "I got %d objects" % (len(objs))

  def output(self):
    s = self.change_objs_decl()
    
    # append programs
    for obj in self.objs:
      s += obj.write_programs()
    return s

  def change_objs_decl(self):
    '''
    return the completed C source code
    '''
    objs = self.objs
    if objs == []: return ""
    src = self.src
    p0 = P(0,0)
    s = p0.gettext(src, objs[0].begin) # to the beginning of the first object
    nobjs = len(objs)
    for i in range(nobjs):
      obj  = objs[i]
      s += obj.write_decl()
      pnext = objs[i+1].begin if i < nobjs - 1 else P(len(src), 0)
      s += obj.end.gettext(src, pnext) # this object's end to the next object's beginning
    return s
   

def handle(file, template = ""):
  ''' generate a C file according to a template '''
  # guess the template name if it's empty
  if template == "":
    pt = os.path.splitext(file)
    template = pt[0] + '.0' + pt[1]

  # read in the template
  src = open(template, 'r').readlines()
  # construct a parser
  psr = Parser(src)
  # write the processed file
  open(file, 'w').write(psr.output())

def main():
  handle("test2.c")

if __name__ == "__main__":
  main()

