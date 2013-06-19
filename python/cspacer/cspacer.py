#!/usr/bin/env python

'''
add spaces for C source code
Copyright (c) 2010-2013 Cheng Zhang
the main function is addspacef()

Example 1:
     if(a>b&&a>3)c=d;
-->  if (a > b && a > 3) c = d;

Example 2:
     int foo(int a,int b){
        return (a>b)?a:b;
     }
---> int foo(int a, int b)
     {
        return (a > b) ? a : b;
     }
'''

import os, sys, shutil, getopt, re

# module attributes
use_rule_add      = 0
use_rule_comma    = 1
use_rule_assign   = 1
use_rule_knr      = 1
use_rule_paren0   = 1
use_rule_paren1   = 1
use_rule_paren2   = 0
use_rule_cmp      = 1
use_rule_ter      = 1
use_rule_bitws    = 1
use_rule_else     = 1
use_rule_scolon   = 1
use_rule_spb4quo  = 1
use_rule_nocppcmt = 0

assign  = r"([\+\-\*\/\|\%\&\^]?=)"    # =, +=, -=, *=, ...
assign2 = r"(\>\>|\<\<)="             # >>=  or <<=
cmps    = r"([<>])"
cmp2    = r"([\!<=>]=)"

op1    = r"([\w0-9\'\"\)\]])"
op1b   = r"([\w0-9\)\]])"

op2    = r"([\w0-9\-\'\"\(])"
op2b   = r"([\w0-9\-\(\{])"
op2c   = r"([\*\&\!\+\-][\w0-9\-\(])"

spa    = r"(\s)"
bitws  = r"(\|\||\&\&)"  # || or &&

# Rules
# pattern, replacement, description
rules_basic = [];
rule_else = [
    (r"\Else{",                   "else {",       r"\belse{"),
    (r"}else\b",                   "} else",       r"}else\b"),
    ]

rule_scolon = [
    (";}",                      "; }",          ";}"),
    (r"(;)" + op2b,             r"\1 \2",       ";a"),    # i=0;i<n --> i=0; i<n
    ];

# the additional comma rule
rule_comma = [
    (r"(,)" + op2b,             r"\1 \2",       ",a"),    # a,i     --> a, i
    (r"(,)" + op2c,             r"\1 \2",       ",*a"),   # a,*p    --> a, *p
    ]

rule_assign = [
    (op1b    + assign  + op2,   r"\1 \2 \3",    "a=b"),   # a=b     --> a = b
    (op1b    + assign  + op2c,  r"\1 \2 \3",    "a=*b"),  # a=*b    --> a = *b
    (op1b    + assign2 + op2,   r"\1 \2= \3",   "a>>=b"), # a>>=b   --> a >>= b
    (op1b    + assign2  + op2c, r"\1 \2= \3",   "a>>=*b"),  # a>>=*b  --> a >>= *b
    ]

rule_cmp = [
    (op1     + cmps   + op2c,   r"\1 \2 \3",    "a<*b"),  # a<*b    --> a < *b
    (op1     + cmps   + op2,    r"\1 \2 \3",    "a<b"),   # a<b     --> a < b
    (spa     + cmps   + op2,    r"\1\2 \3",     "a <b"),  # a <b    --> a < b
    (op1     + cmps   + spa,    r"\1 \2\3",     "a> b"),  # a> b    --> a > b
    (op1     + cmp2   + op2c,   r"\1 \2 \3",    "a==*b"), # a==*b   --> a == *b
    (op1     + cmp2   + op2,    r"\1 \2 \3",    "a==b"),  # a==b    --> a == b
    (spa     + cmp2   + op2,    r"\1\2 \3",     "a ==b"), # a <=b   --> a <= b
    (op1     + cmp2   + spa,    r"\1 \2\3",     "a== b"), # a>= b   --> a >= b
    ]

rule_ter = [
    (op1 +"\?"+ op2,            r"\1 ? \2",     "a?b"),   # a?b     --> a ? b
    (op1 +"\:"+ op2,            r"\1 : \2",     "b:c"),   # b:c     --> b : c
    ("\)([\?\:])\(",            r") \1 (",       ")?("),  # (a)?(b) --> (a) ? (b)
    ("\)([\?\:])([\w$])",       r") \1 \2",     ")?a"),   # (a)?b   --> (a) ?b
    ("([\w])([\?\:])\(",        r"\1 \2 (",     "a?("),   # a?(b)   --> a ? (b)
    ]

rule_bitws = [
    (op1b    + bitws  + op2b,   r"\1 \2 \3",    "a||b"),  # a||b    --> a || b
    (spa     + bitws  + op2b,   r"\1 \2 \3",    "a ||b"), # a ||b   --> a || b
    (op1b    + bitws  + spa,    r"\1 \2 \3",    "a|| b"), # a|| b   --> a || b
    ]

rule_paren0 = [
    (r"\bif\(",                 "if (",         r"\bif("),
    (r"\bfor\(",                "for (",        r"\bfor("),
    (r"\bswitch\(",             "switch (",     r"\bswitch("),
    ]

rule_paren1 = [
    (r"\)(\w)",                 r") \1",        ")a"),    # )a      --> ) a
    ("\)\{",                    ") {",          "){"),
    (r"(\)\s*{)(\w)",           r"\1 \2",       "){a"),   # ){a     --> ){ a
    ]

rule_paren2 = [
    (r"\bif\( ",                "if (",         r"\bif( "),   # if( a)  --> if (a)
    (" \)\{",                   ") {",          " ){"),       # a ){    --> a) {
    ]

# an additional rule
# a+b, a-b, a&b, a|b
op1d   = r"([\w\)\]])";
op1e   = r"([a-df-zA-DF-Z0-9_\)\]])";  # exclude e for 1e-5
op2d   = r"([\w\(])";
plsmin = r"([\+\-])";
bandor = r"([\|\&])"
bor =    r"(\|)"

# serveral special cases
#  1e-5
#  -1.0
#  +3.0
#  &abc
rule_add = [
    (op1d    + bandor + op2d,   r"\1 \2 \3",    "a&b"),   # a&b     --> a & b
    (op1e    + plsmin + op2d,   r"\1 \2 \3",    "a+b"),   # a+b     --> a + b
    (spa     + bor    + op2d,   r"\1\2 \3",     "a |b"),  # a |b    --> a | b
    (op1d    + bandor + spa,    r"\1 \2\3",     "a& b"),  # a+ b    --> a + b
    (op1d    + plsmin + spa,    r"\1 \2\3",     "a+ b"),  # a+ b    --> a + b
    ]


comments = [
    ("/*", "*/"),
    ("'''", "'''"), # python __doc__ string
    ('"""', '"""'),
    ("//", ""),
    ]

strings = [
    ('"',  '"'),
    ("'",  "'"),
    ]

# comment or literals
cmtstr = [("BEGINTEXT", "ENDTEXT"), # dummy pair, unused
    ] + comments + strings



def c_cs_start(s):
  ''' find the beginning of the first comment/string starter '''
  mtype = 0 # normal code
  minpos = 1000000
  n = len(cmtstr)
  # loop over comment starters
  for tp in range(1, n):
    (sym0, sym1) = cmtstr[tp]
    pos = s.find(sym0, 0)
    if pos >= 0 and pos < minpos: # if it starts first
      mtype = tp
      minpos = pos
  return (mtype, minpos)



def c_cs_end(s, cstype):
  ''' find the end of a comment/string started as cstype '''
  if cstype <= 0:
    raise Exception

  sym0, sym1 = cmtstr[cstype]
  pos0 = len(sym0) if s.startswith(sym0) else 0

  if (sym0 == '/*' or                  # block comments
      sym0 == "'''" or sym0 == '"""'): # docstrings
    pos1 = s.find(sym1, pos0) # search the end from the line
  elif sym0 == '//':
    pos1 = len(s)   # till the end of the line
  elif sym0 == '"' or sym0 == "'": # tricky
    pos1 = pos0
    while True:
      pos1 = s.find(sym1, pos1)
      if pos1 >= 0 and s[pos1 - 1] == '\\':  # fake, \", keep going
        pos1 += 1 # skip over it
      else:
        break
  else:
    print "unknown sym0: [%s] cstype=%d" % (sym0, cstype)
    raise Exception

  if pos1 >= 0: # return to normal code
    return (0, pos1 + len(sym1))
  else:         # remain in this type
    return (cstype, -1)



def c_parse_line(s, cstype = 0):
  '''
  parse a line like
    foo(a, "haha /* cmt */", /* this is a number */ 0, b); // line end
  to an array of code and comments (or strings):
    `foo(a, '
    `"haha /* cmt */"'
    `, '
    `/* this is a number */'
    ` 0, b); '
    `// line end'
  return a list of (text, comment_type)
  the `comment_type' of regular code is 0

  `cstype' is the current cstype
  '''
  lst = [(s, cstype)] # construct an initial list
  id  = 0

  while id < len(lst):
    s, cstype = lst[id]
    if cstype == 0:  # is normal code
      cstype, pos = c_cs_start(s)
      if cstype <= 0: break # no more cmtstr
      # split it into code + comment
      sp = s[:pos]
      s  = s[pos:]
      lst[id] = (sp, 0) # change the current item to normal code
      lst += [(s, cstype)] # push the current symbol
    else: # current in comment/string of type `cstype', look for the end
      cstypep = cstype
      cstype, pos = c_cs_end(s, cstype)
      if cstype > 0: # current comment is not finished yet
        break  # no need to add anything else
      sp = s[:pos]
      s  = s[pos:]
      lst[id] = (sp, cstypep)
      if len(s) > 0: # the leftover of `// ...'
        lst += [(s, 0)]
    id += 1

  # return a list of parsed items
  return lst, cstype



def addspace(ilines, verbose = 0):
  ''' add spaces to `ilines', an array of lines,
      return the output file '''
  nlines = len(ilines)
  olines = []

  nchanges = 0
  inpp = False
  cstype  = 0  # it is currently code

  # I. install rules
  rules = rules_basic
  if use_rule_paren0:
    rules += rule_paren0
  if use_rule_paren1:
    rules += rule_paren1
  if use_rule_paren2:
    rules += rule_paren2
  if use_rule_else:
    rules += rule_else
  if use_rule_cmp:
    rules += rule_cmp
  if use_rule_ter:
    rules += rule_ter
  if use_rule_assign:
    rules += rule_assign
  if use_rule_bitws:
    rules += rule_bitws
  if use_rule_comma:
    rules += rule_comma
  if use_rule_scolon:
    rules += rule_scolon
  if use_rule_add:
    rules += rule_add
  if verbose >= 4:
    i = 1
    for r in rules:
      print "%3d: '%s'" % (i, r[2])
      i += 1
    raw_input()

  # II. add spaces line by line
  for i in range(nlines):
    iline = ilines[i]
    oline = iline.rstrip()

    # apply various rules
    changed = []

    # check if we are inside a preprocessor line
    if cstype == 0 and (
        inpp or iline.lstrip().startswith("#")):  # preprocessor
      if verbose >= 5:
        print "line %6d: %s\npreprocessor" % (i+1, iline.rstrip())
        raw_input()
      ilin = iline.rstrip()
      if len(ilin) > 0 and ilin[-1] == '\\': # line of preprocessor
        inpp = True   # continue the preprocessor line
      else:
        inpp = False  # start code again
    else:
      # split a line to code and comment and strings
      lst, cstype = c_parse_line(oline, cstype)
      nlast = len(lst) - 1
      item, cstype1 = lst[nlast]
      sym0, sym1 = cmtstr[cstype1]

      if verbose >= 5:  # show the parser
        print "line %6d: %s" % (i+1, iline.rstrip())
        print "cstype %6d: %s" % (cstype, lst)
        raw_input()

      ''' function starting brace to K & R style
         int foo(){

      -->int foo()
         {
      '''
      if (use_rule_knr and 0 == cstype and nlast == 0
          and not item[0:1].isspace()):
        pat = r"(\)\s*\{$)";
        repl = ")\n{"
        if re.search(pat, item):
          item = re.sub(pat, repl, item)
          lst[-1] = (item, 0)
          changed += ["K&R function"]
          if verbose > 0:
            print "converting to K & R function"

      # convert `// comments', to `/* comments */'
      if (use_rule_nocppcmt and sym0 == "//"):
        item = "/* " + item[2:].strip() + " */"
        lst[-1] = (item, lst[-1][1])
        if verbose > 0:
          print "convert C++ comment"

      if use_rule_spb4quo:
        for k in range(1, len(lst)):
          cstp = lst[k][1]
          if (cmtstr[cstp][0] in "\'\""
              and len(lst[k-1][0]) > 0
              and lst[k-1][0][-1] in ",;"):
            lst[k-1] = (lst[k-1][0] + " ", lst[k-1][1])
            changed += ['[space before "]']

      # apply rules
      oline = ""  # "%d: " % cstype
      for k in range(len(lst)):
        (item, cstp) = lst[k]
        if cstp == 0: # only apply to code
          for pat,repl,desc in rules:
            # we apply each rule multiple times
            # until it no longer applies
            # Note some pattern may emerge after replacement
            while re.search(pat, item):
              item = re.sub(pat, repl, item)
              changed += ['[' + desc + ']']
        oline += item

      # change cstype for the next line
      sym0, sym1 = cmtstr[ lst[-1][1] ]
      if sym0 in ("'", "//", '"'):
        cstype = 0  # return to code mode

    # print out the change
    if oline.rstrip() != iline.rstrip():
      nchanges += 1
      if verbose >= 2:
        print "%6d INPUT : %s" % (i+1, iline),
        print "%6d OUTPUT: %s" % (i+1, oline.rstrip())
        print "Rules:", ',  '.join(changed)
      if verbose >= 3:
        raw_input()

    # strip away trailing spaces before outputting
    olines += [ oline + '\n', ]

  # re split lines
  olines = ''.join(olines).splitlines(True)
  olines = [s.rstrip() + '\n' for s in olines]

  return olines, nchanges



def addspacef(fninp, fnout = "", overwrite = False, verbose = 0):
  """ add space to file 'fninp' """

  try:
    ilines = open(fninp, 'r').readlines();
  except IOError:
    print "cannot open", fninp
    return

  if verbose >= 2: print "processing", fninp
  olines, nchanges = addspace(ilines, verbose)

  # print the number of changes
  if not nchanges:
    if verbose:
      print "keep", fninp
    return
  else:
    nlines = len(ilines)
    print ("about to make %d changes (%.2f%%) to %s"
        % (nchanges, 100.0*nchanges/(nlines + 1e-6), fninp) )
    if (nchanges >= 20 or nchanges >= 0.1*nlines) and verbose >= 1:
      print "the code looks nasty, did you write it?"
      raw_input()

  if overwrite:  # overwrite mode
    try: # backup
      import zcom
      zcom.safebackup(fninp, ".orig")
      fnout = fninp  # write on the original input
    except ImportError: pass
  else:  #
    if not fnout:
      fnout = fninp + ".spr"
      if verbose > 0:
        print "assume the output is", fnout
  open(fnout, 'w').writelines(olines)



def usage():
  ''' print usage and die '''

  print sys.argv[0], "[OPTIONS] input"
  print """
  add spaces to C source code

  OPTIONS:

   -R, --recursive        recursively apply to subdirectories
                          if `input' is a wildcard pattern like *.c
                          the pattern must be quoted as '*.c'
   -L, --nolinks          skip symbolic links
   -w, --overwrite        overwrite the original file
   -a, --add              add space around +, -, &, |
   --paren2               convert if(_ to if_(, and _){ to )_{
   --noknr                allow { to hang after ) for functions
   --noparen0             don't convert if( to if (
   --noparen1             don't convert ){ to ) {
   --noelse               don't convert }else{ to } else {
   --noassign             don't add space around =
   --nocmp                don't add spaces around == or <
   --nobitws              don't add space around ||
   --noter                don't add spaces around ? :
   --nocomma              don't add space after ,
   --noscolon             don't add space after ;
   --nospb4quo            allow no space before the leading quote
   -c, --conservative     be conservative
   --cppcmt               convert C++ style comments // to /* */
   -v                     be verbose
   --verbose=[0-9]        specify verbose level
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options
      results saved to module attributes '''

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvbwacRL",
         ["help", "verbose=", "backup", "overwrite",
          "add", "conservative", "nocomma", "noscolon",
          "noassign", "nocmp", "noter",
          "noparen0", "noparen1", "paren2",
          "noelse", "noknr", "nobitws",
          "nospb4quo",
          "cppcmt",
          "--recursive", "--nolinks",
         ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global use_rule_nocppcmt
  global use_rule_paren0
  global use_rule_paren1
  global use_rule_paren2
  global use_rule_add
  global use_rule_cmp
  global use_rule_bitws
  global use_rule_ter
  global use_rule_assign
  global use_rule_comma
  global use_rule_scolon
  global use_rule_else
  global use_rule_knr
  global use_rule_spb4quo

  overwrite = False
  verbose = 0
  recur = False
  links = True

  for o, a in opts:
    if o in ("-R", "--recursive",):
      recur = True
    elif o in ("-L", "--nolinks",):
      links = False
    elif o in ("-b", "-w", "--backup", "--overwrite"):
      overwrite = True
    elif o in ("--cppcmt",):
      use_rule_nocppcmt = 1
      print "will convert C++ comments to C ones"
    elif o in ("-a", "--add",):
      use_rule_add = 1
      print "enable the additional add rule"
    elif o in ("--paren2",):
      use_rule_paren2 = 1
      print "enable the additional parentheses rule"
    elif o in ("-c", "--conservative",):
      use_rule_add    = 0
      use_rule_paren2 = 0
      #
      use_rule_paren0 = 0
      use_rule_comma  = 0
      use_rule_bitws  = 0
      use_rule_ter    = 0
      use_rule_cmp    = 0
      use_rule_assign = 0
      use_rule_knr    = 0
      use_rule_spb4quo = 0
      '''
      use_rule_scolon = 0
      use_rule_paren1 = 0
      use_rule_else   = 0
      '''
      print "use basic rules only"
    elif o in ("--noparen0",):
      use_rule_paren0 = 0
      print "disable the parentheses start rule"
    elif o in ("--noparen1",):
      use_rule_paren1 = 0
      print "disable the parentheses end rule"
    elif o in ("--noelse",):
      use_rule_else = 0
      print "disable the else rule"
    elif o in ("--nocmp",):
      use_rule_cmp = 0
      print "disable the comparison rule"
    elif o in ("--noter",):
      use_rule_ter = 0
      print "disable the ternary rule"
    elif o in ("--nobitws",):
      use_rule_bitws = 0
      print "disable the bitws rule"
    elif o in ("--noassign",):
      use_rule_assign = 0
      print "disable the assignment rule"
    elif o in ("--nocomma",):
      use_rule_comma = 0
      print "disable the comma rule"
    elif o in ("--noscolon",):
      use_rule_scolon = 0
      print "disable the semicolon rule"
    elif o in ("--nospb4quo",):
      print "allow space before quotes"
    elif o in ("--noknr",):
      use_rule_knr = 0
      print "disable the K&R rule"
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  ls = args
  try: # limit the dependence on argsglob
    import zcom
    ls = zcom.argsglob(args, "*.c *.cpp *.h *.hpp *.java", recur = recur, links = links)
  except ImportError: pass
  return ls, overwrite, verbose



if __name__ == "__main__":
  ls, overwrite, verbose = doargs()
  for fn in ls:
    addspacef(fn, None, overwrite, verbose)

