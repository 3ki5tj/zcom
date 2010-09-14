#!/usr/bin/env python
import os, sys, re 
from copy import copy

'''
C comment
'''

cmt0 = "/*"
cmt1 = "*/"
lncmt = "//"

class CComment:
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



