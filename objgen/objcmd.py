#!/usr/bin/env python
import os, sys, re 
from copy import copy

class Commands:
  '''
  commands imbed in the comments
  '''
  def __init__(self, s):
    self.raw = s
    self.cmds = {}
    self.persist = {}
    self.parse_commands()

  # make the object looks like a dictionary
  def __getitem__(self, key):
    return self.cmds[key]
  def __contains__(self, item):
    return self.cmds.has_key(item)

  def parse_commands(self):
    '''
    parse raw text in comment to commands, save results to .cmds
    cmds['desc'] being the description
    the process
    '''
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
      '''
      command line:
        $cmd op args;
      op can be one of ":", "=", ":=", "::" or "" (nothing)
        : or = means set the current variable only
        := or :: means a persistent command that apply to 
                 all following items
                 if args is empty, the command is removed
                 
      
      note: $$ or \$ means a literal $
      '''
      pattern = r"[^\$\\]*(\$)(\w+)\s*(:?=|:?:)\s*(.*?)\;"
      m = re.search(pattern, s, re.MULTILINE | re.DOTALL)
      if m:
        '''
        group 2 is the cmd
        group 3 is the operator
        group 4 is the argument
        '''
        cmd  = s[m.start(2) : m.end(2)]
        args = s[m.start(4) : m.end(4)]
        if m.end(3) - m.start(3) >= 2: # the operator
          if len(args.strip()):
            # print "found an persistent command, raw=[%s]" % s; raw_input()
            self.persist[cmd] = 1
          else:  # $cmd := ; to unset a persistent command
            self.persist[cmd] = -1
      else:
        if s.find("$") < 0: break
        # print "look for a lazy command, $cmd with no ; s = [%s]" % s
        pattern = r"[^\$\\]*(\$)(\w+)\;?"
        m = re.search(pattern, s)
        if m == None: break
        # print "found a lazy command, $cmd with no ; s = %s" % s
        cmd = s[m.start(2) : m.end(2)]
        args = ""
      self.cmds[cmd] = args # add to dictionary
      s = s[:m.start(1)] + s[m.end(0):] # remove the command from comment
      #print "a pattern is found [%s]: [%s], rest: %s" % (
      #    cmd, param, s)
      #raw_input()

    # merge the remaining string to description
    if not "desc" in self.cmds:
      s = re.sub(r"[\$\\]\$", "$", s) # literal $
      sa = [a.strip() for a in s.splitlines()] # split to lines
      s = ' '.join(sa).strip()
      if len(s): self.cmds["desc"] = s 
      #print "the remain string is [%s]" % self.cmds["desc"]

