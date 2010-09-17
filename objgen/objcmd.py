#!/usr/bin/env python
import os, sys, re 
from copy import copy

class Commands:
  '''
  commands imbed in the comments
  '''
  def __init__(self, s):
    # do symbol substitution first
    self.raw = self.subst_symbols(s)
    self.cmds = {}
    self.persist = {}
    self.parse_commands()

  # make the object looks like a dictionary
  def __getitem__(self, key):
    return self.cmds[key]
  def __setitem__(self, key, value):
    self.cmds[key] = value
  def __iter__(self):
    return self.cmds.__iter__()
  def __contains__(self, key):
    return self.cmds.has_key(key)
  def __len__(self):
    return len(self.cmds)
  def __str__(self):
    return str(self.cmds)
  def __repr__(self):
    return repr(self.cmds)

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
    s = '\n'.join(sa)
    while 1:
      '''
      command line:
        $cmd op args;
      op can be one of ":", "=", ":=", "::" or "" (nothing)
        : or = means set the current variable only
        := or :: means a persistent command that apply to 
                 all following items
                 if args is empty, the command is removed
      '''
      pattern = r"[^\$\\]*(\$)([\w\-]+)\s*([\+\-\:]?=|:?:)\s*(.*?)((?<!\\)\;|$)"
      m = re.search(pattern, s, re.MULTILINE | re.DOTALL)
      if m:
        '''
        group 2 is the cmd
        group 3 is the operator
        group 4 is the argument
        '''
        cmd  = m.group(2)
        args = m.group(4)
        if m.group(3) in (":=", "::"):
          if args.strip() != "$0":
            #print "set an persistent command, raw=[%s] args=[%s]" % (s, args); raw_input()
            self.persist[cmd] = 1
          else:  # $cmd := $0; to unset a persistent command
            #print "turning off an persistent command, raw=[%s]" % s; raw_input()
            self.persist[cmd] = -1
      else:
        if s.find("$") < 0: break
        # print "look for a lazy command, $cmd with no ; s = [%s]" % s
        pattern = r"[^\$\\]*(\$)(\w+)\;?"
        m = re.search(pattern, s)
        if m: 
          # print "found a lazy command, $cmd with no ; s = %s" % s
          cmd = m.group(2)
          args = "on"  # it means turn on a switch
        else:
          # search for comment command
          pattern = r"[^\$\\]*(\$#)(.*?)([^\\]\;|$)"
          m = re.search(pattern, s)
          if m == None: 
            print "possibly an unknown command [%s]" % s
            break
          cmd = "#"
          args = m.group(2)
          #print "a comment is found [%s]" % m.group(0)
     
      cmd = re.sub(r"\-", "_", cmd)  # map - to _
      self.cmds[cmd] = args # add to dictionary
      s = s[:m.start(1)] + s[m.end(0):] # remove the command from comment
      #print "a pattern is found [%s]: [%s], rest: %s" % (
      #    cmd, param, s)
      #raw_input()

    # join multiple line description
    if not "desc" in self.cmds: 
      sa = [a.strip() for a in s.splitlines()] # split to lines
      s = ' '.join(sa).strip()
      if len(s): self.cmds["desc"] = s 

  def on(self, cmd):
    if (cmd not in self.cmds
        or self.cmds[cmd] not in 
        ("on", "ON", "yes", "YES", "y", "Y", "TRUE", "true")):
      return 0
    else: return 1

  def subst_symbols(self, s, do_at = 0):
    ''' 
    $$ or \$  =>  $ 
    \;        =>  ;
    '''
    # merge the remaining string to description
    s = re.sub(r"[\$\\]\$", "$", s) # literal $
    s = re.sub(r"\\;", ";", s) # literal ;
    return s
