#!/usr/bin/env python

import os, subprocess, sys

'''
backup files/folders in tigger
'''

sysname=""

def run_cmd(input, capture=False, old_fashion=False, silent=False):
  if capture:
    pipe=subprocess.PIPE
  else:
    pipe=None

  # detect if the input is a string or not
  if type(input) == type(""):
    cmdstr=input
    cmd=cmdstr.split()
  else:
    cmd=input
    cmdstr=''.join([s+' '  for s in cmd])

  if not silent:
    print sysname+"CMD:",cmdstr

  if old_fashion:
    retcode=os.system(cmdstr)
    oe=["", ""]  # no stdout or stderr
  else:
    p=subprocess.Popen(cmd, stdout=pipe, stderr=pipe)
    oe=p.communicate()
    retcode=p.returncode

  return (retcode, oe[0], oe[1])

def backupdir(input, remove_after=False):
  if type(input)==type([]):
    names=input
  else:
    names=[input]

  for name in names:
    afrom="tigger.rice.edu"
    cmd="rsync -avz "+afrom+":"+name+" ."
    ret=run_cmd(cmd,False,True)
    if 0 != ret[0]:
      print "abort an error occured during rsync"
      return ret[0]

    # we should not remove the file unless rsync is good
    if remove_after:
      print "about to remove remote file/folder "+ name
      ret=run_cmd(["ssh", afrom, "rm", "-r", name])
      if 0 != ret[0]: return ret[0]

  return 0

def help_msg():
  print ("Usage:\n"
      +"  "+myname+" [-D] file_or_dir_to_backup\n\n"
      +"-D to remove the file/dir after the backup\n")
  exit(1)

sysname = run_cmd("hostname", True, False, True)[1].strip()+" "

myname=sys.argv[0]
if len(sys.argv)<=1: help_msg()

filelist=[]
remove_after=False
for f in sys.argv[1:]:
  if f[0] != '-':
    filelist += [f]
  elif f == '-D':
    remove_after=True

if len(filelist)==0: help_msg()

print "About to backup these file/dirs: ", filelist
if remove_after: print "Warning: I will remove these files after the backup!"

# backing up
backupdir(filelist, remove_after)

