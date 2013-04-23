#!/usr/bin/env python

import os, sys, re, shutil, getopt

verbose = 0
delsrc = 0 # delete remote data, e.g., sugar biou
checkit = 0 # check integrity after updating

def localcmd(cmd):
  ''' run cmd, capture output '''
  out = "ACCIO.OUT"
  os.system(cmd + " > " + out)
  output = open(out).read()
  os.remove(out)
  return output

def remotecmd(addr, pw, cmd, capture = 1):
  ''' run cmd at addr, return output '''
  out = "ACCIO.OUT"
  #print "addr: %s; password: %s; cmd: %s" % (addr, pw, cmd)
  if pw:
    sshcmd = "/home/czhang/app/bin/pwssh %s %s %s" % (
        pw, addr, cmd)
  else: # no password
    sshcmd = "ssh %s %s" % (addr, cmd)
  if capture:
    sshcmd += " > " + out
  print "SSH: %s" % sshcmd
  os.system(sshcmd)
  if capture:
    output = open(out).read()
    os.remove(out)
    if pw:
      arr = output.split('\n')
      if arr[0].strip().endswith("assword:"):
        output = '\n'.join(arr[1:])
    return output
  else:
    return None


def remotecopy(pw, scpargs, verbose = 0):
  ''' like scp '''
  if pw:
    scpargs = scpargs.strip()
    if scpargs.startswith("-r "):
      scpargs = scpargs[3:]
    cmd = "/home/czhang/app/bin/pwscp %s %s" % (pw, scpargs)
  else: # no password
    cmd = "scp %s" % scpargs
  if verbose: print "CMD: %s" % cmd
  os.system(cmd)

def getremote(server):
  ''' dc12_sugar --> dc12@sugar.rice.edu '''
  if server.find("_") >= 0:
    ar = server.split('_', 1)
    usr = ar[0]
    server = ar[1]
    if usr == "dc12":
      pw = "v3pIrqe"
    elif usr == "yg7":
      pw = "f703n1ks"
    elif usr == "cz5":
      pw == "w3r3wlf"
    else:
      print "don't know password of %s" % usr
      raise Exception
  else:
    usr = "cz1"
    pw = None
  return usr, server + ".rice.edu", pw



def findlocalprj(root, code, trjid):
  out = localcmd("ls --color=none -d %s/[a-z]*/%s/trj/%s" % (root, code, trjid)).strip()
  arr = out.split('\n')
  if len(arr) != 1:
    print "cannot find a root for %s %s %s" % (root, code, trjid)
    print "arr = %s" % str(arr)
    raise Exception
  arr = arr[0].split("/")
  if len(arr) >= 4:
    return arr[-4]
  else:
    return None

def update(target, hostname, delmode):
  ''' handle target, e.g. ~/disk2/biou/a3m2/m2c '''

  arr = target.rsplit("/", 3)
  if len(arr) < 3:
    print "invalid local dir [%s]" % target
    raise Exception
  server = arr[-3] # "biou" or "dc12_sugar"
  path = arr[-2]+"/"+arr[-1]
  tiggersrc = "/mnt/disk2/czhang"
  if hostname == "tigger":
    usr, addr, pw = getremote(server) # translate to address
    remotepath = "/projects/jpma/%s/%s"  % (usr, path)
    localdir = "%s/%s/%s" % (tiggersrc, server, path)
  else:
    usr = "czhang"
    addr = "tigger.rice.edu"
    pw = None
    remotepath = "%s/%s/%s" % (tiggersrc, server, path)
    localroot = "/media/Extra/work"
    localprj = findlocalprj(localroot, arr[-2], arr[-1])  # e.g., a3d, a3w
    if localprj:
      localdir = "%s/%s/%s/trj/%s" % (localroot, localprj, arr[-2], arr[-1])
    else: # for pegda
      localdir = "%s/%s/%s" % (localroot, arr[-2], arr[-1])


  if not os.path.exists(localdir):
    print "local dir %s doesn't exist" % localdir
    raise Exception

  # list remote content
  lscmd = "ls --color=none %s/data*/md*.xtc" % remotepath
  remotedirs = []
  out = remotecmd(usr+"@"+addr, pw, lscmd).strip()
  if out.find("No such file") >= 0:
    list = []
  else:
    list = out.split()
  for s in list:
    pt = s.rsplit("/", 2)
    if len(pt) != 3:
      print "invalid entry [%s|%s] during [%s]" % (s, out, localdir)
      raise Exception
    datadir = pt[1] # "data33"
    if pt[2].find("part") >= 0: # gromacs 4.0
      id = datadir[4:]
      if pt[2].find(id) < 0:
        print "corruption detected %s" % s
        raise Exception
    remotedirs += [datadir, ]
  print 'remote: ' + ' '.join(remotedirs)

  # list local content
  lscmd = "ls --color=none -d %s/data*/md*.xtc" % localdir
  localdirs = []
  for s in localcmd(lscmd).split():
    pt = s.rsplit("/", 2)
    datadir = pt[1]
    localdirs += [datadir, ]
  print 'local: ' + ' '.join(localdirs)

  if delmode:  # do folder deletion
    updated = 0
    for dir in remotedirs:
      if dir not in localdirs: # skip
        print "local %s has no correspondence" % dir
        raise Exception
      cmd = "rm %s/%s/*" % (remotepath, dir)
      remotecmd(usr+"@"+addr, pw, cmd, capture = 0)
      updated = 1
    print "complete deletion for %s" % remotepath
  else:
    updated = 0
    for dir in remotedirs:
      if dir in localdirs: # skip
        continue
      # copy from the remote dir
      src = "%s@%s:%s/%s" % (usr, addr, remotepath, dir)
      dest = localdir + "/"
      args = "-r %s %s" % (src, dest)
      remotecopy(pw, args, verbose = 1)
      updated = 1

  return localdir, updated

def handling(localdir):
  ''' do rmsbb stuff '''
  if localdir.find("pegda") >= 0: return
  os.chdir(localdir)
  print "handling %s ..." % localdir
  if os.path.exists("data1/ZE"):
    # list local content
    lscmd = "ls --color=none data*/ZE"
    localdirs = []
    max = -1
    for s in localcmd(lscmd).split():
      pt = s.rsplit("/")
      try:
        id = int(pt[0][4:])
      except Exception:
        print "pt = %s" % pt
        exit(1)

      if id > max: max = id
    os.system("ln -sf data%s/ZE ." % max)
  os.system("../rmsbb")

def check_integrity(dir):
  ''' check trajectory is continuous '''
  os.chdir(dir)

  # list all data folders
  lscmd = "ls --color=none -d data*"
  max = -1
  for s in localcmd(lscmd).split():
    try:
      id = int(s[4:])
    except Exception:
      print "s = %s" % s
      exit(1)
    if id > max: max = id

  # if all data folders exists
  cntze = 0
  for i in range(1, max+1):
    datadir = "data%s" % i
    if not os.path.exists(datadir):
      print "corruption: %s/%s doesn't exist" % (dir, datadir)
      raise Exception
    fnze = datadir + "/ZE"
    if os.path.exists(fnze):
      cntze += 1
  if cntze != 0 and cntze != max:
    print "corruption: ZE cnt %d is not %s" % (cntze, max)
    raise Exception

  # check trajectory continuity
  for i in range(1, max):
    ip = i+1
    print "%s check continuity from data%s to data%s" % (localdir, i, ip)
    def dcmp(a, b):
      cmd = "cmp data%s/%s data%s/%s" % (i, a, ip, b)
      if 0 != os.system(cmd):
        print "corruption in comparison %s\n[%s]" % (dir, cmd)
        raise Exception

    dcmp("state.cpt", "state0.cpt")
    if cntze:
      dcmp("mb.av", "mb.av0")
      dcmp("hist.bin", "hist0.bin")
      dcmp("MTSEED", "MTSEED0")


def update_tigger_script():
  ''' update the tigger_grab.py '''
  me = os.path.abspath(__file__)
  cmd = "rsync -avz %s czhang@tigger.rice.edu:app/bin/" % me
  print "updating %s on tigger" % me
  os.system(cmd)

def usage():
  ''' print usage and die '''
  print sys.argv[0], "[Options]"
  print '''lazy switch: %s -d -c''' % sys.argv[0]
  print '''Options:
  -d: delete remote src
  -c: check integrity
  -C: force checking integrity
  -h: help
  '''
  exit(1)

def doargs(hostname):
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hdcCv",
         ["help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose, delsrc, checkit

  for o, a in opts:
    if o == "-v":
      verbose = 1
    elif o in ("-d",):
      delsrc = 1
    elif o in ("-c",):
      checkit = 1
    elif o in ("-C",):
      checkit = 2
    elif o in ("-h", "--help"):
      usage()

  # force checking integrity
  if delsrc and checkit < 1: checkit = 1

  if len(args) > 0:
    return args
  else:
    if hostname == "tigger":
      dir = localcmd("pwd -P .").strip()
      return [dir,]
    else: # do everything
      return [
        'biou/a3m2/m2a', 'biou/a3m2/m2b', 'biou/a3m2/m2c',
        'biou/leuq/lq1', 'biou/leuq/lq2', 'biou/leuq/lq3',
        'biou/lach/ch1', 'biou/lach/ch2', 'biou/lach/ch3',
        'biou/lala/aa1', 'biou/lala/aa2',
        'biou/lalq/aq1', 'biou/lalq/aq2',
        'dc12_sugar/a3wt/wt2', 'dc12_sugar/a3wt/wt4',
        'biou/a3wt/wt3',
        ]

if __name__ == "__main__":
  hostname = localcmd("hostname").strip()
  print "hostname: [%s]" % hostname
  if not hostname in ("prime", "tigger"):
    print "invalid hostname [%s]" % hostname
    raise Exception

  if hostname == "prime":
    update_tigger_script()

  list = doargs(hostname)
  for t in list:
    if hostname == "tigger":  # tigger version
      localdir, updated = update(t, hostname, delsrc)
      if delsrc: continue
      if checkit >= 2 or (checkit >= 1 and updated):
        check_integrity(localdir)
    else:  # local version
      # call tigger_grab.py
      cmd0 = "ssh czhang@tigger.rice.edu /home/czhang/app/bin/tigger_grab.py"
      cmd = cmd0 + " " + t
      print "Call remote: %s" % cmd
      if 0 != os.system(cmd): raise Exception

      localdir, updated = update(t, hostname, 0)
      if updated:
        handling(localdir)

      if checkit >= 2 or (checkit >= 1 and updated):
        check_integrity(localdir)

      if delsrc: # tell tigger to remove remote
        cmd = cmd0 + " -d " + t
        if 0 != os.system(cmd): raise Exception
