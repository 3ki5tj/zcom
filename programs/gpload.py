#!/usr/bin/env python
''' load Gnuplot script(s) from command line '''

import Gnuplot, os, sys, glob

def main():
  # determine the script to plot
  if len(sys.argv) > 1:
    scripts = sys.argv[1:]
  else:
    scripts = glob.glob("*.gp")
    if len(scripts) < 1:
      usage()

  g = Gnuplot.Gnuplot(persist = 1)
  for fn in scripts:
    print "loading script %s..." % fn
    #g('load "%s"' % fn)
    g.load(fn)
    viewfile(fn)
    g.reset()


def usage():
  ''' print usage and die '''
  print "%s your.gp" % sys.argv[0]

def viewfile(fn):
  ''' search the output and call gnome-open to open it '''
  for s in open(fn).read().split('\n'):
    s = s.strip()
    if s.startswith("set output"):
      i = s.find('"')
      if i < 0: continue
      file = s[i+1:-1]
      print "output is %s" % file
      #os.system("gnome-open %s" % file)


if __name__ == "__main__":
  main()
