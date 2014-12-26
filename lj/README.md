# Lennard Jones system #


## Overview ##

* supports 2d/3d versions, 2d functions are generated
  from the 3d ones by the script add2d.py
* supports square-well and hard-sphere potential
* supports switched potential


## Files ##

  Program             | Description
----------------------|------------------------------------------------------------
  lj.0.h/lj.h         | main source code, latter generate from add2d.py
  ljmc.0.h/ljmc.h     | Monte Carlo algorithms
  ljmd.h              | molecular dynamics algorithms
  ljeos.h             | standard equations of state
  ljrdf.h             | routines to compute RDF
  add2d.py            | generate a version lj.h with 2d support from lj.0.h similarly ljmc.0.h --> ljmc.h

### ljeos ###

  Program             | Description
----------------------|------------------------------------------------------------
  ljeos.gp            | standard Johnson et al. (JZG) equation of states, gnuplot script
  ljeosKN.gp          | Kolafa & Nezbeda PVE-hBH equation of state
  ljeosMBWRKN.gp      | Kolafa & Nezbeda MBWR equation of states, gnuplot script
  KN.f                | Kolafa & Nezbeda PVE-hBH equation of state
  ljeosKN.f           | driver of KN.f, compile this with KN.f
  ljeos.html          | JavaScript equipped webpage for Lennard-Jones equation of states
  dygraph-combined.js | dygraph used in ljeos.html

### Test programs ###

  Program             | Description
----------------------|------------------------------------------------------------
  testmc0.c           | basic Monte Carlo simulation in the NVT ensemble
  testmc.c            | Monte Carlo simulation in the isothermal-isobaric ensemble
  testmd.c            | molecular dynamics with thermostat and barostat
  testctr.c           | prototype of the perturbation temperature
  testref.c           | reference equation of states
  testsw.c            | print out the switch potential
  testeos.c           | test the equation of state


## Notes ##

Monte Carlo pair version: old/ljpr.c

### add2d.py ###

The script uses a generic ccom.py, which is located in the $(ROOT)/python.
But python can find ccom.py without symbolic links because of the
special import statement in add2d.py

```
sys.path.insert(0, "..")
from python import ccom
```

and the existence __init__.py, which treats the entire directory tree
as a python package, such that `python' in the above import statement
means "../python".



### An old speed test ###

testmd.c is somehow slower than ../md/testlj.c
because the array were allocated on the heap instead of the stack

test2.c is a better illustration
Using -DSTACK
         | `icc -O3`    | `gcc -O3 -mfpmath=sse -march=pentium4 -fstrict-aliasing -funroll-all-loops`
---------|--------------|-----------------------------------------------------------------------------
 (HEAP)  | 14.6s        | 15.5s
-DSTACK  | 12.0s        | 12.8s

The speed up of using a stack is > 10%.

Interestingly, swapping xnew(x) and xnew(f) alleviates the problem


### compiling Fortran code ###

```
gfortran ljeos.f KN.f
```
