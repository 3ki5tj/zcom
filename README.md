Overview
---------

* This project collects commonly-used routines into a single header file zcom.h.

* Each individual module has a separate directory, except
  zcom.h, which is a combination of all frequently-used ones.

* Code should be concise and small, module names should be short.

* Use symbolic links to define inter-modular dependence,
Use `_' to define the preceding module.


Files
------

File          | Description
--------------|----------------------------------
zcom.0.h:     | template that creates zcom.h
zcom.h        | generated from zcom.0.h, combined header, by assemble.py
assemble.py   | the python script that assembles modules to a single header ./assemble.py -a
Makefile      | Makefile for common tasks, such as `usb', `usball', etc
gen.mak       | generic Makefile to be used for new modules
license.txt   | GPL less license, version 2
zcompick.py   | end-user tool to remove certain modules from zcom.h



Directories
-----------

Directory     | Description
--------------|------------------------------------
python        | helper python scripts, such as add2d.py, rmdbg.py
programs      | although the project is mainly a library, several useful programs are collected `programs' directory
cspacer.py    | code beautifiers: cspacer.py, rstrip.py


Modules
-------

Module        | Description
--------------|------------------------------------
util:         | useful but unclassified routines
bits:         | bit operations
blkmem:       | block memory
ss:           | smart string
endn:         | handle endianness
bio:          | macros for read/write binary files
osys:         | operating system related
rng:          | random number generator
rc:           | real complex number
rv2:          | real 2-vectors and 2x2 matrices
rv3:          | real 3-vectors and 3x3 matrices
rvn:          | real D-vectors and DxD matrices
eig:          | eigenvalues and eigenvectors
cholesky:     | Cholesky decomposition
lu:           | LU decomposition
svd:          | single value decomposition
savgol:       | compute Savitzky-Golay coefficients
specfunc:     | special functions
opt:          | options
argopt:       | command line arguments and options
cfg:          | configuration file
log:          | log file
av:           | simple mean/deviation
hist:         | simple histograms (1d/2d)
mds:          | multidimensional scaling
pdb:          | read/write pdb
md:           | common molecular dynamics routines
hmc:          | hybrid Monte Carlo/molecular dynamics
ising2:       | two dimensional ising model
potts2:       | two dimensional potts model
lj:           | Lennard Jones system
abpro:        | AB beads protein model
cago:         | C-alpha Go model
glez:         | common uses for GLUT

modules with a leading underscore are experimental, and not added to zcom.h


Notes:
-------
* Try to apply cspacer.py/reindent.py under cspacer to check the source here from time to time

* zcom.h should not have a routine like zcom_finish()
  because in case of multiple inclusion, zcom_finish is defined in the first instance,
  it has no knowledge of modules included later.



TOOLS:
-------
* Makefile
  common tasks
  make: call assemble.py to assemble zcom.h from modules
  make zcom.zip: make a package

* assemble.py -a
  assemble zcom.h from modules

* zcompick.py
  client tool to extract a few modules from zcom.h to a smaller version

