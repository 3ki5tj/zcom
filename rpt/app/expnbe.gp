#!/usr/bin/env gnuplot
reset
set terminal pngcairo enhanced font 'times, 24' linewidth 2 round size 1024, 768
set output "expnbe.png"

# Line style for axes
set style line 80 lt rgb "#404040"
set border 3 back linestyle 80  # Remove border on top and right.

# Line style for grid
set style line 81 lt 0 lw 1.5 # dotted
set style line 81 lt rgb "#a0a0a0" # light gray
set grid mxtics xtics mytics ytics back linestyle 81

set tics

set style line 1 lt rgb "#b02010" lw 3 pt 6 ps 3  # dark red
set style line 2 lt rgb "#e0c040" lw 3 pt 12 ps 3  # yellow to orange
set style line 3 lt rgb "#FF2020" lw 3 pt 6 ps 2.5 # bright red
set style line 4 lt rgb "#60b040" lw 3 pt 8 ps 3 # dark green
set style line 5 lt rgb "#305080" lw 3 pt 4 ps 3 # navy blue
set style line 9 lt rgb "#000000" lw 3 pt 12 ps 2.5      #

set xtics .2
set mxtics 2
set xlabel "{/Symbol b}" offset 0, 0.5

set ytics .2
set mytics 2
set ylabel "{/Symbol \341}exp({/Symbol - b e}){/Symbol \361}" offset 1.0, 0
# try {/Symbol-Oblique b} in postscript


set key left top Left opaque reverse spacing 1.2

# {/Symbol \341}, left angle, <
# {/Symbol \361}, right angle, >

plot [0.0:1.1][0:1.7] \
  1  w l ls 9 not, \
  "ehmclM0.02.eb"      u 1:2 w l  ls 1 not, \
  ""                   u 1:2 every 5 w p ls 1 not, \
  -1                   w lp ls 1 t "LJ, single", \
  "ehmcgM0.002.eb"     u 1:2 w l ls 2 not, \
  ""                   u 1:2 every 5 w p ls 2 not, \
  -1                   w lp ls 2 t "LJ, all", \
  "ehsqmclM0.02.eb"    u 1:2 w l ls 4 not, \
  ""                   u 1:2 every 5 w p ls 4 not, \
  -1                   w lp ls 4 t "Square-well, single", \
  "ehismcl.eb"         u 1:2 w l ls 5 not, \
  ""                   u 1:2 every 5 w p ls 5 not, \
  -1                   w lp ls 5 t "Ising model, single"
  
#  "" u 1:2 w l t "Ising model, single"

unset output
set terminal wxt
reset
