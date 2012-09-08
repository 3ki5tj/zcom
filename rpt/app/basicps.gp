#!/usr/bin/env gnuplot
reset
unset multiplot
set terminal postscript enhanced font 'Arial, 14' linewidth 1
set output "basic.ps"

set multiplot

dx = 0.005
dy = 0.02
set label 1 "(a)" at screen dx, 1-dy font ", 17"
set label 2 "(b)" at screen 0.5+dx, 1-dy font ", 17" 
set label 3 "(c)" at screen dx, 0.5-dy font ", 17" 
set label 4 "(d)" at screen 0.5+dx, 0.5-dy font ", 17"

# Line style for axes
set style line 80 lt 1
set style line 80 lt rgb "#555555" lw 2.0
set border 15 back linestyle 80  # Remove border on top and right.

# Line style for grid
set style line 81 lt 0 lw 1.5 # dotted
set style line 81 lt rgb "#a0a0a0" # light gray
set grid mxtics xtics mytics ytics back linestyle 81

set tics font ", 12"
set key font ", 12" spacing 1.5
set tmargin 1.0



set size 0.5, 0.5
set origin 0.0, 0.0

set style line 9 lt 1
set style line 9 lt rgb "#000000" lw 3 pt 1   ps 2.0  # black for reference
set style line 1 lt rgb "#c04020" lw 2 pt 10  ps 1.5  # dark red
set style line 2 lt rgb "#e0c040" lw 2 pt 12  ps 1.5  # yellow
set style line 3 lt rgb "#e030a0" lw 2 pt 6   ps 1.5  # bright magenta
set style line 4 lt rgb "#80d060" lw 2 pt 8   ps 1.8  # dark green
set style line 5 lt rgb "#305080" lw 2 pt 4   ps 1.5  # navy blue

set logscale x
set xtics offset 0, 0.1
set xlabel "{/Arial-Italic u}_{max}" offset 0, 1.0

set mytics 2
set ylabel "{/Symbol-Oblique b}" offset 2.0, 0

set key left bottom Left reverse width -9

# {/Symbol \341}, left angle, <
# {/Symbol \361}, right angle, >

plot [0.005:0.1][0:] \
  "sizemcl.txt" u 1:14 w lp ls 9 t "Eq. (4)", \
            ""  u 1:2  w lp ls 1 t "Eq. (5), {/Symbol b} = 2{/Symbol \341}{/Symbol e}{/Symbol \361}/{/Symbol \341}{/Symbol e}^2{/Symbol \361}", \
            ""  u 1:3  w lp ls 2 t "Eq. (6'), {/Symbol b} = 2{/Symbol \341}{/Symbol e}{/Symbol \361}/{/Symbol \341}{/Symbol De^2}{/Symbol \361}", \
            ""  u 1:4  w lp ls 3 t "Eq. (6), {/Symbol \341}exp( {/Symbol - b e}){/Symbol \361}= 1", \
            ""  u 1:6  w lp ls 4 t "Eq. (11), k = 1", \
            ""  u 1:8  w lp ls 5 t "Eq. (12), k = 1"

unset logscale x








set origin 0.5, 0.0

set xtics 0.2
set mxtics 2
set xlabel "{/Arial-Italic r_c}" offset 2.0, 1.0

set ytics 0.1
set mytics 4
set ylabel "{/Symbol-Oblique b}" offset 2., 0

set key right top Left reverse 

plot [1.2:2.5][0.9:1.2] \
  "rcmcl.txt" u 1:14 w lp ls 9  t "Eq. (4)", \
          ""  u 1:4  w lp ls 3  t "Eq. (6), {/Symbol \341}exp( {/Symbol - b e}){/Symbol \361}= 1", \
          ""  u 1:6  w lp ls 4  t "Eq. (11), k = 1", \
          ""  u 1:8  w lp ls 5  t "Eq. (12), k = 1"


set origin 0., 0.5

set xtics 2 offset 0, 0.0
set mxtics 2
set xlabel "{/Symbol-Oblique e}" offset 1.0, 0.5

set ytics 0.2
set mytics 2
set ylabel "{/Arial-Italic p}({/Symbol-Oblique e})" offset 1.0, 0.0

set key right top Left reverse width -6 

plot [-2:8][0:1.3] \
  "ehmdlM0.05.dat"   u 1:2 every 40::108  with p  ls 1 ps 1.0   not, \
  ""                 u 1:2                with l  ls 1          not, \
  -1                                      with lp ls 1 ps 1.0   t "u_{max} = 0.05, single, MD", \
  "ehmclM0.05.dat"   u 1:2 every 40::156  with p  ls 4 ps 1.0   not, \
  ""                 u 1:2                with l  ls 4          not, \
  -1                                      with lp ls 4 ps 1.0   t "u_{max} = 0.05, single, MC", \
  "ehmdlM0.1.dat"    u 1:2 every 20::8    with p  ls 3 ps 1.0   not, \
  ""                 u 1:2                with l  ls 3          not, \
  -1                                      with lp ls 3 ps 1.0   t "u_{max} = 0.1, single, MC", \
  "ehmcgM0.01.dat"   u 1:2 every 20::75   with p  ls 5 ps 1.0   not, \
  ""                 u 1:2                with l  ls 5          not, \
  -1                                      with lp ls 5 ps 1.0   t "u_{max} = 0.01, all, MC", \
  -1  not


set origin 0.5, 0.5


set xtics .2
set mxtics 2
set xlabel "{/Symbol-Oblique b}" offset 0, 0.5

set ytics .2
set mytics 2
set ylabel "{/Symbol \341}exp( -{/Symbol-Oblique b e} ){/Symbol \361}" offset 2., 0
# try {/Symbol-Oblique b} in postscript

set arrow 1 from 1, 0.94 to 1, 0.05 ls 9

set key left top Left opaque reverse width -6

# {/Symbol \341}, left angle, <
# {/Symbol \361}, right angle, >

plot [0.0:1.1][0:2.0] \
  1                                   w l   ls 9 not, \
  "ehmclM0.1.eb"        u 1:2         w l   ls 1 not, \
  ""                    u 1:2 every 5 w p   ls 1 not, \
  -1                                  w lp  ls 1 t "LJ, single", \
  "ehmcgM0.01.eb"       u 1:2         w l   ls 2 not, \
  ""                    u 1:2 every 5 w p   ls 2 not, \
  -1                                  w lp  ls 2 t "LJ, all", \
  "ehsqmclM0.1.eb"      u 1:2         w l   ls 4 not, \
  ""                    u 1:2 every 5 w p   ls 4 not, \
  -1                                  w lp  ls 4 t "Square-well, single", \
  "ehismcl.eb"          u 1:2         w l   ls 5 not, \
  ""                    u 1:2 every 5 w p   ls 5 ps 1.5 not, \
  -1                                  w lp  ls 5 ps 1.5 t "Ising model, single"

unset arrow 1



unset multiplot
unset output
set terminal wxt
reset
