#!/usr/bin/env gnuplot
reset
unset multiplot
set terminal pngcairo enhanced font 'Arial, 27' linewidth 2 round size 2048, 1536
set output "basic.png"

set multiplot

dx = 0.005
dy = 0.02
set label 1 "(a)" at screen dx, 1-dy font "Arial, 32"
set label 2 "(b)" at screen 0.5+dx, 1-dy font "Arial, 32"
set label 3 "(c)" at screen dx, 0.5-dy font "Arial, 32"
set label 4 "(d)" at screen 0.5+dx, 0.5-dy font "Arial, 32"


# Line style for axes
set style line 80 lt rgb "#404040" lw 2.0
set border 15 back linestyle 80  # Remove border on top and right.

# Line style for grid
set style line 81 lt 0 lw 1.5 # dotted
set style line 81 lt rgb "#a0a0a0" # light gray
#set grid mxtics xtics mytics ytics back linestyle 81

set tics




set size 0.5, 0.5
set origin 0.0, 0.0

set style line 1 lt rgb "#c04020" lw 2 pt 10  ps 4  # dark red
set style line 2 lt rgb "#e0c040" lw 2 pt 12  ps 4  # yellow
set style line 3 lt rgb "#e030a0" lw 2 pt 6   ps 5  # bright magenta
set style line 4 lt rgb "#80d060" lw 2 pt 8   ps 6  # dark green
set style line 5 lt rgb "#305080" lw 2 pt 4   ps 4  # navy blue
set style line 9 lt rgb "#000000" lw 3 pt 1   ps 4  # black for reference

set logscale x
set xtics offset 0, 0.1
set xlabel "u_{max}" offset 0, 1.0

set mytics 2
set ylabel "{/Symbol b}" offset 1.0, 0
# try {/Symbol-Oblique b} in postscript

set key left bottom Left reverse spacing 1.2

# {/Symbol \341}, left angle, <
# {/Symbol \361}, right angle, >

plot [0.005:0.1][0:] \
  "sizemcl.txt" u 1:14 w lp ls 9 t "Eq. (1)", \
            ""  u 1:2  w lp ls 1 t "Eq. (2), {/Symbol b} = 2{/Symbol \341}{/Symbol e}{/Symbol \361}/{/Symbol \341}{/Symbol e}^2{/Symbol \361}", \
            ""  u 1:3  w lp ls 2 t "Eq. (3'), {/Symbol b} = 2{/Symbol \341}{/Symbol e}{/Symbol \361}/{/Symbol \341}{/Symbol De^2}{/Symbol \361}", \
            ""  u 1:4  w lp ls 3 t "Eq. (3), {/Symbol \341}exp( {/Symbol - b e}){/Symbol \361}= 1", \
            ""  u 1:6  w lp ls 4 t "Eq. (8), k = 1", \
            ""  u 1:8  w lp ls 5 t "Eq. (9), k = 1"

unset logscale x








set origin 0.5, 0.0

set xtics 0.2
set mxtics 2
set xlabel "r_c" offset 2.0, 1.0

set ytics 0.1
set mytics 4
set ylabel "{/Symbol b}" offset 1.0, 0
# try {/Symbol-Oblique b} in postscript

set key right top Left reverse spacing 1.2

plot [1.2:2.5][0.9:1.2] \
  "rcmcl.txt" u 1:14 w lp ls 9  t "Eq. (1)", \
          ""  u 1:4  w lp ls 3  t "Eq. (3), {/Symbol \341}exp( {/Symbol - b e}){/Symbol \361}= 1", \
          ""  u 1:6  w lp ls 4  t "Eq. (8), k = 1", \
          ""  u 1:8  w lp ls 5  t "Eq. (9), k = 1"







set origin 0., 0.5

set xtics 2 offset 0, 0.0
set mxtics 2
set xlabel "{/Symbol e}" offset 1.0, 0.5

set ytics 0.2
set mytics 2
set ylabel "p({/Symbol e})" offset 1.0, 0.0

set style line 11 lt rgb "#c04020" lw 2 pt 10 ps 3  # dark red
set style line 12 lt rgb "#305080" lw 2 pt 8  ps 3  # navy blue
set style line 13 lt rgb "#40ff20" lw 0.5 pt 2 ps 2 # magenta
set style line 14 lt rgb "#c020a0" lw 1 pt 4 ps 2   # green

set key right top Left reverse spacing 2.0

#set arrow 1 from 1.40,1.01 to 0.1, 1 ls 12 lw 2.0 front 
#set arrow 2 from 1.40,0.86 to 0.23, 0.83 ls 11 front
#set arrow 3 from 1.40,0.72 to 0.24, 0.438 ls 14 lw 2.0 front 
#set arrow 4 from 2.40,0.54 to 1.30, 0.257 ls 13 lw 2.0 front 

bin(x,w)=w*(floor(x/w)+.5)

#plot [-3:12][0:1.1] \
#  "ehmclM0.05.dat"    u (bin($1, 0.1)):($2*0.05) smooth freq with boxes fs transparent solid 0.5 noborder ls 12 t "u_{max} = 0.01, single, MC", \
#  "ehmdlM0.05.dat"    u (bin($1, 0.1)):($2*0.05) smooth freq with boxes fs empty border ls 11 t "u_{max} = 0.01, single, MD", \
#  "ehmclM0.1.dat"    u (bin($1, 0.05)):($2*0.2) smooth freq with boxes fs transparent solid 0.8 noborder ls 14 t "u_{max} = 0.02, single, MC", \
#  "ehmcgM0.01.dat"   u (bin($1, 0.05)):($2*0.2) smooth freq with boxes fs transparent solid 0.4 noborder ls 13 t "u_{max} = 0.002, all, MC"

plot [-2:7][0:1.3] \
  "ehmdlM0.05.dat"   u 1:2 every 40::108  with p  ls 1  ps 3   not, \
  ""                 u 1:2                with l  ls 1         not, \
  -1                                      with lp ls 1  ps 3   t "u_{max} = 0.05, single, MD", \
  "ehmclM0.05.dat"   u 1:2 every 40::156  with p  ls 4  ps 3   not, \
  ""                 u 1:2                with l  ls 4         not, \
  -1                                      with lp ls 4  ps 3   t "u_{max} = 0.05, single, MC", \
  "ehmdlM0.1.dat"    u 1:2 every 20::8    with p  ls 3  ps 3   not, \
  ""                 u 1:2                with l  ls 3         not, \
  -1                                      with lp ls 3  ps 3   t "u_{max} = 0.1, single, MC", \
  "ehmcgM0.01.dat"   u 1:2 every 20::75   with p  ls 5  ps 2.5 not, \
  ""                 u 1:2                with l  ls 5         not, \
  -1                                      with lp ls 5  ps 2.5 t "u_{max} = 0.01, all, MC", \
  -1  not

#unset arrow 1
#unset arrow 2
#unset arrow 3
#unset arrow 4


set origin 0.5, 0.5


set xtics .2
set mxtics 2
set xlabel "{/Symbol b}" offset 0, 0.5

set ytics .2
set mytics 2
set ylabel "{/Symbol \341}exp({/Symbol - b e}){/Symbol \361}" offset 1.0, 0
# try {/Symbol-Oblique b} in postscript

set arrow 1 from 1, 0.94 to 1, 0.05 ls 9

set key left top Left opaque reverse spacing 1.2

# {/Symbol \341}, left angle, <
# {/Symbol \361}, right angle, >

plot [0.0:1.1][0:1.7] \
  1  w l ls 9 not, \
  "ehmclM0.1.eb"      u 1:2 w l  ls 1 not, \
  ""                   u 1:2 every 5 w p ls 1 not, \
  -1                   w lp ls 1 t "LJ, single", \
  "ehmcgM0.01.eb"     u 1:2 w l ls 2 not, \
  ""                   u 1:2 every 5 w p ls 2 not, \
  -1                   w lp ls 2 t "LJ, all", \
  "ehsqmclM0.1.eb"    u 1:2 w l ls 4 not, \
  ""                   u 1:2 every 5 w p ls 4 ps 4 not, \
  -1                   w lp ls 4 ps 4 t "Square-well, single", \
  "ehismcl.eb"         u 1:2 w l ls 5 not, \
  ""                   u 1:2 every 5 w p ls 5 ps 4 not, \
  -1                   w lp ls 5 ps 4 t "Ising model, single"

unset arrow 1



unset multiplot
unset output
set terminal wxt
reset
