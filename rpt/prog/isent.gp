#!/usr/bin/env gnuplot
unset multiplot
reset
set terminal pngcairo enhanced background rgb "#ffffff" font 'Arial, 37' linewidth 2 round size 2048, 1536
set output "isent.png"
set multiplot


# system size
N=1024
N16=256

set size 1.0, 1.0
set origin 0.0, 0.0

# Line style for axes
set style line 80 lt rgb "#404040" lw 2
set border 15 back linestyle 80  # Remove border on top and right.

# Line style for grid
set style line 81 lt 0 lw 2 # dotted
set style line 81 lt rgb "#a0a0a0" # light gray
#set grid mxtics xtics mytics ytics back linestyle 81

set tics

set style line 1 lt rgb "#c04020" lw 3 pt 10  ps 6  # dark red
set style line 2 lt rgb "#e0c040" lw 3 pt 12  ps 6  # yellow to orange
set style line 3 lt rgb "#e030a0" lw 3 pt 6   ps 7  # bright magenta
set style line 4 lt rgb "#80d060" lw 3 pt 8   ps 7  # dark green
set style line 5 lt rgb "#305080" lw 3 pt 4   ps 5  # navy blue
set style line 6 lt rgb "#30e0e0" lw 3 pt 14  ps 5  # cyan

set style line 9 lt rgb "#000000" lw 2 pt 1   ps 4  # black line
set style line 8 lt rgb "#808080" lw 3 pt 1   ps 4  # gray line


set mxtics 5
set xtics 1.0 offset 0, 0.0
set xlabel "U/N" offset 0, 0.0

set mytics 5
set ytics 0.5  offset 0, 0.0
set ylabel "{/Symbol b}_U(U)" offset 1.0, 0

set key right top Left reverse spacing 1.5

plot [-2:2][-1:1] 0 w l ls 8 not, \
  "profis32.dat" u ($1/N):10 w l ls 9 t "Reference", \
  ""             u ($1/N):3  w l ls 1 not, \
  ""             u ($1/N):3  every 32 w p ls 1 not, \
  -10                        w lp ls 1 t "Eq. (2)", \
  ""             u ($1/N):4  w l ls 2 not, \
  ""             u ($1/N):4  every 32 w p ls 2 not, \
  -10                        w lp ls 2 t "Eq. (3')", \
  ""             u ($1/N):5  w l ls 3 not, \
  ""             u ($1/N):5  every 32 w p ls 3 not, \
  -10                        w lp ls 3 t "Eq. (3)", \
  ""             u ($1/N):7  w l ls 4 not, \
  ""             u ($1/N):7  every 32 w p ls 4 not, \
  -10                        w lp ls 4 t "Eq. (8), k = 1", \
  ""             u ($1/N):9  w l ls 5 not, \
  ""             u ($1/N):9  every 32 w p ls 5 not, \
  -10                        w lp ls 5 t "Eq. (9), k = 1"

#  "profis16.dat" u ($1/N16):5  w l ls 6 not, \
#  ""             u ($1/N16):5  every 32 w p ls 6 not, \
#  -10                        w lp ls 6 t "Eq. (6), 16x16"

insetx0 = 0.1
insety0 = 0.1
insetw = 0.52
inseth = 0.40

set size insetw, inseth
set origin insetx0, insety0

# erase the background 
set object 1 rectangle from graph 0,0 to graph 1,1 behind fc rgb "#ffffff"
#set object 1 rectangle from screen insetx0,insety0 to screen insetx0+insetw,insety0+inseth behind fc rgb "#ffffff"

insetfont="Arial,27"
set mxtics 5
set xtics 1.0 offset 0, 0.5 font insetfont
unset xlabel 
# set xlabel "U/N" offset 0, 1.5 font insetfont

set mytics 5
set ytics 0.5  offset 0.5, 0.0  font insetfont
set ylabel "{/=27 {/Symbol e}[log g_U(U)] = {/Symbol=32 \362}@_{/*0.5 0}^{/*.5 U} {/Symbol e}[{/Symbol b}_U(U')] dU' }" offset 2.5, 0 font insetfont
unset key

plot [-2:2][-1:0.5] 0 w l ls 8 not, \
  "profis32.dat" u ($1/N):11  w l ls 3 not, \
  ""             u ($1/N):11  every 32 w p ls 3 ps 4 not, \
  -10                         w lp ls 3 not, \
  ""             u ($1/N):13  w l ls 4 not, \
  ""             u ($1/N):13  every 32 w p ls 4 ps 4 not, \
  -10                         w lp ls 4 not, \
  ""             u ($1/N):15  w l ls 5 not, \
  ""             u ($1/N):15  every 32 w p ls 5 ps 4 not, \
  -10                         w lp ls 5 not
  
  
#  "profis16.dat" u ($1/N16):11  w l ls 6 not, \
#  ""             u ($1/N16):11  every 32 w p ls 6 not, \
#  -10                         w lp ls 6 not

unset multiplot
unset output
set terminal wxt
reset
