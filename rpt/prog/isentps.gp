#!/usr/bin/env gnuplot
unset multiplot
reset

set encoding cp1250 # make minus sign longer
set terminal postscript enhanced font 'Arial, 11' size 10, 3.5
set output "isent.ps"
set multiplot


# system size
N=1024
N16=256


dx = 0.005
dy = 0.03

set label "(a)" at screen dx, 1 - dy
set label "(b)" at screen 0.5 + dx, 1 - dy

set size 0.5, 1.0
set origin 0.0, 0.0

# Line style for axes
set style line 80 lt 1
set style line 80 lt rgb "#555555" lw 2
set border 15 back linestyle 80  # Remove border on top and right.

# Line style for grid
set style line 81 lt 0 lw 1 # dotted
set style line 81 lt rgb "#cccccc" # light gray
#set grid mxtics xtics mytics ytics back linestyle 81

#set tics font "Arial, 18"

set style line 1 lt rgb "#c04020" lw 2 pt 10  ps 1.5  # dark red
set style line 2 lt rgb "#e0c040" lw 2 pt 12  ps 1.5  # yellow to orange
set style line 3 lt rgb "#e030a0" lw 2 pt 6   ps 1.5  # bright magenta
set style line 4 lt rgb "#80d060" lw 2 pt 8   ps 1.5  # dark green
set style line 5 lt rgb "#4060cc" lw 2 pt 4   ps 1.  # navy blue
set style line 6 lt rgb "#30e0e0" lw 2 pt 14  ps 1.  # cyan

set style line 9 lt rgb "#000000" lw 1 pt 1   ps 1  # black line
set style line 8 lt rgb "#808080" lw 2 pt 1   ps 1  # gray line

set style line 13 lt rgb "#a01060" lw 3 pt 7   ps 1.5  # bright magenta
set style line 14 lt rgb "#20a020" lw 3 pt 9   ps 1.5  # dark green
set style line 15 lt rgb "#305080" lw 3 pt 5   ps 1.  # navy blue

#set style line 13 lt 1
#set style line 14 lt 1
#set style line 15 lt 1

set mxtics 5
set xtics 1.0 offset 0, 0.0
set xlabel "{/Arial-Italic U}/{/Arial-Italic N}" offset 0, 0.0

set mytics 5
set ytics 0.5  offset 0, 0.0
set ylabel "{/Symbol-Oblique b}_{/Arial-Italic U} ({/Arial-Italic U} )" offset 2.0, 0

set key right top Left reverse spacing 2.0

plot [-2:2][-1:1] -100 w l ls 8 not, \
  "profis32.dat" u ($1/N):10 w l ls 9 t "Reference", \
  ""             u ($1/N):3  w l ls 1 not, \
  ""             u ($1/N):3  every 32 w p ls 1 not, \
  -10                        w lp ls 1 t "Eq. (2)", \
  ""             u ($1/N):4  w l ls 2 not, \
  ""             u ($1/N):4  every 32 w p ls 2 not, \
  -10                        w lp ls 2 t "Eq. (3')", \
  ""             u ($1/N):5  w l ls 3 not, \
  ""             u ($1/N):5  every 32 w p ls 3 not, \
  -10                        w lp ls 3 t "Eq. (4), {/Arial-Italic k} = 0", \
  ""             u ($1/N):7  w l ls 4 not, \
  ""             u ($1/N):7  every 32 w p ls 4 not, \
  -10                        w lp ls 4 t "Eq. (8), {/Arial-Italic k} = 1", \
  ""             u ($1/N):9  w l ls 5 not, \
  ""             u ($1/N):9  every 32 w p ls 5 not, \
  -10                        w lp ls 5 t "Eq. (9), {/Arial-Italic k} = 1"

#  "profis16.dat" u ($1/N16):5  w l ls 6 not, \
#  ""             u ($1/N16):5  every 32 w p ls 6 not, \
#  -10                        w lp ls 6 t "Eq. (3), 16x16"

#insetx0 = 0.1
#insety0 = 0.085
#insetw = 0.51
#inseth = 0.42

set size 0.5, 1.0
set origin 0.5, 0.0

# erase the background 
#set object 1 rectangle from graph 0,0 to graph 1,1 behind fc rgb "#ffffff"
#set object 1 rectangle from screen insetx0,insety0 to screen insetx0+insetw,insety0+inseth behind fc rgb "#ffffff"

insetfont=", 11"
set mxtics 5
set xtics 1.0 offset 0, 0.0 font insetfont
# set xlabel "U/N" offset 0, 1.5 font insetfont

set mytics 5
set ytics 0.5  offset 0.0, 0.0  font insetfont
# {/Symbol \362} is the integral sign
# make it a subscript but with larger font
# &{i} is a thin space
set ylabel "{/=11 {/Symbol-Oblique D}&{i}[&{i}log&{i}{/Arial-Italic g}&{i}({/Arial-Italic U}&{i})] = &{i}_{/*2.0 {/Symbol-Oblique \362}}@_{/*0.8 &{i}0}^{/*.8 &{n}{/Arial-Italic U}} {/Symbol-Oblique D}&{i}[&{i}{/Symbol-Oblique b}_{/Arial-Italic U}&{i}({/Arial-Italic U'}&{i})] d{/Arial-Italic U'} }" offset 1.0, 0 font insetfont

set key right bottom width -10

delx = 0.01

plot [-2.0:0][-1.5:0.5] \
  "profis32.dat" u ($1/N):11  w l ls 3 not, \
  ""             u ($1/N):11  every 32 w p ls 3 ps 1.2 not, \
  -10                         w lp ls 3 t "Eq. (4), {/Arial-Italic k} = 0", \
  "is32x.dat"    u ($1/N):((abs($1/N) > delx) ? ($7+$9+$10+7.38) : 0.) w l ls 13 not, \
  ""             u ($1/N):((abs($1/N) > delx) ? ($7+$9+$10+7.38) : 0.) every 32 w p ls 13 ps 1.2 not, \
  -10                         w lp ls 13 t "Eq. (10), {/Arial-Italic k} = 0, [Eq. (4) corrected]", \
  "profis32.dat" u ($1/N):13  w l ls 4 not, \
  ""             u ($1/N):13  every 32 w p ls 4 ps 1.5 not, \
  -10                         w lp ls 4 t "Eq. (8), {/Arial-Italic k} = 1", \
  "is32x.dat"    u ($1/N):((abs($1/N) > delx) ? ($15+$17+$18-2.04) : 0.)  w l ls 14 not, \
  ""             u ($1/N):((abs($1/N) > delx) ? ($15+$17+$18-2.04) : 0.) every 32 w p ls 14 not, \
  -10                         w lp ls 14 t "Eq. (D4), {/Arial-Italic k} = 1 [Eq. (8) corrected]", \
  "profis32.dat" u ($1/N):15  w l ls 5 not, \
  ""             u ($1/N):15  every 32 w p ls 5 ps 1.2 not, \
  -10                         w lp ls 5 t "Eq. (9), {/Arial-Italic k} = 1", \
  "is32x.dat"    u ($1/N):((abs($1/N) > delx) ? ($22+$23) : 1/0) w l ls 15 not, \
  ""             u ($1/N):((abs($1/N) > delx) ? ($22+$23) : 1/0) every 32 w p ls 15 not, \
  -10                         w lp ls 15 t "Eq. (D5), {/Arial-Italic k} = 1 [Eq. (9) corrected]"
  
  
#  "profis16.dat" u ($1/N16):11  w l ls 6 not, \
#  ""             u ($1/N16):11  every 32 w p ls 6 not, \
#  -10                         w lp ls 6 not

unset multiplot
unset output
set terminal wxt
reset
