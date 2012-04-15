reset
set terminal push
set terminal postscript landscape enhanced "Helvetica" 20


set output "vol.ps"

binh = 0.5
nfreq = 50
ps4 = 1.5
ps6 = 1.5
tcfont = "Helvetica, 14"
tcsfont = "Helvetica, 12"
keyfont = "Helvetica, 16"
zcut = 1e-9
darkgray = "#404040"
lightgray = "#aaaaaa"
nicegray = "#606060"
lightgreen = "#ccffcc"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set multiplot

# (a) distribution

set size 1.0, 1.0
set origin 0.0, 0.0
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "V" offset 0, 0.5 font keyfont
set xtics 500 font tcfont
set mxtics 5
set logscale y
set format y "10^{%L}"
set ytics 10 font tcfont
set mytics 10
set ylabel '~{/Symbol r}{.3\^}(V) ({/Symbol \264} 10^{-3})' offset 0.0, 0 font keyfont
set key font keyfont
pres = 0.115
tp = 1.24
bp = pres / tp
pcorrhis = 0.0905044
plot [400:2100][1e-5:3e-3] \
     "ds.dat"  u 1:($6) w l lt 1 lw 3 t "reference", \
     "dsb.dat" u ($1 + binh):($2/1e5*(bp + $9)/pcorrhis) w l lt 4 lc rgb lightgray t "histogram", \
     "dsb.dat" u 1:($6) every nfreq w p pt 6 ps ps6 lt 1 lw 1.5 lc rgb niceblue not, \
     "dsb.dat" u 1:($6) w l lt 2 lw 1.5 lc rgb niceblue not, \
     1e8 w lp lt 2 lw 1.5 lc rgb niceblue pt 6 ps ps6 t "fractional identity", \
     "dsb.dat" u 1:(exp($7)) every nfreq w p pt 4 ps ps4 lt 1 lw 1.5 lc rgb nicered not, \
     "dsb.dat" u 1:(exp($7)) w l lt 5 lw 1.5 lc rgb nicered not, \
     1e8 w lp lt 5 lw 1.5 lc rgb nicered pt 4 ps ps4 t "m.f. integration", \
     0 w l lt 1 lw 0.5 notitle
#
# about the first line
# we attempt to correct the histogram 
# $9 is the mean force, the difference between beta (pvir - p)
#   so $9 + bp = beta * pvir,
# pcorrhis is the averaged beta*pressure, which serves as a normalization,
#   this can be found from the output of volcorr as
#
# volcorr test, sample size: 100000, p corr 0.0905044
#                                           ~~~~~~~~~
#

unset logscale y

# (a) inset: window size
set size 0.36, 0.36
set origin 0.3, 0.24
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set xlabel "V" offset 0, 1. font tcfont
set xtics 500 offset 0, 0.5 font tcsfont
set mxtics 10
set ylabel "window size" offset 2, 0 font tcfont
unset logscale y
set format "%g"
set ytics 20 font tcsfont
set mytics 2
#set ytics 10 offset .5, 0 font tcsfont 
#set mytics 10
set tics front
plot [400:2200][5:120] \
  "dsb.dat" u 1:($12-$11) every nfreq w p pt 6 ps ps6 lt 1 lw 1 lc rgb niceblue notitle, \
  "dsb.dat" u 1:($12-$11) w l lt 2 lw 1 lc rgb niceblue notitle
unset logscale 

unset multiplot
unset output

set terminal pop 
