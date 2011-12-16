reset
set terminal push
set terminal postscript portrait enhanced "Helvetica" 16


set output "vol.ps"

binh = 0.5
tcfont = "Helvetica, 12"
tcsfont = "Helvetica, 10"
keyfont = "Helvetica, 14"
zcut = 1e-9
lightgray = "#cccccc"
nicegray = "#606060"
lightgreen = "#ccffcc"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set multiplot

set label "(a)" at screen 0.02, 0.98
set label "(b)" at screen 0.02, 0.48

# (a) distribution

set size 1.0, 0.5
set origin 0.0, 0.5
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "V" offset 0, 1 font keyfont
set xtics 500 font tcfont
set mxtics 5
set logscale y
set format y "10^{%L}"
set ytics 10 font tcfont
set mytics 10
set ylabel "{/Symbol r}(V) ({/Symbol \264} 10^{-3})" offset 0.0, 0 font keyfont
set key font keyfont
plot [400:2100][1e-5:3e-3] \
     "dsb.dat" u ($1 + binh):($2/1e5) w l lt 4 lc rgb "#C0C0C0" t "histogram", \
     "dsb.dat" u 1:(exp($7)) w l lt 2 lw 3 lc rgb nicered t "m.f. integration", \
     "dsb.dat" u 1:($6) w l lt 1 lw 4 lc rgb niceblue t "fractional identity", \
     "ds.dat"  u 1:($6) w l lt 1 lw 1 t "reference", \
     0 w l lt 1 lw 0.5 notitle
unset logscale y

# (b) inset: window size
set size 0.38, 0.19
set origin 0.3, 0.62
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set xlabel "V" offset 0, 1. font tcfont
set xtics 500 offset 0, 0.5 font tcsfont
set mxtics 10
set ylabel "window size" offset 2, 0 font tcfont
set logscale y
set format "%g"
set ytics 10 offset .5, 0 font tcsfont 
set mytics 10
set tics front
plot [400:2200][5:400] \
  "dsb.dat" u 1:($12-$11) w l lt 1 lw 4 lc rgb niceblue notitle
unset logscale 

# (b) volume vs pressure
set size 1.0, 0.5
set origin 0.0, 0.0
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "Pressure (p)" font keyfont
set xtics 0.01 offset 0, 0.3 font tcfont
set ylabel "Volume (V)" offset 1, 0 font keyfont
set ytics 500 font tcfont
set format y "%g"
pmin = .095
pmax = .18
vmin = 400
vmax = 2100
set arrow 7 from 0.115, vmin to 0.115, vmax nohead
plot [pmin:pmax][vmin:vmax] \
  "fe.dat" u 1:3 w l lt 4 lw 4 lc rgb nicegray title "histogram", \
  ""       u 1:5 w l lt 2 lw 3 lc rgb nicered title "m.f. integration", \
  ""       u 1:7 w l lt 1 lw 4 lc rgb niceblue title "fractional identity", \
  ""       u 1:9 w l lt 1 lw 1 lc rgb "#000000" title "reference"

# (b) inset: relative error of the volume
set size 0.4, 0.2
set origin 0.5, 0.17
set bmargin 0
set tmargin 0
set lmargin 0
set rmargin 0
unset xlabel
set xtics 0.02 font tcsfont 
set mxtics 10
set logscale y
set ylabel "|{/Symbol e}(V)|/V" offset 2., 0 font tcfont
set format y "10^{%L}"
set ytics 10 offset .5, 0 font tcsfont
set mytics 10
errmin = 1e-5
errmax = 0.1
set arrow 7 from 0.115, errmin to 0.115, errmax nohead
plot [pmin:pmax][errmin:errmax] \
  "fe.dat" u 1:(abs($3/$9 - 1)) w l lt 4 lw 4 lc rgb nicegray notitle, \
  ""       u 1:(abs($5/$9 - 1)) w l lt 2 lw 3 lc rgb nicered notitle, \
  ""       u 1:(abs($7/$9 - 1)) w l lt 1 lw 4 lc rgb niceblue notitle, \
  0 w l lt 1 lw 1 notitle
 
unset multiplot
unset output

set terminal pop 
