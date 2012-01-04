reset
set terminal push
set terminal postscript landscape enhanced "Helvetica" 16


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

# (a) distribution

set size 1.0, 1.0
set origin 0.0, 0.0
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

# (a) inset: window size
set size 0.38, 0.38
set origin 0.3, 0.24
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

unset multiplot
unset output

set terminal pop 
