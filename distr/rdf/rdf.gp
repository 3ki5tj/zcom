reset
set terminal push
set terminal postscript portrait enhanced "Helvetica" 16


set output "rdf.ps"

bin = 0.002
binh = bin * 0.5
nfreq = 20
ps4 = 0.9
ps6 = 0.9
# normalization factor for rho = 0.7, and 5 samples
hisnorm(x) = 1.0/(.5*0.7*255*bin*4*pi*(x+binh)**2)/5

tcfont = "Helvetica, 12"
tcsfont = "Helvetica, 10"
keyfont = "Helvetica, 14"
zcut = 1e-9
darkgray = "#404040"
lightgray = "#aaaaaa"
nicegray = "#cccccc"
lightgreen = "#ccffcc"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set multiplot

set label "(a)" at screen 0.02, 0.98
set label "(b)" at screen 0.02, 0.48

rmin = 0.7
rmax = 3.5

# (a) distribution at T = 0.85

set size 1.0, 0.5
set origin 0.0, 0.5
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "r" offset 0, 0 font keyfont
set xtics 1 font tcfont
set mxtics 5
set format y "%g"
set ytics 1 font tcfont
set mytics 4
set ylabel "Radial distribution function g(r), T_1 = 0.85" offset 0.0, 0 font keyfont
set key font keyfont
plot [rmin:rmax][0:3] \
     "dsb0.85.dat" u ($1 + binh):($2*hisnorm($1)) every nfreq w p pt 8 lt 1 lw 1 lc rgb darkgray t "histogram", \
     "dsb0.85.dat" u 1:($6) every nfreq w p pt 6 ps ps6 lt 1 lw 1.5 lc rgb niceblue t "fractional identity", \
     "dsbaj0.85.dat" u 1:($6) every nfreq w p pt 4 ps ps4 lt 1 lw 1.5 lc rgb nicered t "AJ identity", \
     "ds0.85.dat"  u 1:($6) w l lt 1 lw 3 t "reference", \
     1 w l lt 1 lw 0.5 lc rgb nicegray notitle
unset logscale y

#     "dsb0.85.dat" u ($1 + binh):($2*hisnorm($1)) w l lw 1 lt 4 lc rgb lightgray t "histogram",
#     "dsbaj0.85.dat" u 1:($6) w l lt 1 lw 1 lc rgb nicered notitle,
#     "dsb0.85.dat" u 1:($6) w l lt 1 lw 1 lc rgb niceblue notitle,


# "dsaj0.85.dat" u 1:($6) w l lt 6 lw 3 lc rgb nicegreen t "AJ, larger data set", \

# (b) distribution at T = 0.4
set size 1.0, 0.5
set origin 0.0, 0.0
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "r" offset 0, 0 font keyfont
set xtics 1 font tcfont
set mxtics 5
set format y "%g"
set ytics 1 font tcfont
set mytics 4
set ylabel "Radial distribution function g(r), T_2 = 0.4" offset 0.0, 0 font keyfont
set key font keyfont
plot [rmin:rmax][:5.4] \
     "dsb0.4.dat" u ($1 + binh):($2*hisnorm($1)) every nfreq w p pt 8 lt 1 lw 1 lc rgb darkgray t "histogram", \
     "dsb0.4.dat" u 1:($6) every nfreq w p pt 6 ps ps6 lt 1 lw 1.5 lc rgb niceblue t "fractional identity", \
     "dsbaj0.4.dat" u 1:($6) every nfreq w p pt 4 ps ps4 lt 1 lw 1.5 lc rgb nicered t "AJ identity", \
     "ds0.4.dat"  u 1:($6) w l lt 1 lw 3 t "reference", \
     1 w l lt 1 lw 0.5 lc rgb nicegray notitle
unset logscale y

#     "dsb0.4.dat" u ($1 + binh):($2*hisnorm($1)) w l lw 1 lt 4 lc rgb lightgray t "histogram",
#     "dsbaj0.4.dat" u 1:($6) w l lt 1 lw 1 lc rgb nicered notitle,
#     "dsbaj0.4.dat" u 1:($6) every 25 w p pt 4 ps 0.8 lt 1 lw 1.5 lc rgb nicered t "AJ identity",
#     "dsb0.4.dat" u 1:($6) w l lt 1 lw 1 lc rgb niceblue notitle,
#     "dsb0.4.dat" u 1:($6) every 25 w p pt 6 ps 1.0 lt 1 lw 1.5 lc rgb niceblue t "fractional identity",
#     "ds0.4.dat"  u 1:($6) w l lt 1 lw 1 t "reference",
# "dsaj0.4.dat" u 1:($6) w l lt 6 lw 3 lc rgb nicegreen t "AJ, larger data set",
 
unset multiplot
unset output

set terminal pop 
