reset
set terminal push
set terminal postscript portrait enhanced "Helvetica" 16


set output "rdf.ps"

binh = 0.005
# normalization factor for rho = 0.7, and 5 samples
hisnorm(x) = 1.0/(.5*0.7*255*.01*4*pi*(x+binh)**2)/5

tcfont = "Helvetica, 12"
tcsfont = "Helvetica, 10"
keyfont = "Helvetica, 14"
zcut = 1e-9
lightgray = "#606060"
nicegray = "#aaaaaa"
lightgreen = "#ccffcc"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set multiplot

set label "(a)" at screen 0.02, 0.98
set label "(b)" at screen 0.02, 0.48

rmin = 0.7
rmax = 3.5

# (a) distribution

set size 1.0, 0.5
set origin 0.0, 0.5
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "r" offset 0, 1 font keyfont
set xtics 1 font tcfont
set mxtics 5
set format y "%g"
set ytics 1 font tcfont
set mytics 4
set ylabel "Radial distribution function g(r), T_1 = 0.85" offset 0.0, 0 font keyfont
set key font keyfont
plot [rmin:rmax][0:3] \
     "dsb0.85.dat" u ($1 + binh):($2*hisnorm($1)) w l lw 4 lt 4 lc rgb lightgray t "histogram", \
     "dsbaj0.85.dat" u 1:($6) w l lt 2 lw 3 lc rgb nicered t "AJ identity", \
     "dsb0.85.dat" u 1:($6) w l lt 1 lw 4 lc rgb niceblue t "fractional identity", \
     "dsaj0.85.dat" u 1:($6) w l lt 6 lw 3 lc rgb nicegreen t "AJ, larger data set", \
     "ds0.85.dat"  u 1:($6) w l lt 1 lw 1 t "reference", \
     1 w l lt 1 lw 0.5 lc rgb nicegray notitle
unset logscale y


# (b) volume vs pressure
set size 1.0, 0.5
set origin 0.0, 0.0
set rmargin 3.0
set lmargin 8.0
set tmargin 1.0
set bmargin 3.0
set xlabel "r" offset 0, 1 font keyfont
set xtics 1 font tcfont
set mxtics 5
set format y "%g"
set ytics 1 font tcfont
set mytics 4
set ylabel "Radial distribution function g(r), T_2 = 0.4" offset 0.0, 0 font keyfont
set key font keyfont
plot [rmin:rmax][:5.4] \
     "dsb0.4.dat" u ($1 + binh):($2*hisnorm($1)) w l lw 4 lt 4 lc rgb lightgray t "histogram", \
     "dsbaj0.4.dat" u 1:($6) w l lt 2 lw 3 lc rgb nicered t "AJ identity", \
     "dsb0.4.dat" u 1:($6) w l lt 1 lw 4 lc rgb niceblue t "fractional identity", \
     "dsaj0.4.dat" u 1:($6) w l lt 6 lw 3 lc rgb nicegreen t "AJ, larger data set", \
     "ds0.4.dat"  u 1:($6) w l lt 1 lw 1 t "reference", \
     1 w l lt 1 lw 0.5 lc rgb nicegray notitle
unset logscale y

 
unset multiplot
unset output

set terminal pop 
