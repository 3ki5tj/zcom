# gray scale

unset multiplot
reset
set terminal push
set terminal postscript landscape enhanced "Helvetica" 14 

set output "ref.ps"

binh = 0.05  # half bin width
nfreq = 20
ps4 = 0.8
ps6 = 0.8
tcfont = "Helvetica, 9"
tcsfont = "Helvetica, 7"
keyfont = "Helvetica, 10"
zcut = 1e-14
darkgray = "#404040"
lightgray = "#cccccc"
nicegray = "#aaaaaa"
lightgreen = "#99ffaa"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set xrange [-1410:-1200]

set multiplot

set label 1 "(a)" at screen 0.01, 0.98
set label 2 "(b)" at screen 0.01, 0.48
set label 3 "(c)" at screen 0.51, 0.98
set label 4 "(d)" at screen 0.51, 0.48
#########################################################################################
# (a) energy distribution

set size 0.5, 0.5
set origin 0.0, 0.5
set tmargin 1.0
set bmargin 2.0
set lmargin 6.0
set rmargin 3.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
set xlabel "U" offset 0, 1 font keyfont
set ytics font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(U)" offset 2, 0 font keyfont
set key font keyfont
set logscale y
plot [-1380:-1240][1e-7:1] \
     "../ene/data/ds_c_ii.dat" u ($1):($6 < zcut ? zcut : $6)  w l lt 4 lw 2 lc rgb darkgray   t "weighted histogram", \
     "../ene/data/ds_c_ii.dat" u 1:6 w l lt 1 lw 2 lc rgb niceblue t "fractional identity", \
     "../ene/data/ds_c_aj.dat" u 1:6 w l lt 2 lw 2 lc rgb nicered t "AJ identity"


# using the WHAM reference
#plot [][1e-10:1] \
#     "../ene/data/dsmerge0.dat" u ($1):($6 < zcut ? zcut : $6)  w l lt 4 lw 2 lc rgb darkgray   t "weighted histogram", \
#     "../ene/data/dsmerge.dat"   u 1:6 w l lt 1 lw 2 lc rgb niceblue t "fractional identity"
unset logscale
#########################################################################################
# volume distribution

set origin 0.0, 0.0
set xlabel "V" offset 0, 1.0 font keyfont
set xtics 500 font tcfont
set mxtics 5
set logscale y
set format y "10^{%L}"
set ytics 10 font tcfont
set mytics 10
set ylabel '~{/Symbol r}{.3\^}(V) ({/Symbol \264} 10^{-3})' offset 1.0, 0 font keyfont
set key font keyfont
pres = 0.115
tp = 1.24
bp = pres / tp
sampsiz = 1e8
pcorrhis = 0.0901898
plot [400:2100][1e-5:5e-3] \
     "../vol/data/ds.dat" u ($1 + binh):($2/sampsiz*(bp + $9)/pcorrhis) w l lw 2.0 lt 4 lc rgb darkgray t "histogram", \
     "../vol/data/ds.dat" u 1:($6) w l lt 1 lw 2 lc rgb niceblue t "fractional identity", \
     "../vol/data/ds.dat" u 1:(exp($7)) w l lt 2 lw 2 lc rgb nicered t "m.f. integration", \
     0 w l lt 1 lw 0.5 notitle



unset logscale
#########################################################################################
# RDF
rho = 0.7
nsamp = 5000
hisnorm(x) = 1.0/(.5*rho*255*bin*4*pi*(x+binh)**2)/nsamp
bin = 0.002
binh = bin * 0.5
rmin = 0.7
rmax = 3.5

# (c) RDF at T = 0.85

set origin 0.5, 0.5
set rmargin 3.0
set tmargin 1.0
set xlabel "r" offset 0, 1 font keyfont
set xtics 1 font tcfont
set mxtics 5
set format y "%g"
set ytics 1 font tcfont
set mytics 4
set ylabel "g(r), T_1 = 0.85" offset 0.0, 0 font keyfont
set key font keyfont
plot [rmin:rmax][0:3] \
     "../rdf/data/ds0.85.dat" u ($1 + binh):($2*hisnorm($1)) w l lt 4 lw 2 lc rgb darkgray t "histogram", \
     "../rdf/data/ds0.85.dat" u 1:($6) w l lt 1 lw 2 lc rgb niceblue t "fractional identity", \
     "../rdf/data/dsaj0.85.dat" u 1:($6) w l lt 2 lw 2 lc rgb nicered t "AJ identity", \
     1 w l lt 1 lw 0.5 lc rgb nicegray notitle
unset logscale y


# (d) RDF at T = 0.4
set origin 0.5, 0.0
set rmargin 3.0
set tmargin 1.0
set xtics 1 font tcfont
set mxtics 5
set format y "%g"
set ytics 1 font tcfont
set mytics 4
set ylabel "g(r), T_2 = 0.4" offset 0.0, 0 font keyfont
set key font keyfont
plot [rmin:rmax][:6] \
     "../rdf/data/ds0.4.dat" u ($1 + binh):($2*hisnorm($1)) w l lt 4 lw 2 lc rgb darkgray t "histogram", \
     "../rdf/data/ds0.4.dat" u 1:($6) w l lt 1 lw 2 lc rgb niceblue t "fractional identity", \
     "../rdf/data/dsaj0.4.dat" u 1:($6) w l lt 2 lw 2 lc rgb nicered t "AJ identity", \
     1 w l lt 1 lw 0.5 lc rgb nicegray notitle
unset logscale y



unset multiplot
unset output
reset
set terminal pop

