unset multiplot
reset
set terminal push
set terminal pngcairo enhanced font "Helvetica, 12" lw 2 size 1024, 512

set output "bb.png"
set multiplot

r2d = 180/pi
r2d2 = r2d*r2d
binh = .5 

#############################################################################################
reset
set size .55, .55
set size square
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0
set pm3d map
set xlabel "{/Symbol f}" offset 0, 1 font "Helvetica, 16"
set xtics -180,120,180 offset 0, 0.5 
set mxtics 6
set ylabel "{/Symbol y}" font "Helvetica, 16"
set ytics -180,120,180
set mytics 6
unset zlabel
unset ztics
set colorbox
#set contour
#set cntrparam levels 3
#set palette defined (0 1 1 1, 1 .3 .3 .3, 2 .1 .1 .1, 3 0 0 0) 

set palette rgbformulae 30,31,32 
tot = 3.68733e+07
dx = pi/180
den = tot*dx*dx


set label "log[{/Symbol r}({/Symbol f}, {/Symbol y})]" \
  at screen 0.45, 0.95 font "Helvetica, 16" enhanced

#set view 0, 0, 1.4, 1.4
set size 0.5, 1.0
set origin 0, .0
set title "(a) histogram" font "Helvetica, 16"
#set format cb "1.t{/Symbol \264}10^{%T}"
logmin = -20
logmax = -8
set cbtics 2 offset -.5, 0 font "Helvetica, 12"
set mcbtics 2
set cbrange [logmin:logmax]
splot [-180:180][-180:180] \
        "bbii.txt" u (r2d*$1 + binh):(r2d*$2 + binh):(log($3/den/r2d2)) w pm3d notitle
unset cbrange

set origin .5, 0.0
set title "(b) fractional identity" font "Helvetica, 16"
splot [-180:180][-180:180] \
        "bbii.txt" u (r2d*$1):(r2d*$2):(log($9/r2d2)) w pm3d notitle


unset multiplot
unset output
set terminal pop
reset

