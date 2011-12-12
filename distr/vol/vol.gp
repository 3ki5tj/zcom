reset
set terminal push
# large size 10,14
set terminal postscript landscape enhanced "Helvetica" 20


set output "vol.ps"

binh = 0.5

#set logscale y
set multiplot

set size 1.0, 1.0
set origin 0.0, 0.0
set rmargin 3.0
set lmargin 8.0
unset tmargin
unset bmargin
set xlabel "V"
set xtics 500
set mxtics 5
set ytics .5
set mytics 5
set ylabel "{/Symbol r}(V) ({/Symbol \264} 10^{-3})" offset 1.0,0
plot [300:2050][:] \
     "ds.dat" u ($1 + binh):($2*1e3) w l lt 4 lc rgb "#C0C0C0" t "histogram", \
     "ds.dat" u 1:($6*1e3) w l lt 1 lw 2 t "fractional identity", \
     0 w l lt 1 lw 0.5 notitle

unset multiplot
unset output

set terminal pop 
