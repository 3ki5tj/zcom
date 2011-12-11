# gray scale

reset
set terminal push
# large size 10,14
set terminal postscript landscape enhanced "Helvetica" 14 

set output "eh.ps"

binh = 0.05  # half bin width
tcfont = "Helvetica, 9"
tcsfont = "Helvetica, 7"
keyfont = "Helvetica, 10"
zcut = 1e-9
lightgray = "#cccccc"
lightgreen = "#ccffcc"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set xrange [-1380:-1240]

set multiplot

set label 1 "(a)" at screen 0.01, 0.98
set label 2 "(b)" at screen 0.51, 0.98
set label 3 "(c)" at screen 0.01, 0.48
set label 4 "(d)" at screen 0.51, 0.48
#########################################################################################
# (a) canonical ensemble, basic comparison

set size 0.5, 0.5
set origin 0.0, 0.5
set tmargin 1.0
set bmargin 1.5
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
unset xlabel
set ytics font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(E)" offset 1, 0 font keyfont
set key font keyfont
set logscale y
plot [][1e-6:.5] \
     "dsb_c_ii.dat"     u ($1+binh):($2/1e3) w l lt 4 lc rgb lightgray t "histogram", \
     "dsb_c_ii.dat"     u 1:6 w l lt 1 lw 3 lc rgb niceblue t "fractional (this work)", \
     "dsb_c_aj.dat"     u 1:($6 > zcut ? $6 : zcut) w l lt 2 lw 3 lc rgb nicered t "AJ identity", \
     "ds_c_ix.dat"      u 1:6 w l lt 1 lw 1  t "reference"

#########################################################################################
# (b) canonical ensemble, improving mean force

set size 0.5, 0.5
set origin 0.5, 0.5
set tmargin 1.0
set bmargin 1.5
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
unset xlabel
set ytics font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(E)" offset 1, 0 font keyfont
set key font keyfont
set logscale y
plot [][1e-6:.5] \
     "dsb_c_ii.dat"       u ($1+binh):6 w l lt 1 lw 4 lc rgb lightgreen t "mf. from a single bin", \
     "dsb_c_ii_mfii.dat"  u ($1+binh):6 w l lt 1 lw 3 lc rgb niceblue   t "unbiased mf. identity", \
     "dsb_c_ii_mfav.dat"  u ($1+binh):6 w l lt 2 lw 3 lc rgb nicered    t "plain mf. average", \
     "ds_c_ix.dat"   u ($1+binh):6 w l lt 1 lw 1  t "reference"

# mean force key for (b)
set size 0.18, 0.18
set origin 0.68, 0.59
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set format x "%g"
unset xlabel
set xtics 50 offset 0,.5 font tcsfont
set mxtics 5
unset logscale y
set format y "%g"
set ytics .2 offset .5,0 font tcsfont
set mytics 4
set ylabel "d ln{/Symbol r} / dE" offset 4., 0 font tcfont
unset key
plot [][-0.4:0.4] \
     "dsb_c_ii.dat"       u ($1+binh):9 w l lt 1 lw 1 lc rgb lightgreen t "single bin", \
     "dsb_c_ii_mfii.dat"  u ($1+binh):9 w l lt 1 lw 3 lc rgb niceblue   t "unbiased identity", \
     "dsb_c_ii_mfav.dat"  u ($1+binh):9 w l lt 2 lw 3 lc rgb nicered    t "plain average", \
     "ds_c_ix_mfii.dat"        u ($1+binh):9 w l lt 1 lw 1  t "reference", \
      0 w l lt 1 lw 1 lc rgb "#000000" notitle

#########################################################################################
# (c) canonical ensemble, variable window width

set size 0.5, 0.5
set origin 0.0, 0.0
set tmargin 1.0
set bmargin 1.5
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
unset xlabel
set logscale y
set ytics 10 font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(E)" offset 1, 0 font keyfont
set key font keyfont
plot [][1e-6:.15] \
     "dsb_c_ii.dat"     u 1:6 w l lt 2 lw 3 lc rgb niceblue  t "fixed window", \
     "dsb_c_ix.dat"     u 1:6 w l lt 1 lw 3 lc rgb nicegreen t "variable window", \
     "ds_c_ix.dat"      u 1:6 w l lt 1 lw 1  t "reference"

# mean force key for (c)
set size 0.18, 0.18
set origin 0.19, 0.09
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set format x "%g"
unset xlabel
set xtics 50 offset 0,.5 font tcsfont
set mxtics 5
unset logscale y
set ytics 50 offset .5,0 font tcsfont
set mytics 5
set format y "%g"
set ylabel "window boundaries" offset 4, 0 font tcfont
plot [][-1380:-1240] \
     "dsb_c_ix.dat" u 1:12:11 w filledcurves lt 1 lw 0 lc rgb lightgreen notitle, \
     "dsb_c_ix.dat" u 1:11 w l lt 1 lw 1 lc rgb nicegreen notitle, \
     "dsb_c_ix.dat" u 1:12 w l lt 1 lw 1 lc rgb nicegreen notitle, \
      x + 7 w l lt 2 lw 1 lc rgb niceblue notitle, \
      x - 7 w l lt 2 lw 1 lc rgb niceblue notitle, \
      x w l lt 1 lw 1 lc rgb "#000000" notitle


#########################################################################################
# microcanonical ensemble

set size 0.5, 0.5
set origin 0.5, 0.
set tmargin 1.0
set bmargin 1.5
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
unset xlabel
set logscale y
set ytics 10 font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(E)" offset 1, 0 font keyfont
set key font keyfont
plot [][1e-6:.15] \
     "ds_c_ix.dat"      u 1:6 w l lt 1 lw 3 lc rgb niceblue t "canonical", \
     "ds_m_ix.dat"      u 1:6 w l lt 2 lw 3 lc rgb nicered t "microcanonical"


# mean force key for (d)
set size 0.16, 0.16
set origin 0.69, 0.09
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set format x "%g"
unset xlabel
set xtics 50 offset 0,.5 font tcsfont
set mxtics 5
unset logscale y
set format y "%g"
set ytics .2 offset .5,0 font tcsfont
set mytics 4
set ylabel "d ln{/Symbol r} / dE" offset 4., 0 font tcfont
plot [][-0.4:0.5] \
     "ds_c_ix.dat" u 1:9 w l lt 1 lw 2 lc rgb niceblue notitle, \
     "ds_m_ix.dat" u 1:9 w l lt 2 lw 2 lc rgb nicered notitle, \
      0 w l lt 1 lw 1 lc rgb "#000000" notitle

unset multiplot
unset output
reset


set terminal pop

