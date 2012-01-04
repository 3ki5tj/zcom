# gray scale

reset
set terminal push
# large size 10,14
set terminal postscript landscape enhanced "Helvetica" 14 

set output "ene.ps"

binh = 0.05  # half bin width
tcfont = "Helvetica, 9"
tcsfont = "Helvetica, 7"
keyfont = "Helvetica, 10"
zcut = 1e-9
lightgray = "#cccccc"
lightgreen = "#99ffaa"
nicered = "#ff4040"
niceblue = "#40b0ff"
nicegreen = "#00cc00"

set xrange [-1380:-1240]

set multiplot

set label 1 "(a)" at screen 0.01, 0.98
set label 2 "(b)" at screen 0.51, 0.98
set label 3 "(c)" at screen 0.76, 0.98
set label 4 "(d)" at screen 0.01, 0.48
set label 5 "(e)" at screen 0.51, 0.48
#########################################################################################
# (a) canonical ensemble, basic comparison

set size 0.5, 0.5
set origin 0.0, 0.5
set tmargin 1.0
set bmargin 2.5
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
set xlabel "U" offset 0, 1 font keyfont
set ytics font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(U)" offset 1, 0 font keyfont
set key font keyfont
set logscale y
plot [][1e-7:1] \
     "dsb_c_ii.dat"     u ($1+binh):($2/1e3) w l lt 4 lc rgb lightgray t "histogram", \
     "dsb_c_ii.dat"     u 1:6 w l lt 1 lw 3 lc rgb niceblue t "fractional, fixed window", \
     "dsb_c_aj.dat"     u 1:($6 > zcut ? $6 : zcut) w l lt 2 lw 3 lc rgb nicered t "AJ identity", \
     "dsmerge.dat"      u 1:6 w l lt 1 lw 1  t "reference"

#########################################################################################
# (b)/(c) errors

# (b) KS difference (left)
set size 0.27, 0.5
set origin 0.5, 0.5
set tmargin 1.0
set bmargin 2.0
set lmargin 6.0
set rmargin 2.0
unset logscale
set xtics 10 offset 0, 0.4 font tcfont
set mxtics 5
set format x "%g"
set xlabel "window size, {/Symbol D}U" offset 0, 1 font keyfont 
set ytics 0.2 offset .5, 0 font tcsfont
set mytics 2
set ylabel "KS difference" offset 3, 0 font keyfont
set format y "%g"
set key font keyfont
#unset logscale y
plot [0:40][0:1.15] \
  "err.dat" u ($1*.2):2 w l lt 4 lw 3 lc rgb "#000000" t "histogram", \
  ""        u ($1*.2):5 w l lt 1 lw 3 lc rgb niceblue t "fractional", \
  ""        u ($1*.2):8 w l lt 2 lw 3 lc rgb nicered t "AJ"

# (c) entropic distance 
set size 0.23, 0.5
set origin 0.77, 0.5
set tmargin 1.0
set bmargin 2.0
set lmargin 3.0
set rmargin 2.0
unset logscale 
set xtics 10 offset 0, 0.4 font tcfont
set mxtics 5
set format x "%g"
set logscale y
set ytics 10 offset .5, 0 font tcsfont
set mytics 10
set ylabel "entropic distance" offset 3, 0 font keyfont
set format y "10^{%L}"
set key font keyfont
plot [0:40][1e-4:0.5] \
  "err.dat" u ($1*0.2):4  w l lt 4 lw 3 lc rgb "#000000" t "histogram", \
  ""        u ($1*0.2):7  w l lt 1 lw 3 lc rgb niceblue t "fractional", \
  ""        u ($1*0.2):10 w l lt 2 lw 3 lc rgb nicered t "AJ"
unset logscale y


#########################################################################################
# (d) canonical ensemble, improving mean force

set size 0.5, 0.5
set origin 0.0, 0.0
set tmargin 1.0
set bmargin 2.0
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
set xlabel "U" offset 0, 1 font keyfont
set ytics 100 font tcfont
set mytics 100
set format y "10^{%L}"
set ylabel "{/Symbol r}(U)" offset 2, 0 font keyfont
set key font keyfont
set logscale y
plot [-1410:-1200][1e-10:1e1] \
     "../datal/dsb_c.dat" u ($1+binh):($2/1000) w l lt 4 lw 1 lc rgb lightgray t "histograms", \
     "../data/dsb_c.dat"  u ($1+binh):($2/1000) w l lt 4 lw 1 lc rgb lightgray notitle, \
     "../datah/dsb_c.dat" u ($1+binh):($2/1000) w l lt 4 lw 1 lc rgb lightgray notitle, \
     "dsbmerge0.dat" u ($1):6 w l lt 2 lw 1 lc rgb nicered   t "weighted histogram", \
     "dsbmerge.dat"  u ($1):6 w l lt 1 lw 3 lc rgb niceblue  t "fractional identity", \
     "dsmerge.dat"   u ($1):6 w l lt 1 lw 1 lc rgb "#000000" t "reference"


#########################################################################################
# (e) comparison with microcanonical ensemble

set size 0.5, 0.5
set origin 0.5, 0.
set tmargin 1.0
set bmargin 2.0
set lmargin 6.0
set rmargin 2.0
set xtics 50 offset 0, 0.4 font tcfont
set mxtics 5
set xlabel "U" offset 0, 1 font keyfont
set logscale y
set ytics 10 font tcfont
set mytics 10
set format y "10^{%L}"
set ylabel "{/Symbol r}(U)" offset 1, 0 font keyfont
set key font keyfont
plot [][1e-7:.15] \
     "dsmerge.dat"      u 1:6 w l lt 1 lw 3 lc rgb niceblue t "canonical", \
     "ds_m_ix.dat"      u 1:6 w l lt 2 lw 3 lc rgb nicered t "microcanonical"


# lower inset: mean force for (d)
set size 0.16, 0.12
set origin 0.69, 0.09
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set format x "%g"
unset xlabel
set xrange [-1360:-1270]
set xtics 50 offset 0,.5 font tcsfont
set mxtics 5
unset logscale y
set format y "%g"
set ytics .2 offset .5,0 font tcsfont
set mytics 4
set ylabel "d ln{/Symbol r} / dE" offset 4., 0 font tcfont
plot [][-0.36:0.38] \
     "ds_c_ix.dat" u ($1+binh):9 w l lt 1 lw 2 lc rgb niceblue notitle, \
     "ds_m_ix.dat" u ($1+binh):9 w l lt 2 lw 2 lc rgb nicered notitle, \
     "ds_c_ix.dat" u ($1+binh):($9+1) w l lt 1 lw 2 lc rgb niceblue notitle, \
     "ds_m_ix.dat" u ($1+binh):($9+(762/2 - 1)/(-931.5-$1-binh)) w l lt 2 lw 2 lc rgb nicered notitle, \
      0 w l lt 1 lw 1 lc rgb "#000000" notitle

# upper inset: configurational temperature for (d)
set size 0.16, 0.08
set origin 0.69, 0.21
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set format x ""
unset xlabel
set xtics 50 offset 0,.5 font tcsfont
set mxtics 5
unset logscale y
set format y "%g"
set ytics .2 offset .5,0 font tcsfont
set mytics 4
set ylabel "{/Symbol \341}{/Symbol \321}\267v{/Symbol \361}" offset 3.0, -0.0 font tcfont
plot [][0.78:1.25] \
     "ds_c_ix.dat" u ($1+binh):($9+1) w l lt 1 lw 2 lc rgb niceblue notitle, \
     "ds_m_ix.dat" u ($1+binh):($9+(762/2 - 1)/(-931.5-$1-binh)) w l lt 2 lw 2 lc rgb nicered notitle, \
      1 w l lt 1 lw 1 lc rgb "#000000" notitle


unset multiplot
unset output
reset


set terminal pop

