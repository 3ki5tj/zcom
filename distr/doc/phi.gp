unset multiplot
reset
set terminal push

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript portrait enhanced "Helvetica" 18

set output "phi.ps"

set multiplot
set size 1.0, 0.5
set origin 0.0, 0.5

set parametric

set lmargin 5
set bmargin 2
set label 12 "(a)" at -0.2, 1.4

af = .25
xaf = exp(-af)

set arrow 3 from 1,0 to 1, 1.1 lt 1 lw 2.0 nohead front 

set label 5 "{/Symbol r}(x*)" at 1.13, 0.6 front
set arrow 5 from 1.1, 0 to 1.1, 1.11 heads front 
set arrow 15 from 1.18, 1.11 to 1, 1.11 nohead front

set arrow 17 from 0,1.1 to 0.5, 0.5 lt 1 head front
set label 17 "{/Symbol \362@_{/Helvetica=8 x_{-}}^{/Helvetica=8 x_{+}}} {/=14 {/Symbol r}(x) dx}" \
  at -0.2, 1.2 front

set format x ""
set xtics add ("x_{-}" 0, "x*" 1, "x_{+}" 2)
unset ytics

f(x) = exp(-af*x**2 + 0.2*sin(1.5*pi*(x+1)/2)**2)

plot [t=0:1][-0.3:2.3][0:1.49] \
    2*t, f(2*t-1) notitle w filledcurves x1 lt 1 lc rgb "#E0E0E0", \
    2.4*t-0.2, f(2.4*t-1.2) t "distribution {/Symbol r}(x)" w l lt 1 lw 4 lc rgb "#808080"

set size 1.0,0.5
set origin 0.0,0.0


unset label
unset arrow
set tmargin .5
unset bmargin

bmin=0.0
bmax=2.0
bj=1.0
bjm=2.0-bj
am=0.6
ap=0.4

set xtics 1 format ""
set mxtics 5
unset ylabel
set ytics 0.5
set mytics 5

set label 11 "(b)" at -0.2, 1.1

sig(x, e, c) = (1+c)*x**e/(c + x**e)
der(x, e, c) = (1+c)*c*e*x**(e-1)/(c + x**e)**2
em = 2
cm = 5.0
ep = 3
cp = .4
sigm(x) = sig(x, em, cm)
sigp(x) = sig(x, ep, cp)
derm(x) = der(x, em, cm)
derp(x) = der(x, ep, cp)

set arrow 1 from bj+0.1, -ap to bj+0.1, am heads lw .5
set arrow 2 from bj, -ap to bj+.22, -ap nohead lw .5
set arrow 3 from bj,  am to bj+.22,  am nohead lw .5
set label 1 "x*" at bj-0.12,-0.08
set label 2 "x_{-}" at 0,-0.08
set label 3 "x_{+}" at 2,-0.08
set label 4 "1.0" at bj + 0.15, (am-ap)*.5 font "Helvetica, 16"

set label 5 "{/Symbol f}(x*^{-})" at bj - 0.3, am + .18
set arrow 5 from bj - 0.15, am + 0.1 to bj - 0.01, am lw .5
set label 6 "{/Symbol f}(x*^{+})" at bj + 0.05, -ap - .15
set arrow 6 from bj + 0.1, -ap - 0.1 to bj + 0.01, -ap lw .5

plot [t=0:1][-0.3:2.3][-0.8:1.2] \
   bj*t, am*sigm(t) t "{/Symbol f}({/Symbol b})" w l lt 1 lw 4, \
   2-bjm*t, -ap*sigp(t) notitle w l lt 1 lw 4, \
   bj, bj*am/bj*(1-t)-ap*t notitle w l lt 1 lw 4, \
   bj*t, am*derm(t) t "{/Symbol f}{/Symbol \242}({/Symbol b})" w l lt 2, \
   2-bjm*t, ap*derp(t) notitle w l lt 2, \
   5*t-1,0 notitle w l lt 1, \
   1,1*t notitle w l lt 1 


unset multiplot
unset output

#set terminal windows enhanced
set terminal pop
reset



