% mpost this_file.mp

prologues := 3;

beginfig(1)
numeric u;
u = .4in;

pair org, i, j, k, l, jp, jpp, m, n, mb, nb, x, y, z;
org = (2*u, 2*u);
i = (2.2*u, -.3*u) shifted org; % slightly shrink from jp
j = (0, 0) shifted org;
k = (-2,-2)*u shifted org;
l = (-1,.5)*u shifted org;
jp = (.5,-2.5)*u shifted org;
jpp = (.8,2.2)*u shifted org; % slightly shrink from l
n = .5[j, l];
m = .5[i, k];
mb = m shifted (0*u, -1.5*u);
nb = n shifted (1.2*u, -0.9*u);
x = (3.5, 0)*u shifted org;
y = (0, 3)*u shifted org;
z = (-3., -3.)*u shifted org;


% the i, j, k plane
fill i--j--k--jp--cycle withcolor 0.7 white;

% the j, k, l plane
fill j--k--l--jpp--cycle withcolor 0.85 white;

% xij
%drawarrow (j--i);
dotlabel.rt(btex $i$ etex, i);
dotlabel.urt(btex $j$ etex, j);

% xkj
%drawarrow (j--k);
dotlabel.bot(btex $k$ etex, k);

% xkl
%drawarrow (l--k);
dotlabel.top(btex $l$ etex, l);

% m
drawarrow (m--mb) dashed evenly;
label.bot(btex $m$ etex, mb);

% n
drawarrow (n--nb);
label.rt(btex $n$ etex, nb);

drawarrow ((.9*u, -.15*u)..(.8*u, .3*u)..(.3*u, .8*u)) shifted org;
label.urt(btex $\phi$ etex, (.8*u, .3*u) shifted org);

% x
drawarrow (j--x);
label.bot(btex $x$ etex, (x+(-.2*u, -.1*u)));

% y
drawarrow (j--y);
label.rt(btex $y$ etex, (y+(.1*u, -.2*u)));

% z
drawarrow (j--z);
label.top(btex $z$ etex, (z+(.1*u, .2*u)));

endfig;
end.

