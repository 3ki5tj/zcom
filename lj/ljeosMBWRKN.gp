# Lennard-Jones equation of states
# J. Kolafa et al.
# The Lennard-Jones fluid: An accurate analytic
# and theoretically-based equation of state
# Fluid Phase Equilibria (1994) Vol. 100, 1-34
# regressed from  data with T <= 6
# equivalent to the lj_eos3dK() in lj.0.c
x1  =     0.86230851;
x2  =     2.97621877;
x3  =    -8.40223012;
x4  =     0.10541366;
x5  =    -0.85645838;
x6  =     1.39945300;
x7  =    -0.20682219;
x8  =     2.66555449;
x9  =  1205.90355811;
x10 =     0.24414200;
x11 =     6.17927577;
x12 =   -41.33848427;
x13 =    15.14482295;
x14 =    88.90243729;
x15 = -2425.74868591;
x16 =  -148.52651854;
x17 =    68.73779789;
x18 =  2698.26346845;
x19 = -1216.87158315;
x20 = -1199.67930914;
x21 =    -7.28265251;
x22 = -4942.58001124;
x23 =    24.87520514;
x24 = -6246.96241113;
x25 =  -235.12327760;
x26 = -7241.61133138;
x27 =  -111.27706706;
x28 = -2800.52326352;
x29 =  1109.71518240;
x30 =  1455.47321956;
x31 = -2577.25311109;
x32 =   476.67051504;
gam =     4.52000000;

a1(T) = x1*T + x2*sqrt(T) + x3 + x4/T + x5/(T*T);
a2(T) = x6*T + x7 + x8/T + x9/(T*T);
a3(T) = x10*T + x11 + x12/T;
a4(T) = x13;
a5(T) = x14/T + x15/(T*T);
a6(T) = x16/T;
a7(T) = x17/T + x18/(T*T);
a8(T) = x19/(T*T);
b1(T) = (x20 + x21/T)/(T*T);
b2(T) = (x22 + x23/(T*T))/(T*T);
b3(T) = (x24 + x25/T)/(T*T);
b4(T) = (x26 + x27/(T*T))/(T*T);
b5(T) = (x28 + x29/T)/(T*T);
b6(T) = (x30 + x31/T + x32/(T*T))/(T*T);
c1(T) = x2*sqrt(T)/2 + x3 + 2*x4/T + 3*x5/(T*T);
c2(T) = x7 + 2*x8/T + 3*x9/(T*T);
c3(T) = x11 + 2*x12/T;
c4(T) = x13;
c5(T) = 2*x14/T + 3*x15/(T*T);
c6(T) = 2*x16/T;
c7(T) = 2*x17/T + 3*x18/(T*T);
c8(T) = 3*x19/(T*T);
d1(T) = (3*x20 + 4*x21/T)/(T*T);
d2(T) = (3*x22 + 5*x23/(T*T))/(T*T);
d3(T) = (3*x24 + 4*x25/T)/(T*T);
d4(T) = (3*x26 + 5*x27/(T*T))/(T*T);
d5(T) = (3*x28 + 4*x29/T)/(T*T);
d6(T) = (3*x30 + 4*x31/T + 5*x32/(T*T))/(T*T);
F(rho) = exp(-gam*rho*rho);
G1(rho) = (1 - F(rho))/(2*gam);
G2(rho) = -(F(rho)*rho**2 - 2*G1(rho))/(2*gam);
G3(rho) = -(F(rho)*rho**4 - 4*G2(rho))/(2*gam);
G4(rho) = -(F(rho)*rho**6 - 6*G3(rho))/(2*gam);
G5(rho) = -(F(rho)*rho**8 - 8*G4(rho))/(2*gam);
G6(rho) = -(F(rho)*rho**10 - 10*G5(rho))/(2*gam);

# the potential part of the free energy of the Lennard-Jones fluid
Fex(rho, T) = rho*(a1(T) + rho*(a2(T)/2 + rho*(a3(T)/3 + rho*(a4(T)/4 + rho*(a5(T)/5 + rho*(a6(T)/6 + rho*(a7(T)/7 + rho*a8(T)/8))))))) + b1(T)*G1(rho) + b2(T)*G2(rho) + b3(T)*G3(rho) + b4(T)*G4(rho) + b5(T)*G5(rho) + b6(T)*G6(rho);

# pressure
P(rho, T) = rho*T + rho*rho*(a1(T) + rho*(a2(T) + rho*(a3(T) + rho*(a4(T) + rho*(a5(T) + rho*(a6(T) + rho*(a7(T) + rho*a8(T)))))))) + F(rho)*rho*rho*rho*(b1(T) + rho*rho*(b2(T) + rho*rho*(b3(T) + rho*rho*(b4(T) + rho*rho*(b5(T) + rho*rho*b6(T))))));

# Gibbs free energy
muex(rho, T) = Fex(rho, T) + P(rho,T)/rho - T;

# internal energy
U(rho, T) = rho*(c1(T) + rho*(c2(T)/2 + rho*(c3(T)/3 + rho*(c4(T)/4 + rho*(c5(T)/5 + rho*(c6(T)/6 + rho*(c7(T)/7 + rho*c8(T)/8))))))) + d1(T)*G1(rho) + d2(T)*G2(rho) + d3(T)*G3(rho) + d4(T)*G4(rho) + d5(T)*G5(rho) + d6(T)*G6(rho);

