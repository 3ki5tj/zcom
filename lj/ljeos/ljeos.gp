# Lennard-Jones equation of states
# J. Karl Johnson et al.
# The Lennard-Jones equation of states revisited,
# Molecular Physics (1993) Vol. 78, No 3, 591-618
# equivalent to the lj_eos3d() in lj.0.c
x1  =  0.8623085097507421;
x2  =  2.976218765822098;
x3  = -8.402230115796038;
x4  =  0.1054136629203555;
x5  = -0.8564583828174598;
x6  =  1.582759470107601;
x7  =  0.7639421948305453;
x8  =  1.753173414312048;
x9  =  2.798291772190376e+03;
x10 = -4.8394220260857657e-02;
x11 =  0.9963265197721935;
x12 = -3.698000291272493e+01;
x13 =  2.084012299434647e+01;
x14 =  8.305402124717285e+01;
x15 = -9.574799715203068e+02;
x16 = -1.477746229234994e+02;
x17 =  6.398607852471505e+01;
x18 =  1.603993673294834e+01;
x19 =  6.805916615864377e+01;
x20 = -2.791293578795945e+03;
x21 = -6.245128304568454;
x22 = -8.116836104958410e+03;
x23 =  1.488735559561229e+01;
x24 = -1.059346754655084e+04;
x25 = -1.131607632802822e+02;
x26 = -8.867771540418822e+03;
x27 = -3.986982844450543e+01;
x28 = -4.689270299917261e+03;
x29 =  2.593535277438717e+02;
x30 = -2.694523589434903e+03;
x31 = -7.218487631550215e+02;
x32 =  1.721802063863269e+02;
gam = 3.0;

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
muex(rho, T) = Fex(rho, T) + P(rho, T)/rho - T;

# internal energy
U(rho, T) = rho*(c1(T) + rho*(c2(T)/2 + rho*(c3(T)/3 + rho*(c4(T)/4 + rho*(c5(T)/5 + rho*(c6(T)/6 + rho*(c7(T)/7 + rho*c8(T)/8))))))) + d1(T)*G1(rho) + d2(T)*G2(rho) + d3(T)*G3(rho) + d4(T)*G4(rho) + d5(T)*G5(rho) + d6(T)*G6(rho);

