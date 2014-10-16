c===================================================================
C      Package supplying the thermodynamic properties of the
C      LENNARD-JONES fluid
c
c      J. Kolafa, I. Nezbeda, Fluid Phase Equil. 100 (1994), 1
c
c      ALJ(rho, T)...Helmholtz free energy (including the ideal term)
c      PLJ(rho, T)...Pressure
c      ULJ(rho, T)...Internal energy
c===================================================================

C     Helmholtz free energy (including the ideal-gas term)
      DOUBLE PRECISION FUNCTION ALJ(rho, T)
        DOUBLE PRECISION rho, T, pfac, dC, BC, betaAHS, gammaBH, DALJ
        ALJ = dlog(rho)*T + betaAHS(pfac(rho, T))*T
     +        + rho*BC(T)/dexp(gammaBH()*rho**2)*T
     +        + DALJ(rho, T)
        RETURN
      END

C     Helmholtz free energy (without the ideal-gas term)
      DOUBLE PRECISION FUNCTION ALJres(rho, T)
        DOUBLE PRECISION rho, T, pfac, dC, BC, betaAHS, gammaBH, DALJ
        ALJres = betaAHS(pfac(rho, T))*T
     +           + rho*BC(T)*dexp(-gammaBH()*rho**2)*T
     +           + DALJ(rho, T)
        RETURN
      END

C pressure
      DOUBLE PRECISION FUNCTION PLJ(rho, T)
        DOUBLE PRECISION rho, T, pfac, dC, BC, zHS, gammaBH, sum
        sum = ((2.01546797d0*2+rho*(
     +        (-28.17881636d0)*3+rho*(
     +        28.28313847d0*4+rho*
     +        (-10.42402873d0)*5)))
     +        +((-19.58371655d0)*2+rho*(
     +        +75.62340289d0*3+rho*(
     +        (-120.70586598d0)*4+rho*(
     +        +93.92740328d0*5+rho*
     +        (-27.37737354d0)*6))))/dsqrt(T)
     +        + ((29.34470520d0*2+rho*(
     +        (-112.35356937d0)*3+rho*(
     +        +170.64908980d0*4+rho*(
     +        (-123.06669187d0)*5+rho*
     +        34.42288969d0*6))))+
     +        ((-13.37031968d0)*2+rho*(
     +        65.38059570d0*3+rho*(
     +        (-115.09233113d0)*4+rho*(
     +        88.91973082d0*5+rho*
     +        (-25.62099890d0)*6))))/T)/T)*rho**2
        PLJ = ((zHS(pfac(rho, T))
     +    + BC(T)*dexp(-gammaBH()*rho**2)
     +    *rho*(1-2*gammaBH()*rho**2))*T
     +    +sum )*rho
        RETURN
      END

C internal energy
      DOUBLE PRECISION FUNCTION ULJ(rho, T)
        DOUBLE PRECISION rho, T, dCdT, BCdT, dC, pfac, zHS, gammaBH
        ULJ = 3*(zHS(pfac(rho, T))-1)*dCdT(T)/dC(T)
     +        + rho*BCdT(T)*dexp(-gammaBH()*rho**2)
     +        + ((2.01546797d0+rho*(
     +        (-28.17881636d0)+rho*(
     +        +28.28313847d0+rho*
     +        (-10.42402873d0))))
     +        + (-19.58371655d0*1.5d0+rho*(
     +        75.62340289d0*1.5d0+rho*(
     +        (-120.70586598d0)*1.5d0+rho*(
     +        93.92740328d0*1.5d0+rho*
     +        (-27.37737354d0)*1.5d0))))/dsqrt(T)
     +        + ((29.34470520d0*2+rho*(
     +        -112.35356937d0*2+rho*(
     +         170.64908980d0*2+rho*(
     +        -123.06669187d0*2+rho*
     +        34.42288969d0*2)))) +
     +        (-13.37031968d0*3+rho*(
     +         65.38059570d0*3+rho*(
     +         -115.09233113d0*3+rho*(
     +        88.91973082d0*3+rho*
     +        (-25.62099890d0)*3))))/T)/T) *rho*rho
        RETURN
      END

C packing fraction
      DOUBLE PRECISION FUNCTION pfac(rho, T)
        DOUBLE PRECISION rho, T, dC
        DATA PI /3.141592653589793238d0/
        pfac = PI/6*rho*(dC(T))*3
        RETURN
      END

      DOUBLE PRECISION FUNCTION zHS(eta)
        DOUBLE PRECISION eta
        zHS = (1+eta*(1+eta*(1-eta/1.5*(1+eta)))) / (1-eta)**3
        RETURN
      END

      DOUBLE PRECISION FUNCTION betaAHS(eta)
        DOUBLE PRECISION eta
        betaAHS = dlog(1-eta)*5/3
     +            +eta*((4*eta - 33)*eta + 34)/6/(1-eta)**2
        RETURN
      END

C hBH diameter
      DOUBLE PRECISION FUNCTION dLJ(T)
        DOUBLE PRECISION T, isT
        isT = 1/dsqrt(T)
        dLJ = ((( 0.011117524191338 *isT-0.076383859168060)
     +        *isT)*isT+0.000693129033539)/isT+1.080142247540047
     +        +0.127841935018828*dlog(isT)
        RETURN
      END

      DOUBLE PRECISION FUNCTION dC(T)
        DOUBLE PRECISION T, sT
        sT = dsqrt(T)
        dC = -0.063920968d0*dlog(T)+0.011117524d0/T
     +       -0.076383859d0/sT+1.080142248d0+0.000693129d0*sT
        RETURN
      END

      DOUBLE PRECISION FUNCTION dCdT(T)
        DOUBLE PRECISION T, sT
        sT = dsqrt(T)
        dCdT = 0.063920968d0*T+0.011117524d0+(-0.5d0*0.076383859d0
     +         -0.5d0*0.000693129d0*T)*sT
        RETURN
      END

      DOUBLE PRECISION FUNCTION BC(T)
        DOUBLE PRECISION T, isT
        isT = 1/dsqrt(T)
        BC = (((((-0.58544978d0*isT+0.43102052d0)*isT
     +       +.87361369d0)*isT-4.13749995d0)*isT+2.90616279d0)*isT
     +       -7.02181962d0)/T+0.02459877d0
        RETURN
      END

      DOUBLE PRECISION FUNCTION BCdT(T)
        DOUBLE PRECISION T, iST
        isT = 1/dsqrt(T)
        BCdT = ((((-0.58544978d0*3.5*isT+0.43102052d0*3)*isT
     +    +0.87361369d0*2.5d0)*isT-4.13749995d0*2)*isT
     +    +2.90616279d0*1.5d0)*isT-7.02181962d0
        RETURN
      END

      DOUBLE PRECISION FUNCTION gammaBH()
        gammaBH = 1.92907278d0
        RETURN
      END

      DOUBLE PRECISION FUNCTION DALJ(rho, T)
        DOUBLE PRECISION rho, T
        DALJ = ((+2.01546797d0+rho*(-28.17881636d0
     +       +rho*(+28.28313847d0+rho*(-10.42402873d0))))
     +       +(-19.58371655d0+rho*(75.62340289d0+rho*((-120.70586598d0)
     +       +rho*(93.92740328d0+rho*(-27.37737354d0)))))/dsqrt(T)
     +       + ( (29.34470520d0+rho*((-112.35356937d0)
     +       +rho*(+170.64908980d0+rho*((-123.06669187d0)
     +       +rho*34.42288969d0))))
     +       +(-13.37031968d0+rho*(65.38059570d0+
     +       rho*((-115.09233113d0)+rho*(88.91973082d0
     +       +rho* (-25.62099890d0)))))/T)/T) *rho*rho
        RETURN
      END

