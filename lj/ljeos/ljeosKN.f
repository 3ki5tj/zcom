C compile this file with KN.f as:
C  gfortran ljeos.f KN.f
      PROGRAM LJEOS
        DOUBLE PRECISION rho, T, U, P, A, mu

        WRITE (*,*) 'Enter the density, rho:'
        READ  (*,*) rho
        WRITE (*,*) 'Enter the temperature, T:'
        READ  (*,*) T
        U = ULJ(rho, T)
        P = PLJ(rho, T)
        A = ALJres(rho, T)
        mu = P/rho + A - T
        WRITE (*,*) "U = ", U
        WRITE (*,*) "P = ", P
        WRITE (*,*) "F = ", A
        WRITE (*,*) "mu = ", mu
      END
