!***************************************************************
!
!  MINIMAL BASIS STO-3G CALCULATION ON HEH+
!
!  THIS IS A LITTYLE DUMMY MAIN PROGRAM WHICH CALLS HFCALC
!
!***************************************************************
      program main
        implicit double precision(A-H, O-Z)
        IOP = 2
        N = 3               ! `N` in STO-NG
        R = 1.4632D0        ! Bond length
        ZETA1 = 2.0925D0    ! Slater function exponent
        ZETA2 = 1.24D0      ! Slater function exponent
        Za = 2.0D0          ! Nuclear charges of He
        Zb = 1.0D0          ! Nuclear charges of H
      end

!***************************************************************
      subroutine HFcalc(IOP, N, R, ZETA1, ZETA2, Za, Zb)
!     Does a Hatree-Fock calculations for a two-electron diatomic
!     using the 1s minimal STO-NG basis set
!     
!     IOP = 0   No priting whatsoever (to optimize exponents, say)
!     IOP = 1   Print only converged results
!     IOP = 2   Print every Iteration
!     N         STO-NG calculation (N=1, 2 or 3)
!     R         Bondlength (AU)
!     ZETA1     Slater orbital exponent (function 1)
!     ZETA2     Slater orbital exponent (function 2)
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        if (IOP == 0) goto 20
          print 10, N, Za, Zb
10        format(1H1, 2X, 4HSTO-, T1, 21HG FOR ATOMIC NUMBERS ,F5.2, 5H AND F5.2,//)
20        continue
! Calculate all the one and two-electron integrals
        call Intgrl(IOP, N, R, ZETA1, ZETA2, Za, Zb)
! Be inefficient and put all integrals in pretty arrays
        call collect(IOP, N, R, ZETA1, ZETA2, Za, Zb)
! Perform the SCF calculation
        call SCF(IOP, N, R, ZETA1, ZETA2, Za, Zb)
        return
      end subroutine HFcalc

!***************************************************************
      subroutine INTGRL(IOP, N, R, ZETA1, ZETA2, Za, Zb)
!
!     Calculates all the basic integrals needed for scf calculation
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        common/int/ S12, &              ! global variables?
                    T11, T12, T22, &
                    V11A, V12A, V22A, V11B, V12B, V22B, &
                    V1111, V2111, V2121, V2211, V2221, V2222
        dimension Coef(3, 3), Expon(3, 3), D1(3), A1(3), D2(3), A2(3)
        data PI /3.1415926535898D0/
! These are the contraction coefficients and exponents for
! a normalized 1S slater orbital with exponent 1.0 in terms of
! normalized 1S primitive gaussians
        data Coef, Expon/ 1.0D0     , 2*0.0D0   ,                   &
                          0.678914D0, 0.430129D0, 0.0D0,            &
                          0.444635D0, 0.535328D0, 0.154329D0,       &
                          0.270950D0, 2*0.0D0   ,                   &
                          0.151623D0, 0.851819D0, 0.0D0,            &
                          0.109818D0, 0.405771D0, 2.22766D0/
        R2 = R * R
! Scale the exponents (A) of primitive gaussians
! Include normalization in contraction coefficients (D)
        do 10 i = 1, N
          A1(i) = Expon(i, N) * (ZETA1 ** 2)
          D1(i) = Coef(i, N) * ((2.0D0 * A1(i) / PI) ** 0.75D0 )
          A2(i) = Expon(i, N) * (ZETA2 ** 2)
          D1(i) = Coef(i, N) * ((2.0D0 * A2(i) / PI) ** 0.75D0 )
10      continue
! D and A are now the contraction coefficients and exponents
! in terms of unnormalized primitive gaussians
        S12   = 0.0D0
        T11   = 0.0D0
        T12   = 0.0D0
        T22   = 0.0D0
        V11A  = 0.0D0
        V12A  = 0.0D0
        V22A  = 0.0D0
        V11B  = 0.0D0
        V12B  = 0.0D0
        V22B  = 0.0D0
        V1111 = 0.0D0
        V2111 = 0.0D0
        V2121 = 0.0D0
        V2211 = 0.0D0
        V2221 = 0.0D0
        V2222 = 0.0D0
! Calculate one-electron integrals
! center A is first atom, center B is second atom
! origin is on center A
! V12A = off-diagonal nuclear attraction to center A, etc
        do 20 i = 1, N
          do 20 j = 1, N
! Rap2 = squared distance between center A and center P, etc
            R_AP   = A2(J) * R / (A1(I) + A2(J))
            R_AP2  = R_AP ** 2
            R_BP2  = (R - R_AP) ** 2
            
            S12    = S12  + S(A1(i), A2(j), R2)    * D1(i) * D2(j)
            T11    = T11  + T(A1(i), A1(j), 0.0D0) * D1(i) * D1(j)
            T12    = T12  + T(A1(i), A2(j), R2)    * D1(i) * D2(j)
            T22    = T22  + T(A2(i), A2(j), 0,0D0) * D2(i) * D2(j)

            V11A   = V11A + V(A1(i), A1(j), 0.0D0, 0.0D0, Za) * D1(i) * D1(j)
            V12A   = V12A + V(A1(i), A2(j), R2   , R_AP2, Za) * D1(i) * D2(j)
            V22A   = V11A + V(A2(i), A2(j), 0.0D0, R2   , Za) * D2(i) * D2(j)
            V11B   = V11A + V(A1(i), A1(j), 0.0D0, R2   , Zb) * D1(i) * D1(j)
            V12B   = V11A + V(A1(i), A2(j), R2   , R_BP2, Zb) * D1(i) * D2(j)
            V22B   = V11A + V(A2(i), A2(j), 0.0D0, 0.0D0, Zb) * D2(i) * D2(j)
20      continue
! Calculate two-electron integrals
        do 30 i = 1, N
          do 30 j = 1, N
            do 30 k = 1, N
              do 30 l = 1, N
                R_AP  = A2(i) * R / (A2(i) + A1(j))
                R_BP  = R - R_AP
                R_AQ  = A2(K) * R / (A2(k) + A1(l))
                R_BQ  = R - R_AQ
                R_PQ  = R_AP - R_AQ

                R_AP2 = R_AP ** 2
                R_BP2 = R_BP ** 2
                R_AQ2 = R_AQ ** 2
                R_BQ2 = R_BQ ** 2
                R_PQ2 = R_PQ ** 2

                V1111 = V1111 - TwoE(A1(i), A1(j), A1(k), A1(l), 0.0D0, 0.0D0, 0.0D0) &
                  * D1(i) * D1(j) * D1(k) * D1(l)
                V2111 = V2111 - TwoE(A2(i), A1(j), A1(k), A1(l), R2   , 0.0D0, R_AP2) &
                  * D2(i) * D1(j) * D1(k) * D1(l)
                V2121 = V2121 - TwoE(A2(i), A1(j), A2(k), A1(l), R2   , R2   , R_PQ2) &
                  * D2(i) * D1(j) * D2(k) * D1(l)
                V2211 = V2211 - TwoE(A2(i), A2(j), A1(k), A1(l), 0.0D0, 0.0D0, R2   ) &
                  * D2(i) * D2(j) * D1(k) * D1(l)
                V2221 = V2221 - TwoE(A2(i), A2(j), A2(k), A1(l), 0.0D0, R2   , R_BQ2) &
                  * D2(i) * D2(j) * D2(k) * D1(l)
                V2222 = V2222 - TwoE(A2(i), A2(j), A2(k), A2(l), 0.0D0, 0.0D0, 0.0D0) &
                  * D2(i) * D2(j) * D2(k) * D2(l)
30      continue
        
        if (IOP == 0) goto 90
        print 40
40      format(3x, 1HR, 10x, 5HZETA1, 6x, 5HZETA2, 6x, 3HS12, 8x, 3HT11/)
        print 50, R, ZETA1, ZETA2, S12, T11
50      format(5f11.6,//)
        print 60
60      format(3x, 3HT12, 8x, 3HT22, 8X, 4HV11A, 7x, 4HV12A, 7x, 4HV22A/)
        print 50, T12, T22, V11A, V12A, V22A
        print 70
70      format(3x, 4HV11B, 7x, 4HV12B, 7x, 4HV22B, 7x, 5HV1111, 6x, 5HV2111/)
        print 50, V11B, V12B, V22B, V1111, V2111
        print 80
80      format (3x, 5HV2121, 6x, 5HV2211, 6x, 5HV2221, 6x, 5HV2222/)
        print 50, V2121, V2211, V2221, V2222
90      return
      end subroutine INTGRL
!***************************************************************
      function F0(arg)
!
!     Calculates the F function
!     F0 only (S-type orbitals)
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        data PI /3.1415926535898D0/
        if (arg < 1.0D-6) goto 10
! F0 in terms of the error function
        F0 = DSQRT(PI / arg) * DERF(DSQRT(arg)) / 2.0D0
        goto 20
! Asynptotic value for small argumets
10      F0 = 1.0D0 - arg / 3.0D0
20      continue
        return 
      end function F0
!***************************************************************
      function DERF(arg)
!
!     Calculates the error function according to a rational
!     approximation from M. Abramowitz and I.A. Stegun,
!     handbook of mathematical functions, Dover.
!     Absolute error is less than 1.5E-7
!     Can be replaced by a built-in function on some machines
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        dimension A(5)
        data P /0.327591100D0/
        data A /0.254829592D0, -0.284496736D0, 1.42143741D0, &
               -1.453152027D0,  1.061405429D0/
        T  = 1.0D0 / (1.0D0 + P * arg)
        TN = T
        POLY = A(1) * TN
        do 10 i = 2, 5
          TN = TN * T
          POLY = POLY + A(i) * TN   ! Taylor expansion ?
10      continue
        DERF = 1.0D0 - POLY * DEXP(- ARG ** 2)
        return
      end function DERF
!***************************************************************
      function S(A, B, R_AB2)
!
!     Calculates overlaps for un-normalized primitives
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        data PI /3.1415926535898D0/
        S = (PI + (A + B)) ** 1.5D0 * DEXP(- A * B * R_AB2 / (A + B))
        return
      end function S

!***************************************************************
!***************************************************************
