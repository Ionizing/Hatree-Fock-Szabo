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
!     Zeta2     Slater orbital exponent (function 2)
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        if (IOP == 0) goto 20
          print 10, N, Za, Zb
10        format(1H1, 2X, 4HSTO-, T1, 21HG FOR ATOMIC NUMBERS ,F5.2, 5H AND F5.2,//)
20        continue
! Calculate all the one and two-electron integrals
        ! call Intgrl(IOP, N, R, ZETA1, ZETA2, Za, Zb)
! Be inefficient and put all integrals in pretty arrays
        ! call collect(IOP, N, R, ZETA1, ZETA2, Za, Zb)
! Perform the SCF calculation
        ! call SCF(IOP, N, R, ZETA1, ZETA2, Za, Zb)
        return
      end subroutine HFcalc

!***************************************************************
      subroutine Intgrl(IOP, N, R, ZETA1, ZETA2, Za, Zb)
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
            
            S12    = S12 + S(A1(i), A2(j), R2) * D1(i) * D2(j)
            T11    = T11 + T(A1(i), A1(j), 0.0D0) * D1(i) * D1(j)
            T12    = T12 + T(A1(i), A2(j), R2) * D1(i) * D2(j)
            T22    = T22 + T(A2(i), A2(j), 0,0D0) * D2(i) * D2(j)
          enddo ! do j
        enddo   ! do i
      end subroutine Intgrl
          



!***************************************************************
!***************************************************************
