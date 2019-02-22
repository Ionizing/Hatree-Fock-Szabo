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
10        format(1H1, 2X, 4HSTO-, T1, 21HG FOR ATOMIC NUMBERS ,F5.2, 5H AND $F5.2,//)
20        continue
      end subroutine HFcalc





!***************************************************************
!***************************************************************
