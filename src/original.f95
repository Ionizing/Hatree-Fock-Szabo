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
        call HFcalc(IOP, N, R, ZETA1, ZETA2, Za, Zb)
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
10        format(1H1, 2X, 4HSTO-, I1, 21HG FOR ATOMIC NUMBERS ,F5.2, 5H AND F5.2,//)
20        continue
! Calculate all the one and two-electron integrals
        call Intgrl(IOP, N, R, ZETA1, ZETA2, Za, Zb)
! Be inefficient and put all integrals in pretty arrays
        call Colect(IOP, N, R, ZETA1, ZETA2, Za, Zb)
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
            T22    = T22  + T(A2(i), A2(j), 0.0D0) * D2(i) * D2(j)

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
      function T(A, B, R_AB2)
!
!     Calculates kinetic energy integrals for un-normalized primitives
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        data PI /3.1415926535898D0/
        T = A * B / (A + B) * (3.0D0 - 2.0D0 * A * B * R_AB2 / (A + B)) * (PI / (A + B)) ** 1.5D0 &
          * DEXP(-A * B * R_AB2 / (A + B))
        return
      end function T
!***************************************************************
      function V(A, B, R_AB2, R_CP2, ZC)
!
!     Calculates un-normalized nuclear attraction integrals
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        data PI/3.1415926535898D0/
        V = 2.0D0 * PI / (A + B) * F0((A + B) * R_CP2) * DEXP(-A * B * R_AB2 / (A + B))
        V = -V * ZC
        return
      end function V
!***************************************************************
      function TwoE(A, B, C, D, R_AB2, R_CD2, R_PQ2)
! 
!     Calculates two-electron integrals for un-normalized primitives
!     A, B, C, D are the exponents alpha, beta, etc.
!     R_AB2 equals squared distance between center A and center B
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        data PI /3.1415926535898D0/
        TwoE = 2.0D0 * (PI ** 2.5D0) / ((A + B) * (C + D) * DSQRT(A + B + C + D)) &
          * F0((A + B) * (C + D) * R_PQ2 / (A + B + C + D)) &
          * DEXP(-A * B * R_AB2 / (A + B) - C * D * R_CD2 / (C + D) )
        return
      end function TwoE
!***************************************************************
      subroutine Colect(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
!  
!     This takes the basic integrals from common and assembles the
!     relevant matrices, that is S, H, X, XT, and two-electron integrals
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        common /Matrix/ S(2, 2), X(2, 2), XT(2, 2), H(2, 2), F(2, 2), G(2, 2), C(2, 2), &
          FPrime(2, 2), CPrime(2, 2), P(2, 2), OldP(2, 2), TT(2, 2, 2, 2), E(2, 2)
        common/int/ S12, &              ! global variables?
                    T11, T12, T22, &
                    V11A, V12A, V22A, V11B, V12B, V22B, &
                    V1111, V2111, V2121, V2211, V2221, V2222
! From core hamiltonian
        H(1, 1) = T11 + V11A + V11B
        H(1, 2) = T12 + V12A + V12B
        H(2, 1) = H(1, 2)             ! Hermitian operator
        H(2, 2) = T22 + V22A + V22B
! From overlap matrix
        S(1, 1) = 1.0D0
        S(1, 2) = S12
        S(2, 1) = S(1, 2)
        S(2, 2) = 1.0D0
! Use canonical orthogonalization
        X(1, 1) = 1.0D0 / DSQRT(2.0D0 * (1.0D0 + S12))
        X(2, 1) = X(1, 1)
        X(1, 2) = 1.0D0 / DSQRT(2.0D0 * (1.0D0 - S12))
        X(2, 2) = X(1, 2)
! Transpose of transformation matrix
        XT(1, 1) = X(1, 1)
        XT(1, 2) = X(2, 1)
        XT(2, 1) = X(1, 2)
        XT(2, 2) = X(2, 2)
! Matrix of two-electron integrals
        TT(1, 1, 1, 1) = V1111
        TT(2, 1, 1, 1) = V2111
        TT(1, 2, 1, 1) = V2111
        TT(1, 1, 2, 1) = V2111
        TT(1, 1, 1, 2) = V2111
        TT(2, 1, 2, 1) = V2121
        TT(1, 2, 2, 1) = V2121
        TT(2, 1, 1, 2) = V2121
        TT(1, 2, 1, 2) = V2121
        TT(2, 2, 1, 1) = V2211
        TT(1, 1, 2, 2) = V2211
        TT(2, 2, 2, 1) = V2221
        TT(2, 2, 1, 2) = V2221
        TT(2, 1, 2, 2) = V2221
        TT(1, 2, 2, 2) = V2221
        TT(2, 2, 2, 2) = V2222

        if (IOP == 0) goto 40
        call MatOut(S, 2, 2, 2, 2, 4HS    )
        call MatOut(X, 2, 2, 2, 2, 4HX    )
        call MatOut(H, 2, 2, 2, 2, 4HH    )
        print 10
10      format(//)
        do 30 i = 1, 2
          do 30 j = 1, 2
            do 30 k = 1, 2
              do 30 l = 1, 2
                print 20, i, j, k, l, TT(i, j, k, l)
20      format(3X, 1H(,4I2,2H ), F10.6)
30      continue
40      return
      end subroutine Colect
!***************************************************************
      subroutine SCF(IOP, N, R, ZETA1, ZETA2, Za, Zb)
!
!     Performs the SCF iterations
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        common /Matrix/ S(2, 2), X(2, 2), XT(2, 2), H(2, 2), F(2, 2), G(2, 2), C(2, 2), &
          FPrime(2, 2), CPrime(2, 2), P(2, 2), OldP(2, 2), TT(2, 2, 2, 2), E(2, 2)
        data PI /3.1415926535898D0/
! convergence criterion for density matrix
        data crit /1.0D-4/
! Maximum number of iterations
        data MaxIt /25/
! Iteration number
        iter = 0
! Use core-hamiltonian for initial guess at F. i.e. (P = 0)
        do 10 i = 1, 2
          do 10 j = 1, 2
10          P(i, j) = 0.0D0
        if (IOP < 2) goto 20
        call MatOut(P, 2, 2, 2, 2, 4HP    )
! start of iteration loop
20      iter = iter + 1
        if (IOP < 2) goto 40
        print 30, iter
30      format (/, 4x, 28HSTART OF ITERATION NUMBER = , I2)
40      continue
! From two-elextron part of fock matrix from P
        call FormG
        if (IOP < 2) goto 50
        call MatOut(G, 2, 2, 2, 2, 4HG    )
50      continue
! Add core hamiltonian to get fock matrix
        do 60 i = 1, 2
          do 60 j = 1, 2
            F(i, j) = H(i, j) + G(i, j)
60      continue
! Calculate electronic energy
        en = 0.0D0
        do 70 i = 1, 2
          do 70 j = 1, 2
            en = en + 0.5D0 * P(i, j) * (H(i, j) + F(i, j))
70      continue
        if (IOP < 2) goto 90
        call MatOut(F, 2, 2, 2, 2, 4HF    )
        print 80, en
80      format (///, 4x, 20HELECTRONIC ENERGY = ,D20.12)
90      continue
! Transform fock matrix using G for temporary storage
        call Mult(F, X, G, 2, 2)
        call Mult(XT, G, FPrime, 2, 2)
! Diagonalize transformed fock matrix
        call Diag(FPrime, CPrime, E)
! Transform eigenvectors to get matrix C
        call mult(X, CPrime, C, 2, 2)
! Form new density matrix
        do 100 i = 1, 2
          do 100 j = 1, 2
! save present density matrix
! before creating new one
            OldP(i, j) = P(i, j)
            P(i, j) = 0.0D0
            do 100 k = 1, 1
              P(i, j) = P(i, j) + 2.0D0 * C(i,k) * C(j, k)  ! eval product ?
100     continue
        if (IOP < 2) goto 110
        call MatOut(FPrime, 2, 2, 2, 2, 4HF'   )
        call MatOut(CPrime, 2, 2, 2, 2, 4HC'   )
        call MatOut(E, 2, 2, 2, 2, 4HE    )
        call MatOut(C, 2, 2, 2, 2, 4HC    )
        call MatOut(P, 2, 2, 2, 2, 4HP    )
110     continue
! Calculate delta
        Delta = 0.0D0
        do 120 i = 1, 2
          do 120 j = 1, 2
            Delta = Delta + (P(i, j) - OldP(i, j)) ** 2
120     continue
        Delta = DSQRT(Delta / 4.0D0)
        if (IOP == 0) goto 140
        print 130, Delta
130     format(/, 4x, 39HDELTA(CONVERGENCE OF DENSITY MATRIX) = &
          F10.6, /)
140     continue
! Check for convergence
        if (Delta < crit) goto 160
! Not yet converged
! Test for maximum nunber of iterations
! If maximum number not yet reached then
! go back for another iteration
        if (iter < MaxIt) goto 20
! Something wrong here
        print 150
150     format (4x, 21HNO CONVERGENCE IN SCF)
        stop
160     continue
! Calculation converged if it got here
! Add nuclear repulsion to get total energy
        ENT = en + Za * Zb / R
        if (IOP == 1) goto 180
        print 170, en, ENT
170     format(//, 4x, 21HCALCULATION CONVERGED, //, &
          4x, 20HELECTRONIC ENERGY = , D20.12, // &
          4x, 20HTOTAL ENERGY =      , D20.12)
180     continue
        if (IOP /= 1) goto 190
! Print out the final results if
! have not done so already
        call MatOut (G, 2, 2, 2, 2, 4HG    )
        call MatOut (F, 2, 2, 2, 2, 4HF    )
        call MatOut (E, 2, 2, 2, 2, 4HE    )
        call MatOut (C, 2, 2, 2, 2, 4HC    )
        call MatOut (P, 2, 2, 2, 2, 4HP    )
190     continue
! PS matrix has mulliken populations
        call mult(P, S, OldP, 2, 2)
        if (IOP == 0) goto 200
        call MatOut(OldP, 2, 2, 2, 2, 4HPS    )
200     continue
        return
      end subroutine SCF
!***************************************************************
      subroutine FormG
!
!     Calculates the G matrix from the density matrix
!     and two-electron integrals
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        common /Matrix/ S(2, 2), X(2, 2), XT(2, 2), H(2, 2), F(2, 2), G(2, 2), C(2, 2), &
          FPrime(2, 2), CPrime(2, 2), P(2, 2), OldP(2, 2), TT(2, 2, 2, 2), E(2, 2)
        do 10 i = 1, 2
          do 10 j = 1, 2
            G(i, j) = 0.0D0
            do 10 k = 1, 2
              do 10 l = 1, 2
                G(i, j) = G(i, j) + P(k, l) * (TT(i, j, k, l) - 0.5D0 * TT(i, l, k, j))
10      continue
        return
      end subroutine FormG
!***************************************************************
      subroutine Diag(F, C, E)
! Diagonalizes F to give eigenvectors in C and eigenvalues in E
! Theta is the angle describing solution
!***************************************************************
        implicit double precision(A-H, O-Z)
        dimension F(2, 2), C(2, 2), E(2, 2)
        data PI/3.1415926535898D0/
        if (DABS(F(1, 1) - F(2, 2)) > 1.0D-20) goto 10
! Here is symmetry determined solution (Homonuclear diatomic)
        Theta = PI / 4.0D0
        goto 20
10      continue
! Solution for heteronuclear diatomic
        Theta = 0.5D0 * DATAN(2.0D0 * F(1, 2) / (F(1, 1) - F(2, 2)))
20      continue
        C(1, 1) =  DCOS(Theta)
        C(2, 1) =  DSIN(Theta)
        C(1, 2) =  DSIN(Theta)
        C(2, 2) = -DCOS(Theta)
        E(1, 1) = F(1, 1) * DCOS(Theta) ** 2 + F(2, 2) * DSIN(Theta) ** 2 &
          + F(1, 2) * DSIN(2.0D0 * Theta)
        E(2, 2) = F(2, 2) * DCOS(Theta) ** 2 + F(1, 1) * DSIN(Theta) ** 2 &
          - F(1, 2) * DSIN(2.0D0 * Theta)
        E(2, 1) = 0.0D0
        E(1, 2) = 0.0D0
! Order eigenvalues and eigenvectors
        if (E(2, 2) > E(1, 1)) goto 30
        temp = E(2, 2)
        E(2, 2) = E(1, 1)
        E(1, 1) = temp
        temp = C(1, 2)
        C(1, 2) = C(1, 1)
        C(1, 1) = temp
        temp = C(2, 2)
        C(2, 2) = C(2, 1)
        C(2, 1) = temp
30      return
      end subroutine Diag
!***************************************************************
      subroutine mult(A, B, C, IM, M)
!
!     Multiplies two square matrices A and B to get C
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        dimension A(IM, IM), B(IM, IM), C(IM, IM)
        do 10 i = 1, M
          do 10 j = 1, M
            C(i, j) = 0.0D0
            do 10 K = 1, M
10      C(i, j) = C(i, j) + A(i, k) * B(k, j)
        return
      end subroutine mult

!***************************************************************
      subroutine MatOut(A, IM, IN, M, N, LABEL)
!
!     Print Matrices of site M by N
!
!***************************************************************
        implicit double precision(A-H, O-Z)
        dimension A(IM, IN)
        ihigh  = 0
10      low    = ihigh + 1
        ihigh  = ihigh + 8
        ihigh  = MIN0(ihigh, N)
        print 20, LABEL, (i, i=low, ihigh)
20      format(///, 3x, 5H THE , a4, 6H ARRAY, /, 15x, 5(10x, I3, 6X)//)
        do 30 i = 1, M
30      print 40, I, (A(i, j), j=low, ihigh)
40      format (I10, 5x, 5(1x, d18.10))
        if (N - ihigh) 50, 50, 10
50      return
      end subroutine MatOut
