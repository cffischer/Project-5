                                                                  
SUBROUTINE SIXJ(I,J,K,L,M,N,ITIK,SI)
!*******************************************************************  
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | I/2  J/2  K/2 |                                            *
!     | L/2  M/2  N/2 |          (29.1A) [J.B.77]                  *
!                                                                  *
!   This restructured package uses the CONTAINS keyword to group   *
!   related subroutines together into self-contained blocks,       *
!   rather than relying on calls to subroutines in separate files. *
!   The ultimate aim is to reduce the size and complexity of the   *
!   "librang90" library.                                           *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Reststructured by A. Senchuk                    May      2019  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE vast_kind_param, ONLY: DOUBLE 
  USE CONS_C,          ONLY: ZERO, ONE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
  USE ixjtik_I 
  USE dracah_I 
  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER             :: I, J, K, L, M, N
  INTEGER, INTENT(IN) :: ITIK 
  REAL(DOUBLE)        :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: IFA 
  REAL(DOUBLE), DIMENSION(0:4,0:4,0:4,0:4,0:4,0:4) :: RACA 
  REAL(DOUBLE) :: UNDEF, A 
  LOGICAL :: SAVE 
!-----------------------------------------------
  DATA RACA/ 15625*1.D-20/  
  DATA UNDEF/ 1.D-20/  
  SI = ZERO 
  IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
     IF (IXJTIK(I,J,K,L,M,N) == 0) RETURN  
  ENDIF
  SAVE = .FALSE. 
  IF (MAX0(I,J,K,L,M,N) <= 4) THEN 
     SI = RACA(I,J,K,L,M,N) 
     IF (SI == UNDEF) THEN 
        SAVE = .TRUE. 
     ELSE 
        RETURN  
     ENDIF
  ENDIF
!
!     CALCULATED IN CASE WHEN ONE OF PERAMETERS EQUAL 0
!
  IF (I*J*K*L*M*N == 0) THEN 
     IF (I == 0) THEN 
        A = DBLE((M + 1)*(K + 1)) 
        IFA = L + M + K 
     ELSE IF (J == 0) THEN 
        A = DBLE((L + 1)*(K + 1)) 
        IFA = I + M + N 
     ELSE IF (K == 0) THEN 
        A = DBLE((I + 1)*(L + 1)) 
        IFA = I + M + N 
     ELSE IF (L == 0) THEN 
        A = DBLE((J + 1)*(K + 1)) 
        IFA = I + J + K 
     ELSE IF (M == 0) THEN 
        A = DBLE((I + 1)*(K + 1)) 
        IFA = I + J + K 
     ELSE 
        A = DBLE((I + 1)*(J + 1)) 
        IFA = I + J + K 
     ENDIF
     SI = ONE/DSQRT(A) 
     IF (MOD(IFA,4) /= 0) SI = -SI 
!
!     THE CASE 1/2
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 1) THEN 
     IF (I == 1) THEN 
        CALL SIXJ5 (M, K, L, J, N, 0, SI) 
     ELSE IF (J == 1) THEN 
        CALL SIXJ5 (I, N, M, L, K, 0, SI) 
     ELSE IF (K == 1) THEN 
        CALL SIXJ5 (I, M, N, L, J, 0, SI) 
     ELSE IF (L == 1) THEN 
        CALL SIXJ5 (J, K, I, M, N, 0, SI) 
     ELSE IF (M == 1) THEN 
        CALL SIXJ5 (I, K, J, L, N, 0, SI) 
     ELSE 
        CALL SIXJ5 (I, J, K, L, M, 0, SI) 
     ENDIF
!
!     THE CASE 1
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 2) THEN 
     IF (I == 2) THEN 
        CALL SIXJ1 (M, K, L, J, N, 0, SI) 
     ELSE IF (J == 2) THEN 
        CALL SIXJ1 (I, N, M, L, K, 0, SI) 
     ELSE IF (K == 2) THEN 
        CALL SIXJ1 (I, M, N, L, J, 0, SI) 
     ELSE IF (L == 2) THEN 
        CALL SIXJ1 (J, K, I, M, N, 0, SI) 
     ELSE IF (M == 2) THEN 
        CALL SIXJ1 (I, K, J, L, N, 0, SI) 
     ELSE 
        CALL SIXJ1 (I, J, K, L, M, 0, SI) 
     ENDIF
!
!     THE CASE 3/2
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 3) THEN 
     IF (I == 3) THEN 
        CALL SIXJ35 (M, K, L, J, N, 0, SI) 
     ELSE IF (J == 3) THEN 
        CALL SIXJ35 (I, N, M, L, K, 0, SI) 
     ELSE IF (K == 3) THEN 
        CALL SIXJ35 (I, M, N, L, J, 0, SI) 
     ELSE IF (L == 3) THEN 
        CALL SIXJ35 (J, K, I, M, N, 0, SI) 
     ELSE IF (M == 3) THEN 
        CALL SIXJ35 (I, K, J, L, N, 0, SI) 
     ELSE 
        CALL SIXJ35 (I, J, K, L, M, 0, SI) 
     ENDIF
!
!     THE CASE 2
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 4) THEN 
     IF (I == 4) THEN 
        CALL SIXJ2 (M, K, L, J, N, 0, SI) 
     ELSE IF (J == 4) THEN 
        CALL SIXJ2 (I, N, M, L, K, 0, SI) 
     ELSE IF (K == 4) THEN 
        CALL SIXJ2 (I, M, N, L, J, 0, SI) 
     ELSE IF (L == 4) THEN 
        CALL SIXJ2 (J, K, I, M, N, 0, SI) 
     ELSE IF (M == 4) THEN 
        CALL SIXJ2 (I, K, J, L, N, 0, SI) 
     ELSE 
        CALL SIXJ2 (I, J, K, L, M, 0, SI) 
     ENDIF
!
!     THE CASE 5/2
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 5) THEN 
     CALL DRACAH (I, J, M, L, K, N, SI) 
     IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
!
!     CASES 3
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 6) THEN 
     IF (I == 6) THEN 
        CALL SIXJ3 (M, K, L, J, N, 0, SI) 
     ELSE IF (J == 6) THEN 
        CALL SIXJ3 (I, N, M, L, K, 0, SI) 
     ELSE IF (K == 6) THEN 
        CALL SIXJ3 (I, M, N, L, J, 0, SI) 
     ELSE IF (L == 6) THEN 
        CALL SIXJ3 (J, K, I, M, N, 0, SI) 
     ELSE IF (M == 6) THEN 
        CALL SIXJ3 (I, K, J, L, N, 0, SI) 
     ELSE 
        CALL SIXJ3 (I, J, K, L, M, 0, SI) 
     ENDIF
!
!     THE CASE 7/2
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 7) THEN 
     CALL DRACAH (I, J, M, L, K, N, SI) 
     IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
!
!     CASES 4
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 8) THEN 
     IF (I == 8) THEN 
        CALL SIXJ4 (M, K, L, J, N, 0, SI) 
     ELSE IF (J == 8) THEN 
        CALL SIXJ4 (I, N, M, L, K, 0, SI) 
     ELSE IF (K == 8) THEN 
        CALL SIXJ4 (I, M, N, L, J, 0, SI) 
     ELSE IF (L == 8) THEN 
        CALL SIXJ4 (J, K, I, M, N, 0, SI) 
     ELSE IF (M == 8) THEN 
        CALL SIXJ4 (I, K, J, L, N, 0, SI) 
     ELSE 
        CALL SIXJ4 (I, J, K, L, M, 0, SI) 
     ENDIF
!
!     THE CASE 9/2
!
  ELSE IF (MIN0(I,J,K,L,M,N) == 9) THEN 
     CALL DRACAH (I, J, M, L, K, N, SI) 
     IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
!
!     CALCULATED OTHER CASES
!
  ELSE 
     CALL DRACAH(I,J,M,L,K,N,SI)
     IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
  ENDIF
  IF (SAVE) RACA(I,J,K,L,M,N) = SI 
  RETURN  

CONTAINS 

  SUBROUTINE SIXJ5(J,K,L,M,N,ITIK,SI) 
!*******************************************************************
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2  1/2 |             [B.M.X. 75]                    *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             March    1995  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructed by C. Froese Fischer                 May      2019  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER                   :: J, K, L, M, N
    INTEGER, INTENT(IN)       :: ITIK 
    REAL(DOUBLE), INTENT(OUT) :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: I1 
    REAL(DOUBLE) :: AS, A, B, C, AKA 
!-----------------------------------------------
    SI = ZERO 
    IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
       IF (IXJTIK(J,K,L,M,N,1) == 0) RETURN  
    ENDIF
    I1 = (J + K + L)/2 
    AS = DBLE(I1) 
    A = DBLE(L) 
    B = DBLE(K) 
    C = DBLE(J) 
    AKA = ONE 
    IF (MOD(I1,2) /= 0) AKA = -AKA 
    IF (K < M) THEN 
       IF (J < N) THEN
!              M > K,  J < N.
          SI = -AKA*DSQRT((AS + TWO)*(AS - A + ONE)/((B + ONE)*(B + TWO)*(C    &
               + ONE)*(C + TWO))) 
       ELSE IF (J > N) THEN 
!              M > K,  J > N.
          SI = AKA*DSQRT((AS - C + ONE)*(AS - B)/((B + ONE)*(B + TWO)*C*(C +   &
               ONE))) 
       ENDIF
    ELSE IF (K > M) THEN 
       IF (J < N) THEN 
!             M < K,  J < N.
          SI = AKA*DSQRT((AS - C)*(AS - B + ONE)/(B*(B + ONE)*(C + ONE)*(C +   &
               TWO))) 
       ELSE IF (J > N) THEN 
!             M < K,  J > N.
          SI = AKA*DSQRT((AS + ONE)*(AS - A)/(B*(B + ONE)*C*(C + ONE))) 
       ENDIF
    ENDIF
    RETURN  
  END SUBROUTINE SIXJ5
!
!##############################################################################
!      
  SUBROUTINE SIXJ1(I,J,K,L,M,ITIK,SI) 
!*******************************************************************
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | I/2  J/2  K/2 |                                            *
!     | L/2  M/2   1  |               [B.M.X.  75].                *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             March    1995  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructed by A. Senchuk                        September 2019 *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER                   :: I, J, K, L, M
    INTEGER, INTENT(IN)       :: ITIK 
    REAL(DOUBLE), INTENT(OUT) :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IFA 
    REAL(DOUBLE) :: AS, AKA, A, B, C 
!-----------------------------------------------
    SI = ZERO 
    IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
       IF (IXJTIK(I,J,K,L,M,2) == 0) RETURN  
    ENDIF
    IFA = (I + J + K)/2 
    AS = DBLE(IFA) 
    AKA = ONE 
    IF (MOD(IFA,2) /= 0) AKA = -AKA 
    A = DBLE(K) 
    B = DBLE(J) 
    C = DBLE(I) 
    IF (I < M) THEN 
       IF (J < L) THEN 
!              M > I,   L > J.
          SI = AKA*DSQRT((AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A + TWO    &
               )/((B + ONE)*(B + TWO)*(B + THREE)*(C + ONE)*(C + TWO)*(C +       &
               THREE))) 
       ELSE IF (J == L) THEN 
!              M > I,  L = J.
          SI = (-AKA)*DSQRT(TWO*(AS + TWO)*(AS - C)*(AS - B + ONE)*(AS - A +     &
               ONE)/(B*(B + ONE)*(B + TWO)*(C + ONE)*(C + TWO)*(C + THREE))) 
       ELSE 
!              M > I,  L < J.
          SI = AKA*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO    &
               )/((B - ONE)*B*(B + ONE)*(C + ONE)*(C + TWO)*(C + THREE))) 
       ENDIF
    ELSE IF (I == M) THEN 
       IF (J < L) THEN 
!              M = L,  L > J.
          SI = (-AKA)*DSQRT((AS + TWO)*(AS - C + ONE)*(AS - B)*(AS - A + ONE)    &
               *TWO/((B + ONE)*(B + TWO)*(B + THREE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (J == L) THEN 
!              M = I,  L = J.
          SI = (-AKA)*((B*B + C*C - A*A)*HALF + B + C - A)/DSQRT(B*(B + ONE)*    &
               (B + TWO)*C*(C + ONE)*(C + TWO)) 
       ELSE 
!              M = I,  L < J.
          SI = AKA*DSQRT((AS + ONE)*(AS - C)*(AS - B + ONE)*(AS - A)*TWO/((B     &
               - ONE)*B*(B + ONE)*C*(C + ONE)*(C + TWO))) 
       ENDIF
    ELSE 
       IF (J < L) THEN 
!              M < I,   L > J.
          SI = AKA*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - B - ONE)*(AS - B    &
               )/((B + ONE)*(B + TWO)*(B + THREE)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (J == L) THEN 
!              M < I,   L = J.
          SI = AKA*DSQRT((AS + ONE)*(AS - C + ONE)*(AS - B)*(AS - A)*TWO/(B*(    &
               B + ONE)*(B + TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE 
!              M < I,   L < J.
          SI = AKA*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B - ONE)*B*(    &
               B + ONE)*(C - ONE)*C*(C + ONE))) 
       ENDIF
    ENDIF
    RETURN  
  END SUBROUTINE SIXJ1
!
!##############################################################################
!      
  SUBROUTINE SIXJ35(J,K,L,M,N,ITIK,SI) 
!*******************************************************************
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2  3/2 |             [B.M.X. 75]                    *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  * 
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!   Restructed by C. Froese Fischer                 May      2019  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER                   :: J, K, L, M, N
    INTEGER, INTENT(IN)       :: ITIK 
    REAL(DOUBLE), INTENT(OUT) :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: I1 
    REAL(DOUBLE) :: AS, A, B, C, AKA 
!-----------------------------------------------
    SI = ZERO 
    IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
       IF (IXJTIK(J,K,L,M,N,3) == 0) RETURN  
    ENDIF
    I1 = (J + K + L)/2 
    AS = DBLE(I1) 
    A = DBLE(L) 
    B = DBLE(J) 
    C = DBLE(K) 
    AKA = ONE 
    IF (MOD(I1,2) /= 0) AKA = -AKA 
    IF (J - N == 3) THEN 
! -3
       IF (K - M == 3) THEN 
!  I                      -3/2  -3/2
          SI = AKA*DSQRT((AS - ONE)*AS*(AS + ONE)*(AS - A - TWO)*(AS - A -       &
               ONE)*(AS - A)/((B - TWO)*(B - ONE)*B*(B + ONE)*(C - TWO)*(C -     &
               ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 3) THEN 
!  IV  P(12)              3/2   -3/2
          SI = AKA*DSQRT((AS - C - TWO)*(AS - C - ONE)*(AS - C)*(AS - B + ONE    &
               )*(AS - B + TWO)*(AS - B + THREE)/((C + 1)*(C + TWO)*(C + THREE)  &
               *(C + FOUR)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
       ELSE IF (K - M == 1) THEN 
!  II  P(12)             -1/2   -3/2
          SI = AKA*DSQRT(THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)*(AS - C)    &
               *(AS - B + ONE)/((C - ONE)*C*(C + ONE)*(C + TWO)*(B - TWO)*(B -   &
               ONE)*B*(B + ONE))) 
       ELSE IF (M - K == 1) THEN 
!  III P(12)              1/2   -3/2
          SI = AKA*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - C - ONE)*(AS - C)*(     &
               AS - B + ONE)*(AS - B + TWO)/(C*(C + ONE)*(C + TWO)*(C + THREE)*  &
               (B - TWO)*(B - ONE)*B*(B + ONE))) 
       ENDIF
    ELSE IF (N - J == 3) THEN 
!  3
       IF (K - M == 3) THEN 
!  IV                     -3/2   3/2
          SI = AKA*DSQRT((AS - B - TWO)*(AS - B - ONE)*(AS - B)*(AS - C + ONE    &
               )*(AS - C + TWO)*(AS - C + THREE)/((B + ONE)*(B + TWO)*(B +       &
               THREE)*(B + FOUR)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 3) THEN 
!  2       pataisyta               3/2   3/2
          SI = -AKA*DSQRT((AS + TWO)*(AS + THREE)*(AS + FOUR)*(AS - A + ONE)*    &
               (AS - A + TWO)*(AS - A + THREE)/((B + ONE)*(B + TWO)*(B + THREE)  &
               *(B + FOUR)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
       ELSE IF (K - M == 1) THEN 
!  1   P(12)   pataisytas          -1/2    3/2
          SI = -AKA*DSQRT(THREE*(AS + TWO)*(AS - A + ONE)*(AS - C + ONE)*(AS     &
               - C + TWO)*(AS - B - ONE)*(AS - B)/((C - ONE)*C*(C + ONE)*(C +    &
               TWO)*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))) 
       ELSE IF (M - K == 1) THEN 
!  3  P(12)     taisyta           1/2    3/2
          SI = AKA*DSQRT(THREE*(AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A    &
               + TWO)*(AS - B)*(AS - C + ONE)/(C*(C + ONE)*(C + TWO)*(C +        &
               THREE)*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))) 
       ENDIF
! -1
    ELSE IF (J - N == 1) THEN 
       IF (K - M == 3) THEN 
!  II                   -3/2   -1/2
          SI = AKA*DSQRT((THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)*(AS - B    &
               )*(AS - C + ONE))/((B - ONE)*B*(B + ONE)*(B + TWO)*(C - TWO)*(C   &
               - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 3) THEN 
!  1                     3/2   -1/2
          SI = -AKA*DSQRT(THREE*(AS + TWO)*(AS - A + ONE)*(AS - B + ONE)*(AS     &
               - B + TWO)*(AS - C - ONE)*(AS - C)/((B - ONE)*B*(B + ONE)*(B +    &
               TWO)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
       ELSE IF (K - M == 1) THEN 
!  V                    -1/2   -1/2
          SI = AKA*(TWO*(AS - B)*(AS - C) - (AS + TWO)*(AS - A - ONE))*DSQRT(    &
               (AS + ONE)*(AS - A)/((B - ONE)*B*(B + ONE)*(B + TWO)*(C - ONE)*C  &
               *(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 1) THEN 
!  VI P(12)              1/2   -1/2
          SI = AKA*((AS - B + TWO)*(AS - C + ONE) - TWO*(AS - A + ONE)*(AS +     &
               ONE))*DSQRT((AS - C)*(AS - B + ONE)/(C*(C + ONE)*(C + TWO)*(C +   &
               THREE)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
       ENDIF
! 1
    ELSE IF (N - J == 1) THEN 
       IF (K - M == 3) THEN 
!  III                  -3/2    1/2
          SI = AKA*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - B - ONE)*(AS - B)*(     &
               AS - C + ONE)*(AS - C + TWO)/(B*(B + ONE)*(B + TWO)*(B + THREE)*  &
               (C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 3) THEN 
!  3              pataisyta       3/2    1/2
          SI = AKA*DSQRT(THREE*(AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A    & 
               + TWO)*(AS - B + ONE)*(AS - C)/(B*(B + ONE)*(B + TWO)*(B +        &
               THREE)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
       ELSE IF (K - M == 1) THEN 
!  VI                   -1/2    1/2
          SI = AKA*((AS - C + TWO)*(AS - B + ONE) - TWO*(AS - A + ONE)*(AS +     &
               ONE))*DSQRT((AS - B)*(AS - C + ONE)/(B*(B + ONE)*(B + TWO)*(B +   &
               THREE)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 1) THEN 
!  4      pataisyta               1/2    1/2
          SI = -AKA*(TWO*(AS - B)*(AS - C) - (AS + THREE)*(AS - A))*DSQRT((AS    &
               + TWO)*(AS - A + ONE)/(B*(B + ONE)*(B + TWO)*(B + THREE)*C*(C     &
               + ONE)*(C + TWO)*(C + THREE))) 
       ENDIF
    ENDIF
    RETURN  
  END SUBROUTINE SIXJ35
!
!################################################################################
!      
  SUBROUTINE SIXJ2(J,K,L,M,N,ITIK,SI) 
!*******************************************************************
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2   2  |             [B.M.X. 75]                    *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructed by A. Senchuk                        September 2019 *
!                                                                  *
!*******************************************************************
      
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER                   :: J, K, L, M, N
    INTEGER, INTENT(IN)       :: ITIK 
    REAL(DOUBLE), INTENT(OUT) :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: I1 
    REAL(DOUBLE) :: AS, A, B, C, AKA, X1, X2, X3 
!-----------------------------------------------
    SI = ZERO 
    IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
       IF (IXJTIK(J,K,L,M,N,4) == 0) RETURN  
    ENDIF
    I1 = (J + K + L)/2 
    AS = DBLE(I1) 
    A = DBLE(L) 
    B = DBLE(J) 
    C = DBLE(K) 
    AKA = ONE 
    IF (MOD(I1,2) /= 0) AKA = -AKA 
    IF (J - N == 4) THEN 
! -2
       IF (K - M == 4) THEN 
!  I                      -2  -2
          SI = AKA*DSQRT((AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - THREE)*(B      &
               - TWO)*(B - ONE)*B*(B + ONE))) 
          SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS -      &
               A)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 4) THEN 
!  V   P(12)               2  -2
          SI = AKA*DSQRT((AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS       &
               - C)/((C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)*(C + TWO +        &
               THREE))) 
          SI = SI*DSQRT((AS - B + ONE)*(AS - B + TWO)*(AS - B + THREE)*(AS -      &
               B + FOUR)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
       ELSE IF (K - M == 2) THEN 
!  II   P(12)              -1   -2
          SI = AKA*TWO*DSQRT((AS - ONE)*AS*(AS + ONE)/((C - TWO)*(C - ONE)*C*     &
               (C + ONE)*(C + TWO))) 
          SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)*(AS - C)*(AS       &
               - B + ONE)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
       ELSE IF (M - K == 2) THEN 
!  IV   P(12)               1   -2
          SI = AKA*TWO*DSQRT((AS + ONE)*(AS - A)*(AS - C - TWO)*(AS - C - ONE     &
               )*(AS - C)/(C*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
          SI = SI*DSQRT((AS - B + ONE)*(AS - B + TWO)*(AS - B + THREE)/((B -      &
               THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
       ELSE IF (K - M == 0) THEN 
!  III  P(12)               0   -2
          SI = AKA*DSQRT(TWO*THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((C      &
               - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
          SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO)     &
               /((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
       ENDIF
!  2
    ELSE IF (N - J == 4) THEN 
       IF (K - M == 4) THEN 
!  V                      -2   2
          SI = AKA*DSQRT((AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS       &
               - B)/((B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)*(B + TWO +        &
               THREE))) 
          SI = SI*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - C + THREE)*(AS -      &
               C + FOUR)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 4) THEN 
!  1                       2   2
          SI = AKA*DSQRT((AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS      &
               - A + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*    &
               (B + ONE))) 
          SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO     &
               )/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))   &
               ) 
       ELSE IF (K - M == 2) THEN 
!  3                      -1   2
          SI = -AKA*DSQRT((AS - A + ONE)*(AS + TWO)*(AS - B - TWO)*(AS - B -      &
               ONE)*(AS - B)/((B + TWO + THREE)*(B + FOUR)*(B + THREE)*(B + TWO   &
               )*(B + ONE))) 
          SI = SI*TWO*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((     &
               C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 2) THEN 
!  2                       1   2
          SI = -AKA*DSQRT((AS - B)*(AS - C + ONE)*(AS - A + THREE)*(AS - A +      &
               TWO)*(AS - A + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B   &
               + TWO)*(B + ONE))) 
          SI = SI*TWO*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + FOUR)*(     &
               C + THREE)*(C + TWO)*(C + ONE)*C)) 
       ELSE IF (K - M == 0) THEN 
!  5                        0   2
          SI = AKA*DSQRT(THREE*TWO*(AS - B)*(AS - B - ONE)*(AS - C + TWO)*(AS     &
               - C + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*    &
               (B + ONE))) 
          SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)     &
               /((C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE))) 
       ENDIF
    ELSE IF (J - N == 2) THEN 
! -1
       IF (K - M == 4) THEN 
!  II   P(12)              -2  -1
          SI = AKA*TWO*DSQRT((AS - ONE)*AS*(AS + ONE)/((B - TWO)*(B - ONE)*B*     &
               (B + ONE)*(B + TWO))) 
          SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)*(AS - B)*(AS       &
               - C + ONE)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 4) THEN 
!  3   P(12)                2  -1
          SI = -AKA*DSQRT((AS - A + ONE)*(AS + TWO)*(AS - C - TWO)*(AS - C -      &
               ONE)*(AS - C)/((C + TWO + THREE)*(C + FOUR)*(C + THREE)*(C + TWO   &
               )*(C + ONE))) 
          SI = SI*TWO*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((     &
               B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
       ELSE IF (K - M == 2) THEN 
!  VI                       -1  -1
          SI = AKA*((A + B)*(A - B + TWO) - (C - TWO)*(C - B + TWO))/DSQRT((B     &
               - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)) 
          SI = SI*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((C - TWO)*(C       &
               - ONE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 2) THEN 
!  VIII P(12)              1  -1
          SI = AKA*((A + C + FOUR)*(A - C - TWO) - (B - TWO)*(B + C + FOUR))/     &
               DSQRT(C*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)) 
          SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO)     &
               /((B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
       ELSE IF (K - M == 0) THEN
!  VII  P(12)              0  -1
          SI = AKA*HALF*((A + C + TWO)*(A - C) - B*B + FOUR)/DSQRT((C - ONE)*     &
               C*(C + ONE)*(C + TWO)*(C + THREE)) 
          SI = SI*DSQRT(THREE*TWO*(AS + ONE)*(AS - A)*(AS - C)*(AS - B + ONE)     &
               /((B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
       ENDIF
    ELSE IF (N - J == 2) THEN 
!  1
       IF (K - M == 4) THEN 
!  IV                     -2   1
          SI = AKA*TWO*DSQRT((AS + ONE)*(AS - A)*(AS - B - TWO)*(AS - B - ONE     &
               )*(AS - B)/(B*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))) 
          SI = SI*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - C + THREE)/((C -      &
               THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 4) THEN 
!  2 P(12)                 2   1
          SI = -AKA*DSQRT((AS - C)*(AS - B + ONE)*(AS - A + THREE)*(AS - A +      &
               TWO)*(AS - A + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C   &
               + TWO)*(C + ONE))) 
          SI = SI*TWO*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((B + FOUR)*(     &
               B + THREE)*(B + TWO)*(B + ONE)*B)) 
       ELSE IF (K - M == 2) THEN 
!  VIII                   -1   1
          SI = AKA*((A + B + FOUR)*(A - B - TWO) - (C - TWO)*(B + C + FOUR))/     &
               DSQRT(B*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)) 
          SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + ONE)*(AS - C + TWO)     &
               /((C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 2) THEN 
!  4                       1   1
          SI = AKA*(THREE*(AS - B)*(AS - C) - (AS - A)*(AS + FOUR))/DSQRT((B      &
               + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B) 
          SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)     &
               /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
       ELSE IF (K - M == 0) THEN 
!  6 P(12)                 0   1
          SI = -AKA*((AS - B - ONE)*(AS - C) - (AS - A)*(AS + THREE))/DSQRT((     &
               B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B) 
          SI = SI*DSQRT(THREE*TWO*(AS - B)*(AS - C + ONE)*(AS - A + ONE)*(AS      &
               + TWO)/((C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE))) 
       ENDIF
    ELSE IF (N - J == 0) THEN 
! 0
       IF (K - M == 4) THEN 
!  III                     -2   0
          SI = AKA*DSQRT(THREE*TWO*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B      &
               - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE))) 
          SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + ONE)*(AS - C + TWO)     &
               /((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 4) THEN 
!  5                        2   0
          SI = AKA*DSQRT(THREE*TWO*(AS - C)*(AS - C - ONE)*(AS - B + TWO)*(AS     &
               - B + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*    &
               (C + ONE))) 
          SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)     &
               /((B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE))) 
       ELSE IF (K - M == 2) THEN 
!  VII                     -1   0
          SI = AKA*HALF*((A + B + TWO)*(A - B) - C*C + FOUR)/DSQRT((B - ONE)*&
               B*(B + ONE)*(B + TWO)*(B + THREE)) 
          SI = SI*DSQRT(THREE*TWO*(AS + ONE)*(AS - A)*(AS - B)*(AS - C + ONE)&
               /((C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 2) THEN 
!  6                        1   0
          SI = -AKA*((AS - C - ONE)*(AS - B) - (AS - A)*(AS + THREE))/DSQRT((&
               C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C) 
          SI = SI*DSQRT(THREE*TWO*(AS - C)*(AS - B + ONE)*(AS - A + ONE)*(AS&
               + TWO)/((B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE))) 
       ELSE IF (K - M == 0) THEN 
!  IX                       0   0
          X1 = (AS - B)*(AS - B - ONE)*(AS - C)*(AS - C - ONE) 
          X2 = FOUR*(AS - B)*(AS - C)*(AS - A)*(AS + TWO) 
          X3 = (AS - A)*(AS - A - ONE)*(AS + THREE)*(AS + TWO) 
          SI = AKA*(X1 - X2 + X3)/DSQRT((B - ONE)*B*(B + ONE)*(B + TWO)*(B + &
               THREE)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)) 
       ENDIF
    ENDIF
    RETURN  
  END SUBROUTINE SIXJ2
!
!################################################################################
!      
  SUBROUTINE SIXJ3(J, K, L, M, N, ITIK, SI) 
!*******************************************************************
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2   2  |             [B.M.X. 75]                    *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructed by A. Senchuk                        September 2019 *
!                                                                  *
!*******************************************************************

    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER  :: J, K, L, M, N
    INTEGER, INTENT(IN) :: ITIK 
    REAL(DOUBLE), INTENT(OUT) :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: I1 
    REAL(DOUBLE) :: AS, A, B, C, AKA 
!-----------------------------------------------
    SI = ZERO 
    IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
       IF (IXJTIK(J,K,L,M,N,6) == 0) RETURN  
    ENDIF
    I1 = (J + K + L)/2 
    AS = DBLE(I1) 
    A = DBLE(L) 
    B = DBLE(J) 
    C = DBLE(K) 
    AKA = ONE 
    IF (MOD(I1,2) /= 0) AKA = -AKA 
    IF (J - N == 6) THEN 
! -3
       IF (K - M == 6) THEN 
!                        -3  -3
          SI = AKA*DSQRT((AS + ONE)*AS*(AS - ONE)*(AS - TWO)*(AS - THREE)*(AS&
               - FOUR)/((B + ONE)*B*(B - ONE)*(B - TWO)*(B - THREE)*(B - FOUR)&
               *(B - THREE - TWO))) 
          SI = SI*DSQRT((AS - A - THREE - TWO)*(AS - A - FOUR)*(AS - A - &
               THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - A)/((C + ONE)*C*(C - &
               ONE)*(C - TWO)*(C - THREE)*(C - FOUR)*(C - THREE - TWO))) 
       ELSE IF (M - K == 6) THEN 
!                        3  -3
          SI = AKA*DSQRT((AS - C)*(AS - C - ONE)*(AS - C - TWO)*(AS - C - &
               THREE)*(AS - C - FOUR)*(AS - C - THREE - TWO)/((B + ONE)*B*(B - &
               ONE)*(B - TWO)*(B - THREE)*(B - FOUR)*(B - THREE - TWO))) 
          SI = SI*DSQRT((AS - B + THREE + THREE)*(AS - B + THREE + TWO)*(AS&
               - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C&
               + FOUR + THREE)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR&
               )*(C + THREE)*(C + TWO)*(C + ONE))) 
       ELSE IF (K - M == 4) THEN 
!                       -2  -3
          SI = AKA*DSQRT(TWO*THREE*(AS - C)*(AS - B + ONE)*(AS - THREE)*(AS&
               - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - THREE - TWO)*(B - FOUR)*(&
               B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
          SI = SI*DSQRT((AS - A - FOUR)*(AS - A - THREE)*(AS - A - TWO)*(AS&
               - A - ONE)*(AS - A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)&
               *C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 4) THEN 
!                        2  -3
          SI = AKA*DSQRT(TWO*THREE*(AS + ONE)*(AS - A)*(AS - C - FOUR)*(AS - &
               C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((C + THREE + &
               THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
               ONE)*C)) 
          SI = SI*DSQRT((AS - B + THREE + TWO)*(AS - B + FOUR)*(AS - B + &
               THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - THREE - TWO)*(B - &
               FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
       ELSE IF (K - M == 2) THEN 
!                       -1   -3
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B&
               + TWO)*(AS - B + ONE)*(AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B&
               - THREE - TWO)*(B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
               + ONE))) 
          SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
               A)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + &
               THREE))) 
       ELSE IF (M - K == 2) THEN 
!                        1   -3
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - ONE)*(AS&
               - A)*(AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((&
               B - THREE - TWO)*(B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
               + ONE))) 
          SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
               - B + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*&
               (C + ONE)*C*(C - ONE))) 
       ELSE IF (M - K == 0) THEN 
!                        0   -3
          SI = AKA*TWO*DSQRT((THREE + TWO)*(AS + ONE)*AS*(AS - ONE)*(AS - A&
               - TWO)*(AS - A - ONE)*(AS - A)/((B - THREE - TWO)*(B - FOUR)*(B&
               - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
          SI = SI*DSQRT((AS - C - TWO)*(AS - C - ONE)*(AS - C)*(AS - B + &
               THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + FOUR)*(C + THREE)*(C&
               + TWO)*(C + ONE)*C*(C - ONE)*(C - TWO))) 
       ENDIF
    ELSE IF (N - J == 6) THEN 
!  3
       IF (K - M == 6) THEN 
!                        -3  3
          SI = AKA*DSQRT((AS - B)*(AS - B - ONE)*(AS - B - TWO)*(AS - B - &
               THREE)*(AS - B - FOUR)*(AS - B - THREE - TWO)/((C + ONE)*C*(C - &
               ONE)*(C - TWO)*(C - THREE)*(C - FOUR)*(C - THREE - TWO))) 
          SI = SI*DSQRT((AS - C + THREE + THREE)*(AS - C + THREE + TWO)*(AS&
               - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B&
               + FOUR + THREE)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR&
               )*(B + THREE)*(B + TWO)*(B + ONE))) 
       ELSE IF (M - K == 6) THEN 
!                        3   3
          SI = AKA*DSQRT((AS - A + THREE + THREE)*(AS - A + THREE + TWO)*(AS&
               - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)/((B&
               + SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
               THREE)*(B + TWO)*(B + ONE))) 
          SI = SI*DSQRT((AS + SEVEN)*(AS + THREE + THREE)*(AS + THREE + TWO)*&
               (AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + SEVEN)*(C + THREE + &
               THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
               ONE))) 
       ELSE IF (K - M == 4) THEN 
!                       -2   3
          SI = -AKA*DSQRT(TWO*THREE*(AS - A + ONE)*(AS + TWO)*(AS - B - FOUR)&
               *(AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((B + &
               SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
               THREE)*(B + TWO)*(B + ONE))) 
          SI = SI*DSQRT((AS - C + THREE + TWO)*(AS - C + FOUR)*(AS - C + &
               THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - FOUR)*(C - THREE)*(C&
               - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
       ELSE IF (M - K == 4) THEN 
!                        2   3
          SI = -AKA*DSQRT(TWO*THREE*(AS - B)*(AS - C + ONE)*(AS - A + THREE&
               + TWO)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A&
               + ONE)/((B + SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + &
               FOUR)*(B + THREE)*(B + TWO)*(B + ONE))) 
          SI = SI*DSQRT((AS + THREE + THREE)*(AS + THREE + TWO)*(AS + FOUR)*(&
               AS + THREE)*(AS + TWO)/((C + THREE + THREE)*(C + THREE + TWO)*(C&
               + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
       ELSE IF (K - M == 2) THEN 
!                       -1   3
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - A + ONE)*(AS - A + TWO)*(&
               AS + THREE)*(AS + TWO)*(AS - B - THREE)*(AS - B - TWO)*(AS - B&
               - ONE)*(AS - B)/((B + SEVEN)*(B + THREE + THREE)*(B + THREE + &
               TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE))) 
          SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
               - C + ONE)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + &
               TWO)*(C + THREE))) 
       ELSE IF (M - K == 2) THEN 
!                        1   3
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C&
               + TWO)*(AS - C + ONE)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A&
               + TWO)*(AS - A + ONE)/((B + SEVEN)*(B + THREE + THREE)*(B + &
               THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE))) 
          SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
               )/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*&
               C*(C - ONE))) 
       ELSE IF (M - K == 0) THEN 
!                        0   3
          SI = -AKA*TWO*DSQRT((THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - B&
               - TWO)*(AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + &
               SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
               THREE)*(B + TWO)*(B + ONE))) 
          SI = SI*DSQRT((AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)*(AS + &
               FOUR)*(AS + THREE)*(AS + TWO)/((C + FOUR)*(C + THREE)*(C + TWO)*&
               (C + ONE)*C*(C - ONE)*(C - TWO))) 
       ENDIF
    ELSE IF (J - N == 4) THEN 
! -2
       IF (K - M == 6) THEN 
!                       -3  -2
          SI = AKA*DSQRT(TWO*THREE*(AS - B)*(AS - C + ONE)*(AS - THREE)*(AS&
               - TWO)*(AS - ONE)*AS*(AS + ONE)/((C - THREE - TWO)*(C - FOUR)*(&
               C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
          SI = SI*DSQRT((AS - A - FOUR)*(AS - A - THREE)*(AS - A - TWO)*(AS&
               - A - ONE)*(AS - A)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)&
               *B*(B + ONE)*(B + TWO))) 
       ELSE IF (M - K == 6) THEN 
!                       3   -2
          SI = -AKA*DSQRT(TWO*THREE*(AS - A + ONE)*(AS + TWO)*(AS - C - FOUR)&
               *(AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((C + &
               SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + &
               THREE)*(C + TWO)*(C + ONE))) 
          SI = SI*DSQRT((AS - B + THREE + TWO)*(AS - B + FOUR)*(AS - B + &
               THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - FOUR)*(B - THREE)*(B&
               - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
       ELSE IF (K - M == 4) THEN 
!                      -2  -2
          SI = AKA*((TWO + THREE)*(AS - C)*(AS - B) - (AS + TWO)*(AS - A - &
               FOUR))*DSQRT((AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - FOUR)*(B&
               - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
          SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
               A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + &
               TWO))) 
       ELSE IF (M - K == 4) THEN 
!                       2  -2
          SI = -AKA*((TWO + THREE)*(AS + ONE)*(AS - A + ONE) - (AS - C + ONE)&
               *(AS - B + THREE + TWO))*DSQRT((AS - C - THREE)*(AS - C - TWO)*(&
               AS - C - ONE)*(AS - C)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - &
               ONE)*B*(B + ONE)*(B + TWO))) 
          SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
               - B + ONE)/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C&
               + THREE)*(C + TWO)*(C + ONE)*C)) 
       ELSE IF (K - M == 2) THEN 
!                       -1  -2
          SI = AKA*(TWO*(AS - C - ONE)*(AS - B) - (AS + TWO)*(AS - A - THREE)&
               )*DSQRT(TWO*(THREE + TWO)*(AS - C)*(AS - B + ONE)*(AS - ONE)*AS*&
               (AS + ONE)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B + &
               ONE)*(B + TWO))) 
          SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)/((C - THREE)*(&
               C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
       ELSE IF (M - K == 2) THEN 
!                        1  -2
          SI = -AKA*(TWO*AS*(AS - A + ONE) - (AS - C + ONE)*(AS - B + FOUR))*&
               DSQRT(TWO*(THREE + TWO)*(AS + ONE)*(AS - A)*(AS - C - TWO)*(AS&
               - C - ONE)*(AS - C)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)&
               *B*(B + ONE)*(B + TWO))) 
          SI = SI*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + &
               THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - &
               ONE))) 
       ELSE IF (M - K == 0) THEN 
!                        0  -2
          SI = -AKA*((AS - ONE)*(AS - A + ONE) - (AS - C + ONE)*(AS - B + &
               THREE))*DSQRT(TWO*THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - &
               ONE)*(AS - A)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
               + ONE)*(B + TWO))) 
          SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + TWO)*(AS - B + ONE)&
               /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - &
               TWO))) 
       ENDIF
    ELSE IF (N - J == 4) THEN 
!  2
       IF (K - M == 6) THEN 
!                        -3  2
          SI = AKA*DSQRT(TWO*THREE*(AS + ONE)*(AS - A)*(AS - B - FOUR)*(AS - &
               B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((B + THREE + &
               THREE)*(B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + &
               ONE)*B)) 
          SI = SI*DSQRT((AS - C + THREE + TWO)*(AS - C + FOUR)*(AS - C + &
               THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - THREE - TWO)*(C - &
               FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
       ELSE IF (M - K == 6) THEN 
!                        3  2
          SI = -AKA*DSQRT(TWO*THREE*(AS - C)*(AS - B + ONE)*(AS - A + THREE&
               + TWO)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A&
               + ONE)/((C + SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + &
               FOUR)*(C + THREE)*(C + TWO)*(C + ONE))) 
          SI = SI*DSQRT((AS + THREE + THREE)*(AS + THREE + TWO)*(AS + FOUR)*(&
               AS + THREE)*(AS + TWO)/((B + THREE + THREE)*(B + THREE + TWO)*(B&
               + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)) 
       ELSE IF (K - M == 4) THEN 
!                       -2  2
          SI = -AKA*((TWO + THREE)*(AS + ONE)*(AS - A + ONE) - (AS - B + ONE)&
               *(AS - C + THREE + TWO))*DSQRT((AS - B - THREE)*(AS - B - TWO)*(&
               AS - B - ONE)*(AS - B)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - &
               ONE)*C*(C + ONE)*(C + TWO))) 
          SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
               - C + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B&
               + THREE)*(B + TWO)*(B + ONE)*B)) 
       ELSE IF (M - K == 4) THEN 
!                        2  2
          SI = AKA*((TWO + THREE)*(AS - B)*(AS - C) - (AS - A)*(AS + THREE + &
               THREE))*DSQRT((AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(&
               AS - A + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*&
               (B + THREE)*(B + TWO)*(B + ONE)*B)) 
          SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
               )/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*&
               (C + TWO)*(C + ONE)*C)) 
       ELSE IF (K - M == 2) THEN 
!                       -1  2
          SI = AKA*(TWO*(AS - A + TWO)*(AS + ONE) - (AS - B + ONE)*(AS - C +   &
               FOUR))*DSQRT(TWO*(THREE + TWO)*(AS - A + ONE)*(AS + TWO)*(AS - B&
               - TWO)*(AS - B - ONE)*(AS - B)/((B + THREE + THREE)*(B + THREE  &
               + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)) 
          SI = SI*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((C -   &
               THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
       ELSE IF (M - K == 2) THEN 
!                        1  2
          SI = -AKA*(TWO*(AS - B - ONE)*(AS - C) - (AS - A)*(AS + THREE + TWO  &
               ))*DSQRT(TWO*(THREE + TWO)*(AS - B)*(AS - C + ONE)*(AS - A +    &
               THREE)*(AS - A + TWO)*(AS - A + ONE)/((B + THREE + THREE)*(B +  &
               THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)) 
          SI = SI*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + THREE + TWO  &
               )*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE))) 
       ELSE IF (M - K == 0) THEN 
!                        0  2
          SI = AKA*((AS - B - TWO)*(AS - C) - (AS - A)*(AS + FOUR))*DSQRT(TWO  &
               *THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C + TWO)*(AS &
               - C + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B &
               + THREE)*(B + TWO)*(B + ONE)*B)) 
          SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)  &
               /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C -   &
               TWO))) 
       ENDIF
    ELSE IF (J - N == 2) THEN 
! - 1
       IF (K - M == 6) THEN 
!                       -3   -1
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C   &
               + TWO)*(AS - C + ONE)*(AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((C   &
               - THREE - TWO)*(C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C  &
               + ONE))) 
          SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS -   &
               A)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B +  &
               THREE))) 
       ELSE IF (M - K == 6) THEN 
!                       3    -1
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - A + ONE)*(AS - A + TWO)*(   &
               AS + THREE)*(AS + TWO)*(AS - C - THREE)*(AS - C - TWO)*(AS - C  &
               - ONE)*(AS - C)/((C + SEVEN)*(C + THREE + THREE)*(C + THREE +   &
               TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))) 
          SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS    &
               - B + ONE)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B +    &
               TWO)*(B + THREE))) 
       ELSE IF (K - M == 4) THEN 
!                       -2   -1
          SI = AKA*(TWO*(AS - B - ONE)*(AS - C) - (AS + TWO)*(AS - A - THREE)  &
               )*DSQRT(TWO*(THREE + TWO)*(AS - B)*(AS - C + ONE)*(AS - ONE)*AS*&
               (AS + ONE)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C +   &
               ONE)*(C + TWO))) 
          SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)/((B - THREE)*(  &
               B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE))) 
       ELSE IF (M - K == 4) THEN 
!                       2    -1
          SI = AKA*(TWO*(AS - A + TWO)*(AS + ONE) - (AS - C + ONE)*(AS - B +   &
               FOUR))*DSQRT(TWO*(THREE + TWO)*(AS - A + ONE)*(AS + TWO)*(AS - C&
               - TWO)*(AS - C - ONE)*(AS - C)/((C + THREE + THREE)*(C + THREE  &
               + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
          SI = SI*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((B -   &
               THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE))) 
       ELSE IF (K - M == 2) THEN 
!                        -1   -1
          SI = AKA*(TWO*THREE*(AS - C)*(AS - C - ONE)*(AS - B - ONE)*(AS - B)  &
               - TWO*FOUR*(AS - C)*(AS - B)*(AS + TWO)*(AS - A - TWO) + (AS +  &
               THREE)*(AS + TWO)*(AS - A - THREE)*(AS - A - TWO)) 
          SI = SI*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B - THREE)*(B  &
               - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)*(C - THREE)* &
               (C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
       ELSE IF (M - K == 2) THEN 
!                         1   -1
          SI = AKA*(TWO*THREE*(AS + ONE)*AS*(AS - A + ONE)*(AS - A + TWO) -    &
               TWO*FOUR*(AS + ONE)*(AS - A + ONE)*(AS - C + ONE)*(AS - B +     &
               THREE) + (AS - C + ONE)*(AS - C + TWO)*(AS - B + FOUR)*(AS - B  &
               + THREE)) 
          SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + TWO)*(AS - B + ONE)  &
               /((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B +    &
               THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C +  &
               ONE)*C*(C - ONE))) 
       ELSE IF (M - K == 0) THEN 
!                         0   -1
          SI = AKA*TWO*(AS*(AS - ONE)*(AS - A + ONE)*(AS - A + TWO) - THREE*   &
               AS*(AS - A + ONE)*(AS - C + ONE)*(AS - B + TWO) + (AS - C + ONE)&
               *(AS - C + TWO)*(AS - B + THREE)*(AS - B + TWO)) 
          SI = SI*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - C)*(AS - B + ONE)/((B  &
               - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)* &
               (C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - TWO)&
               )) 
       ENDIF
    ELSE IF (N - J == 2) THEN 
! - 1
       IF (K - M == 6) THEN 
!                        -3   1
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - ONE)*(AS  &
               - A)*(AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/(( &
               C - THREE - TWO)*(C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C&
               + ONE))) 
          SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS    &
               - C + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)* &
               (B + ONE)*B*(B - ONE))) 
       ELSE IF (M - K == 6) THEN 
!                        3    1
          SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B   &
               + TWO)*(AS - B + ONE)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A  &
               + TWO)*(AS - A + ONE)/((C + SEVEN)*(C + THREE + THREE)*(C +     &
               THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))) 
          SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO  &
               )/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*&
               B*(B - ONE))) 
       ELSE IF (K - M == 4) THEN 
!                        -2   1
          SI = -AKA*(TWO*AS*(AS - A + ONE) - (AS - B + ONE)*(AS - C + FOUR))*  &
               DSQRT(TWO*(THREE + TWO)*(AS + ONE)*(AS - A)*(AS - B - TWO)*(AS  &
               - B - ONE)*(AS - B)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE) &
               *C*(C + ONE)*(C + TWO))) 
          SI = SI*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B +   &
               THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B -  &
               ONE))) 
       ELSE IF (M - K == 4) THEN 
!                         2   1
          SI = -AKA*(TWO*(AS - C - ONE)*(AS - B) - (AS - A)*(AS + THREE + TWO  &
               ))*DSQRT(TWO*(THREE + TWO)*(AS - C)*(AS - B + ONE)*(AS - A +    &
               THREE)*(AS - A + TWO)*(AS - A + ONE)/((C + THREE + THREE)*(C +  &
               THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
          SI = SI*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((B + THREE + TWO  &
               )*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE))) 
       ELSE IF (M - K == 2) THEN 
!                         1   1
          SI = AKA*(TWO*THREE*(AS - B)*(AS - B - ONE)*(AS - C)*(AS - C - ONE)  &
               - TWO*FOUR*(AS - B)*(AS - C)*(AS - A)*(AS + FOUR) + (AS - A)*(  &
               AS - A - ONE)*(AS + THREE + TWO)*(AS + FOUR)) 
          SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)  &
               /((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B&
               *(B - ONE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C&
               + ONE)*C*(C - ONE))) 
       ELSE IF (K - M == 2) THEN 
!                        -1   1
          SI = AKA*(TWO*THREE*(AS + ONE)*AS*(AS - A + ONE)*(AS - A + TWO) -    &
               TWO*FOUR*(AS + ONE)*(AS - A + ONE)*(AS - B + ONE)*(AS - C +     &
               THREE) + (AS - B + ONE)*(AS - B + TWO)*(AS - C + FOUR)*(AS - C  &
               + THREE)) 
          SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + TWO)*(AS - C + ONE)  &
               /((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C +    &
               THREE)*(B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B +  &
               ONE)*B*(B - ONE))) 
       ELSE IF (M - K == 0) THEN 
!                         0   1
          SI = -AKA*TWO*((AS - B - ONE)*(AS - B - TWO)*(AS - C)*(AS - C - ONE  &
               ) - THREE*(AS - B - ONE)*(AS - C)*(AS - A)*(AS + THREE) + (AS - &
               A)*(AS - A - ONE)*(AS + FOUR)*(AS + THREE)) 
          SI = SI*DSQRT(THREE*(AS - B)*(AS - C + ONE)*(AS - A + ONE)*(AS +     &
               TWO)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B +   &
               ONE)*B*(B - ONE)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C&
               - ONE)*(C - TWO))) 
       ENDIF
    ELSE IF (J - N == 0) THEN 
!  0
       IF (K - M == 6) THEN 
!                       -3    0
          SI = AKA*TWO*DSQRT((THREE + TWO)*(AS + ONE)*AS*(AS - ONE)*(AS - A    &
               - TWO)*(AS - A - ONE)*(AS - A)/((C - THREE - TWO)*(C - FOUR)*(C &
               - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
          SI = SI*DSQRT((AS - B - TWO)*(AS - B - ONE)*(AS - B)*(AS - C +       &
               THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + FOUR)*(B + THREE)*(B &
               + TWO)*(B + ONE)*B*(B - ONE)*(B - TWO))) 
       ELSE IF (M - K == 6) THEN 
!                        3    0
          SI = -AKA*TWO*DSQRT((THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - C    &
               - TWO)*(AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C +     &
               SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C +    &
               THREE)*(C + TWO)*(C + ONE))) 
          SI = SI*DSQRT((AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)*(AS +   &
               FOUR)*(AS + THREE)*(AS + TWO)/((B + FOUR)*(B + THREE)*(B + TWO)*&
               (B + ONE)*B*(B - ONE)*(B - TWO))) 
       ELSE IF (K - M == 4) THEN 
!                       -2    0
          SI = -AKA*((AS - ONE)*(AS - A + ONE) - (AS - B + ONE)*(AS - C +      &
               THREE))*DSQRT(TWO*THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A -   &
               ONE)*(AS - A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C  &
               + ONE)*(C + TWO))) 
          SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + TWO)*(AS - C + ONE)  &
               /((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B -   &
               TWO))) 
       ELSE IF (M - K == 4) THEN 
!                        2    0
          SI = AKA*((AS - C - TWO)*(AS - B) - (AS - A)*(AS + FOUR))*DSQRT(TWO  &
               *THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B + TWO)*(AS &
               - B + ONE)/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C &
               + THREE)*(C + TWO)*(C + ONE)*C)) 
          SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)  &
               /((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B -   &
               TWO))) 
       ELSE IF (K - M == 2) THEN 
!                       -1    0
          SI = AKA*TWO*(AS*(AS - ONE)*(AS - A + ONE)*(AS - A + TWO) - THREE*   &
               AS*(AS - A + ONE)*(AS - B + ONE)*(AS - C + TWO) + (AS - B + ONE)&
               *(AS - B + TWO)*(AS - C + THREE)*(AS - C + TWO)) 
          SI = SI*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - B)*(AS - C + ONE)/((C  &
               - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)* &
               (B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B - TWO)&
               )) 
       ELSE IF (M - K == 2) THEN 
!                        1    0
          SI = -AKA*TWO*((AS - C - ONE)*(AS - C - TWO)*(AS - B)*(AS - B - ONE  &
               ) - THREE*(AS - C - ONE)*(AS - B)*(AS - A)*(AS + THREE) + (AS - &
               A)*(AS - A - ONE)*(AS + FOUR)*(AS + THREE)) 
          SI = SI*DSQRT(THREE*(AS - C)*(AS - B + ONE)*(AS - A + ONE)*(AS +     &
               TWO)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C +   &
               ONE)*C*(C - ONE)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B&
               - ONE)*(B - TWO))) 
       ELSE IF (K - M == 0) THEN 
!                        0    0
          SI = AKA*((AS - B)*(AS - B - ONE)*(AS - B - TWO)*(AS - C)*(AS - C    &
               - ONE)*(AS - C - TWO) - THREE*THREE*(AS - B)*(AS - B - ONE)*(AS &
               - C)*(AS - C - ONE)*(AS - A)*(AS + TWO) + THREE*THREE*(AS - B)* &
               (AS - C)*(AS - A)*(AS - A - ONE)*(AS + THREE)*(AS + TWO) - (AS  &
               - A)*(AS - A - ONE)*(AS - A - TWO)*(AS + FOUR)*(AS + THREE)*(AS &
               + TWO)) 
          SI = SI/DSQRT((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE  &
               )*(B - TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C -   &
               ONE)*(C - TWO)) 
       ENDIF
    ENDIF
    RETURN  
  END SUBROUTINE SIXJ3
!
!#################################################################################
!  
  SUBROUTINE SIXJ4(JC,JE,JD,JB,JF,ITIK,SI) 
!*******************************************************************
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | JC/2  JE/2  JD/2 |                                         *
!     | JB/2  JF/2    4  |                                         *
!                                                                  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructed by A. Senchuk                        September 2019 *
!                                                                  *
!*******************************************************************
      
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER             :: JC, JE, JD, JB, JF
    INTEGER, INTENT(IN) :: ITIK 
    REAL(DOUBLE)        :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    REAL(DOUBLE) :: A, C, E, D, B, F, X1, X2, X3, S2, S3 
!-----------------------------------------------
    SI = ZERO 
    IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
       IF (IXJTIK(JC,JE,JD,JB,JF,8) == 0) RETURN  
    ENDIF
    IF (IXJTIK(JC,JE,JD,JB,JF,6) == 0) THEN 
       CALL DRACAH (JC, JE, JF, JB, JD, 8, SI) 
       IF (MOD(JC + JE + JF + JB,4) /= 0) SI = -SI 
    ELSE 
       A = THREE 
       C = DBLE(JC)*HALF 
       E = DBLE(JE)*HALF 
       D = DBLE(JD)*HALF 
       B = DBLE(JB)*HALF 
       F = DBLE(JF)*HALF 
       X1 = A*DSQRT((A+B+E+TWO)*(A-B+E+ONE)*(A+B-E+ONE)*((-           &
            A)+B+E)*(A+C+F+TWO)*(A-C+F+ONE)*(A+C-F+ONE)*(             &
            (-A)+C+F)) 
       X2 = (A+ONE)*DSQRT((A+B+E+ONE)*(A-B+E)*(A+B-E)*((-A)           &
            + B+E+ONE)*(A+C+F+ONE)*(A-C+F)*(A+C-F)*((-A)+C            &
            + F+ONE)) 
       X3 = (TWO*A+ONE)*(TWO*(A*(A+ONE)*D*(D+ONE)-B*(B+ONE)*C*(C+     &
            ONE)-E*(E+ONE)*F*(F+ONE))+(A*(A+ONE)-B*(B+ONE)-E*(E       &
            +ONE))*(A*(A+ONE)-C*(C+ONE)-F*(F+ONE))) 
       IF (DABS(X2) < EPS) THEN 
          S2 = ZERO 
       ELSE 
          CALL SIXJ2 (JC, JE, JD, JB, JF, 0, S2) 
       ENDIF
       IF (DABS(X3) < EPS) THEN 
          S3 = ZERO 
       ELSE 
          CALL SIXJ3 (JC, JE, JD, JB, JF, 0, S3) 
       ENDIF
       SI = (X3*S3 - X2*S2)/X1 
    ENDIF
    RETURN  
  END SUBROUTINE SIXJ4
  
END SUBROUTINE SIXJ
