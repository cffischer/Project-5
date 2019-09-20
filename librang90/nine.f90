SUBROUTINE NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,INN,AA) 
!*******************************************************************
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  L3/2 |                                        *
!     |  K1/2  K2/2  K3/2 |                                        *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             March    1995  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructured by A. Senchuk                     September 2019  *  
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE vast_kind_param 
  USE CONS_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
  USE ittk_I  
  USE sixj_I
  
  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER  :: J1, J2, J3, L1, L2, L3, K1, K2, K3
  INTEGER, INTENT(IN)  :: I 
  INTEGER, INTENT(OUT) :: INN
  REAL(DOUBLE)  :: AA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: N1, N2, N3, N4, N5, N6, MAX_, MIN_, IX 
  REAL(DOUBLE) :: S1, S2, S3, X 
!-----------------------------------------------
  IF (I == 1) THEN 
     INN = 0 
     IF (ITTK(J1,J2,J3) == 0) RETURN  
     IF (ITTK(L1,L2,L3) == 0) RETURN  
     IF (ITTK(K1,K2,K3) == 0) RETURN  
     IF (ITTK(J1,L1,K1) == 0) RETURN  
     IF (ITTK(J2,L2,K2) == 0) RETURN  
     IF (ITTK(J3,L3,K3) == 0) RETURN  
     INN = 1 
     RETURN  
  ENDIF
  IF (J1*J2*J3*L1*L2*L3*K1*K2*K3 == 0) THEN 
     INN = 1 
     CALL NINE0 (J1, J2, J3, L1, L2, L3, K1, K2, K3, AA) 
  ELSE 
     N1 = IABS(J1 - K3) 
     N2 = IABS(L3 - J2) 
     N3 = IABS(L1 - K2) 
     N4 = IABS(J2 - L3) 
     N5 = IABS(K2 - L1) 
     N6 = IABS(J1 - K3) 
     MAX_ = MAX0(N1,N2,N3,N4,N5,N6) 
     N1 = J1 + K3 
     N2 = L3 + J2 
     N3 = J2 + L3 
     N4 = K2 + L1 
     N5 = J1 + K3 
     N6 = L1 + K2 
     MIN_ = MIN0(N1,N2,N3,N4,N5,N6) 
     INN = 1 
     AA = ZERO 
     DO IX = MAX_, MIN_, 2 
        CALL SIXJ (J1, J2, J3, L3, K3, IX, 0, S1) 
        CALL SIXJ (L1, L2, L3, J2, IX, K2, 0, S2) 
        CALL SIXJ (K1, K2, K3, IX, J1, L1, 0, S3) 
        X = S1*S2*S3*DBLE(IX + 1) 
        IF (MOD(IX,2) /= 0) X = -X 
        AA = X + AA 
     END DO
  ENDIF
  RETURN

CONTAINS

  SUBROUTINE NINE0(J1,J2,J3,L1,L2,L3,K1,K2,K3,AA) 
!*******************************************************************
!    THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT          *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  L3/2 |                                        *
!     |  K1/2  K2/2    0  |                                        *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER  :: J1, J2, J3, L1, L2, L3, K1, K2, K3
    REAL(DOUBLE), INTENT(OUT) :: AA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IFA 
    REAL(DOUBLE) :: A, B 
!-----------------------------------------------
    IF (J1 == 0) THEN 
       CALL SIXJ (L2, K2, J2, K3, L3, L1, 0, A) 
       B = DBLE((J2 + 1)*(L1 + 1)) 
       IFA = K2 + J2 + L3 + L1 
    ELSE IF (J2 == 0) THEN 
       CALL SIXJ (L3, K3, J3, K1, L1, L2, 0, A) 
       B = DBLE((J3 + 1)*(L2 + 1)) 
       IFA = K3 + J3 + L1 + L2 
    ELSE IF (J3 == 0) THEN 
       CALL SIXJ (L1, K1, J1, K2, L2, L3, 0, A) 
       B = DBLE((J1 + 1)*(L3 + 1)) 
       IFA = K1 + J1 + L2 + L3 
    ELSE IF (L1 == 0) THEN 
       CALL SIXJ (K2, J2, L2, J3, K3, K1, 0, A) 
       B = DBLE((L2 + 1)*(K1 + 1)) 
       IFA = J2 + L2 + K3 + K1 
    ELSE IF (L2 == 0) THEN 
       CALL SIXJ (K3, J3, L3, J1, K1, K2, 0, A) 
       B = DBLE((L3 + 1)*(K2 + 1)) 
       IFA = J3 + L3 + K1 + K2 
    ELSE IF (L3 == 0) THEN 
       CALL SIXJ (K1, J1, L1, J2, K2, K3, 0, A) 
       B = DBLE((L1 + 1)*(K3 + 1)) 
       IFA = J1 + L1 + K2 + K3 
    ELSE IF (K1 == 0) THEN 
       CALL SIXJ (J2, J3, J1, L3, L2, K2, 0, A) 
       B = DBLE((J1 + 1)*(K2 + 1)) 
       IFA = J3 + J1 + L2 + K2 
    ELSE IF (K2 == 0) THEN 
       CALL SIXJ (J3, J1, J2, L1, L3, K3, 0, A) 
       B = DBLE((J2 + 1)*(K3 + 1)) 
       IFA = J1 + J2 + L3 + K3 
    ELSE IF (K3 == 0) THEN 
       CALL SIXJ (J1, J2, J3, L2, L1, K1, 0, A) 
       B = DBLE((J3 + 1)*(K1 + 1)) 
       IFA = J2 + J3 + L1 + K1 
    ELSE 
       A = ZERO 
       B = ONE 
    ENDIF
    AA = A/DSQRT(B) 
    IF (MOD(IFA,4) /= 0) AA = -AA 
    RETURN  
  END SUBROUTINE NINE0
  
END SUBROUTINE NINE
