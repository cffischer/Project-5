SUBROUTINE REC3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  RECO3                                    *
!                                                                  *
!   Written by  G. Gaigalas                                        *
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
  USE m_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
  USE diaga1_I
  USE diaga2_I
  USE diaga3_I
  USE diaga4_I
  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN)       :: JA1,JA2,JA3,K1,K2,K3,IRE
  INTEGER, INTENT(OUT)      :: IAT
  REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: IFAZ
!-----------------------------------------------
  IF((JA3 > JA1).AND.(JA3 > JA2)) THEN
     IF(JA1-JA2 < 0) THEN
        CALL RECO3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,RECC)
     ELSE IF(JA1-JA2 > 0) THEN
        CALL RECO3(JA2,JA1,JA3,K2,K1,K3,IRE,IAT,RECC)
        IFAZ=K1+K2-K3
        IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
     ELSE
        GO TO 10
     END IF
  ELSE IF((JA3 < JA1).AND.(JA3 < JA2)) THEN
     IF(JA1-JA2 < 0) THEN
        CALL RECO3(JA3,JA1,JA2,K3,K1,K2,IRE,IAT,RECC)
        IF((K3/2)*2 /= K3)RECC=-RECC
     ELSE IF(JA1-JA2 > 0) THEN
        CALL RECO3(JA3,JA2,JA1,K3,K2,K1,IRE,IAT,RECC)
        IFAZ=K1+K2+K3
        IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
     ELSE
        GO TO 10
     END IF
  ELSE
     IF(JA1-JA2 < 0)THEN
        CALL RECO3(JA1,JA3,JA2,K1,K3,K2,IRE,IAT,RECC)
        IFAZ=K1-K2-K3
        IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
     ELSE IF(JA1-JA2 > 0) THEN
        CALL RECO3(JA2,JA3,JA1,K2,K3,K1,IRE,IAT,RECC)
        IF((K1/2)*2 /= K1)RECC=-RECC
     ELSE
        GO TO 10
     END IF
  END IF
  RETURN
10 WRITE(99,100)
100 FORMAT(5X,'ERRO IN REC')
  STOP

CONTAINS

  SUBROUTINE RECO3(JA1,JA2,JA3,K1,K2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3,DIAGA4              *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1,JA2,JA3,K1,K2,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IA3, IB3, IJ1, IJ2, IJ3, ISKR
    REAL(DOUBLE) :: S, S1, S2, S3, RE
!-----------------------------------------------
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IJ3=JLIST(JA3)
    S1=JJQ1(3,IJ1)
    S2=JJQ1(3,IJ2)
    S3=JJQ1(3,IJ3)
    S=S1*S2*S3
    RECC=ONE/DSQRT(S)
    IA3=JJQ1(3,IJ3)-1
    IB3=JJQ2(3,IJ3)-1
    RECC=RECC*DSQRT(DBLE(IA3+1))/DSQRT(DBLE((KA+1)*(IB3+1)))
!
    IAT=0
    ISKR=JA3-JA2
    IF(ISKR > 1) THEN
       CALL DIAGA3(JA2,JA3,KA,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
    END IF
!
    IAT=0
    CALL DIAGA2(JA1,JA3,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    IAT=0
    CALL DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    IF(JA1 == 1.AND.JA2 == 2)RETURN
!
    IAT=0
    CALL DIAGA1(JA1,K1,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    ISKR=JA2-JA1
    IF(JA1 == 1)ISKR=JA2-1-JA1
    IF(ISKR <= 1)RETURN
    IAT=0
    CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECO3
      
END SUBROUTINE REC3
    
