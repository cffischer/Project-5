MODULE RECO_M

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ONE
      USE m_C,             ONLY: JLIST, JJQ1, JJQ2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE diaga

CONTAINS

  SUBROUTINE RECO2(JA1,JA2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1,JA2,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IA1,IA2,IB1,IB2,IJ1,IJ2,ISKR
    REAL(DOUBLE) :: S, SS, RE
!-----------------------------------------------
    IAT=0
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    S=DBLE(JJQ1(3,IJ1))
    SS=DBLE(JJQ1(3,IJ2))
    SS=S*SS
    RECC=ONE/DSQRT(SS)
    IF(IRE == 0) THEN
       IAT=0
    ELSE IF(KA /= 0) THEN
       IAT=0
    ELSE
       IAT=1
       RETURN
    END IF
    IA1=JJQ1(3,IJ1)-1
    IA2=JJQ1(3,IJ2)-1
    IB1=JJQ2(3,IJ1)-1
    IB2=JJQ2(3,IJ2)-1
!
    CALL DIAGA2(JA1,JA2,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC*DSQRT(DBLE(IA2+1))/DSQRT(DBLE((KA+1)*(IB2+1)))
    IF(JA1 == 1.AND.JA2 == 2)RETURN
!
    IAT=0
    CALL DIAGA1(JA1,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    ISKR=JA2-JA1
    IF(JA1 == 1)ISKR=JA2-1-JA1
    IF(ISKR <= 1)RETURN
!
    IAT=0
    CALL DIAGA3(JA1,JA2,KA,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECO2
  !
  !##############################################################
  !
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
10  WRITE(99,100)
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
  !
  !######################################################################
  !
  SUBROUTINE RECO4(JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 09  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3,DIAGA4              *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)  :: JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,IRE
    INTEGER, INTENT(OUT) :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IA4, IB4, IJ1, IJ2, IJ3, IJ4, ISKR
    REAL(DOUBLE) :: S, S1, S2, S3, S4, RE
!-----------------------------------------------
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IJ3=JLIST(JA3)
    IJ4=JLIST(JA4)
    S1=JJQ1(3,IJ1)
    S2=JJQ1(3,IJ2)
    S3=JJQ1(3,IJ3)
    S4=JJQ1(3,IJ4)
    S=S1*S2*S3*S4
    RECC=ONE/DSQRT(S)
    IA4=JJQ1(3,IJ4)-1
    IB4=JJQ2(3,IJ4)-1
    RECC=RECC*DSQRT(DBLE(IA4+1))/DSQRT(DBLE((K4+1)*(IB4+1)))
!
    ISKR=JA3-JA2
    IF(ISKR > 1) THEN
       IAT=0
       CALL DIAGA3(JA2,JA3,KA,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
    END IF
!
    ISKR=JA4-JA3
    IF(ISKR > 1) THEN
       IAT=0
       CALL DIAGA3(JA3,JA4,K4,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
    END IF
!
    IAT=0
    CALL DIAGA2(JA1,JA4,K4,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    IAT=0
    CALL DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    IAT=0
    CALL DIAGA4(JA2,JA3,KA,K3,K4,IRE,IAT,RE)
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
  END SUBROUTINE RECO4

  
END MODULE RECO_M
