MODULE RECOP
!*******************************************************************
!                                                                  *
!   Subroutines written by  G. Gaigalas                            *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Module created by A. Senchuk                   September 2019  *  
!                                                                  *
!*******************************************************************

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
  USE vast_kind_param
  USE CONS_C
  USE m_C
  USE diaga
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
  USE ittk_I

CONTAINS

  SUBROUTINE RECOP00(NS,JA1,JA2,KA,IAT)
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)  :: NS, JA1, JA2, KA
    INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: I, ISKR, IJ1, J, JI, JJ, NPEELGG, KK, LK, LD
!-----------------------------------------------
    IAT=1
    IF(NPEEL == 1 .AND. NS == -1)RETURN
    IAT=0
    IF(NS == -1) THEN
       NPEELGG = NPEEL
    ELSE
       NPEELGG = NS
    END IF
    IF(NPEELGG > 1) THEN
       LK=JJC1(NPEELGG-1)-1
       LD=JJC2(NPEELGG-1)-1
    ELSE IF(NPEELGG == 1) THEN
       LK=JJC1(NPEELGG)-1
       LD=JJC2(NPEELGG)-1
    ELSE
       PRINT*, "ERROR in RECOP00"
       STOP
    END IF
    IF(ITTK(LK,LD,2*KA) == 0)RETURN
    IAT=1
    IF(NPEEL == 1)RETURN
    DO I=1,NPEEL
       IF(JA1 /= I) THEN
          IF(JA2 /= I) THEN
             IJ1=JLIST(I)
             IF(JJQ1(1,IJ1) /= JJQ2(1,IJ1))IAT=0
             IF(JJQ1(2,IJ1) /= JJQ2(2,IJ1))IAT=0
             IF(JJQ1(3,IJ1) /= JJQ2(3,IJ1))IAT=0
          ENDIF
       ENDIF
    END DO
    IF(IAT == 0)RETURN
    IF(NPEELGG <= 2)RETURN      
    IF(JA1 <= 2)RETURN
    DO J=3,JA1
       JJ=J-2
       IF(JJC1(JJ) /= JJC2(JJ))IAT=0
       IF(IAT == 0)RETURN
    END DO
    ISKR=NPEELGG-JA2
    IF(ISKR > 0) THEN
       DO JI=1,ISKR
          KK=JA2-2+JI
          LK=JJC1(KK)-1
          LD=JJC2(KK)-1
          IF(ITTK(LK,LD,2*KA) == 0)IAT=0
          IF(IAT == 0)RETURN
       END DO
    ENDIF
    RETURN
  END SUBROUTINE RECOP00
!
!###################################################################
!
  SUBROUTINE RECOP1(NS,JA1,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA3,DIAGA5                     *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: NS, JA1, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IJ1, ISKR, NPEELGG
    REAL(DOUBLE) :: S, RE
!-----------------------------------------------
    IAT=1
    IJ1=JLIST(JA1)
    S=DBLE(JJQ1(3,IJ1))
    RECC=ONE/DSQRT(S)
    IF(NPEEL == 1 .AND. NS == -1)RETURN
    IF(NS == -1) THEN
       NPEELGG = NPEEL
    ELSE
       NPEELGG = NS
    END IF
    IAT=0
    IF(IRE /= 0) THEN
       IF(KA == 0) THEN
          IAT=1
          RETURN
       END IF
    END IF
    IAT=1
    IF(NPEELGG == 1) RETURN
    IAT=0
    IF(NPEELGG /= 2) THEN
       CALL DIAGA5(NPEELGG,JA1,2*KA,IRE,IAT,RE)
       RECC=RE*RECC
       IF(IAT == 0) RETURN
       IF(JA1 == NPEELGG) RETURN
       IAT=0
    END IF
    CALL DIAGA1(JA1,2*KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    IF(NPEELGG == 2) RETURN
    ISKR=NPEELGG-JA1
    IF(JA1 == 1)ISKR=NPEELGG-1-JA1
    IF(ISKR <= 1)RETURN
    IAT=0
    CALL DIAGA3(JA1,NPEELGG,2*KA,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECOP1
!
!###################################################################
!
  SUBROUTINE RECOP1(NS,JA1,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA3,DIAGA5                     *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: NS, JA1, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IJ1, ISKR, NPEELGG
    REAL(DOUBLE) :: S, RE
!-----------------------------------------------
    IAT=1
    IJ1=JLIST(JA1)
    S=DBLE(JJQ1(3,IJ1))
    RECC=ONE/DSQRT(S)
    IF(NPEEL == 1 .AND. NS == -1)RETURN
    IF(NS == -1) THEN
       NPEELGG = NPEEL
    ELSE
       NPEELGG = NS
    END IF
    IAT=0
    IF(IRE /= 0) THEN
       IF(KA == 0) THEN
          IAT=1
          RETURN
       END IF
    END IF
    IAT=1
    IF(NPEELGG == 1) RETURN
    IAT=0
    IF(NPEELGG /= 2) THEN
       CALL DIAGA5(NPEELGG,JA1,2*KA,IRE,IAT,RE)
       RECC=RE*RECC
       IF(IAT == 0) RETURN
       IF(JA1 == NPEELGG) RETURN
       IAT=0
    END IF
    CALL DIAGA1(JA1,2*KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    IF(NPEELGG == 2) RETURN
    ISKR=NPEELGG-JA1
    IF(JA1 == 1)ISKR=NPEELGG-1-JA1
    IF(ISKR <= 1)RETURN
    IAT=0
    CALL DIAGA3(JA1,NPEELGG,2*KA,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECOP1
!
!###################################################################
!
  SUBROUTINE RECOP2(NS,JA1,JA2,K1,K2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA3,DIAGA4,DIAGA5              *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: NS, JA1, JA2, K1, K2, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ, ISKR, NPEELGG
    REAL(DOUBLE) :: S, SS, RE
!-----------------------------------------------
    IAT=1
    IJ=JLIST(JA1)
    S=DBLE(JJQ1(3,IJ))
    IJ=JLIST(JA2)
    SS=DBLE(JJQ1(3,IJ))
    SS=S*SS
    RECC=ONE/DSQRT(SS)
    IF(NPEEL == 1 .AND. NS == -1)RETURN
    IF(NS == -1) THEN
       NPEELGG = NPEEL
    ELSE
       NPEELGG = NS
    END IF
    IAT=0
    ISKR=NPEELGG-JA2
    IF(ISKR > 1) THEN
       CALL DIAGA3(JA2,NPEELGG,2*KA,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
       IAT=0
    END IF
    IF(JA2 /= NPEELGG) THEN
       CALL DIAGA5(NPEELGG,JA2,2*KA,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
       IAT=0
    ENDIF
    CALL DIAGA4(JA1,JA2,K1,K2,2*KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    IF(JA1 == 1.AND.JA2 == 2)RETURN
    IAT=0
    CALL DIAGA1(JA1,K1,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    ISKR=JA2-JA1
    IF(JA1 == 1)ISKR=JA2-1-JA1
    IF(ISKR <= 1)RETURN
    IAT=0
    CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECOP2

END MODULE RECOP
