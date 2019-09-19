MODULE diaga
!*******************************************************************
!                                                                  *  
!   This MODULE groups together the subroutines in the series      *
!   diaga1 - diaga5, rather than relying on calls to the           *
!   subroutines in separate files.                                 *
!                                                                  *   
!                                                                  *
!   SUBROUTINES written by  G. Gaigalas                            *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   MODULE created by A. Senchuk                   September 2019  *  
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
  USE ixjtik_I
  USE sixj_I
  USE nine_I
      
  IMPLICIT NONE

CONTAINS

  SUBROUTINE DIAGA1(JA1,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ1,IA1,IB1,K1,J1,IT1,IT1S,IFAZ,LL1
    REAL(DOUBLE) :: A1
!-----------------------------------------------
    RECC = ZERO
    IJ1=JLIST(JA1)
    IA1=JJQ1(3,IJ1)-1
    IB1=JJQ2(3,IJ1)-1
    IF(JA1 <= 2) THEN
       LL1=JLIST(1)
       IF(JA1 == 1)LL1=JLIST(2)
       J1=JJQ1(3,LL1)-1
       IT1=JJC1(1)-1
       IT1S=JJC2(1)-1
    ELSE
       K1=JA1-2
       J1=JJC1(K1)-1
       K1=JA1-1
       IT1=JJC1(K1)-1
       IT1S=JJC2(K1)-1
    END IF
    IF(IRE /= 0) THEN
       CALL SIXJ(KA,IB1,IA1,J1,IT1,IT1S,0,A1)
       A1=A1*DSQRT(DBLE((IA1+1)*(IT1S+1)))
       IFAZ=J1+IT1+IB1+KA
       IF((IFAZ/4)*4 /= IFAZ)A1=-A1
       RECC=A1
       IAT=1
       IF(JA1 /= 1)RETURN
       IFAZ=IA1+IB1+2*J1-IT1-IT1S
       IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
    ELSE
       IF(IXJTIK(KA,IB1,IA1,J1,IT1,IT1S) == 0)RETURN
       IAT=1
    END IF
    RETURN
  END SUBROUTINE DIAGA1
!
!###################################################################
!
  SUBROUTINE DIAGA2(JA1,JA2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1, JA2, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ1,IJ2,IA1,IA2,IB1,IB2,N1,N2,J2,IT2,IT2S,IFAZ
    REAL(DOUBLE) :: A2
!-----------------------------------------------
    RECC = ZERO
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IA1=JJQ1(3,IJ1)-1
    IA2=JJQ1(3,IJ2)-1
    IB1=JJQ2(3,IJ1)-1
    IB2=JJQ2(3,IJ2)-1
    IF(JA1 == 1.AND.JA2 == 2) THEN
       IT2=IA1
       IT2S=IB1
       J2=JJC1(1)-1
    ELSE
       N1=JA2-1
       J2=JJC1(N1)-1
       N2=JA2-2
       IT2=JJC1(N2)-1
       IT2S=JJC2(N2)-1
    END IF
    IF(IRE /= 0) THEN
       CALL SIXJ(KA,IB2,IA2,J2,IT2,IT2S,0,A2)
       RECC=A2*DSQRT(DBLE((IB2+1)*(IT2+1)))
       IFAZ=J2+IT2S+IA2+KA
       IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
       IAT=1
    ELSE
       IF(IXJTIK(KA,IB2,IA2,J2,IT2,IT2S) == 0)RETURN
       IAT=1
    END IF
    RETURN
  END SUBROUTINE DIAGA2
!
!###################################################################
!
  SUBROUTINE DIAGA3(JA1,JA2,KA,IRE,IAT,REC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)  :: JA1, JA2, KA, IRE
    INTEGER, INTENT(OUT) :: IAT
    REAL(DOUBLE), INTENT(OUT) :: REC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: I,LL1,JI,KK1,KK2,ITI,ITI1,ITIS,ITI1S,IFAZ
    REAL(DOUBLE) :: AA, A3
!-----------------------------------------------
    REC = ZERO
    AA=ONE
    I=JA1+1
    IF(JA1 == 1)I=I+1
    IF(I >= JA2) THEN
       REC=AA
       IAT=1
    ELSE
1      LL1=JLIST(I)
       JI=JJQ1(3,LL1)-1
       KK2=I-2
       ITI=JJC1(KK2)-1
       ITIS=JJC2(KK2)-1
       KK1=I-1
       ITI1=JJC1(KK1)-1
       ITI1S=JJC2(KK1)-1
       IF(IRE /= 0) THEN
          CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
          A3=A3*SQRT(DBLE((ITI+1)*(ITI1S+1)))
          IFAZ=KA+JI+ITI+ITI1S
          IF((IFAZ/4)*4 /= IFAZ)A3=-A3
          AA=AA*A3
       ELSE
          IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) == 0)RETURN
       END IF
       I=I+1
       IF(I == JA2) THEN
          REC=AA
          IAT=1
       ELSE
          GOTO 1
       END IF
    END IF
    RETURN
  END SUBROUTINE DIAGA3
! 
!###################################################################
!
  SUBROUTINE DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  NINE                                     *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1,JA2,K1,K2,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ1,IJ2,IA1,IA2,IB1,IB2,N1,N2,J2,J2S,IT2,IT2S
    REAL(DOUBLE) :: A2
!-----------------------------------------------
    RECC = ZERO
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IA1=JJQ1(3,IJ1)-1
    IA2=JJQ1(3,IJ2)-1
    IB1=JJQ2(3,IJ1)-1
    IB2=JJQ2(3,IJ2)-1
    IF(JA1 == 1.AND.JA2 == 2) THEN
       IT2=IA1
       IT2S=IB1
       J2=JJC1(1)-1
       J2S=JJC2(1)-1
    ELSE
       N1=JA2-1
       J2=JJC1(N1)-1
       J2S=JJC2(N1)-1
       N2=JA2-2
       IT2=JJC1(N2)-1
       IT2S=JJC2(N2)-1
    END IF
    IF(IRE /= 0) THEN
       CALL NINE(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,0,IAT,A2)
       RECC=A2*DSQRT(DBLE((IT2+1)*(KA+1)*(IA2+1)*(J2S+1)))
    ELSE
       CALL NINE(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,1,IAT,A2)
    END IF
    RETURN
  END SUBROUTINE DIAGA4
!
!###################################################################
!
  SUBROUTINE DIAGA5(NPEELGG,JA1,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: NPEELGG,JA1,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: ITI1,ITI1S,IJ1,ITI,ITIS,JI
    REAL(DOUBLE) :: A3
!-----------------------------------------------
    RECC = ZERO
    ITI1=JJC1(NPEELGG-1)-1
    ITI1S=JJC2(NPEELGG-1)-1
    IJ1=JLIST(NPEELGG)
    IF(JA1 == NPEELGG) THEN
       ITI=JJQ1(3,IJ1)-1
       ITIS=JJQ2(3,IJ1)-1
       JI=JJC1(NPEELGG-2)-1
    ELSE
       JI=JJQ1(3,IJ1)-1
       ITI=JJC1(NPEELGG-2)-1
       ITIS=JJC2(NPEELGG-2)-1
    END IF
    IF(IRE == 0) THEN
       IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) /= 0) IAT=1
    ELSE
       CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
       RECC=A3*DSQRT(DBLE((ITI+1)*(ITI1S+1)))
       IF(MOD(KA+JI+ITIS+ITI1,4) /= 0)RECC=-RECC
       IAT=1
       IF(JA1 == NPEELGG)RETURN
       IF(MOD(ITI+ITIS-ITI1S-ITI1+2*JI,4) /= 0)RECC=-RECC
    END IF
    RETURN
  END SUBROUTINE DIAGA5

END MODULE diaga
