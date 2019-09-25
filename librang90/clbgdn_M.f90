MODULE clebsch_gordan
!*******************************************************************
!                                                                  *  
!   This MODULE groups together the subroutines to calculate the   *
!   the Clebsch-Gordan coefficients in various circumstances,      *
!   rather than relying on calls to the subroutines in separate    *
!   files.                                                         *
!                                                                  *   
!                                                                  *
!   SUBROUTINES written by  G. Gaigalas             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   MODULE created by A. Senchuk                   September 2019  *  
!                                                                  *
!*******************************************************************

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
  USE vast_kind_param
  USE CONS_C

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
  USE ittk_I
      
  IMPLICIT NONE


CONTAINS
  
  SUBROUTINE C0T5S(Q,QM,SM,C,CM,A)
!*******************************************************************
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                ---          ---  *
!                                                I  Q  1/2  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:              I              I  *
!                                                I  QM  SM  CM  I  *
!                                                ---          ---  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    REAL(DOUBLE), INTENT(IN)  :: Q, QM, SM, C, CM
    REAL(DOUBLE), INTENT(OUT) :: A
!      DIMENSION GC(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER                    :: IIQ, IIC, IE
    REAL(DOUBLE), DIMENSION(2) :: GC
!-----------------------------------------------
    GC(1)=ONE
    GC(2)=-ONE
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IF(ITTK(IIQ,IIC,1) == 0)RETURN
    IF(DABS(QM+SM-CM) > EPS)RETURN
    IF((HALF+TENTH) < DABS(SM))RETURN
    IF((Q+TENTH) < DABS(QM))RETURN
    IF((C+TENTH) < DABS(CM))RETURN
    IE=DABS(HALF-SM)+ONE+TENTH
    IF(DABS(Q+HALF-C) < EPS) THEN
       A=DSQRT((C+GC(IE)*CM)/(TWO*C))
    ELSE
       IF(DABS(Q-HALF-C) > EPS)RETURN
       A=-GC(IE)*DSQRT((C-GC(IE)*CM+ONE)/(TWO*C+TWO))
    ENDIF
    RETURN
  END SUBROUTINE C0T5S
!
!####################################################################
!
  SUBROUTINE CLE0SM(Q,QM,S,C,CM,A)
!*******************************************************************
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   S  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    REAL(DOUBLE), INTENT(IN)  :: Q, QM, S, C, CM
    REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IIQ, IIC, IIS
!-----------------------------------------------
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IIS=TWO*S+TENTH
    IF(ITTK(IIQ,IIC,IIS).EQ.0)RETURN
    IF(S.LT.EPS) THEN
       IF((Q+TENTH).LT.DABS(QM))RETURN
       IF((C+TENTH).LT.DABS(CM))RETURN
       IF(DABS(Q-C).GT.EPS)RETURN
       IF(DABS(QM-CM).GT.EPS)RETURN
       A=ONE
    ELSE
       CALL C1E0SM(Q,QM,C,CM,A)
    END IF
    RETURN
  END SUBROUTINE CLE0SM
!
!###################################################################
!
  SUBROUTINE C1E0SM(Q,QM,C,CM,A)
!*******************************************************************
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    REAL(DOUBLE), INTENT(IN)  :: Q, QM, C, CM
    REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IIQ, IIC, IS, IG
!-----------------------------------------------
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IF(ITTK(IIQ,IIC,2) == 0)RETURN
    IF(DABS(QM-CM) > EPS) RETURN 
    IF((Q+TENTH) < DABS(QM)) RETURN
    IF((C+TENTH) < DABS(CM)) RETURN
    IF(DABS(QM) <= EPS) THEN
       IS=Q+C+ONE+TENTH
       IF((IS/2)*2 /= IS) RETURN
    END IF
    IG=Q-C+TWO+TENTH
    IF(IG <= 0) RETURN
    IF(IG > 3) RETURN
    IF (IG == 1) THEN
       A=DSQRT(((C+CM)*(C-CM))/((TWO*C-ONE)*C))
    ELSE IF (IG == 2) THEN
       A=CM/DSQRT(C*(C+ONE))
    ELSE IF (IG == 3) THEN
       A=-DSQRT(((C+CM+ONE)*(C-CM+ONE))/((C+ONE)*(TWO*C+THREE)))
    END IF
    RETURN
  END SUBROUTINE C1E0SM
!
!###################################################################
!
  SUBROUTINE C1E1SM(Q,QM,SM,C,CM,A)
!*******************************************************************
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  1  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    REAL(DOUBLE), INTENT(IN)  :: Q, QM, SM, C, CM
    REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER                    :: IE, IIQ, IIC
    REAL(DOUBLE), DIMENSION(2) :: GC
!-----------------------------------------------
    GC(1)=ONE
    GC(2)=-ONE
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IF(ITTK(IIQ,IIC,2).EQ.0)RETURN
    IF(DABS(QM+SM-CM).GT.EPS)RETURN
    IF((Q+TENTH).LT.DABS(QM))RETURN
    IF((C+TENTH).LT.DABS(CM))RETURN
    IE=0
    IF(DABS(SM-ONE).LT.EPS)IE=1
    IF(DABS(SM+ONE).LT.EPS)IE=2
    IF(IE.EQ.0)RETURN
    IF(DABS(Q+ONE-C).LT.EPS) THEN
       A=DSQRT((C+GC(IE)*CM-ONE)*(C+GC(IE)*CM)/ &
            ((TWO*C-ONE)*TWO*C))
    ELSE IF(DABS(Q-C).LT.EPS) THEN
       A=-GC(IE)*DSQRT((C+GC(IE)*CM)*(C-GC(IE)*CM+ONE)/ &
            ((C+ONE)*TWO*C))
    ELSE IF(DABS(Q-ONE-C).GT.EPS) THEN
       RETURN
    ELSE
       A=DSQRT((C-GC(IE)*CM+ONE)*(C-GC(IE)*CM+TWO)/ &
            ((TWO*C+TWO)*(TWO*C+THREE)))
    END IF
    RETURN
  END SUBROUTINE C1E1SM

END MODULE clebsch_gordan
