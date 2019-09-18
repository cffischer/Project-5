SUBROUTINE PERKO2(JA1,JA2,JA3,JA4,I)
!*******************************************************************
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (GENERAL CASE)     *
!                                                                  *
!     SUBROUTINE CALLED: PERKO1                                    *
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
  USE trk_C
  USE CONS_C
  USE m_C
  USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
  USE nmtejj_I
  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN) :: JA1, JA2, JA3, JA4, I
!-----------------------------------------------
  CALL PERKO1(JA1,BK1,IK1,BD1,ID1)
  IF(I == 1)RETURN
  CALL PERKO1(JA2,BK2,IK2,BD2,ID2)
  IF(I == 2)RETURN
  CALL PERKO1(JA3,BK3,IK3,BD3,ID3)
  IF(I == 3)RETURN
  CALL PERKO1(JA4,BK4,IK4,BD4,ID4)
  RETURN

CONTAINS

  SUBROUTINE PERKO1(JA,BK,IK,BD,ID)
!*******************************************************************
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (FOR ONE SHELL)    *
!     CALLS TO: NMTEJJ                                             *    
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER,      INTENT(IN)                :: JA
    INTEGER,      INTENT(OUT), DIMENSION(7) :: IK, ID
    REAL(DOUBLE), INTENT(OUT), DIMENSION(3) :: BK, BD
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ
!-----------------------------------------------
    IJ=JLIST(JA)
    IK(2)=NP(IJ)
    ID(2)=IK(2)
    IK(3)=(IABS(NAK(IJ))*2)-1
    ID(3)=IK(3)
    IK(4)=NQ1(IJ)
    ID(4)=NQ2(IJ)
    IK(5)=(IK(3)+NAK(IJ)/IABS(NAK(IJ)))/2
    ID(5)=IK(5)
    IK(6)=JJQ1(3,IJ)-1
    ID(6)=JJQ2(3,IJ)-1
    IK(7)=IABS(NAK(IJ))-JJQ1(1,IJ)
    ID(7)=IABS(NAK(IJ))-JJQ2(1,IJ)
    BK(1)=HALF*DBLE(IK(7))
    BD(1)=HALF*DBLE(ID(7))
    BK(2)=HALF*DBLE(IK(6))
    BD(2)=HALF*DBLE(ID(6))
    BK(3)=-HALF*DBLE(IABS(NAK(IJ))-IK(4))
    BD(3)=-HALF*DBLE(IABS(NAK(IJ))-ID(4))
    IK(1)=NMTEJJ(IK(7),IK(6),IK(3),ID(4),IK(4))
    ID(1)=NMTEJJ(ID(7),ID(6),ID(3),ID(4),IK(4))
    RETURN
  END SUBROUTINE PERKO1
  
END SUBROUTINE PERKO2
