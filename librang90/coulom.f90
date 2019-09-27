SUBROUTINE COULOM(J1,J2,J3,J4,L1,L2,L3,L4,K,AA)
!*******************************************************************
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
!                                                                  *
!                          k   k+1  (k) (k)                        *
!     (n l j T  n l j T ::r  / r  ( C   C )::n l j T  n l j T )    *
!       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
!                                                                  *
!     SUBROUTINE CALLED:  CRE                                      *
!                                                                  *
!   Written by G. Gaigalas,                                        *
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
  USE clrx_I 
  
  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN)        :: J1,J2,J3,J4,L1,L2,L3,L4,K
  REAL(DOUBLE), INTENT(OUT)  :: AA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: I, IFAZ
!-----------------------------------------------
  AA=ZERO
  IF(ITTK(L1,L3,K) == 0)RETURN
  IF(ITTK(L2,L4,K) == 0)RETURN
  I=(2*K+1)/2
  AA=CRE (J1,I,J3)
  IF(DABS(AA) < EPS)RETURN
  AA=AA*CRE (J2,I,J4)
  IF(DABS(AA) < EPS)RETURN
  IFAZ=L3-2*K-L1+L4-L2
  IF((IFAZ/4)*4 /= IFAZ)AA=-AA
  RETURN

CONTAINS

  REAL(KIND(0.0D0)) FUNCTION CRE (KAP1, K, KAP2) 
!***********************************************************************
!                                                                      *
!   Computes the relativistic reduced matrix element                   *
!                                                                      *
!                         (j1 || C(K) || j2),                          *
!                                                                      *
!   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
!   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
!   conditions are tested by the routine CLRX.                         *
!                                                                      *
!   Call(s) to: [LIB92] CLRX.                                          *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:10   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER  :: KAP1 
    INTEGER  :: K 
    INTEGER  :: KAP2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: K1 
    REAL(DOUBLE) :: DK1K2 
!-----------------------------------------------
!
    K1 = ABS(KAP1) 
    DK1K2 = DBLE(4*K1*IABS(KAP2)) 
    CRE = SQRT(DK1K2)*CLRX(KAP1,K,KAP2) 
    IF (MOD(K1,2) == 1) CRE = -CRE 
!
    RETURN  

  END FUNCTION CRE
!
!#######################################################################
!
!!$  REAL(KIND(0.0D0)) FUNCTION CLRX (KAPPAA, K, KAPPAB) 
!!$!***********************************************************************
!!$!                                                                      *
!!$!   The value of CLRX is the 3-j symbol:                               *
!!$!                                                                      *
!!$!                    ( JA        K        JB  )                        *
!!$!                    ( 1/2       0       -1/2 )                        *
!!$!                                                                      *
!!$!   The  K'S are kappa angular quantum numbers. The formula is taken   *
!!$!   from D M Brink and G R Satchler, <Angular Momentum>, second edi-   *
!!$!   tion (Oxford: Clarendon press, 1968), p 138.   The logarithms of   *
!!$!   the first  MFACT  factorials must be available in  COMMON/FACTS/   *
!!$!   for this program to function correctly. Note that  N!  is stored   *
!!$!   in FACT(N+1)                                                       *
!!$!                                                                      *
!!$!   No subroutines called.                                             *
!!$!                                                                      *
!!$!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!!$!                                                                      *
!!$!***********************************************************************
!!$!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:52   2/14/04  
!!$!...Modified by Charlotte Froese Fischer 
!!$!                     Gediminas Gaigalas  10/05/17
!!$!-----------------------------------------------
!!$!   M o d u l e s 
!!$!-----------------------------------------------
!!$    USE vast_kind_param, ONLY:  DOUBLE 
!!$    USE FACTS_C 
!!$    IMPLICIT NONE
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$    INTEGER, INTENT(IN) :: KAPPAA, K, KAPPAB 
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$    INTEGER :: KA, KB, KAPKB, KABKP, KAMKB, KBMKA 
!!$    REAL(DOUBLE) :: EXPTRM 
!!$!-----------------------------------------------
!!$!
!!$!
!!$!   Determine the absolute values of the kappas
!!$!
!!$    KA = ABS(KAPPAA) 
!!$    KB = ABS(KAPPAB) 
!!$!
!!$!   Perform the triangularity check
!!$!
!!$    IF (ABS(KA - KB)<=K .AND. KA+KB-1>=K) THEN 
!!$!
!!$!   Triangularity satisfied; compute the 3j coefficient
!!$!
!!$!   Begin with the logarithm of the square of the leading term
!!$!
!!$       EXPTRM = -LOG(DBLE(KA*KB)) 
!!$!
!!$!   Compute the logarithm of the square root of the leading term
!!$!   and the factorial part that doesn't depend on the parity of
!!$!   KA+KB+K (the delta factor)
!!$!
!!$       KAPKB = KA + KB 
!!$       KABKP = KAPKB + K 
!!$       KAMKB = KA - KB 
!!$       KBMKA = KB - KA 
!!$       EXPTRM = 0.5D00*(EXPTRM + GAM(KAPKB-K)+GAM(KAMKB+K+1)+GAM(KBMKA+K+1)-&
!!$            GAM(KABKP+1)) 
!!$!
!!$!   The remainder depends on the parity of KA+KB+K
!!$!
!!$       IF (MOD(KABKP,2) == 0) THEN 
!!$!
!!$!   Computation for even parity case
!!$!
!!$!   Include the phase factor: a minus sign if necessary
!!$!
!!$          IF (MOD(3*KABKP/2,2) == 0) THEN 
!!$             CLRX = 1.0D00 
!!$          ELSE 
!!$             CLRX = -1.0D00 
!!$          ENDIF
!!$!
!!$!   Include the contribution from the factorials
!!$!
!!$          EXPTRM = EXPTRM + GAM((KABKP+2)/2) - GAM((KAPKB-K)/2) - GAM((KAMKB+&
!!$               K+2)/2) - GAM((KBMKA+K+2)/2) 
!!$!
!!$       ELSE 
!!$!
!!$!   Computation for odd parity case
!!$!
!!$!   Include the phase factor: a minus sign if necessary
!!$!
!!$          IF (MOD((3*KABKP - 1)/2,2) == 0) THEN 
!!$             CLRX = 1.0D00 
!!$          ELSE 
!!$             CLRX = -1.0D00 
!!$          ENDIF
!!$!
!!$!   Include the contribution from the factorials
!!$!
!!$          EXPTRM = EXPTRM + GAM((KABKP+1)/2) - GAM((KAPKB-K+1)/2) - GAM((&
!!$               KAMKB+K+1)/2) - GAM((KBMKA+K+1)/2) 
!!$!
!!$       ENDIF
!!$!
!!$!   Final assembly
!!$!
!!$       CLRX = CLRX*EXP(EXPTRM) 
!!$!
!!$    ELSE 
!!$!
!!$!   Triangularity violated; set the coefficient to zero
!!$!
!!$       CLRX = 0.0D00 
!!$!
!!$    ENDIF
!!$!
!!$    RETURN  
!!$!
!!$  END FUNCTION CLRX

END SUBROUTINE COULOM
