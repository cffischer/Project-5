
SUBROUTINE RKCO_GG (JA,JB,CORD,INCOR,ICOLBREI)
!***********************************************************************
!   Configurations JA, JB. Analyse the tables of quantum numbers set   *
!   in the COMMON  blocks M0 , M1, M2, M3  to determine all possible   *
!   sets of interacting  orbitals which give a non-vanishing Coulomb   *
!   matrix element,  and  initiates the calculation of coefficients.   *
!   The following conventions are in force: (1) labels 1, 2 refer to   *
!   left, right sides of matrix element respectively;   (2) pointers   *
!   JA1, JB1, JA2, JB2 point to the JLIST array of active  orbitals;   *
!   IA1, IB1, IA2, IB2 point to the complete list of orbitals.         *
!                                                                      *
!   Call(s) to: [LIB92]: COR, CORD, ISPAR, ITJPO, SETQNA, VIJOUT.      *
!                                                                      *
!   Rewrite by G. Gaigalas                                             *
!   Transform to fortran 90/95 by G. Gaigalas           December 2012  *
!   The last modification made by G. Gaigalas           October  2017  *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
  USE vast_kind_param, ONLY:  DOUBLE
  USE m_C
  USE orb_C
  USE dumx_C
  USE CONS_C,          ONLY: ZERO, HALF, EPS
  USE trk_C
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE cord_I
  USE itjpo_I
  USE ispar_I
  USE ichkq2_I
  USE setqna_I
  USE reco_I
  USE reco2_I
  USE perko2_I
  USE itrig_I
  USE itrexg_I
  USE ixjtik_I
  USE snrc_I
  USE speak_I
  USE coulom_I
  USE ww1_I
  USE sixj_I
  USE cxk_I
  USE talk_I
  USE gg1122_I
  USE gg1112_I
  USE eile_I
  USE rec3_I
  USE jfaze_I
  USE gg1233_I
  USE reco4_I
  USE gg1234_I
  
  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      EXTERNAL CORD
  INTEGER, INTENT(IN) :: JA,JB,INCOR,ICOLBREI
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: I,IB1,IB2,IDQ,IDQG,II,I1,IJ,IJW,IM,IPCA,IT1,IT2, &
       JA1,JB1,JA2,JB2,J,JW,JT1,JT2,JT3,NDQ, &
       KLAST,KW,KWA,K1,KW1,KW2,NPEELM
!-----------------------------------------------
!
!   The Hamiltonian is an even scalar operator
  IF (ITJPO (JA) .NE. ITJPO (JB)) RETURN
  IF (ISPAR (JA) .NE. ISPAR (JB)) RETURN
  IF (ICHKQ2(JA,JB) .EQ. 0) RETURN
  CALL SETQNA (JA,JB)
!
!   1.0 Analyse peel shell interactions
!
!   1.1 Analyse electron distribution in peel. (The full procedure is
!       needed only if the number of peel orbitals NPEEL .GE. 2)
  IF (NW .LT. 1) THEN
     PRINT *, 'RKCO_GG: No subshells.'
     STOP
  ENDIF
  IF (NPEEL .EQ. 0) GOTO 48
  IF (NPEEL .EQ. 1) GOTO 43
!
!   Find differences in occupations, NDQ, for each peel orbital in
!   turn and use to set up labels of active orbitals maintaining the
!   convention JA1 .LE. JB1, JA2 .LE. JB2.
  IDQ = 0
  JA1 = 0
  JB1 = 0
  JA2 = 0
  JB2 = 0
  IDQG = 0
  DO JW = 1,NPEEL
     J = JLIST(JW)
     NDQ = NQ1(J) - NQ2(J)
     IF (IABS (NDQ) .GT. 2) RETURN
     IF (NDQ .LT. 0) THEN
        IF (NDQ+1 .GT. 0) THEN
           CYCLE
        ELSE IF (NDQ+1 .EQ. 0) THEN
           IF (JA2 .GT. 0) THEN
              JB2 = JW
           ELSE
              JA2 = JW
           END IF
           IDQ = IDQ+1
           IDQG=1+IDQG
        ELSE IF (NDQ+1 .LT. 0) THEN
           JA2 = JW
           IDQ = IDQ+2
           IDQG=20+IDQG
        END IF
     ELSE
        IF (NDQ-1 .GT. 0) THEN
           JA1 = JW
           IDQ = IDQ+2
           IDQG=20+IDQG
           CYCLE
        ELSE IF (NDQ-1 .EQ. 0) THEN
           IF (JA1 .GT. 0) THEN
              JB1 = JW
           ELSE
              JA1 = JW
           END IF
           IDQ = IDQ+1
           IDQG=1+IDQG
           CYCLE
        ELSE
           CYCLE
        END IF
     END IF
  END DO
!
!   1.2 Calculate coefficients for all possible sets of active shells.
!
!   There are 4 cases, depending on the value of IDQ, the sum of the
!   absolute differences NDQ:
!
!   1.2.1 IDQ .GT. 4: matrix element null
  IF (IDQ .GT. 4) RETURN
  IF (IDQ .EQ. 4) THEN
!
!   1.2.2 IDQ .EQ. 4: matrix element uniquely defined
     IF (JB1 .EQ. 0) THEN
        JB1 = JA1
     END IF
     IF (JB2 .EQ. 0) THEN
        JB2 = JA2
     END IF
     CONTINUE
     IF(IDQG.NE.40) THEN
        IF(IDQG.NE.22) THEN
!
!         TARP KONFIGURACIJU
!                        KAI      N'1=N1+-1
!                                 N'2=N2+-1
!                        KAI      N'3=N3+-1
!                                 N'4=N4+-1
           CALL EL5(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
        ELSE
           CALL EL4(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
        END IF
     ELSE
!
!       TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
        CALL EL2(JA,JB,JA1,JA2,ICOLBREI)
     END IF
     RETURN
  END IF
!
  IF (IDQ .NE. 2) THEN
     IF (IDQ .NE. 0) THEN
!
!         3.0 Diagnostic print - NW .LT. 1
        WRITE (*,300)
        STOP
     END IF
  ELSE
     KLAST = 1
     GOTO 16
  END IF
!
  IF (JA .EQ. JB) GOTO 43
  KLAST = NPEEL
!
!   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
!                     possible spectators.
!
!   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
!   only. Must sum over all pairs of orbitals excluding core-core
!   terms
16 DO KWA = 1,KLAST
     IF (IDQ .NE. 2) THEN
        JA1 = KWA
        JA2 = KWA
     END IF
     JT1 = JA1
     JT2 = JA2
     IT1 = JLIST(JA1)
     IT2 = JLIST(JA2)
     DO KW = KWA,NPEEL
        K1 = JLIST(KW)
        IF (NQ1(K1)*NQ2(K1) .EQ. 0) CYCLE
        JB1 = KW
        JB2 = KW
        JA1 = JT1
        JA2 = JT2
!
!   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
        IF (JA1-JB1 .GT. 0) THEN
           JT3 = JB1
           JB1 = JA1
           JA1 = JT3
        ELSE IF (JA1-JB1 .EQ. 0) THEN
           IB1 = JLIST(JB1)
           IF (NQ1(IB1) .LE. 1) CYCLE
        END IF
        IF (JA2-JB2 .GT. 0) THEN
           JT3 = JB2
           JB2 = JA2
           JA2 = JT3
        ELSE IF (JA2-JB2 .EQ. 0) THEN
           IB2 = JLIST(JB2)
           IF (NQ2(IB2) .LE. 1) CYCLE
        END IF
        IF(IDQ.NE.0) THEN
!
!     TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
           CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
        ELSE
!
!     TARP TU PACIU KONFIGURACIJU
           CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
        END IF
     END DO
     IF ((IDQ .EQ. 0) .AND. (NCORE .EQ. 0)) CYCLE
     IF ((NCORE .EQ. 0) .OR. (NAK(IT1) .NE. NAK(IT2))) RETURN
!
!   This section calculates the terms arising from active electrons
!   which are in closed shells
     NPEELM = NPEEL-1
     DO I = 1,NPEEL
        JLIS(I) = JLIST(I)
     END DO
     DO I = 1,NPEELM
        JC1S(I) = JJC1(I)
        JC2S(I) = JJC2(I)
     END DO
     DO KW = 1,NCORE
        IJW = KLIST(KW)
        DO I = 1,NPEEL
           IJ = JLIST(I)
           IF (IJW .LT. IJ) GOTO 29
        END DO
        I = NPEEL+1
        GOTO 31
29      IM = NPEEL-I+1
        DO II = 1,IM
           JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
           IF (NPEEL .EQ. II) GOTO 31
           JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
           JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
        END DO
31      CONTINUE
        IF (I .LT. 3) THEN
           I1 = JLIST(1)
           JJC1(1) = JJQ1(3,I1)
           JJC2(1) = JJQ2(3,I1)
        ELSE
           JJC1(I-1) = JJC1(I-2)
           JJC2(I-1) = JJC2(I-2)
        END IF
        JLIST(I) = IJW
        JA1 = JT1
        IF (JT1 .GE. I) JA1 = JA1+1
        JB1 = I
        JA2 = JT2
        IF (JT2 .GE. I) JA2 = JA2+1
        JB2 = I
        IF (JA1-JB1 .GT. 0) THEN
           JT3 = JB1
           JB1 = JA1
           JA1 = JT3
        END IF
        IF (JA2-JB2 .GT. 0) THEN
           JT3 = JB2
           JB2 = JA2
           JA2 = JT3
        END IF
        NPEEL = NPEEL+1
        IF(IDQ.NE.0) THEN
           IF(IDQG.NE.40) THEN
              IF(IDQG.NE.2) THEN
                 WRITE(99,995)
                 RETURN
              ELSE
!
!     TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
                 CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
              END IF
           ELSE
              WRITE(99,994)
              RETURN
           END IF
        ELSE
!
!     TARP TU PACIU KONFIGURACIJU
           CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
        END IF
        NPEEL = NPEEL-1
        NPEELM = NPEEL-1
        DO I = 1,NPEEL
           JLIST(I) = JLIS(I)
        END DO
        DO I = 1,NPEELM
           JJC1(I)  = JC1S(I)
           JJC2(I)  = JC2S(I)
        END DO
     END DO
  END DO
  RETURN
!
!   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
!         JA1 = JA2, JB1 = JB2.
43 DO KW1 = 1,NPEEL
     K1 = JLIST(KW1)
     JB1 = KW1
     JB2 = KW1
     DO KW2 = 1,KW1
        JA1 = KW2
        IF (JA1 .EQ. JB1) THEN
           IF (NQ1(K1) .LE. 1) CYCLE
        END IF
        JA2 = JA1
        IF(JA.NE.JB) THEN
           IF(IDQG.NE.2) THEN
              WRITE(99,996)
              WRITE(99,*)"JA,JB,JA1,JB1",JA,JB,JA1,JB1
              WRITE(99,*)"IDQG IDQ",IDQG,IDQ
              RETURN
           ELSE
!
!                TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
              CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
           END IF
        ELSE
!
!             TARP TU PACIU BUSENU
           CALL EL1(JA,JB,JA1,JB1,0,ICOLBREI)
        END IF
     END DO
  END DO
48 IF (INCOR .LT. 1) RETURN
  IF (NCORE .EQ. 0) RETURN
!
!   2.0 The diagonal case. deal with contributions from core orbitals
!       if INCOR .EQ. 1.
  IF(JA.NE.JB) RETURN
  DO KW1 = 1,NCORE
     JB1 = KW1
     JB2 = KW1
!
!   2.1 Calculate contribution from core/core terms
     IPCA = 2
     DO KW2 = 1,KW1
        JA1 = KW2
        JA2 = KW2
        CALL CORD (JA,JB,JA1,IPCA,JB1)
     END DO
!
!   2.2 Calculate contribution from peel/core terms
     IF (NPEEL .EQ. 0) CYCLE
     IPCA = 1
     DO KW2 = 1,NPEEL
        JA1 = KW2
        JA2 = KW2
        CALL CORD (JA,JB,JA1,IPCA,JB1)
     END DO
  END DO
  RETURN
300 FORMAT ('RKCO_GG: Error.')
994 FORMAT('   rie zymes 38?? atv N=N-N  !!!!!!')
995 FORMAT('   rie zymes 38 atv N=N-N  !!!!!!')
996 FORMAT('   rie zymes 45 atv N=N-N  !!!!!!')

CONTAINS

  SUBROUTINE EL1(JJA,JJB,JA,JB,IIRE,ICOLBREI)
!*******************************************************************
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!      SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,      *
!                         RECO,RECO2,SIXJ,SPEAK,WW1                *
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
    INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,IIRE,ICOLBREI 
!      DIMENSION CONE(7,20),S(12),IS(4),KAPS(4),KS(4)
!      DIMENSION PMGG(30),RAGG(30),J(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: II,IA,IB,IAT,IP1,IP2,IP3,IG1,IG2,IG3,IKK,I1,I2,I3,I4,& 
         IFAZ,J12,IBRD,IBRE,KRA,KRA1,L1,L2,MU,N,NU,ND1,ND2,   &
         NE1,NE2,NUP1
    INTEGER, DIMENSION(2) :: J
    INTEGER, DIMENSION(4) :: IS,KAPS,KS
    REAL(DOUBLE)          :: QM1,QM2,QM3,QM4,AA,AB,A1,BB,SI,RECC,RAG
    REAL(DOUBLE), DIMENSION(12)   :: S
    REAL(DOUBLE), DIMENSION(30)   :: PMGG,RAGG
    REAL(DOUBLE), DIMENSION(7,20) :: CONE 
!-----------------------------------------------
    IF(JA /= JB)GO TO 9
!
    IF(JJA /= JJB) RETURN
!
!     THE CASE 1111   + + - -
!
    IF(IIRE /= 0) THEN
       CALL RECO(JA,JA,JA,JA,0,IAT)
       IF(IAT /= 0)RETURN
    END IF
    CALL PERKO2(JA,JA,JA,JA,1)
    QM1=HALF
    QM2=HALF
    QM3=-HALF
    QM4=-HALF
    IA=JLIST(JA)
    J(1)=ID1(3)
    IP2=ITREXG(J(1),J(1),J(1),J(1),IKK)+1
    IF(IKK <= 0) RETURN
    IG2=IP2+IKK-1
    L1=(J(1)+1)/2
    IP1=IP2
    IG1=IG2
    IF (ICOLBREI == 2) THEN
       IS(1)=IA
       IS(2)=IA
       IS(3)=IA
       IS(4)=IA
       KAPS(1)=2*NAK(IS(1))
       KAPS(2)=2*NAK(IS(2))
       KAPS(3)=2*NAK(IS(3))
       KAPS(4)=2*NAK(IS(4))
       KS(1)=IABS(KAPS(1))
       KS(2)=IABS(KAPS(2))
       KS(3)=IABS(KAPS(3))
       KS(4)=IABS(KAPS(4))
       CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
       IF(IBRD <= 0)RETURN 
    END IF
    DO I2=IP1,IG1,2
       KRA=(I2-1)/2
       IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L1,L1,ID1(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS)CYCLE
          A1=-A1*HALF
       END IF
       AB=ZERO
       DO I3=IP2,IG2,2
          J12=(I3-1)/2
          IF(IXJTIK(J(1),J(1),KRA*2,J(1),J(1),J12*2) == 0)CYCLE
          CALL WW1(IK1,BK1,ID1,BD1,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS)CYCLE
          CALL SIXJ(J(1),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=IK1(3)+J12+KRA
          IF((IFAZ/2)*2 /= IFAZ)AA=-AA
          AB=AB+AA
       END DO
!
!     RECOUPLING COEFFICIENTS
!
       IF (ICOLBREI == 1) THEN
          BB=AB*A1
          BB=BB/DSQRT(DBLE(IK1(6)+1))
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IA,IA,KRA,BB)
       ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
             CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
             IF(DABS(S(1)) > EPS) THEN
                BB =-HALF*S(1)*AB/DSQRT(DBLE(IK1(6)+1))
                IF(DABS(BB) > EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IA,IA,4,BB)
             END IF
          END IF
       END IF
    END DO
    RETURN
!  ............................................................
9   IF(NPEEL <= 1)RETURN
    IF(IIRE /= 0) THEN
       CALL RECO(JA,JB,JB,JB,1,IAT)
       IF(IAT == 0)RETURN
    END IF
    IA=JLIST(JA)
    IB=JLIST(JB)
    QM1=HALF
    QM2=-HALF
    QM3=HALF
    QM4=-HALF
    CALL PERKO2(JA,JB,JA,JA,2)
    J(1)=ID1(3)
    J(2)=ID2(3)
    L1=(J(1)+1)/2
    L2=(J(2)+1)/2
    IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
    IF(IKK <= 0)RETURN
    IG1=IP1+IKK-1
    IP3=IP1
    IG3=IG1
    DO I4=IP1,IG1,2
       KRA=(I4-1)/2
       KRA1=KRA+1
       IF(KRA1 > 30)GO TO 10
       RAGG(KRA1)=ZERO
       PMGG(KRA1)=ZERO
       CALL RECO2(JA,JB,KRA*2,0,IAT,RECC)
       IF(IAT == 0) CYCLE
       CALL GG1122(KRA,KRA,QM1,QM2,QM3,QM4,RAG)
       IF(DABS(RAG) < EPS) CYCLE
       RAGG(KRA1)=RAG
       CALL RECO2(JA,JB,KRA*2,1,IAT,RECC)
       PMGG(KRA1)=RECC
    END DO
! * * *                      * * *                      * * *
!     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
!           2121                                1122
!
    IF (ICOLBREI == 2) THEN
       IS(1)=IA
       IS(2)=IB
       IS(3)=IA
       IS(4)=IB
       KAPS(1)=2*NAK(IS(1))
       KAPS(2)=2*NAK(IS(2))
       KAPS(3)=2*NAK(IS(3))
       KAPS(4)=2*NAK(IS(4))
       KS(1)=IABS(KAPS(1))
       KS(2)=IABS(KAPS(2))
       KS(3)=IABS(KAPS(3))
       KS(4)=IABS(KAPS(4))
       CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
       DO II=1,20
          CONE(1,II)=ZERO
          CONE(2,II)=ZERO
          CONE(3,II)=ZERO
          CONE(4,II)=ZERO
          CONE(5,II)=ZERO
          CONE(6,II)=ZERO
          CONE(7,II)=ZERO
       END DO
       IF(IBRD == 0 .AND. IBRE == 0)RETURN 
    END IF
    DO I1=IP1,IG1,2
       KRA=(I1-1)/2
       KRA1=KRA+1
       IF(KRA1 > 30)GO TO 10
       IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L2,L1,L2,ID1(5),ID2(5),ID1(5),ID2(5),KRA,AA)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA/DSQRT(DBLE(I1))
          IF(DABS(AA) > EPS) CALL SPEAK(JJA,JJB,IA,IB,IA,IB,KRA,AA)
       ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
             CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
             IF(DABS(S(1)) > EPS) THEN
                BB=S(1)*PMGG(KRA1)*RAGG(KRA1)/DSQRT(DBLE(I1))
                IF(DABS(BB) > EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IB,IB,4,BB)
             END IF
          END IF
       END IF
    END DO
! * * *                      * * *                      * * *
!     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
!           2112                                1122
!
    IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
    IF(IKK <= 0) RETURN
    IG2=IP2+IKK-1
    DO I2=IP2,IG2,2
       KRA=(I2-1)/2
       IF(KRA > 30)GO TO 10
       IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L2,L2,L1,ID1(5),ID2(5),ID2(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
       END IF
       AB=ZERO
       DO I3=IP3,IG3,2
          J12=(I3-1)/2
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) == 0)CYCLE
          CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          AB=AB+AA
       END DO
       IF (ICOLBREI == 1) THEN
          BB=A1*AB
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IB,IB,IA,KRA,BB)
       ELSE IF (ICOLBREI == 2) THEN
          NU=KRA 
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
             IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
                  (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
                IF(NU > 0) THEN
                   N=(NU-NE1)/2+1
                   CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                   DO MU = 1,3
                      CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                   END DO
                END IF
             END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
             IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
                  (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
                IF(NU >= 0) THEN
                   N=(NU-NE1)/2+1
                   IF(N <= NE2) THEN
                      CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                      DO MU = 1,3
                         CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                      END DO
                   END IF
                END IF
             END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
             IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
                  (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
                IF(NU >= 0) THEN
                   N=(NU-NE1)/2+1
                   IF(N < NE2) THEN
                      CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                      DO MU = 1,7
                         CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       END IF
    END DO
    IF (ICOLBREI == 2) THEN
       DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IB,IA,IB,IA,5,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IA,IB,IB,IA,5,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IA,IB,IA,IB,5,CONE(3,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IA,IB,IA,IB,6,CONE(4,N))
          CALL TALK(JJA,JJB,NUP1,IB,IA,IB,IA,6,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IA,IB,IB,IA,6,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IB,IA,IA,IB,6,CONE(7,N))
       END DO
    END IF
    RETURN
10  WRITE(99,100)
100 FORMAT(5X,'ERRO IN EL1  PMGG RAGG')
    STOP
  END SUBROUTINE EL1
!
!#####################################################################
!
  SUBROUTINE EL2(JJA,JJB,JA,JB,ICOLBREI)
!*******************************************************************  
!   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,RECO2,SIXJ,SPEAK                     *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT,IA,IB,IBRD,IBRE,IP1,IP2,IG1,IG2,IKK,II,I2,I3, &
                 IFAZ,IFAZ1,IFAZFRCS,J12,JAA,JBB, &
                 KRA,L1,L2,N,ND1,ND2,NE1,NE2,NUP1,NU,MU
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,A1,AB,BB,QM1,QM2,QM3,QM4,RECC,SI
      REAL(DOUBLE), DIMENSION(12)    :: S
      REAL(DOUBLE), DIMENSION(12,20) :: COND
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
! * * *                      * * *                      * * *
!     THE CASE 1122   + + - -
!
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IB
        IS(4)=IB
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
        END DO
      END IF
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L2,L2,ID1(5),ID1(5),ID2(5),ID2(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-HALF*A1
        END IF
        AB=ZERO
          DO I3=IP1,IG1,2
            J12=(I3-1)/2
            CALL RECO2(JAA,JBB,J12*2,0,IAT,RECC)
            IF(IAT /= 0) THEN
              IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) /= 0) THEN
                CALL GG1122(J12,J12,QM1,QM2,QM3,QM4,AA)
                IF(DABS(AA) > EPS) THEN
                  CALL RECO2(JAA,JBB,J12*2,1,IAT,RECC)
                  AA=AA*RECC
                  CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
                  AA=AA*SI*DSQRT(DBLE(I3))
                  IFAZ=IK1(3)+IK2(3)+KRA*2+J12*2
                  IF((IFAZ/4)*4 /= IFAZ)AA=-AA
                  AB=AB+AA
                END IF
              END IF
            END IF
          END DO
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS=1
        IFAZ1=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ1/4)*4 /= IFAZ1)IFAZFRCS=-IFAZFRCS
!
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IB,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI .EQ. 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
      RETURN
    END SUBROUTINE EL2
!
!######################################################################
!
    SUBROUTINE EL3(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!*******************************************************************
!   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     SUBROUTINE CALLED: EL31,EL32,EL33                            *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JB == JD) THEN
         IF(JA == JB.OR.JC == JB) THEN
            IF(JA == JC)GO TO 10
            IF(JC /= JB) THEN
               CALL EL32(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
            ELSE
               CALL EL31(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
            END IF
         ELSE
            CALL EL33(JJA,JJB,JC,JA,JB,1,JA,JB,JC,JD,ICOLBREI)
         END IF
         RETURN
      ELSE IF(JA == JC) THEN
         IF(JB == JA.OR.JD == JA) THEN
            IF(JB == JD)GO TO 10
            IF(JD /= JA) THEN
               CALL EL32(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
            ELSE
               CALL EL31(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
            END IF
         ELSE
            CALL EL33(JJA,JJB,JD,JB,JA,1,JA,JB,JC,JD,ICOLBREI)
         END IF
         RETURN
      ELSE IF(JA == JD) THEN
         IF(JB == JA.OR.JC == JA) THEN
            IF(JB == JC)GO TO 10
            IF(JC /= JD) THEN
               CALL EL32(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
            ELSE
               CALL EL31(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
            END IF
         ELSE
            CALL EL33(JJA,JJB,JC,JB,JA,2,JA,JB,JD,JC,ICOLBREI)
         END IF
         RETURN
      ELSE IF(JB == JC) THEN
         IF(JA == JB.OR.JD == JB) THEN
            IF(JA == JD)GO TO 10
            IF(JD /= JB) THEN
               CALL EL32(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
            ELSE
               CALL EL31(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
            END IF
         ELSE
            CALL EL33(JJA,JJB,JD,JA,JB,2,JA,JB,JD,JC,ICOLBREI)
         END IF
         RETURN
      END IF
10    WRITE(99,100)
100   FORMAT(5X,'ERRO IN EL3  PMGG RAGG')
      STOP
    END SUBROUTINE EL3
!
!###################################################################
!

      SUBROUTINE EL31(JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI)
!*******************************************************************
!   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1222,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,RECO2,SIXJ,SPEAK                     *
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
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IA,IB,II,I2,I3,IAT,IIA,IIB,IIC,IID,IKK,IP1,IG1, &
                  IBRD,IBRE,IFAZ,IFAZFRCS,INN,JAA,JBB,JB1,J12,L1, &
                  L2,KRA,ND1,ND2,NE1,NE2,N,NN
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,SI,QM1,QM2,QM3,QM4,RECC
      REAL(DOUBLE), DIMENSION(12) :: S
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      CALL RECO2(JAA,JBB,J(2),0,IAT,RECC)
      IF(IAT == 0)RETURN
      IP1=ITREXG(J(1),J(1),J(1),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      CALL RECO2(JAA,JBB,J(2),1,IAT,RECC)
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
      END IF
! * * *                      * * *                      * * *
!     CASES 2111   + + - -        TRANSFORM TO  1112   + - - +
!           1211                                1112
!
      DO I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L2,L1,L1,L1,ID2(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-A1
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(2)-J12+1
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          IF(IXJTIK(J(2),J(1),KRA*2,J(1),J(1),J12*2) == 0)CYCLE
          CALL GG1222(IK2,IK1,BK2,BK1,ID2,ID1,BD2,BD1,J12,  &
                      QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(2),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(1)+KRA*2+J12*2
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        AB=AB*RECC
        IF(DABS(AB) < EPS) CYCLE
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS = 1
        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
        NN=0
        JB1=JBB-1
        DO II=JAA,JB1
          INN=JLIST(II)
          NN=NQ1(INN)+NN
        END DO
        IF((NN/2)*2 == NN)AB=-AB
        IF (ICOLBREI == 1) THEN
           BB=A1*AB*DBLE(IFAZFRCS)
           CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
            CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
            IF(DABS(S(1)) > EPS) THEN
              BB=-S(1)*AB
              IF(DABS(BB) > EPS)                                &
              CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
            END IF
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EL31
      !
      !###
      !
      !*******************************************************************
!                                                                  *
      SUBROUTINE EL32(JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 07  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2221, 2212 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1222,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,RECO2,SIXJ,SPEAK                     *
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
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIA,IIB,IIC,IID,IATIKK,IP1,IG1,II,I2,I3,IBRD,IBRE,    &
                 IA,IB,IAT,IFAZ,IFAZFRCS,IKK,INN,JB1,JAA,JBB,J12,KRA,  &
                 L1,L2,N,NN,ND1,ND2,NE1,NE2
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,SI,RECC
      REAL(DOUBLE), DIMENSION(12) :: S
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      CALL RECO2(JAA,JBB,J(1),0,IAT,RECC)
      IF(IAT == 0)RETURN
      IP1=ITREXG(J(2),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
      END IF
      CALL RECO2(JAA,JBB,J(1),1,IAT,RECC)
! * * *                      * * *                      * * *
!     CASES 2221   + + - -        TRANSFORM TO  1222   - + + -
!           2212                                1222
!
      DO I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L2,L2,L2,L1,ID2(5),ID2(5),ID2(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-A1
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J12+1
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          IF(IXJTIK(J(2),J(2),KRA*2,J(1),J(2),J12*2) == 0) CYCLE
          CALL GG1112(IK2,IK1,BK2,BK1,ID2,ID1,BD2,     &
                      BD1,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(2),J(2),KRA*2,J(1),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(2)+KRA*2+J12*2
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        AB=AB*RECC
        IF(DABS(AB) < EPS) CYCLE
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS=1
        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
        NN=0
        JB1=JBB-1
        DO II=JAA,JB1
          INN=JLIST(II)
          NN=NQ1(INN)+NN
        END DO
        IF((NN/2)*2 == NN)AB=-AB
        IF (ICOLBREI == 1) THEN
          BB=AB*A1*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
            CALL CXK(S,IS,KAPS,KRA,KRA,IBRD,1)
            IF(DABS(S(1)) > EPS) THEN
              BB=-S(1)*AB
              IF(DABS(BB) > EPS)                              &
              CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
            END IF
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EL32
      !
      !##
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL33(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD   &
                                                        ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 08  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2313, 3231, 3213, 2331    *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,EILE,GG1233,ITREXG,IXJTIK,         *
!                      JFAZE,PERKO2,RECO,REC3,SIXJ,SPEAK           *
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
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD &
                                                        ,ICOLBREI
!      DIMENSION PMGG(30),RAGG(30),J(3)
!      DIMENSION CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,IAT,II,IIA,IIB,IIC,IID,IP1,IP2,IG1,IG2,   &
                 IKK,INN,I1,I2,I3,I4,IBRD,IBRE,IFAZ,IFAZP,IFAZFRCS, &
                 JAA,JBB,JCC,JB1,J12,KRA,KRA1,L1,L2,L3,ND1,ND2,NE1, &
                 NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(3) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: A1,AA,AB,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG,RAGG
      REAL(DOUBLE), DIMENSION(12,20) :: CONE
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      QM1=HALF
      QM2=-HALF
      QM3=HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        DO II=1,20
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        RAGG(KRA1)=ZERO
        PMGG(KRA1)=ZERO
        CALL REC3(JB,JA,JC,J(2),J(1),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL GG1233(IK2,IK1,IK3,BK2,BK1,BK3,ID2,ID1,ID3,BD2,  &
                    BD1,BD3,KRA,QM1,QM2,QM3,QM4,RAG)
        IF(DABS(RAG) < EPS) CYCLE
        RAGG(KRA1)=RAG
        CALL REC3(JB,JA,JC,J(2),J(1),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
      IFAZP=JFAZE(JB,JA,JC,JC)
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      NN=0
      JB1=JBB-1
      DO II=JAA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      IF(IREZ == 2)GO TO 5
! * * *                      * * *                      * * *
!     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
!           3231                                2133
!
    6 CONTINUE
      DO I1=IP1,IG1,2
        KRA=(I1-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L2,L3,L1,L3,ID2(5),ID3(5),ID1(5),ID3(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I1))
        AA=AA*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AA*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
          IF(DABS(S(1)) > EPS) THEN
            BB=S(1)*AA
            IF(DABS(BB) > EPS)                    &
          CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
          END IF
        END IF
      END DO
      IF(IREZ == 2) GO TO 7
! * * *                      * * *                      * * *
!     CASES 3213   + + - -        TRANSFORM TO  2133   + - + -
!           2331                                2133
!
    5 IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF(KRA > 30)GO TO 10
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L3,L2,L1,L3,ID3(5),ID2(5),ID1(5),ID3(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IID,IIC,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              IF(NU > 0) THEN
                N=(NU-NE1)/2+1
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N <= NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,4
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >=  0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJJA,JJJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJJA,JJJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      IF(IREZ == 2)GO TO 6
    7 CONTINUE
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL33  PMGG RAGG')
      STOP
      END SUBROUTINE EL33
      !
      !#
      !

      
!*******************************************************************
!                                                                  *
      SUBROUTINE EL4(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 09  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!                                                                  *
!     SUBROUTINE CALLED: EL41                                      *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 2)RETURN
      IF(JA == JB) THEN
        CALL EL41(JJA,JJB,JC,JD,JA,1,JA,JB,JC,JD,ICOLBREI)
      ELSE IF(JC == JD) THEN
        CALL EL41(JJA,JJB,JA,JB,JC,2,JA,JB,JC,JD,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL4 ')
      END SUBROUTINE EL4
      !
      !#
      !
      !*******************************************************************
!                                                                  *
      SUBROUTINE EL41(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD,    &
                                                          ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 10  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 + 1        *
!                                              N'3 = N3 - 2,       *
!     WHEN IREZ = 1   . . . . . . . . . . . . . . . . . . .        *
!                                              N'1 = N1 - 1        *
!                                              N'2 = N2 - 1        *
!                                              N'3 = N3 + 2,       *
!     WHEN IREZ = 2   . . . . . . . . . . . . . . . . . . .        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,EILE,GG1233,ITREXG,IXJTIK,         *
!                        JFAZE,PERKO2,RECO,REC3,SIXJ,SPEAK         *
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
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB, &
                             JJC,JJD,ICOLBREI
!      DIMENSION J(3)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IA,IB,IC,II,IIA,IIB,IIC,IBRD,IBRE,IP1,IG1,IP2,IG2, &
                  IAT,IID,IKK,IFAZ,IFAZP,IFAZFRCS,INN,I2,I3,JB1,JAA, &
                  JBB,JCC,J12,KRA,L1,L2,L3,ND1,ND2,NE1,NE2,N,NN,NU,  &
                  NUP1,MU
      INTEGER, DIMENSION(3) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,SI,RECC
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(12,20) :: COND
!-----------------------------------------------
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      IF(NPEEL <= 1)RETURN
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(IREZ == 1) THEN
        QM1=-HALF
        QM2=-HALF
        QM3=HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=HALF
        QM3=-HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
        END DO
      END IF
      IFAZP=JFAZE(JC,JA,JB,JC)
      IFAZFRCS = 1
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)+ &
      IK3(5)*IK3(4)-ID3(5)*ID3(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      NN=0
      JB1=JBB-1
      DO II=JAA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
!           3321                                1233
!                                                    (IREZ = 1)
!     OR
!     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
!           2133                                1233
!                                                    (IREZ = 2)
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L1,L2,L3,L3,ID1(5),ID2(5),ID3(5),ID3(5),KRA,A1)
          ELSE
            CALL COULOM(L3,L3,L1,L2,ID3(5),ID3(5),ID1(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(2)-J12+1
          IF(IREZ == 2)IFAZ=J(1)-J12+1
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          CALL REC3(JA,JB,JC,J(1),J(2),J12*2,0,IAT,AA)
          IF(IAT == 0) CYCLE
          IF(IXJTIK(J(3),J(1),KRA*2,J(2),J(3),J12*2) == 0) CYCLE
          CALL GG1233(IK1,IK2,IK3,BK1,BK2,BK3,ID1,ID2,ID3,BD1,    &
                    BD2,BD3,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS) CYCLE
          CALL REC3(JA,JB,JC,J(1),J(2),J12*2,1,IAT,RECC)
          AA=AA*RECC
          CALL SIXJ(J(3),J(1),KRA*2,J(2),J(3),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(3)+J(1)+2*J12+2*KRA
          IF(IREZ == 2)IFAZ=J(2)+J(3)+2*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL41  PMGG RAGG')
      STOP
      END SUBROUTINE EL41
      !
      !#
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL5(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!      SUBROUTINE CALLED: EL51,EL52,EL53                           *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      IF(JB < JC) THEN
        CALL EL51(JJA,JJB,JA,JB,JC,JD,1,ICOLBREI)
      ELSE IF(JA > JD.AND.JB > JD) THEN
        CALL EL51(JJA,JJB,JC,JD,JA,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA < JC) THEN
        CALL EL52(JJA,JJB,JA,JC,JB,JD,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA > JC) THEN
        CALL EL52(JJA,JJB,JC,JA,JD,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA < JC) THEN
        CALL EL53(JJA,JJB,JA,JC,JD,JB,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA > JC) THEN
        CALL EL53(JJA,JJB,JC,JA,JB,JD,2,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL5 ')
      END SUBROUTINE EL5
      !
      !#
      !
      !*******************************************************************
!                                                                  *
      SUBROUTINE EL51(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 12  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1234, 2134, 1243, 2134    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 + 1   *
!     AND    3412, 4321, 3421, 4312                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,REC4,SIXJ,SPEAK                      *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG
      REAL(DOUBLE), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=-HALF
        QM3=HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=HALF
        QM3=-HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=IB
          IS(3)=IC
          IS(4)=ID
        ELSE IF(IREZ == 2) THEN
          IS(1)=IC
          IS(2)=ID
          IS(3)=IA
          IS(4)=IB
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,   &
                  ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)  &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1234   + + - -
!           2134                  TRANSFORM TO  1234   + + - -
!                                                    (IREZ = 1)
!     OR
!     CASES 3412   + + - -        TRANSFORM TO  1234   - - + +
!           3421                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(2),J(4),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L3,L4,L1,L2,ID3(5),ID4(5),ID1(5),ID2(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L2,L3,L4,ID1(5),ID2(5),ID3(5),ID4(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)+J(4)-2*J12
          IF(IREZ == 2)IFAZ=J(2)+J(3)-2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2) == 0) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(2)+J(3)+2*J12+2*KRA
          IF(IREZ == 2)IFAZ=J(1)+J(4)-2*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IB,IC,ID,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IC,ID,IA,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
! * * *                      * * *                      * * *
!     CASES 1243   + + - -        TRANSFORM TO  1234   + + - -
!           2134                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 3421   + + - -        TRANSFORM TO  1234   - - + +
!           4321                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(4),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L3,L4,L2,L1,ID3(5),ID4(5),ID2(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L2,L4,L3,ID1(5),ID2(5),ID4(5),ID3(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J(4)
          IF(IREZ == 2)IFAZ=J(3)-J(2)
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(4),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(4),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(2)+J(3)+2*J(4)+2*KRA
          IF(IREZ == 2)IFAZ=J(1)+2*J(2)+J(4)+4*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IB,ID,IC,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IC,ID,IB,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL51  PMGG RAGG')
      STOP
      END SUBROUTINE EL51
      !
      !#
      !
      !*******************************************************************
!                                                                  *
      SUBROUTINE EL52(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 13  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1324, 3142, 1342, 3124    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!     AND    2413, 4231, 2431, 4213                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,REC4,SIXJ,SPEAK                      *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!     DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG
      REAL(DOUBLE), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=HALF
        QM3=-HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=-HALF
        QM3=HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=IC
          IS(3)=IB
          IS(4)=ID
        ELSE IF(IREZ == 2) THEN
          IS(1)=IB
          IS(2)=ID
          IS(3)=IA
          IS(4)=IC
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)  &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1324   + + - -        TRANSFORM TO  1234   + - + -
!           1342                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2413   + + - -        TRANSFORM TO  1234   - + - +
!           4231                                1234
!                                                    (IREZ = 2)
      DO I3=IP1,IG1,2
        KRA=(I3-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L4,L1,L3,ID2(5),ID4(5),ID1(5),ID3(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L3,L2,L4,ID1(5),ID3(5),ID2(5),ID4(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAG
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I3))
        AB=AA*DBLE(IFAZP)
        IF(IREZ == 2) THEN
          IFAZ=J(1)+J(2)+J(4)+J(3)+4*KRA
          IF((IFAZ/4)*4 /= IFAZ)AB=-AB
        END IF
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IC,IB,ID,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,ID,IA,IC,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU >  0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
! * * *                      * * *                      * * *
!     CASES 1342   + + - -        TRANSFORM TO  1234   + - + -
!           3124                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2431   + + - -        TRANSFORM TO  1234   - + - +
!           4213                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(4),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L4,L3,L1,ID2(5),ID4(5),ID3(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L3,L4,L2,ID1(5),ID3(5),ID4(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J(3)+2*J12
          IF(IREZ == 2)IFAZ=J(4)-J(2)-2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(4),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(4),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(3)-4*KRA-4*J12
          IF(IREZ == 2)IFAZ=J(1)+J(2)+J(3)-J(4)-4*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IC,ID,IB,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,ID,IC,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL52  PMGG RAGG')
      STOP
      END SUBROUTINE EL52
      !
      !#
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL53(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 14  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1423, 4132, 1432, 4123    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!     AND    2314, 3241, 2341, 3214                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                      RECO,REC4,SIXJ,SPEAK                        *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG
      REAL(DOUBLE), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=HALF
        QM3=HALF
        QM4=-HALF
      ELSE
        QM1=HALF
        QM2=-HALF
        QM3=-HALF
        QM4=HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=ID
          IS(3)=IB
          IS(4)=IC
        ELSE IF(IREZ == 2) THEN
          IS(1)=IB
          IS(2)=IC
          IS(3)=IA
          IS(4)=ID
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4) &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1423   + + - -        TRANSFORM TO  1234   + - - +
!           4132                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2314   + + - -        TRANSFORM TO  1234   - + + -
!           3241                                1234
!                                                    (IREZ = 2)
      DO I3=IP1,IG1,2
        KRA=(I3-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L3,L1,L4,ID2(5),ID3(5),ID1(5),ID4(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L4,L2,L3,ID1(5),ID4(5),ID2(5),ID3(5),KRA,A1)
          END IF
        IF(DABS(A1) < EPS) CYCLE
        END IF
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAG
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I3))
        IFAZ=J(4)+J(3)-2*KRA+2
        IF(IREZ == 2)IFAZ=J(1)+J(2)-2*KRA+2
        IF((IFAZ/4)*4 /= IFAZ)AA=-AA
        AB=AA*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,ID,IB,IC,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,IC,IA,ID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND. &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
! * * *                      * * *                      * * *
!     CASES 1432   + + - -        TRANSFORM TO  1234   + - - +
!           4132                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2341   + + - -        TRANSFORM TO  1234   - + + -
!           3214                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(4),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L3,L4,L1,ID2(5),ID3(5),ID4(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L4,L3,L2,ID1(5),ID4(5),ID3(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)+J(4)+2*J12
          IF(IREZ == 2)IFAZ=J(2)+J(3)+2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(3)-J(4)-4*KRA+2*J12
          IF(IREZ == 2)IFAZ=J(1)-J(2)-4*KRA-2*J12
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,ID,IC,IB,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,IC,ID,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)KRA1
  100 FORMAT(5X,'ERRO IN EL53  PMGG RAGG KRA1=',I100)
      STOP
      END SUBROUTINE EL53
      
END SUBROUTINE RKCO_GG
