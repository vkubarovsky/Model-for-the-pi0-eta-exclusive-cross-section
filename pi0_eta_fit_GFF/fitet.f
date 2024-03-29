      REAL FUNCTION fitet(XX)
      IMPLICIT NONE
      REAL XX
      REAL MESONMASS
      INTEGER K1,K2
      REAL Q,X,T,A,T0,E
      REAL EPS,FLUX
      INTEGER IND,INDMAX
      REAL EPSLT,XSIGMA_T,XSIGMA_TT,EPSILON
      REAL XST,XSTT,COST2
      REAL HT,ET
      COMMON/PIETA/K
      INTEGER K
      
C      INCLUDE 'pi0eta.par'

      VECTOR _Q2(*)
      VECTOR _XB(*)
      VECTOR _TT(*)

      VECTOR _ST(*)
      VECTOR _SSTT(*)

      VECTOR _ET2(*)
      VECTOR _SET2(*)

      VECTOR _HT2(*)
      VECTOR _SHT2(*)

      VECTOR AMAS(1)
      VECTOR DIM1(1)
C----
      MESONMASS=AMAS(1)
      CALL XSINIT(MESONMASS)

      IF(K.eq.1) INDMAX=96
      IF(K.eq.2) INDMAX=69
      
      IND=XX
      IF(IND.LE.0.OR.IND.GT.INDMAX) THEN
         PRINT *,'ERROR IN FITXS, X=', xx, ' IND=',IND, '  K=',k
         fitet=0.
         READ *,IND
         return
      ENDIF
      E=5.75
      Q=_Q2(IND)
      X=_XB(IND)
      T=_TT(IND)

      FITET =ET(-T,X,Q)**2
C      print *,XX,IND,K1,T,X,Q, 'FUN=',fitet

      RETURN
      END

c      INCLUDE 'dvmpx_paw.F'
