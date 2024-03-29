      REAL FUNCTION fitsf(XX)
      IMPLICIT NONE
      REAL XX
      REAL MESONMASS
      INTEGER K1,K2

      INCLUDE 'pi0eta.par'

      VECTOR _Q2(*)
      VECTOR _XB(*)
      VECTOR _TT(*)

      VECTOR _ST(*)
      VECTOR _SSTT(*)
      
      VECTOR AMAS(1)

      VECTOR DIM1(1)
      
      REAL Q,X,T,A,T0,E
      REAL EPS,FLUX
      INTEGER IND,INDMAX
      REAL EPSLT,XSIGMA_T,XSIGMA_TT,EPSILON
      REAL XST,XSTT,COST2

      MESONMASS=AMAS(1)
      CALL XSINIT(MESONMASS)

      IF(K.eq.1) INDMAX=96*2
      IF(K.eq.2) INDMAX=69*2
      
      IND=XX
      IF(IND.LE.0.OR.IND.GT.INDMAX) THEN
         PRINT *,'ERROR IN FITXS, X=', xx, ' IND=',IND, '  K=',k
         fitSF=0.
         READ *,IND
         return
      ENDIF
      E=5.75
      Q=_Q2(IND)
      X=_XB(IND)
      T=_TT(IND)

      XST =XSIGMA_T (-T,X,Q,E)
      XSTT=XSIGMA_TT(-T,X,Q,E)
      fitsf = 0.
      K1=DIM1(1)
      IF(IND.le.K1) fitsf=XST
      IF(IND.GT.K1) fitsf=XSTT

C      print *,XX,IND,K1,T,X,Q, 'FUN=',fitsf

      RETURN
      END

c      INCLUDE 'dvmpx_paw.F'
