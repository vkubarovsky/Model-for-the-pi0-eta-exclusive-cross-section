      REAL FUNCTION fitxs(xx)
      IMPLICIT NONE
      REAL XX
      PARAMETER (PI=3.1415926536)
      REAL MESONMASS

      VECTOR QV(*)
      VECTOR XV(*)
      VECTOR TV(*)
      VECTOR AV(*)
      VECTOR AMAS(1)
      VECTOR MAXIND(1)

      REAL Q,X,T,A,T0,E
      REAL EPS,FLUX
      INTEGER IND,INDMAX
      REAL EPSLT,XSIGMA_T,XSIGMA_TT,EPSILON
      REAL XST,XSTT,XSLT,COST2,COST
      REAL HELI,DVMPX,XSIGMA_LT
      MESONMASS=AMAS(1)
      CALL XSINIT(MESONMASS)
      IND=XX
      IF(IND.LE.0) THEN
         PRINT *,'ERROR IN FITXS, X=', xx, ' IND=',IND
         fitxs=0.
         READ *,IND
         return
      ENDIF
      E=5.75
      HELI=0.
      Q=QV(IND)
      X=XV(IND)
      T=TV(IND)
      A=AV(IND)

      EPS=EPSILON(x,q,E)
      EPSLT=SQRT(2*EPS*(1+EPS))
      COST2=COS(2.*A)
      COST=COS(A)
c      PRINT 200,IND,-T,X,Q,E,A
 200  format(/'=========='/' IND=',I5,' t=',f6.3,' X=',f5.3,' Q=',f5.3,
     * ' E=',f5.2,
     * ' A=',f7.4)
      XST =XSIGMA_T (-T,X,Q,E)
      XSTT=XSIGMA_TT(-T,X,Q,E)
      XSLT=XSIGMA_LT(-T,X,Q,E)
      fitxs = (XST + EPS*XSTT*COST2 +
     *   SQRT(2.*EPS*(1.+EPS))*XSLT*COST)/2./pi

c      fitxs = DVMPX(-T,X,Q,A,E,HELI,MESONMASS)/(2.*pi)

C      PRINT 100,FITXS,XST,XSTT,XSLT,EPS,A
 100  format(
     *     ' XS,XST,XSTT,XSLT=',4F10.1, ' Eps=',f10.3,' a=',F7.5)
      RETURN
      END


C      INCLUDE 'dvmpx_paw.f'
