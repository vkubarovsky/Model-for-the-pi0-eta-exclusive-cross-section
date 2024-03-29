      REAL FUNCTION fitxs_test(xx)
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
      REAL XST,XSTT,COST2

      MESONMASS=AMAS(1)
      CALL XSINIT(MESONMASS)
      IND=XX
      IF(IND.LE.0) THEN
         PRINT *,'ERROR IN FITXS, X=', xx, ' IND=',IND
         fitxs_test=0.
         READ *,IND
         return
      ENDIF
      E=5.75
      Q=QV(IND)
      X=XV(IND)
      T=TV(IND)
      A=AV(IND)

      EPS=EPSILON(x,q,E)
      EPSLT=SQRT(2*EPS*(1+EPS))
      COST2=COS(2.*A)
c      PRINT 200,IND,-T,X,Q,E,A
 200  format(/'=========='/' IND=',I5,' t=',f6.3,' X=',f5.3,' Q=',f5.3,
     * ' E=',f5.2,
     * ' A=',f7.4)
      XST =XSIGMA_T (-T,X,Q,E)
      XSTT=XSIGMA_TT(-T,X,Q,E)
      fitxs_test = (XST + EPS * XSTT * COST2)/2./pi
c      PRINT 100,FITXS_test,XST,XSTT,EPS,COST2,P
 100  format(
     *     ' XS,XST,XSTT=',3F11.3, ' Eps=',f10.3,' COST2=',F7.5/
     *     ' HT=',5f8.2,'  ET=',5F8.2) 
      RETURN
      END


      INCLUDE 'dvmpx_test.f'
