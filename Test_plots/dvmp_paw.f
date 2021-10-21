      SUBROUTINE DVMP_PAW(MESONMASS)
      IMPLICIT NONE
      REAL MESONMASS
      INTEGER K
      COMMON/PIETA/K
      REAL PAR
      COMMON /PAWPAR/ PAR(12)
      REAL AM(2)
      DATA AM/0.134976,0.547300/

      CALL XSINIT(MESONMASS)

      PRINT 100,AM(K),PAR
 100  FORMAT('M=',F10.5/ '/PAWPAR/',12F8.3)

C      DVMPX_PAW=0.
      RETURN
      END
C--------1--------2---------3---------4---------5---------6---------7--
      REAL function DVMPP(COST,W,Q2, Phi_g,E,heli,MESONMASS)
      IMPLICIT NONE
      REAL    COST,W,del2,xb,Q2, Phi_g,E,mesonmass
      INTEGER heli
      REAL    Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     sigma0, S_T, S_L, S_LT, S_TT,S_LTP
      REAL     SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT,SIGMA_LTP
      REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,SIGMA_T,SIGMA_L, SIGMA_LT, SIGMA_TT
      EXTERNAL SIGMA_LTP
      REAL     JACG,JACR
      LOGICAL  CHECK_KINE,KINE
      REAL S_TOT,Ee,nu

      DVMPP=0.
      KINE=CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)
      PRINT *,'--------DVMPP------- '
      PRINT *,DVMPP,COST,W,Q2,Phi_g,E,heli,mesonmass
      PRINT *, 't,xB ',del2,xb
      PRINT *, 'KINE ',KINE
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR))
     *                   RETURN
      PRINT *,'---- KINE is TRUE'
      
      CALL  DVMPW(COST,W,Q2, Phi_g,E,heli,MESONMASS,
     *            S_T,S_L,S_TT,S_LT,S_LTP)
      PRINT *,'-----After call dvmpw'
      EPS=EPSILON(xb,Q2,E)
      DVMPP =FLUXW(xb,Q2,E)/(2.*PI)*(
     *                       S_T   + EPS*S_L          + 
     * EPS                  *S_TT  * COS(2*PHI_G)       + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * )

      PRINT *,'-----DVMPP ',DVMPP,COST,W,Q2, Phi_g,E
      PRINT *,  FLUXW(xb,Q2,E),EPS
      nu=Q2/(2*Mp*xB)
      Ee=E-nu
      PRINT *,'E electron ',Ee,del2,w,q2,xb,nu
      RETURN
      END


      INCLUDE '../Functions/dvmpw.F'
      INCLUDE '../Functions/dvmpx.F'
