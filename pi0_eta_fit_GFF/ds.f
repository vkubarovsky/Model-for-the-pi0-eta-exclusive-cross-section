      REAL FUNCTION DS(del2,x,Q2,Phi,MASS)
      implicit none
      REAL del2,x,Q2,Phi,MASS

      REAL       DVMP,fluxw, pi,E,phir
      PARAMETER (PI=3.1415926536)

      COMMON /PAWPAR/ P(10)
       REAL P
       REAL AM(2)
       DATA AM/0.134976,0.547300/
       INTEGER K
       COMMON/PIETA/K
c     INCLUDE 'pi0eta.par'
      E=5.75
      CALL XSINIT(MASS)
      phir=phi/180.*pi
      DS=DVMP(del2,x,Q2,Phir,E,0.,MASS)/FLUXW(x,Q2,E)
     *      *(2.*PI)

      RETURN
      END
      

      REAL FUNCTION DVMP(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
C
c  dsigma/dQ2 dX dt dphi for ep-->ep pi0
C
C exc pi0 x-section
c
cinput:
cdel2=t (negative GeV^2)           NEGATIVE !!!!
cxb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c MESONMASS is the mass of th epi0 or eta

      IMPLICIT NONE
      REAL del2,xb,Q2, Phi_g,E,heli,MESONMASS
      real Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     S_T, S_L, S_LT, S_TT,S_LTP
      REAL     XSIGMA_T, XSIGMA_L, XSIGMA_LT, XSIGMA_TT,XSIGMA_LTP
      REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,XSIGMA_T,XSIGMA_L, XSIGMA_LT, XSIGMA_TT
      EXTERNAL XSIGMA_LTP
      REAL     FL
      LOGICAL  XCHECK_KINE
      
      COMMON /PAWPAR/ P(10)
       REAL P
       REAL AM(2)
       DATA AM/0.134976,0.547300/
       INTEGER K
       COMMON/PIETA/K

c      INCLUDE 'pi0eta.par'
c
       K=1
       CALL XSINIT(MESONMASS)
      
      DVMP = 0.0
      
      IF(.NOT.XCHECK_KINE(del2,xb,Q2,E)) RETURN
       EPS  =  EPSILON(XB,q2,e)
       S_T  =  104.6
       S_L  =  0.
       S_LT =  7.6
       S_TT =  -80.1
       S_LTP = 0.

       DVMP =FLUXW(xb,Q2,E)/(2.*PI)*(
     *                       S_T   + EPS*S_L          + 
     * EPS                  *S_TT  * COS(2*PHI_G)       + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * ) 
      if(DVMP.lt.0.) DVMP=0.

      RETURN
      END
C====================================================================================c

      real function fluxw(x,Q2,E)
      implicit none
      real x,Q2,E,y,eps,e1
      real alpha,Mp,PI
      parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)

c
      y=Q2/(2*Mp*x*E)
      e1=(y*x*Mp)**2/Q2
      eps=(1.0-y-e1)/(1-y+y**2/2+e1)
      fluxw=alpha/(2*PI)*y*y/(1-eps)*(1-x)/x/Q2
      return
      end
C====================================================================================c

      real function EPSILON(x,Q2,E)
      implicit none
      real x,Q2,E,y,eps,e1
      real alpha,Mp,PI
      parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)
c
      y=Q2/(2*Mp*x*E)
      e1=(y*x*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
      return
      end
c-------------------------------------------------------------------------

      SUBROUTINE XSINIT(MESONMASS)
C
      IMPLICIT NONE
      REAL MESONMASS
      INCLUDE 'pi0eta.par'

      IF(MESONMASS.GT.0.140) THEN
           K=2
      ELSE
           K=1
      ENDIF

      RETURN
      END

C====================================================================================c

      LOGICAL FUNCTION XCHECK_KINE(del2,xb,Q2,E)
      IMPLICIT NONE
      REAL del2,xb,Q2,E
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
           

      real nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real  xmin1,xmax1
      real y, e1, epsilon

       REAL AM(2)
       DATA AM/0.134976,0.547300/
       INTEGER K
       COMMON/PIETA/K

      MPI0=AM(K)
      XCHECK_KINE=.FALSE.
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2.0*Mp*nu - Q2
      IF(W2.LT.(Mp+Mpi0)**2)                  RETURN
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)

      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2.*W)
      IF(E2cm.LE.Mp)                             RETURN
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)

      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
  
      IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

      XCHECK_KINE=.TRUE.

      RETURN
      END
