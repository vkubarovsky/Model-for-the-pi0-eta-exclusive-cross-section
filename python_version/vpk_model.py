#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import math

####   Constants
Mp=0.938272
Me=0.0005109989461
Mpi0=0.1349766
Meta=0.547862
Pi=3.1415926536
alpha = 0.00729927
HC2 = 389379.36


par_pi0 = np.array([ 6.2054,  0.8020, -0.1066,  1.8364,  0.00, 94.6474,  3.4426, -1.9769, -0.1324,  0.0000, 17.0001,  2.1228])
par_eta = np.array([ 6.8774,  1.7856, -1.3181,  0.7827,  0.00, 17.7821,  1.1451,   0.0053, 1.5816,  0.0000, 21.9122,  4.0069])
par     = par_pi0
Mmes    = Mpi0

def xsinit(mesonmass):
    global Mmes, par
    if mesonmass > 300.:
       par=par_eta
       Mmes = Meta
    else:
       par=par_pi0
       Mmes = Mpi0
       
def dvmpx(del2,xb,Q2,phi_g,E,heli,mesonmass):
#
#  dsigma/dQ2 dX dt dphi for ep-->ep pi0/eta
# input:
# del2=t (negative GeV^2)           NEGATIVE !!!!
# xb,Q2 x and Q^2 (GeV^2)
# Phi_g angle in the photon frame (radians)
# E beam energy of the electron in GeV
# heli electron helicity -1 or +1
# MESONMASS is the mass of the pi0 or eta
#
    global Mmes, par
    if (mesonmass > 0.300):
       par=par_eta
       Mmes = Meta
    else:
       par=par_pi0
       Mmes = Mpi0
    if xcheck_kine(del2,xb,Q2,E) == False: return(0.0)
    eps=epsilon(xb,Q2,E)
    ds = fluxw(xb,Q2,E)/(2*Pi)*(
        xsigma_T(del2,xb,Q2,E) + eps      * xsigma_L  (del2,xb,Q2,E)
        + eps                             * xsigma_TT (del2,xb,Q2,E) * math.cos(2*phi_g)
        + math.sqrt(2*eps*(1+eps))        * xsigma_LT (del2,xb,Q2,E) * math.cos(phi_g)
        + heli * math.sqrt(2*eps*(1-eps)) * xsigma_LTP(del2,xb,Q2,E) * math.sin(2*phi_g))
    if (ds < 0): ds=0.0
    return(ds)

def dvmpw(cost,W,Q2,phi_g,E,heli,mesonmass):
#
#  dsigma/dQ2 dW dCosT dphi for ep-->ep pi0/eta
# input:
# COST = Cos(Theta*), Theta* is the angle (pi0-gamma* or  p-p')in the CM system
# W is W
# Q2 is Q2
# Phi_g angle in the photon frame (radians)
# E energy of the electron in GeV
# heli electron helicity -1 or +1

    global Mmes, par
    if (mesonmass > 0.300):
       par=par_eta
       Mmes = Meta
    else:
       par=par_pi0
       Mmes = Mpi0
    if (check_kine(cost,W,Q2,E) == False): return(0.0)
    eps=epsilon(xB(Q2,W),Q2,E)
    ds = fluxw(xB(Q2,W),Q2,E)*JACG(Q2,W)/(2*Pi)*(
        sigma_T(cost,W,Q2,E) + eps        * sigma_L  (cost,W,Q2,E)
        + eps                             * sigma_TT (cost,W,Q2,E) * math.cos(2*phi_g)
        + math.sqrt(2*eps*(1+eps))        * sigma_LT (cost,W,Q2,E) * math.cos(  phi_g)
        + heli * math.sqrt(2*eps*(1-eps)) * sigma_LTP(cost,W,Q2,E) * math.sin(2*phi_g))
    if (ds<0): ds=0.0
    return(ds)

def printglobal():
    print('Mass= ',Mmes)
    print('Par = ',par )

def HT(del2,xb,Q2):
    return(par[0] * math.exp((par[1]+par[2]*(math.log(xb)-math.log(0.15)))*del2) * math.pow(Q2,par[3]/2))

def ET(del2,xb,Q2):
    return(par[5]*math.exp((par[6]+par[7]*(math.log(xb)-math.log(0.15)))*del2) * math.pow(Q2,par[8]/2))

def HTEBAR(del2,xb,Q2):
    return(par[10] * math.exp(par[11]*del2))

def xsigma_T(del2,xb,Q2,E):
    if xcheck_kine(del2,xb,Q2,E) == False: return(0.0)
    t0=tminq(Q2,xb)
    ksi=xb/(2-xb) * (1+Mp*Mp/Q2)
    return(4*Pi*alpha/2/phase(Q2,xb)/(Q2*Q2)*(
        (1-ksi*ksi)          * HT(del2,xb,Q2)*HT(del2,xb,Q2) +
        (-del2-t0)/8/(Mp*Mp) * ET(del2,xb,Q2)*ET(del2,xb,Q2)) * HC2)

def xsigma_LT(del2,xb,Q2,E):
    if xcheck_kine(del2,xb,Q2,E) == False: return(0.0)
    t0=tminq(Q2,xb)
    ksi=xb/(2-xb) * (1+Mp*Mp/Q2)
    return(4*Pi*alpha/math.sqrt(2.)/phase(Q2,xb)/math.pow(Q2,3./2.)*
        ksi*math.sqrt(1-ksi*ksi) *
        math.sqrt(-del2-t0)/2/Mp *   HTEBAR(del2,xb,Q2)*HTEBAR(del2,xb,Q2) * HC2)

def xsigma_TT(del2,xb,Q2,E):
    if xcheck_kine(del2,xb,Q2,E) == False: return(0.0)
    t0=tminq(Q2,xb)
    return(-4*Pi*alpha/2/phase(Q2,xb)/(Q2*Q2)*
        (-del2-t0)/8/(Mp*Mp) * ET(del2,xb,Q2)*ET(del2,xb,Q2) * HC2)

def xsigma_L(del2,xb,Q2,E): return(0.0)
def xsigma_LTP(del2,xb,Q2,E): return(0.0)
def xsigma_U(del2,xb,Q2,E): return (xsigma_T(del2,xb,Q2,E)+epsilon(xb,Q2,E)*xsigma_L(del2,xb,Q2,E))

def sigma_L  (cost,W,Q2,E): return(0.0)
def sigma_T  (cost,W,Q2,E): return(JACR(Q2,W)*xsigma_T  (Del2(cost,W,Q2),xB(Q2,W),Q2,E))
def sigma_TT (cost,W,Q2,E): return(JACR(Q2,W)*xsigma_TT (Del2(cost,W,Q2),xB(Q2,W),Q2,E))
def sigma_LT (cost,W,Q2,E): return(JACR(Q2,W)*xsigma_LT (Del2(cost,W,Q2),xB(Q2,W),Q2,E))
def sigma_LTP(cost,W,Q2,E): return(JACR(Q2,W)*xsigma_LTP(Del2(cost,W,Q2),xB(Q2,W),Q2,E))

#def sigma_LT (cost,W,Q2,E): return(2*JACR(Q2,W)*xsigma_LT (Del2(cost,W,Q2),xB(Q2,W),Q2,E))
#def sigma_LTP(cost,W,Q2,E): return(2*JACR(Q2,W)*xsigma_LTP(Del2(cost,W,Q2),xB(Q2,W),Q2,E))

def xcheck_kine(del2,xb,Q2,E):
      xmin1 = Q2/(2*Mp*E)
      xmax1 = 1.0
      nu  = Q2/(2*Mp*xb)
      W2  = Mp**2 + 2*Mp*nu - Q2
      if (W2 <= (Mp+Mmes)*(Mp+Mmes)):                       return(False)
      W   = math.sqrt(W2)
      qmod = math.sqrt(nu*nu + Q2)
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp*Mp-Mmes*Mmes)/(2*W)
#      print('E2cm ',E2cm) 
      if (E2cm <  Mp):                                      return(False)
      P2cm = math.sqrt(E2cm*E2cm - Mp*Mp)
      del2max = 2*(Mp*Mp - E1cm*E2cm - P1cm*P2cm)
      del2min = 2*(Mp*Mp - E1cm*E2cm + P1cm*P2cm)
#      print('xb ',xb,xmin1,xmax1) 
      if (xb <= xmin1 or xb>=xmax1):                        return(False)
#      print('del2 ',del2,del2min,del2max); 
      if (del2>=del2min or del2<=del2max):                  return(False)
      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)*(y*xb*Mp) / Q2
      eps=(1.0-y-e1)/(1-y+y*y/2+e1)
      if (eps<= 0. or eps >= 1.):                           return(False)
      ksi=xb/(2-xb)*(1+Mp*Mp/Q2)
      if (ksi >= 1.):                                       return(False)
      return(True)

def   xB(Q2,W):
      return(Q2/(W*W+Q2-Mp*Mp))

def   JACG(Q2,W):
# Jacobian d(cost,W)/d(t,xB)    
      xb=Q2/(W*W+Q2-Mp*Mp)
      return(2*W*xb/(W*W+Q2-Mp*Mp))

def   JACR(Q2,W):
# Jacobian d(Q2,W)/d(Q2.W)
      W2=W*W
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      qmod = math.sqrt(nu*nu + Q2)
      E1cm = Mp*(Mp + nu)/W
      E2cm = (W2 + Mp*Mp-Mmes*Mmes)/(2*W)
      P1cm = Mp*qmod/W
      if (E2cm <= Mp):                       return(0.0)
      P2cm = math.sqrt(E2cm*E2cm - Mp*Mp)
      return(2*P1cm*P2cm)

def   check_kine(cost,W,Q2,E):
      W2 = W*W
      xb = Q2/(W2+Q2-Mp*Mp)
      xmin1 = Q2/(2*Mp*E)
      xmax1 = 1.0
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      if (W <= (Mp+Mmes)):                       return(False)
      qmod = math.sqrt(nu*nu + Q2)
      if(W2              >= -Q2+2*Mp*E+Mp*Mp):    return(False)
      if(Q2/(4*E*(E-nu)) >= 1.0):                 return(False)
      if(xb              <= Q2/(2*Mp*E)):         return(False)
      if(abs(cost)       >= 1.0):                 return(False)
      xmin1 = Q2/(2*Mp*E)
      xmax1 = 1.0
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp*Mp-Mmes*Mmes)/(2*W)
      if (E2cm <  Mp):                                      return(False)
      P2cm = math.sqrt(E2cm*E2cm - Mp*Mp)
      del2max = 2*(Mp*Mp - E1cm*E2cm - P1cm*P2cm)
      del2min = 2*(Mp*Mp - E1cm*E2cm + P1cm*P2cm)
      del2=del2min-2*P1cm*P2cm*(1-cost)
      if (xb <= xmin1 or xb>=xmax1):                        return(False)
      if (del2>=del2min or del2<=del2max):                  return(False)
      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)*(y*xb*Mp) / Q2
      eps=(1.0-y-e1)/(1-y+y*y/2+e1)
      if (eps <= 0. or eps >= 1.):                          return(False)
      ksi=xb/(2-xb)*(1+Mp*Mp/Q2)
      if (ksi >= 1.):                                        return(False)
      return(True)
  
def   Del2(cost,W,Q2):
      W2 = W*W
      xb = Q2/(W2+Q2-Mp*Mp)
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      if (W <= (Mp+Mmes)):                       return(0.0)
      qmod = math.sqrt(nu*nu + Q2)

      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp*Mp-Mmes*Mmes)/(2*W)
      if (E2cm <= Mp):                            return(0.0)
      P2cm = math.sqrt(E2cm*E2cm - Mp*Mp)
      del2min = 2*(Mp*Mp - E1cm*E2cm + P1cm*P2cm)
      return       (del2min-2*P1cm*P2cm*(1-cost))


def fluxw(xb,Q2,E):
      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)*(y*xb*Mp)/Q2
      eps=(1.0-y-e1)/(1-y+y*y/2+e1)
      return(alpha/(2*Pi)*y*y/(1-eps)*(1-xb)/xb/Q2)
  
def epsilon(xb,Q2,E):
      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)*(y*xb*Mp)/Q2
      return( (1.0-y-e1)/(1-y+y*y/2+e1))

def tminq(Q2,xb):
      if(xb<=0. or xb>=1.):      return(0.0)
      W2 = Q2*(1./xb-1.)+Mp**2
      W=math.sqrt(W2)
      if(W<=Mp+Mmes):           return(0.0)

      E1CM=(W2+Q2+Mp*Mp)/(2*W)
      P1CM=math.sqrt(E1CM*E1CM-Mp*Mp)

      E3CM=(W2-Mmes*Mmes+Mp*Mp)/(2*W)
      P3CM=math.sqrt(E3CM*E3CM-Mp*Mp)
      
      return( -((Q2+Mmes*Mmes)*(Q2+Mmes*Mmes)/4./W2-(P1CM-P3CM)*(P1CM-P3CM)) )
    
def phase(Q2,xb):
      if(xb>=1.):             return(0.0)
      W2=Q2*(1/xb-1.)+Mp*Mp
      W4=W2*W2
      MP2=Mp*Mp
      MP4=MP2*MP2
      Q4=Q2*Q2
      LAMBDA=W4+Q4+MP4+2*W2*Q2-2*W2*MP2+2*Q2*MP2
      if(LAMBDA<=0.):          return(0.0)
      return(16.*Pi*(W2-MP2)*math.sqrt(LAMBDA))
  
#------------------------------------------------------------#

def plotsf(t,fn,xb,Q2,E,label,color):
    xs = t*0
    for i in range(np.size(t)):
        xs[i] = fn(-t[i],xb,Q2,E)
    plt.plot(t,xs,color,label=label)

def plotw(c,fn,w,Q2,E,label,color):
    xs = t*0
    for i in range(np.size(c)):
        xs[i] = fn(c[i],w,Q2,E)
    plt.plot(t,xs,color,label=label)
    
fig=plt.figure(figsize=(12,10),dpi=100)
plt.title("Structure Function vs -t", fontsize='16')	
#plt.plot([1, 2, 3, 4], [6,2,8,4])  #plot the points
plt.xlabel("-t",fontsize='13')	    
plt.ylabel("xs,nb",fontsize='13')
plt.legend(('sigma_T'),loc='upper center')	  
#plt.savefig('Y_X.png')	                  #saves the figure in the present directory
plt.grid()                                          

Ebeam=10.6
Ebeam=24.0
q2=1.14; xb=0.131
t = np.linspace(0., 2., 1000)
plt.plot([0.,2.],[0.,0.],'black')

plotsf(t,xsigma_T, xb,q2,Ebeam,"$\sigma_T$",   'red')
plotsf(t,xsigma_TT,xb,q2,Ebeam,"$\sigma_{TT}$",'blue')
plotsf(t,xsigma_LT,xb,q2,Ebeam,"$\sigma_{LT}$",'green')

plt.legend(loc='upper right', prop={'size': 24})
plt.show(block=False)
plt.savefig('ds_dxB.png')

#-----------------
fig=plt.figure(figsize=(12,10),dpi=100)
plt.title("Structure Function vs -t", fontsize='16')
#plt.plot([1, 2, 3, 4], [6,2,8,4])	                #plot the points
plt.xlabel("-t",fontsize='13')	 
plt.ylabel("xs,nb",fontsize='13')
plt.legend(('sigma_T'),loc='upper center')
plt.grid()                                          
cost=0.9; w=2.2
cst = np.linspace(-1., 1., 1000)
plt.xlabel("CosT",fontsize='13')	            
plt.plot([0.,2.],[0.,0.],'black')

plotw(cst,sigma_T, w,q2,Ebeam,"$\sigma_T$",   'red')
plotw(cst,sigma_TT,w,q2,Ebeam,"$\sigma_{TT}$",'blue')
plotw(cst,sigma_LT,w,q2,Ebeam,"$\sigma_{LT}$",'green')

plt.legend(loc='upper center', prop={'size': 24})
plt.show(block=False)
plt.savefig('ds_dW.png')
