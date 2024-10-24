import numpy as np
from .constants import *

def xsinit(mesonmass):
    global Mmes, par
    if mesonmass > 300.:
        par = par_eta
        Mmes = Meta
    else:
        par = par_pi0
        Mmes = Mpi0

def dvmpx_vectorized(del2, xb, Q2, phi_g, E, heli, mesonmass):
    global Mmes, par
    par, Mmes = (par_eta, Meta) if mesonmass > 0.300 else (par_pi0, Mpi0)
    
    valid_kine = np.vectorize(xcheck_kine)(del2, xb, Q2, E)
    
    eps = epsilon(xb, Q2, E)
    flux = fluxw(xb, Q2, E) / (2 * np.pi)
    
    sigma_T = xsigma_T(del2, xb, Q2, E)
    sigma_L = xsigma_L(del2, xb, Q2, E)
    sigma_TT = xsigma_TT(del2, xb, Q2, E)
    sigma_LT = xsigma_LT(del2, xb, Q2, E)
    sigma_LTP = xsigma_LTP(del2, xb, Q2, E)
    
    ds = flux * (
        sigma_T + eps * sigma_L
        + eps * sigma_TT * np.cos(2 * phi_g)
        + np.sqrt(2 * eps * (1 + eps)) * sigma_LT * np.cos(phi_g)
        + heli * np.sqrt(2 * eps * (1 - eps)) * sigma_LTP * np.sin(2 * phi_g)
    )
    
    ds = np.maximum(ds, 0)  # Ensuring non-negative cross-section values
    
    return ds * valid_kine  # Return only valid kinematic values

def dvmpw_vectorized(cost, W, Q2, phi_g, E, heli, mesonmass):
    global Mmes, par
    par, Mmes = (par_eta, Meta) if mesonmass > 0.300 else (par_pi0, Mpi0)
    
    del2 = -Mp**2 + W**2 - Q2 + 2 * Mp * (E - np.sqrt(E**2 - Q2))
    valid_kine = np.vectorize(xcheck_kine)(del2, xb, Q2, E)
    
    eps = epsilon(xb, Q2, E)
    flux = fluxw(xb, Q2, E) / (2 * np.pi)
    
    sigma_T = sigma_T(del2, xb, Q2, E)
    sigma_L = sigma_L(del2, xb, Q2, E)
    sigma_TT = sigma_TT(del2, xb, Q2, E)
    sigma_LT = sigma_LT(del2, xb, Q2, E)
    sigma_LTP = sigma_LTP(del2, xb, Q2, E)
    
    ds = flux * (
        sigma_T + eps * sigma_L
        + eps * sigma_TT * np.cos(2 * phi_g)
        + np.sqrt(2 * eps * (1 + eps)) * sigma_LT * np.cos(phi_g)
        + heli * np.sqrt(2 * eps * (1 - eps)) * sigma_LTP * np.sin(2 * phi_g)
    )
    
    ds = np.maximum(ds, 0)  # Ensuring non-negative cross-section values
    
    return ds * valid_kine  # Return only valid kinematic values

def xcheck_kine(del2, xb, Q2, E):
    W2 = Mp**2 + Q2 * (1 - xb) / xb - del2
    return (W2 >= (Mp + Mpi0)**2) & (W2 <= (Mp + Meta)**2) & (del2 <= 0)

def epsilon(xb, Q2, E):
    nu = Q2 / (2 * Mp * xb)
    return (1 - nu**2 / E**2) / (1 + nu**2 / (2 * E**2))

def fluxw(xb, Q2, E):
    nu = Q2 / (2 * Mp * xb)
    W2 = Mp**2 + Q2 * (1 - xb) / xb
    return (alpha / (2 * np.pi**2)) * (nu / Q2) * (1 - xb) / (1 - epsilon(xb, Q2, E))

def xsigma_T(del2, xb, Q2, E):
    return np.ones_like(del2)

def xsigma_L(del2, xb, Q2, E):
    return np.ones_like(del2)

def xsigma_TT(del2, xb, Q2, E):
    return np.ones_like(del2)

def xsigma_LT(del2, xb, Q2, E):
    return np.ones_like(del2)

def xsigma_LTP(del2, xb, Q2, E):
    return np.ones_like(del2)

def sigma_T(del2, xb, Q2, E):
    return np.ones_like(del2)

def sigma_L(del2, xb, Q2, E):
    return np.ones_like(del2)

def sigma_TT(del2, xb, Q2, E):
    return np.ones_like(del2)

def sigma_LT(del2, xb, Q2, E):
    return np.ones_like(del2)

def sigma_LTP(del2, xb, Q2, E):
    return np.ones_like(del2)
