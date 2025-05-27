import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from scipy.special import erf
import analysis.utils as utils # here's rate calculation functions
import analysis.plotter as plotter
import analysis.config as config
from scipy.special import spherical_jn
import math

N0 = 6.02e26
sigma_pb = 1e-1 #pb
sigma_n = sigma_pb * 1e-40 #m2
rho_D = 7.13e-22 # SI
A = 28
M_T = A*0.932 * 1.78e-27 #SI = GeV * 1.78e-27
M_D = 1 * 1.78e-27
mu_n = M_T*M_D/(M_T+M_D)

def nuclear_form_factor(E_R_KeV, M_T = M_T, A = A):
    E_R_SI = E_R_KeV * 1.6e-16
    q = np.sqrt(2*M_T*E_R_SI)
    r_n = 1.14 * A**(1/3) #fm
    s = 0.9 #fm
    F = 3 * spherical_jn(1, q*r_n)/(q*r_n)*np.exp(-(q*s)**2/2)
    return F

def velo_int(E_R, v_0, v_E, v_esc):
    # get normalization factor k 
    E_R_SI = E_R *1.602e-16
    prefactor = np.pi**(1.5) * v_0**3
    term1 = erf(v_esc / v_0)
    term2 = (2 * v_esc / (np.sqrt(np.pi) * v_0)) * np.exp(-v_esc**2 / v_0**2)
    inverse_k = prefactor * (term1 - term2)
    k =  1 / inverse_k
    
    v_min = np.sqrt(M_T*E_R_SI/2/mu_n**2)
    
    prefactor2 = (np.pi**1.5 * v_0**3 * k) / (2 * v_E)
    v_min = np.asarray(v_min)  # ensure array-like for masking
    eta_val = np.zeros_like(v_min)

    # Case 1: v_min <= v_esc - v_E
    mask1 = v_min <= (v_esc - v_E)
    eta_val[mask1] = (
        erf((v_E - v_min[mask1]) / v_0) +
        erf((v_E + v_min[mask1]) / v_0) -
        (4 * v_E / (np.sqrt(np.pi) * v_0)) *
        np.exp(-v_esc**2 / v_0**2)
    )

    # Case 2: v_esc - v_E < v_min <= v_esc + v_E
    mask2 = ((v_min > (v_esc - v_E)) & (v_min <= (v_esc + v_E)))
    eta_val[mask2] = (
        erf((v_E - v_min[mask2]) / v_0) +
        erf(v_esc / v_0) -
        (2 * (v_E + v_esc - v_min[mask2]) / (np.sqrt(np.pi) * v_0)) *
        np.exp(-v_esc**2 / v_0**2)
    )

    # Case 3: v_min > v_esc + v_E -> already zero

    return prefactor2 * eta_val

def get_rate_ben(E_KeV, M_D = M_D, xs = sigma_n):
    rate_SI = (N0 * xs *rho_D *M_T)/(2*A*mu_n**2*M_D)*A**2*nuclear_form_factor(E_KeV)**2*velo_int(E_KeV, v_0=250e3, v_E=263e3, v_esc = 544e3)
    # rate_SI = (N0 * xs *rho_D *M_T)/(2*A*mu_n**2*M_D)*A**2*nuclear_form_factor(E_KeV)**2*velo_int(E_KeV, v_0=250e3, v_E=263e3, v_esc = 999999e3)
    # rate_SI = (N0 * xs *rho_D *M_T)/(2*A*mu_n**2*M_D)*A**2*nuclear_form_factor(E_KeV)**2*velo_int(E_KeV, v_0=250e3, v_E=0, v_esc = 999999e3)
    # print(N0, xs, rho_D , M_T, A, mu_n, M_D)
    # print(velo_int(E_KeV,v_0=250e3, v_E=263e3, v_esc = 544e3))
    # print(A**2*nuclear_form_factor(E_KeV)**2*velo_int(E_KeV, v_0=250e3, v_E=263e3, v_esc = 544e3))
    seconds_per_day = 24*3600
    keV_to_J = 1.6e-16
    return rate_SI * seconds_per_day * keV_to_J

# print(get_rate_ben(0.1))

def get_rate_13(E_KeV, v_esc = 544, m_D = 1, xs = sigma_pb):
    """
    Calculate the differential rate from numerical method, see derivation in paper"""
    E = E_KeV/1e6 # change to GeV
    R_0 = 540./config.A/m_D*xs*1*(config.v_0/250) 
    E_0 = 1/2 * m_D * (config.v_0/(3*10**5))**2 # in GeV
    m_T = config.A*0.938 
    r = 4*(m_D*m_T)/(m_D+m_T)**2
    
    v_min = np.sqrt(E/E_0/r)*config.v_0
    if v_min > v_esc:
        return 0
    diff1 = R_0/E_0/r*np.pi**(1/2)/4*config.v_0/config.v_E*(math.erf((v_min+config.v_E)/config.v_0)-math.erf((v_min-config.v_E)/config.v_0) )
    x = v_esc/config.v_0
    coef = math.erf(x) - 2/np.pi**0.5 * x * np.exp(-x**2)
    diff2 = 1/coef * (diff1 - R_0/E_0/r*np.exp(-x**2))
    F = np.exp(-1/3*(6.92*10**-3)**2*config.A*((1.14*config.A)**1/3)**2)
    if diff2 < 0:
        print("aaaa negative rate!!")
        return diff2 * F * config.A**2 * 1e-6
        # return 0
    else:
        return diff2 * F * config.A**2 * 1e-6
    
def get_rate_14(E_KeV, A = 28, v_0 = 250, m_D = 1, xs = sigma_pb):
    c1 = 0.751; c2 = 0.561

    E = E_KeV/1e6   # in GeV
    m_T = A*0.938  #[GeV]

    E0 = 1/2 * m_D * (v_0/(3*10**5))**2 # in GeV
    r = 4*m_D*m_T/(m_D+m_T)**2
    R0 = 540/(A*m_D)*(xs)*1*(v_0/250) #[pb]
    k = c2/(E0*r) + (1/3*(6.92*10**-3)**2*A*((1.14*A)**1/3)**2) # k factor in integral, see details in note
    return c1 * R0/E0/r * np.exp(-k*E) *A*A * 1e-6
