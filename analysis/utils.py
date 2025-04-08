import analysis.config as config
import numpy as np
import math

# for rate: consider interaction factor A^2 by default

def get_spectrum_in_GeV_analytical(E):
    """
    Calculate the differential rate from numerical method, only v_esc = inf., see derivation in writeup
    E: recoil energy in GeV
    output = dR/dE in GeV^-1 kg^-1 day^-1
    """
    m_T = config.A*0.938 #[GeV]
    m_D = config.m_D_0
    xs = config.xs_0
    flag_nuclear_form_factor = True

    # most probable incident energy of the dark matter [GeV]
    E0 = 1/2 * m_D * (config.v_0/(3*10**5))**2 # in GeV
    r = 4*m_D*m_T/(m_D+m_T)**2
    R0 = 540/(config.A*m_D)*(xs)*1*(config.v_0/230)
    if flag_nuclear_form_factor == True:
        k = config.c2/(E0*r) + (1/3*(6.92*10**-3)**2*config.A*((1.14*config.A)**1/3)**2) # k factor in integral, see details in note
    else:
        k = config.c2/(E0*r) 
 
    return config.c1 * R0/E0/r * np.exp(-k*E)*config.A**2 #in GeV^-1 kg^-1 day^-1

def get_spectrum_in_GeV_numerical(E_R):
    """
    Calculate the differential rate from numerical method, see derivation in paper"""
    m_D = config.m_D_0
    xs = config.xs_0
    R_0 = 540./config.A/m_D*xs*1*(config.v_0/230) 
    E_0 = 1/2 * m_D * (config.v_0/(3*10**5))**2 # in GeV
    m_T = config.A*0.938 
    r = 4*(m_D*m_T)/(m_D+m_T)**2
    
    v_min = np.sqrt(E_R/E_0/r)*config.v_0
    diff1 = R_0/E_0/r*np.pi**(1/2)/4*config.v_0/config.v_E*(math.erf((v_min+config.v_E)/config.v_0)-math.erf((v_min-config.v_E)/config.v_0) )
    k_0 = np.pi**(3/2)*(config.v_0**3)
    k_1  = k_0 * (math.erf(config.v_esc/config.v_0) - 2/np.pi**0.5 * config.v_esc/config.v_0 * np.exp(-config.v_esc**2/config.v_0**2))
    diff2 = k_0/k_1 * (diff1 - R_0/E_0/r*np.exp(-config.v_esc**2/config.v_0**2))

    F = np.exp(-1/3*(6.92*10**-3)**2*config.A*((1.14*config.A)**1/3)**2)
    if diff2 < 0:
        return 0
    else:
        return diff2 * F * config.A**2