import numpy as np
from scipy.special import erf
from scipy.special import spherical_jn
import math
from dm_params import DMParams


def nuclear_form_factor(E_KeV, dm_params):
    # todo check why it's always so close to !???, compare with F from Levin?
    E_R_SI = E_KeV * 1.6e-16
    q = np.sqrt(2*dm_params.m_T_kg*E_R_SI)
    r_n = 1.14 * dm_params.mass_number_A**(1/3) * 1e-15 #m
    s = 0.9 * 1e-15 #m
    F = 3 * spherical_jn(1, q*r_n)/(q*r_n)*np.exp(-(q*s)**2/2)
    return F

def velo_int(E_R, dm_params):
    # get normalization factor k 
    E_R_SI = E_R *1.602e-16
    prefactor = np.pi**(1.5) * dm_params.v0_mps**3
    term1 = erf(dm_params.vesc_mps / dm_params.v0_mps)
    term2 = (2 * dm_params.vesc_mps / (np.sqrt(np.pi) * dm_params.v0_mps)) * np.exp(-dm_params.vesc_mps**2 / dm_params.v0_mps**2)
    inverse_k = prefactor * (term1 - term2)
    k =  1 / inverse_k
    
    v_min = np.sqrt(dm_params.m_T_kg*E_R_SI/2/dm_params.mu_n**2)
    
    prefactor2 = (np.pi**1.5 * dm_params.v0_mps**3 * k) / (2 * dm_params.vE_mps)
    v_min = np.asarray(v_min)  # ensure array-like for masking
    eta_val = np.zeros_like(v_min)

    # Case 1: v_min <= v_esc - v_E
    mask1 = v_min <= (dm_params.vesc_mps - dm_params.vE_mps)
    eta_val[mask1] = (
        erf((dm_params.vE_mps - v_min[mask1]) / dm_params.v0_mps) +
        erf((dm_params.vE_mps + v_min[mask1]) / dm_params.v0_mps) -
        (4 * dm_params.vE_mps / (np.sqrt(np.pi) * dm_params.v0_mps)) *
        np.exp(-dm_params.vesc_mps**2 / dm_params.v0_mps**2)
    )

    # Case 2: v_esc - v_E < v_min <= v_esc + v_E
    mask2 = ((v_min > (dm_params.vesc_mps - dm_params.vE_mps)) & (v_min <= (dm_params.vesc_mps + dm_params.vE_mps)))
    eta_val[mask2] = (
        erf((dm_params.vE_mps - v_min[mask2]) / dm_params.v0_mps) +
        erf(dm_params.vesc_mps / dm_params.v0_mps) -
        (2 * (dm_params.vE_mps + dm_params.vesc_mps - v_min[mask2]) / (np.sqrt(np.pi) * dm_params.v0_mps)) *
        np.exp(-dm_params.vesc_mps**2 / dm_params.v0_mps**2)
    )

    # Case 3: v_min > v_esc + v_E -> already zero

    return prefactor2 * eta_val
   

def get_rate_ben(E_KeV, dm_params):
    # 'supposedly correct function from ben loer's thesis
    rate_SI = ((dm_params.N0 * dm_params.sigma_m2 *dm_params.rho_SI *dm_params.m_T_kg)/(2*dm_params.mass_number_A*dm_params.mu_n**2*dm_params.m_D_kg)*
                dm_params.mass_number_A**2*nuclear_form_factor(E_KeV, dm_params)**2*velo_int(E_KeV, dm_params))
    # print(dm_params.N0, dm_params.sigma_m2, dm_params.rho_SI, dm_params.m_T_kg, dm_params.mass_number_A, dm_params.mu_n, dm_params.m_D_kg)
    # print((dm_params.N0 * dm_params.sigma_m2 *dm_params.rho_SI *dm_params.m_T_kg)/(2*dm_params.mass_number_A*dm_params.mu_n**2*dm_params.m_D_kg))
    # print(dm_params.mass_number_A**2*nuclear_form_factor(E_KeV, dm_params)**2*velo_int(E_KeV, dm_params))
    # print(rate_SI)
    # print(velo_int(E_KeV, dm_params))
    seconds_per_day = 24*3600
    keV_to_J = 1.6e-16
    return rate_SI * seconds_per_day * keV_to_J #in KeV-1kg-1day-1 unit


# params = DMParams(m_D_GeV = 1,
#                  sigma_pb = 0.1,
#                  rho_GeV_cm3 = 0.4,
#                  v0_kms = 250,
#                  vE_kms = 263,
#                  vesc_kms = 544,
#                  exposure_kgdays = 11,
#                  confidence_level = 0.9,
#                  mass_number_A = 28,
#                  c1 = 0.751, c2 = 0.561,
#                  E_min_KeV= 0.1
#                  )

# print(get_rate_ben(0.1, params))

def get_rate_14(E_KeV, dm_params):
    # from levin paper eq 3.14, (vE = 230 (fixed), vesc = inf)
    A = dm_params.mass_number_A
    m_D = dm_params.m_D_GeV
    v_0 = dm_params.v0_kms
    xs = dm_params.sigma_pb
    c1 = dm_params.c1
    c2 = dm_params.c2
    E_GeV = E_KeV/1e6   # in GeV
    m_T = A*0.938  #[GeV]

    E0 = 1/2 * m_D * (v_0/(3*10**5))**2 # in GeV
    r = 4*m_D*m_T/(m_D+m_T)**2
    R0 = 540/(A*m_D)*(xs)*1*(v_0/230) #[pb]
    k = c2/(E0*r) + (1/3*(6.92*10**-3)**2*A*((1.14*A)**1/3)**2) # k factor in integral, see details in note
    return c1 * R0/E0/r * np.exp(-k*E_GeV) *A*A * 1e-6

def get_rate_12(E_KeV, dm_params):
    # from levin paper eq 3.12
    A = dm_params.mass_number_A
    m_D = dm_params.m_D_GeV
    v_0 = dm_params.v0_kms
    v_E = dm_params.vE_kms
    xs = dm_params.sigma_pb
    m_T = A*0.938  #[GeV]
    E = E_KeV/1e6 # in GeV
    R_0 = 540./A/m_D*xs*1*(v_0/230) 
    E_0 = 1/2 * m_D * (v_0/(3*10**5))**2 # in GeV
    r = 4*(m_D*m_T)/(m_D+m_T)**2
    
    v_min = np.sqrt(E/E_0/r)*v_0
    diff1 = R_0/E_0/r*np.pi**(1/2)/4*v_0/v_E*(math.erf((v_min+v_E)/v_0)-math.erf((v_min-v_E)/v_0) )
    # F = np.exp(-1/3*(6.92*10**-3)**2*A*((1.14*A)**1/3)**2)
    F = nuclear_form_factor(E_KeV, dm_params)**2
    return diff1 * F * A**2 * 1e-6

def get_rate_13(E_KeV, dm_params):
    
    # known incorrect leven eq 3.13
    A = dm_params.mass_number_A
    m_D = dm_params.m_D_GeV
    v_0 = dm_params.v0_kms
    v_E = dm_params.vE_kms
    v_esc = dm_params.vesc_kms
    xs = dm_params.sigma_pb
    m_T = A*0.938  #[GeV]
    R_0 = 540./A/m_D*xs*1*(v_0/230) 
    E_0 = 1/2 * m_D * (v_0/(3*10**5))**2 # in GeV
    r = 4*(m_D*m_T)/(m_D+m_T)**2

    E = E_KeV/1e6 # change to GeV
    R_0 = 540./A/m_D*xs*1*(v_0/230) 
    E_0 = 1/2 * m_D * (v_0/(3*10**5))**2 # in GeV
    m_T = A*0.938 
    r = 4*(m_D*m_T)/(m_D+m_T)**2
    
    v_min = np.sqrt(E/E_0/r)*v_0
    if v_min > v_esc:
        return 0
    diff1 = R_0/E_0/r*np.pi**(1/2)/4*v_0/v_E*(math.erf((v_min+v_E)/v_0)-math.erf((v_min-v_E)/v_0) )
    x = v_esc/v_0
    coef = math.erf(x) - 2/np.pi**0.5 * x * np.exp(-x**2)
    diff2 = 1/coef * (diff1 - R_0/E_0/r*np.exp(-x**2))
    # F = np.exp(-1/3*(6.92*10**-3)**2*A*((1.14*A)**1/3)**2)
    F = nuclear_form_factor(E_KeV, dm_params)**2
    if diff2 < 0:
        print("aaaa negative rate!!")
    return diff2 * F * A**2 * 1e-6
    