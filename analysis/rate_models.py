import numpy as np
from scipy.special import erf
from scipy.special import spherical_jn


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
    # 'supposedly correct function from ben loer's thesis section 1.2.2
    rate_SI = ((dm_params.N0 * dm_params.sigma_m2 *dm_params.rho_SI *dm_params.m_T_kg)/(2*dm_params.mass_number_A*dm_params.mu_n**2*dm_params.m_D_kg)*
                dm_params.mass_number_A**2*nuclear_form_factor(E_KeV, dm_params)**2*velo_int(E_KeV, dm_params))
    seconds_per_day = 24*3600
    keV_to_J = 1.6e-16
    return rate_SI * seconds_per_day * keV_to_J #in KeV-1kg-1day-1 unit