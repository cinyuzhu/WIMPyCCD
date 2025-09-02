from scipy.interpolate import interp1d
import numpy as np

def lindhard(E_KeV):
    """ return the nuclear yield given by lindhard model (analycal)"""
    # lindhard model at k=0.146
    # input: nuclear recoil energy in KeV; atomic number of nuclear Z 
    # see also at Levin 5.1
    # k = 0.146
    Z = 14
    k = 0.15
    eta = 11.5 * E_KeV * Z ** (-7/3)
    g = 3 * eta ** 0.15 + 0.7 * eta ** 0.6 + eta
    nr_yield = k * g /(1 + k * g) # yield
    return nr_yield

def damicm_fit(Enr_KeV):
    """ return the nuclear yield given by damin-m fitting result (analycal)"""
    E_ee_KeV = 0.14* Enr_KeV ** 1.3 - 0.03
    nr_yield = E_ee_KeV/Enr_KeV
    return nr_yield



def get_full_yield(E_KeVnr, filename = ''):
    # input = a array of recoil energy
    # read the E_vr (unit = eVnr) vs. yield files
    # interpolate in measured region
    # E > Emax; yield = lindhard
    # E < Emin: yield = 0
    # return the yield on all energy region
    if filename == '':
        return  lindhard(E_KeVnr, Z = 14)
    
    data = np.loadtxt(filename, delimiter=",")
    E_data = data[:, 0] / 1000 # convert from eV to KeV
    Y_data = data[:, 1]
    interp = interp1d(E_data, Y_data, kind='linear', bounds_error=False, fill_value="extrapolate")
    E_min, E_max = E_data[0], E_data[-1]

    E = np.atleast_1d(E_KeVnr)  # ensure array-like input
    yield_vals = np.zeros_like(E)

    in_range = (E >= E_min) & (E <= E_max)
    above_range = E > E_max

    # Interpolate where valid
    yield_vals[in_range] = interp(E[in_range])

    # Use Lindhard model above range
    yield_vals[above_range] = lindhard(E[above_range])

    # Below E_min, yield remains 0
    return yield_vals if len(yield_vals) > 1 else yield_vals[0]