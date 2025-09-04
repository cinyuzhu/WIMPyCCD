# input = Enr bin, output = Eee bin

from scipy.interpolate import interp1d
import numpy as np

def load_and_interpolate(fname, kind="cubic"):
    # load two columns, skipping lines starting with '#'
    data = np.loadtxt(fname, comments="#")
    x, y = data[:,0], data[:,1]
    f = interp1d(x, y, kind=kind, bounds_error=False, fill_value=np.nan, assume_sorted=False)
    return f, float(np.min(x)), float(np.max(x))


def yield1(Enr_KeV):
    return Enr_KeV


def lindhard(Enr_KeV):
    """ return the nuclear yield given by lindhard model (analycal)"""
    # input: nuclear recoil energy in KeV; atomic number of nuclear Z 
    # see also at Levin 5.1 # k = 0.146
    Z = 14
    k = 0.15
    eta = 11.5 * Enr_KeV * Z ** (-7/3)
    g = 3 * eta ** 0.15 + 0.7 * eta ** 0.6 + eta
    nr_yield = k * g /(1 + k * g) # yield
    return nr_yield * Enr_KeV


def damicm_fit(Enr_KeV):
    """ return the nuclear yield given by damin-m fitting result (analycal)"""
    E_ee_KeV = 0.14* Enr_KeV ** 1.3 - 0.03
    return E_ee_KeV

def damicm_no_fano(Enr_KeV):
    """Inside file range: interpolate; below range: 0; above range: Lindhard."""
    f, xmin, xmax = load_and_interpolate('from_Julian/004-merged-main-SYS_001.txt')
    E = np.asarray(Enr_KeV, dtype=float)
    y = f(E)

    below = E < xmin
    above = E > xmax

    # Inside range -> keep interpolated y; below -> 0; above -> Lindhard
    if np.isscalar(Enr_KeV):
        if below:   return 0.0
        if above:   return float(lindhard(E))
        return float(y)
    else:
        y[below] = 0.0
        y[above] = lindhard(E[above])
        return y
    
def damicm_fano(Enr_KeV):
    """Inside file range: interpolate; below range: 0; above range: Lindhard."""
    f, xmin, xmax = load_and_interpolate('from_Julian/005-merged-main-gamm-overall-1-central-SYS_000.txt')
    E = np.asarray(Enr_KeV, dtype=float)
    y = f(E)

    below = E < xmin
    above = E > xmax

    # Inside range -> keep interpolated y; below -> 0; above -> Lindhard
    if np.isscalar(Enr_KeV):
        if below:   return 0.0
        if above:   return float(lindhard(E))
        return float(y)
    else:
        y[below] = 0.0
        y[above] = lindhard(E[above])
        return y
    


