import numpy as np
from tqdm import tqdm
from analysis.rate_models import get_rate_ben
from analysis.yield_functions import yield1

def get_Nexp_from_CL(CL = 0.90):
    return -np.log(1-CL)

def scan_sensitivity_hist(dm_params, threshhold_KeV = 0.01, exposure_kgdays = 365, yield_function=yield1, lgmass_range = [-1, 2], lgxs_range = [-8, 0], npoints = 10):
    contour_mass = []
    contour_xs = []
    lgmasses  = np.linspace(*lgmass_range, npoints)
    lgxses = np.linspace(*lgxs_range,npoints)

    for m_D in tqdm(10**lgmasses, total=npoints, desc="Scanning m_D", unit="mass"):
        for xs in 10**lgxses:
            dm_params.m_D_GeV = m_D
            dm_params.sigma_pb = xs
            
            ##### start change for hists #####
            E_nr_min, E_nr_max = 0.001, 30.0     # keVnr (example)
            dE_nr = 0.001                      # fine recoil step for smooth remapping
            E_nr = np.arange(E_nr_min, E_nr_max + dE_nr, dE_nr)

            # Compute weights (expected counts per small recoil bin)
            weights = get_rate_ben(E_nr,dm_params) * dE_nr

            # Map to ionized (measured) energy
            E_ion = yield_function(E_nr)   

            # hist of E_ee
            Eion_min, Eion_max = 0, 8
            Eion_width = 0.01
            Eion_edges = np.arange(Eion_min, Eion_max + Eion_width, Eion_width)
            Eion_centers = 0.5 * (Eion_edges[1:] + Eion_edges[:-1])
            hist_Eion, _ = np.histogram(E_ion, bins=Eion_edges, weights=weights)

            # only sum over the bins with energy > threshhold
            mask = Eion_centers > threshhold_KeV
            R_tot = hist_Eion[mask].sum()
            ##### end change for hists #####

            if abs(R_tot * exposure_kgdays - get_Nexp_from_CL(dm_params.confidence_level)) < 0.1:
                contour_mass.append(m_D)
                contour_xs.append(xs)
    return np.array(contour_mass), np.array(contour_xs)
