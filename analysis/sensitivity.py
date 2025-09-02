import numpy as np
from tqdm import tqdm
import scipy.integrate as integrate
import os

def get_Nexp_from_CL(CL = 0.90):
    return -np.log(1-CL)

def scan_sensitivity(dm_params, rate_function, lgmass_range = [0, 1], lgxs_range = [-4, 0], npoints = 10, filetag = '', overwrite = False):
    output_filename = f"data/sensitivity_{filetag}.npz" if filetag else "data/sensitivity_.npz"

    # if overwrite == False:
    #     if os.path.exists(output_filename):
    #         data = np.load(output_filename)
    #         return data["mass"], data["xs"]
    
    contour_mass = []
    contour_xs = []
    lgmasses  = np.linspace(*lgmass_range, npoints)
    lgxses = np.linspace(*lgxs_range,npoints)

    for i, m_D in enumerate(tqdm(10**lgmasses, total=npoints, desc="Scanning m_D", unit="mass")):
        for j, xs in enumerate(10**lgxses):
            dm_params.m_D_GeV = m_D
            dm_params.sigma_pb = xs
            
            R_tot = integrate.quad(lambda E: rate_function(E, dm_params = dm_params), dm_params.E_min_KeV, 40)[0]
            if abs(R_tot * dm_params.exposure_kgdays - get_Nexp_from_CL(dm_params.confidence_level)) < 0.1:
                contour_mass.append(m_D)
                contour_xs.append(xs)

    np.savez(output_filename, mass=np.array(contour_mass), xs=np.array(contour_xs))
    return np.array(contour_mass), np.array(contour_xs)




