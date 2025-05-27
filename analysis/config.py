# put all global constants here
GeV_to_KeV = 1e6 
CL = 0.95
A = 28. # target mass Si
c1 = 0.751
c2 = 0.561
m_T = A*0.938 #[GeV]
v_0 = 230. #[km s-1] d
# v_esc = 999999 #[km s-1]
v_esc = 600. #[km s-1] # as recommended in Lewin and Smith
# v_E = 244. #[km s-1]
v_E = 263.
flag_nuclear_form_factor = True

#base parameter:
m_D_0 = 1 # GeV
xs_0 = 10**-1 #pb

# recoil energy bounds for detector, E1 = min, E2 = max
E1 = 1*10**-6 #in GeV
E2 = 100 #in GeV 

# the length and mass of the experiment data-taking
DAQday = 1
DAQkg = 1


# add functions for dynamically updating constants if needed
def update_constants(**kwargs):
    """dynamically updating constants if needed"""
    globals().update(kwargs)

# usecase:
# config.update_constants(v_0=200.0)  # Modify constants at runtime

def reset_constants():
    """Reset constants to their original values by reloading the module."""
    import importlib
    import sys
    importlib.reload(sys.modules[__name__])

# usecase: 
# config.reset_constants()