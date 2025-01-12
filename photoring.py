"""Compute photoring effect
"""
from geotrans2 import *
import numpy as np

S = RingedSystem(
    system=dict(
        Mstar=1.0*MSUN,
        Rstar=1.0*RSUN,
        Lstar=1.0*LSUN, 
        ap=1.0*AU, #Semi-major axis
        iorb=90.0*DEG, #Orbital inclination
        Rplanet=1.0*RSAT, #Planet radius
        fe=RSAT_ARING/RSAT, #Exterior ring (Rp)
        fi=RSAT_BRING/RSAT, #Interior ring (Rp)
        ir=45.0*DEG, #Ring inclination
        phir=45.0*DEG, #Ring roll angle
        tau=1.0, #Opacity        
    )
)
S.fe=1
S.tau=0
PR = S.PR()
print(np.cos(S.ieff),S.teff*RAD)
print(S.rho_true,S.rho_obs,PR)
