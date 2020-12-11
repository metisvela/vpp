"""
File dove raccogliere tutte le funzioni utilizzate più di una volta nel codice per
facilitarne la modifica
"""
from xfoil import XFoil
import numpy as np

def lift_coefficients(angleOfAttack, aspectRatio):
    angleOfAttack = np.radians(angleOfAttack)
    cL2D = 0.1 # slope of the 2D lift profile
    cL = cL2D / (1 + 2 / aspectRatio) * angleOfAttack

    cDvisc = 0 # For now we don't take it into account
    cDind = cL **2 / (np.pi * aspectRatio)
    cDtot = cDvisc + cDind
    return cL, cDtot

def dynamic_pressure(Sea, boatSpeed):
    return 0.5 * Sea.waterDensity * boatSpeed**2
