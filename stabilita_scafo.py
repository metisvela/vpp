# -*- coding: utf-8 -*-
"""
Script per calcolare la stabilità di forma data da uno scafo parallelepipedo,
per approssimare gli scafi Métis

"""

import numpy as np

def hull_stability(Boat, Sea, boatSpeed):
    h = Boat.freeboard + Boat.canoeDraft
    b = Boat.beamWaterl
    Tc = Boat.canoeDraft
    VCG = h/2 # stima da migliorare

    if np.radians(Boat.rollAngle) < np.arctan(2*Tc/b):
    # caso 1
        theta = np.radians(Boat.rollAngle)
        VCB = Tc/2 + b**2/(24*Tc)*np.tan(theta)**2
        HCB = b**2/(12*Tc) * np.tan(theta)

        GZ = HCB*np.cos(theta) - (VCG - VCB) * np.sin(theta)

    elif (np.radians(Boat.rollAngle) > np.arctan(2*Tc/b)) & (np.radians(Boat.rollAngle) < np.arctan((h-Tc)/b)) :
        #caso 2
        theta = np.radians(Boat.rollAngle)
        VCB = Tc/3 + b/6*np.tan(theta)
        HCB = b/3 - Tc/3 * np.tan(theta)**(-1)

        GZ = HCB*np.cos(theta) - (VCG - VCB) * np.sin(theta)

    else:
        print("Angoli troppo alti per calcolare la stabilità di scafo!")
        GZ = 0
    # caso 3
    # da compilare


    rightingMoment = Boat.maxDisplacement * 1000 * Sea.gravityConstant * GZ # displacement da m3 a kg
    
    try:
        rightingMoment = rightingMoment * np.ones(shape=boatSpeed.size) # adatto a vettore per il plotting
    except AttributeError:
        pass 
        
    return rightingMoment 
