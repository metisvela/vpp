# -*- coding: utf-8 -*-
"""
File dove raccogliere tutte le funzioni utilizzate più di una volta nel codice per
facilitarne la modifica
"""
import xfoil
import numpy as np
from scipy.interpolate import RectBivariateSpline

def lift_coefficients_3D(cl, cd, aspectRatio):
    cl3D = cl / (1 + 2 / aspectRatio)

    cDind = cl3D **2 / (np.pi * aspectRatio)
    cDtot = cd + cDind
    return cl3D, cDtot

def dynamic_pressure(Sea, boatSpeed):
    return 0.5 * Sea.waterDensity * boatSpeed**2

def interpolate_wing_coefficients(xf):
    """
    Precompute cL and cD of a profile using xfoil implementation as functions
    of angle of attack and reynolds number. In this way the xfoil package is
    called only once at the init phase, and it does not slow down following
    modules.
    The interpolating function is created using a spline 2D function provided
    by the scipy package.
    """
    # Please BE AWARE that the interval of reynolds number and angles of attack,
    # as well as the interval's steps, have been carefully selected to guarantee
    # accurate results for the minimum run time.
    # Don't change this arrays if you don't know what you're doing!
    Re = np.linspace(2e4, 1e6, 10) # list of reynolds number that the interpolating
                                   # function will cover
    alfa = np.arange(-5, 5.1, 0.5) # list of angles of attack that the interpolating
                                   # function will cover


    # Initialize lists:
    clalfa = []
    clre   = []
    cdalfa = []
    cdre   = []

    for R in Re:
        xf.Re = R # sets the reynolds number for the analysis
        for a in alfa:
            cl, cd, cm, cp = xf.a(a) # xfoil analysis
            clalfa.append(cl)
            cdalfa.append(cd)
        clre.append(clalfa)
        cdre.append(cdalfa)
        clalfa = []
        cdalfa = []

    x = alfa
    y = Re
    zlift = np.array(clre).T
    zdrag = np.array(cdre).T
    # xfoil analysis often leads to nan values on the results. This seriously affects
    # the interpolating function. To avoid that, nan values are eliminated and
    # missing values are filled in based on nearby values via a linear interpolation.
    mask = np.isnan(zlift)
    zlift[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), zlift[~mask])
    mask = np.isnan(zdrag)
    zdrag[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), zdrag[~mask])

    # Scipy interpolating function
    clFunction = RectBivariateSpline(x,y,zlift)
    cdFunction = RectBivariateSpline(x,y,zdrag)

    return clFunction, cdFunction
