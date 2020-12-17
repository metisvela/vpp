"""
File dove raccogliere tutte le funzioni utilizzate più di una volta nel codice per
facilitarne la modifica
"""
import xfoil
import numpy as np

def lift_coefficients_2D(Foil, angleOfAttack, boatSpeed, Sea):
    """
    Il modo per usare xfoil è il seguente:
    - Salvo le coordinate del profilo in due vettori numpy x,y
    - profilo = xfoil.model.Airfoil(x,y)
    - xf = xfoil.XFoil()
    - xf.airfoil = profilo
    - xf.repanel()
    - xf.Re = xxx
    - xf.n_crit = 3
    - cl, cd, cm, cp = xf.a()
    """
    # calcolo del coefficiente di lift e drag 2D:
    x = np.array([1.     , 0.99893, 0.99572, 0.99039, 0.98296, 0.97347, 0.96194,
       0.94844, 0.93301, 0.91573, 0.89668, 0.87592, 0.85355, 0.82967,
       0.80438, 0.77779, 0.75   , 0.72114, 0.69134, 0.66072, 0.62941,
       0.59755, 0.56526, 0.5327 , 0.5    , 0.4673 , 0.43474, 0.40245,
       0.37059, 0.33928, 0.30866, 0.27886, 0.25   , 0.22221, 0.19562,
       0.17033, 0.14645, 0.12408, 0.10332, 0.08427, 0.06699, 0.05156,
       0.03806, 0.02653, 0.01704, 0.00961, 0.00428, 0.00107, 0.     ,
       0.00107, 0.00428, 0.00961, 0.01704, 0.02653, 0.03806, 0.05156,
       0.06699, 0.08427, 0.10332, 0.12408, 0.14645, 0.17033, 0.19562,
       0.22221, 0.25   , 0.27886, 0.30866, 0.33928, 0.37059, 0.40245,
       0.43474, 0.4673 , 0.5    , 0.5327 , 0.56526, 0.59755, 0.62941,
       0.66072, 0.69134, 0.72114, 0.75   , 0.77779, 0.80438, 0.82967,
       0.85355, 0.87592, 0.89668, 0.91573, 0.93301, 0.94844, 0.96194,
       0.97347, 0.98296, 0.99039, 0.99572, 0.99893, 1.
       ])
    y = np.array([0.     ,  0.00023,  0.00086,  0.00193,  0.00341,  0.00534,
            0.00766,  0.01035,  0.01342,  0.01681,  0.02053,  0.02447,
            0.02864,  0.03298,  0.03747,  0.042  ,  0.04652,  0.05089,
            0.05511,  0.05905,  0.06275,  0.06608,  0.06911,  0.07174,
            0.07409,  0.07596,  0.0775 ,  0.07845,  0.07898,  0.07888,
            0.07838,  0.0772 ,  0.07565,  0.07339,  0.07081,  0.06754,
            0.06404,  0.05989,  0.05569,  0.05086,  0.04609,  0.04056,
            0.03523,  0.02948,  0.02368,  0.01775,  0.01144,  0.00574,
            0.     , -0.002  , -0.00436, -0.00691, -0.0097 , -0.01247,
           -0.01481, -0.0173 , -0.01913, -0.0211 , -0.02246, -0.02377,
           -0.02447, -0.02503, -0.025  , -0.02475, -0.02389, -0.02275,
           -0.021  , -0.01895, -0.01637, -0.01357, -0.01045, -0.00731,
           -0.00405, -0.00092,  0.00217,  0.00496,  0.00748,  0.00951,
            0.01118,  0.01241,  0.01329,  0.01373,  0.01381,  0.01347,
            0.0128 ,  0.01186,  0.01072,  0.00941,  0.00802,  0.00659,
            0.00515,  0.00377,  0.00251,  0.00151,  0.0007 ,  0.00015,
            0.])
    cllist = []
    cdlist = []
    wortmann = xfoil.model.Airfoil(x,y)
    xf = xfoil.XFoil()
    xf.airfoil = wortmann
    xf.repanel()
    for v in boatSpeed:
        xf.Re = v * Foil.chord / Sea.cinematicViscosity
        xf.n_crit = 3
        xf.max_iter = 40
        cl, cd, _, __ = xf.a(angleOfAttack)
        cllist.append(cl)
        cdlist.append(cd)
    return np.array(cllist), np.array(cdlist)

def lift_coefficients_3D(cl, cd, aspectRatio):
    cl3D = cl / (1 + 2 / aspectRatio)

    cDind = cl3D **2 / (np.pi * aspectRatio)
    cDtot = cd + cDind
    return cl3D, cDtot

def dynamic_pressure(Sea, boatSpeed):
    return 0.5 * Sea.waterDensity * boatSpeed**2
