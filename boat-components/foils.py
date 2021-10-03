# -*- coding: utf-8 -*-

from utilities import dynamic_pressure, lift_coefficients_3D, interpolate_wing_coefficients
from scipy.interpolate import interp1d
import pickle
import numpy as np

"""
Script per calcolare il momento raddrizzante dato da un foil a v
Created on Fri Nov 13 17:22:48 2020

@author: giova
"""

def geometry_v_foil(Foil, Boat, Sea, boatSpeed):
    # stabilità foil
    # Parametri di progetto
    B = Boat.foilBeam #larghezza scafo
    Hs = Boat.foilHeight # altezza buco di uscita scassa
    Tc = Boat.canoeDraft # pescaggio
    gamma1 = np.radians(Foil.gamma1) # angolo avambraccio
    gamma2 = np.radians(Foil.gamma2) # angolo braccio
    FB = Boat.freeboard # altezza bordo libero da DWL

    # costanti idrodinamiche
    pressioneDinamica = 0.5 * Sea.waterDensity * boatSpeed**2 #kg/(m2 s) pressione dinamica a 2.5 m/s
    corda = Foil.chord #cm, corda di braccio e avambraccio

    theta = np.radians(Boat.rollAngle)

    lunghezzaMassima = B/np.cos(gamma1) + (Hs - Tc)/np.sin(gamma1)

    immersione = B * np.tan(gamma1) - FB + Hs
    if immersione < 0:
        print("il foil è fuori dall'acqua! Cambia la geometria delle scasse")
        momentoTotale = 0
    else:
        spanAvambraccio = immersione / np.sin(gamma1)
        spanBraccio = immersione / np.sin(gamma2)

        braccioAvambraccio = lunghezzaMassima * (np.cos(gamma1))**2 + spanAvambraccio/2
        braccioBraccio = lunghezzaMassima * np.cos(gamma1) * np.sin(gamma1) + spanBraccio/2

        # variazione con theta delle geometrie dell'avambraccio
        deltaSpanAvambraccio = (lunghezzaMassima * np.cos(gamma1) * np.sin(theta))/(np.sin(gamma1 + theta))
        braccioAvambraccioTheta = braccioAvambraccio - deltaSpanAvambraccio/2
        spanAvambraccioTheta = spanAvambraccio + deltaSpanAvambraccio


        # variazione con theta delle geometrie del braccio
        immersioneTheta = spanAvambraccioTheta * np.sin(gamma1 + theta)

        spanBraccioTheta = 0.5 # immersioneTheta / (np.sin(gamma2 - theta))
        braccioBraccioTheta = lunghezzaMassima * np.cos(gamma1) * np.sin(gamma1) * spanBraccioTheta/2

        geometries = {"strut length" : lunghezzaMassima,
                      "immersion"    : immersioneTheta,
                      "strut span"   : spanAvambraccioTheta,
                      "tip span"     : spanBraccioTheta,
                      "strut lever"  : braccioAvambraccioTheta,
                      "tip lever"    : braccioBraccioTheta,
                      }


    return geometries

"""
Functions to calculate lift, drag, righting moment and scarroccio
Created on Wed Dec 23 10:22:48 2020

@author: giova
"""

def forces_v_foil(Foil, Boat, clStrut, clTip, cdStrut, cdTip, gm, dynPressure):
    """
    Compute righting moment of a v foil
    Input:
    Foil: Foil class
    Boat: Boat class
    clStrut: 3D lift coefficient of the strut
    clTip  : 3D lift coefficient of the tip
    cdStrut: 3D drag coefficient of the strut
    cdTip  : 3D drag coefficient of the tip

    Return:
    rightMoment  : righting moment of the foil in Newton per meter [Nm]
    lift         : lift force of the foil in Newton
    drag         : drag force of the foil in Newton
    leeway       : leeway force of strut and tip combined;
                   positive if it counteracts the wind,
                   negative if it adds to it
    """
    rollAngle  = np.radians(Boat.rollAngle)
    strutSpan  = gm['strut span']
    tipSpan    = gm['tip span']
    strutLever = gm['strut lever']
    tipLever   = gm['tip lever']
    gamma1     = np.radians(Foil.gamma1)
    gamma2     = np.radians(Foil.gamma2)

    forceStrut = dynPressure * Foil.chord * strutSpan * clStrut * np.cos(rollAngle)
    forceTip = dynPressure * Foil.chord * tipSpan * clTip * np.cos(rollAngle)

    # Compute righting moments:
    rightMomStrut = forceStrut * strutLever
    rightMomTip   = forceTip * tipLever
    rightMom      = rightMomStrut + rightMomTip

    # Compute lifts:
    lift = forceStrut*np.cos(gamma1+rollAngle) + forceTip*np.cos(gamma2-rollAngle)

    # Compute resistance:
    dragStrut   = dynPressure * Foil.chord * strutSpan * cdStrut * np.cos(rollAngle)
    dragTip     = dynPressure * Foil.chord * tipSpan   * cdTip   * np.cos(rollAngle)
    drag        = dragStrut + dragTip

    # Compute leeway:
    leeway = forceTip*np.sin(gamma2-rollAngle) - forceStrut*np.sin(gamma1+rollAngle)

    return rightMom, lift, drag, leeway

###############################################################################
# Master Foil class starts below
###############################################################################
class Foil:
    def __init__(self, foilsDict, wingProfiles):
        # See the corresponding dictionary on input_data.py for information
        # on this parameters.
        self.gamma1 = foilsDict["gamma1"]
        self.gamma2 = foilsDict["gamma2"]
        self.chord  = foilsDict["chord"]
        self.camber = foilsDict["keying"] #gradi di camber EFFETTIVO
        self.profile = foilsDict['profile']
        # Load the precomputed interpolation functions from folder
        # "foil-profiles"
        handler = open('foil-profiles/'+self.profile, 'rb')
        f_list = pickle.load(handler)
        self.liftFunction = f_list[0]
        self.dragFunction = f_list[1]
        #self.liftFunction, self.dragFunction = interpolate_wing_coefficients(self.xf)
        return

    def foil_forces(self, Boat, Sea, boatSpeed):
        # Recall the current geometries based on roll angle and foil characteristics
        gm               = geometry_v_foil(self, Boat, Sea, boatSpeed)
        AoATip           = self.camber * np.cos(np.radians(self.gamma2)) \
                            - Boat.leewayAngle*np.sin(np.radians(self.gamma2))# degrees, effective angle of the tip
        AoAStrut         = self.camber * np.cos(np.radians(self.gamma1)) \
                            - Boat.leewayAngle*np.sin(np.radians(self.gamma1)) # degrees, effective angle of the strut
        reynolds         = boatSpeed*self.chord / Sea.cinematicViscosity # np array storing all reynolds number
        aspectRatioStrut = gm['strut span'] / self.chord
        aspectRatioTip   = gm['tip span'] / self.chord

        # Calculate the cl, cd for both strut and tip
        try:
            liftCoeffStrut2D = self.liftFunction(AoAStrut, reynolds.flatten())
            dragCoeffStrut2D = self.dragFunction(AoAStrut, reynolds.flatten())
            liftCoeffTip2D   = self.liftFunction(AoATip  , reynolds.flatten())
            dragCoeffTip2D   = self.dragFunction(AoATip  , reynolds.flatten())
        except AttributeError:
            liftCoeffStrut2D = self.liftFunction(AoAStrut, reynolds)
            dragCoeffStrut2D = self.dragFunction(AoAStrut, reynolds)
            liftCoeffTip2D   = self.liftFunction(AoATip  , reynolds)
            dragCoeffTip2D   = self.dragFunction(AoATip  , reynolds)
        liftCoeffStrut3D, dragCoeffStrut3D = lift_coefficients_3D(liftCoeffStrut2D, dragCoeffStrut2D, aspectRatioStrut)
        liftCoeffTip3D,   dragCoeffTip3D   = lift_coefficients_3D(liftCoeffTip2D,   dragCoeffTip2D,   aspectRatioTip)

        dynPressure = dynamic_pressure(Sea, boatSpeed)

        # Calculate foil forces:
        rightMom, lift, drag, leeway = forces_v_foil(self, Boat, liftCoeffStrut3D, liftCoeffTip3D, dragCoeffStrut3D, dragCoeffTip3D, gm, dynPressure)
        return rightMom, lift, drag, leeway
