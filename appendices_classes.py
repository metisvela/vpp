"""
Questa è una classe per raccogliere tutti i dati e le funzioni necessarie
alle proprietà della deriva
"""
import numpy as np
from xfoil import XFoil
from xfoil.test import naca0012

class Keel:

    def __init__(self):
        self.meanChord       = None
        self.span            = None
        self.area            = self.chord * self.span
        self.aspectRatio     = self.span / self.meanChord
        self.profile         = XFoil()
        self.profile.airfoil = naca0012

        return

    def keel_lift(self, Boat, Sea, leewayAngle, boatSpeed):
        cL2D = np.pi**2 / 10 * leewayAngle # in degrees
        cL = cL2D / (1 + 2/self.aspectRatio) * leewayAngle
        effectiveLeewayAngle = leewayAngle * np.cos(np.radians(Boat.rollAngle))
        lift = 0.5 * Sea.waterDensity * boatSpeed * self.area * cL
        return lift

    def keel_resistance(self, Boat, Sea, leewayAngle, boatSpeed):
        cL2D = np.pi**2 / 10 * leewayAngle # in degrees
        cL = cL2D / (1 + 2/self.aspectRatio) * leewayAngle
        cDvisc = None
        cDind  = cL**2 / (np.pi * self.aspectRatio)
        cD = cDvisc + cDind
        effectiveLeewayAngle = leewayAngle * np.cos(np.radians(Boat.rollAngle))
        keelResistance = 0.5 * Sea.waterDensity * boatSpeed * self.area * cD
        return keelResistance 
