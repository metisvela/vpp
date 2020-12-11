"""
classe per la raccolta dei dati di una barca e delle sue funzioni
di resistenza e stabilità
"""
from hull_resistances import residuary_resistance, frictional_resistance, added_residuary_res
from stabilita_foil import stabilita_foil
from stabilita_scafo import hull_stability
import pandas as pd
import numpy as np
from xfoil import XFoil
from xfoil.test import naca0012

class Boat:
    def __init__(self, boatDict):
        self.lengthWaterl = boatDict["lengthWaterl"]
        self.beamWaterl = boatDict["beamWaterl"]
        self.canoeDraft = boatDict["canoeDraft"]
        self.midshipCoeff = boatDict["midshipCoeff"]
        self.displacement = boatDict["displacement"]
        self.prismCoeff = boatDict["prismCoeff"]
        self.longCentBuoy = boatDict["longCentBuoy"]
        self.longCentFlot = boatDict["longCentFlot"]
        self.wettedArea = boatDict["wettedArea"]
        self.foilBeam = boatDict["foilBeam"] # distanza fra la metà barca e il punto di uscita del foil
        self.foilHeight = boatDict["foilHeight"]# altezza del punto di uscita del foil rispetto al fondo barca
        self.freeboard = boatDict["freeboard"] # altezza del bordo libero rispetto alla DWL
        self.rollAngle = 0
        self.residuaryCoefficients = self.import_DSYHS_coefficients("residuary coefficients")
        self.heelCoefficients = self.import_DSYHS_coefficients("delta residuary coefficients")
        self.deltaWettedAreaCoefficients = self.import_DSYHS_coefficients("wetted area coefficients")
        return

    def hull_resistance(self, boatSpeed, Sea):
        resResistance = residuary_resistance(self, boatSpeed, Sea)
        fricResistance = frictional_resistance(self, boatSpeed, Sea)
        addResResistance = added_residuary_res(self, boatSpeed, Sea)
        return resResistance + fricResistance + addResResistance

    def hull_stability(self, Sea, boatSpeed):
        hullRightMom = hull_stability(self, Sea, boatSpeed)
        return hullRightMom

    def import_DSYHS_coefficients(self,sheet_name):
        coefficientTable = pd.read_excel("coefficientTableDSYHS.xlsx",sheet_name=sheet_name,index_col=0)
        return coefficientTable


class Crew:
    def __init__(self, crewDict):
        self.bowmanWeight = crewDict["bowmanWeight"]
        self.helmsmanWeight = crewDict["helmsmanWeight"]
        self.bowmanHeight   = crewDict["bowmanHeight"]
        self.helmsmanHeight = crewDict["helmsmanHeight"]
        self.crewWeight = self.bowmanWeight + self.helmsmanWeight
        return

    def crew_stability(self, Boat, Sea, boatSpeed):
        bowmRightMom = self.bowmanHeight / 2 * self.bowmanWeight * Sea.gravityConstant * np.cos(np.radians(Boat.rollAngle))
        helmRightMom = self.helmsmanHeight / 2 * self.helmsmanWeight * Sea.gravityConstant * np.cos(np.radians(Boat.rollAngle))
        crewRightMom = bowmRightMom + helmRightMom
        crewRightMom = crewRightMom * np.ones(shape=boatSpeed.size)
        return crewRightMom

class Foil:
    def __init__(self, foilsDict):
        self.gamma1 = foilsDict["gamma1"]
        self.gamma2 = foilsDict["gamma2"]
        self.chord  = foilsDict["chord"]
        self.cL     = foilsDict["cL"]
        return

    def foil_stability(self, Boat, Sea, boatSpeed):
        foilRightingMoment = stabilita_foil(self, Boat, Sea, boatSpeed)
        return foilRightingMoment

    def foil_lift(self):
        print("nothing coded yet")
        return

    def foil_scarroccio(self):
        print("nothing coded yet")
        return

    def foil_resistance(self):
        print("nothing coded yet")
        return

class Sea:
    def __init__(self, seaDict):
        self.waterDensity = seaDict["waterDensity"]
        self.cinematicViscosity = seaDict["cinematicViscosity"]
        self.gravityConstant = seaDict["gravityConstant"]
        return

class Keel:

    def __init__(self,meanChord, span):
        self.meanChord       = meanChord
        self.span            = span
        self.area            = self.meanChord * self.span
        self.aspectRatio     = self.span / self.meanChord
        self.profile         = XFoil()
        self.profile.airfoil = naca0012

        return

    def keel_lift(self, Boat, Sea, leewayAngle, boatSpeed):
        cL2D = 0.1 # slope of the 2D lift profile
        effectiveLeewayAngle = leewayAngle * np.cos(np.radians(Boat.rollAngle))
        cL = cL2D / (1 + 2/self.aspectRatio) * effectiveLeewayAngle
        lift = 0.5 * Sea.waterDensity * boatSpeed * self.area * cL
        return lift

    def keel_resistance(self, Boat, Sea, leewayAngle, boatSpeed):
        cL2D = 0.1 # slope of the 2D lift profile
        effectiveLeewayAngle = leewayAngle * np.cos(np.radians(Boat.rollAngle))
        cL = cL2D / (1 + 2/self.aspectRatio) * effectiveLeewayAngle
        cDvisc = 0
        cDind  = cL**2 / (np.pi * self.aspectRatio)
        cD = cDvisc + cDind
        keelResistance = 0.5 * Sea.waterDensity * boatSpeed * self.area * cD
        return keelResistance
