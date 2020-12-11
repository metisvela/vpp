"""
classe per la raccolta dei dati di una barca e delle sue funzioni
di resistenza e stabilità
"""
from hull_resistances import residuary_resistance, frictional_resistance, added_residuary_res
from stabilita_foil import stabilita_foil
from stabilita_scafo import hull_stability
import pandas as pd
import numpy as np
from utilities import lift_coefficients, dynamic_pressure

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
        self.cD     = foilsDict["cD"]
        return

    def foil_stability(self, Boat, Sea, boatSpeed):
        gm = stabilita_foil(self, Boat, Sea, boatSpeed)
        
        spanAvambraccioTheta = gm["strut span"]
        braccioAvambraccioTheta = gm["strut lever"]
        spanBraccioTheta = gm["tip span"]
        braccioBraccioTheta = gm["tip lever"]
        gamma1 = np.radians(self.gamma1) # angolo avambraccio
        gamma2 = np.radians(self.gamma2) # angolo braccio
        theta = np.radians(Boat.rollAngle)
        corda = self.chord
        pressioneDinamica = dynamic_pressure(Sea, boatSpeed)
        cL = self.cL #lift_coefficients(angleOfAttack, aspectRatio)
        
        # calcolo forze e momento raddrizzante
        forzaAvambraccio = pressioneDinamica * corda * spanAvambraccioTheta * cL *np.cos(gamma1) * np.cos(theta)
        momentoAvambraccio = forzaAvambraccio * braccioAvambraccioTheta
    
        forzaBraccio = pressioneDinamica * corda * spanBraccioTheta * cL *np.cos(gamma2) * np.cos(theta)
        momentoBraccio = forzaBraccio * braccioBraccioTheta
    
        momentoTotale = momentoAvambraccio + momentoBraccio
        return momentoTotale

    def foil_lift(self):
        print("nothing coded yet")
        return

    def foil_scarroccio(self):
        print("nothing coded yet")
        return

    def foil_resistance(self, Boat, Sea, boatSpeed):
        gm = stabilita_foil(self, Boat, Sea, boatSpeed)
        foilImmersedArea = (gm["strut span"] + gm["tip span"]) * self.chord
        dynPressure = dynamic_pressure(Sea, boatSpeed)
        foilResistance = dynPressure * foilImmersedArea * self.cD
        return foilResistance

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
        self.profile         = None

        return

    def keel_lift(self, Boat, Sea, angleOfAttack, boatSpeed):
        effectiveAoA = angleOfAttack * np.cos(np.radians(Boat.rollAngle))
        cL, _ = lift_coefficients(effectiveAoA, self.aspectRatio)
        dynPressure = dynamic_pressure(Sea, boatSpeed)
        lift = dynPressure * self.area * cL
        return lift

    def keel_resistance(self, Boat, Sea, angleOfAttack, boatSpeed):
        effectiveAoA = angleOfAttack * np.cos(np.radians(Boat.rollAngle))
        _, cD = lift_coefficients(effectiveAoA, self.aspectRatio)
        dynPressure = dynamic_pressure(Sea, boatSpeed)
        keelResistance = dynPressure * self.area * cD
        return keelResistance
