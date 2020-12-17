"""
classe per la raccolta dei dati di una barca e delle sue funzioni
di resistenza e stabilità
"""
from hull_resistances import residuary_resistance, frictional_resistance, added_residuary_res
from stabilita_foil import stabilita_foil
from stabilita_scafo import hull_stability
import pandas as pd
import numpy as np
from utilities import dynamic_pressure, lift_coefficients_2D, lift_coefficients_3D, interpolate_wing_coefficients
import xfoil 

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
        bowmRightMom = (self.bowmanHeight / 2 + 1.05) * self.bowmanWeight * Sea.gravityConstant * np.cos(np.radians(Boat.rollAngle))
        helmRightMom = (self.helmsmanHeight / 2 + 1.05) * self.helmsmanWeight * Sea.gravityConstant * np.cos(np.radians(Boat.rollAngle))
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
        self.camber = 2 #gradi di camber
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

        profilo = xfoil.model.Airfoil(x,y)
        self.xf = xfoil.XFoil()
        self.xf.airfoil = profilo
        self.xf.n_crit = 3
        self.xf.repanel()

        self.liftFunction, self.dragFunction = interpolate_wing_coefficients(self.xf)
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

        aspectRatioTip = spanBraccioTheta / self.chord
        aspectRatioStrut = spanAvambraccioTheta / self.chord
        AoATip = self.camber * np.cos(np.radians(self.gamma2)) # in gradi, angolo effettivo della tip
        AoAStrut = self.camber * np.cos(np.radians(self.gamma1)) # in gradi, angolo effettivo della strut
        reynolds = boatSpeed * self.chord / Sea.cinematicViscosity

        liftCoeffStrut2D = self.liftFunction(AoAStrut, reynolds.flatten())
        dragCoeffStrut2D = self.dragFunction(AoAStrut, reynolds.flatten())
        liftCoeffTip2D   = self.liftFunction(AoATip  , reynolds.flatten())
        drafCoeffTip2D   = self.dragFunction(AoATip  , reynolds.flatten())

        liftCoeffStrut3D, dragCoeffStrut3D = lift_coefficients_3D(liftCoeffStrut2D, dragCoeffStrut2D, aspectRatioStrut)
        liftCoeffTip3D,   dragCoeffTip3D   = lift_coefficients_3D(liftCoeffTip2D,   drafCoeffTip2D,   aspectRatioTip)

        pressioneDinamica = dynamic_pressure(Sea, boatSpeed)
        # calcolo forze e momento raddrizzante
        forzaAvambraccio = pressioneDinamica * self.chord * spanAvambraccioTheta * liftCoeffStrut3D * np.cos(theta)
        momentoAvambraccio = forzaAvambraccio * braccioAvambraccioTheta

        forzaBraccio = pressioneDinamica * corda * spanBraccioTheta * liftCoeffTip3D * np.cos(theta)
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

        tipSpan = gm["tip span"]
        strutSpan = gm["strut span"]
        areaTip = tipSpan * self.chord
        areaStrut = strutSpan * self.chord
        aspectRatioTip = tipSpan / self.chord
        aspectRatioStrut = strutSpan / self.chord
        AoATip = self.camber * np.cos(np.radians(self.gamma2)) # in gradi, angolo effettivo della tip
        AoAStrut = self.camber * np.cos(np.radians(self.gamma1)) # in gradi, angolo effettivo della strut

        reynolds = boatSpeed * self.chord / Sea.cinematicViscosity


        dragCoeffStrut2D = self.dragFunction(AoAStrut, reynolds.flatten())

        drafCoeffTip2D   = self.dragFunction(AoATip  , reynolds.flatten())

        _, dragCoeffStrut3D   = lift_coefficients_3D(1, dragCoeffStrut2D, aspectRatioStrut)
        _,   dragCoeffTip3D   = lift_coefficients_3D(1,   drafCoeffTip2D,   aspectRatioTip)

        strutResistance = dynPressure * areaStrut * dragCoeffStrut3D * np.cos(np.radians(Boat.rollAngle))
        tipResistance = dynPressure * areaTip     * dragCoeffTip3D   * np.cos(np.radians(Boat.rollAngle))
        return strutResistance + tipResistance

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
