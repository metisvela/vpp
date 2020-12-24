"""
classe per la raccolta dei dati di una barca e delle sue funzioni
di resistenza e stabilità
"""
from hull_resistances import residuary_resistance, frictional_resistance, added_residuary_res
from foil_geometry import geometry_v_foil
from stabilita_scafo import hull_stability
import pandas as pd
import numpy as np
from utilities import dynamic_pressure, lift_coefficients_2D, lift_coefficients_3D, interpolate_wing_coefficients
from foil_forces_calculation import forces_v_foil
from scipy.interpolate import interp1d
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
        self.camber = 2 #gradi di camber EFFETTIVO
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

    def foil_forces(self, Boat, Sea, boatSpeed):
        # Recall the current geometries based on roll angle and foil characteristics
        gm               = geometry_v_foil(self, Boat, Sea, boatSpeed)
        AoATip           = self.camber * np.cos(np.radians(self.gamma2)) # degrees, effective angle of the tip
        AoAStrut         = self.camber * np.cos(np.radians(self.gamma1)) # degrees, effective angle of the strut
        reynolds         = boatSpeed*self.chord / Sea.cinematicViscosity # np array storing all reynolds number
        aspectRatioStrut = gm['strut span'] / self.chord
        aspectRatioTip   = gm['tip span'] / self.chord

        # Calculate the cl, cd for both strut and tip
        liftCoeffStrut2D = self.liftFunction(AoAStrut, reynolds.flatten())
        dragCoeffStrut2D = self.dragFunction(AoAStrut, reynolds.flatten())
        liftCoeffTip2D   = self.liftFunction(AoATip  , reynolds.flatten())
        dragCoeffTip2D   = self.dragFunction(AoATip  , reynolds.flatten())

        liftCoeffStrut3D, dragCoeffStrut3D = lift_coefficients_3D(liftCoeffStrut2D, dragCoeffStrut2D, aspectRatioStrut)
        liftCoeffTip3D,   dragCoeffTip3D   = lift_coefficients_3D(liftCoeffTip2D,   dragCoeffTip2D,   aspectRatioTip)

        dynPressure = dynamic_pressure(Sea, boatSpeed)

        # Calculate foil forces:
        rightMom, lift, drag, leeway = forces_v_foil(self, Boat, liftCoeffStrut3D, liftCoeffTip3D, dragCoeffStrut3D, dragCoeffTip3D, gm, dynPressure)
        return rightMom, lift, drag, leeway

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

class Sails:
    def __init__(self, sailsDict):
        self.EHM = sailsDict['EHM']
        self.EMDC = sailsDict['EMDC']
        self.BAD = sailsDict['BAD']
        self.BMAX = sailsDict['BMAX']
        self.FA = sailsDict['FA']
        self.H    = sailsDict['H']
        self.F    = sailsDict['F']
        self.Am   = sailsDict['Am']
        self.Aj   = sailsDict['Aj']
        self.Ag   = sailsDict['Ag']
        self.genAngle = sailsDict['genAngle']
        self.An   = self.Am + self.Aj # upwind sail area
        return

    def sails_coefficients(self, AWA):
        """
        Utilizes Hazen model to interpolate cL and cD values provided by ORC
        VPP.
        
        Input:
        AWA : apparent wind angle

        Returns:
        advCoeff : coefficient of forward movement, derived from lift and drag coefficients
        capCoeff : similarly obtained from cL and cD, coefficient of capsizing
        """
        if self.H==1:
            BetaM =	np.array([0, 7, 9, 12, 28, 60, 90, 120, 150, 180])   # Angoli per la randa rispetto al Vento Apparente
            # Main sail coefficients interpolation
            C_LM =	np.array([0.00000, 0.86207, 1.05172, 1.16379, 1.34698, 1.35345,
                              1.26724, 0.93103, 0.38793, -0.11207])
            C_DM =	np.array([0.04310, 0.02586, 0.02328, 0.02328, 0.03259,
                              0.11302, 0.38250, 0.96888, 1.31578, 1.34483])
            # Jib coefficients interpolation
            BetaJ = np.array([7, 15, 20, 27, 50, 60, 100, 150, 180])
            C_LJ = np.array([0.00000, 1.00000, 1.37500, 1.45000, 1.45000, 1.25000,
                             0.40000, 0.00000, -0.10000])
            C_DJ = np.array([0.05000, 0.03200, 0.03100, 0.03700, 0.25000, 0.35000,
                             0.73000, 0.95000, 0.90000])
        elif self.H==2:
        	# Cl e Cd 2D Randa HIGH
            BetaM =	np.array([0, 7, 9, 12, 28, 60, 90, 120, 150, 180])
            C_DM =	np.array([0.03448, 0.01724, 0.01466, 0.01466, 0.02586,
                              0.11302, 0.38250, 0.96888, 1.31578, 1.34483])
            C_LM =	np.array([0.00000, 0.94828, 1.13793, 1.25000, 1.42681,
                              1.38319, 1.26724, 0.93103, 0.38793, -0.11207])

        	# Cl e Cd 2D Fiocco HIGH
            BetaJ =	np.array([7, 15, 20, 27, 50, 60, 100, 150, 180])
            C_DJ =	np.array([0.05000, 0.03200, 0.03100, 0.03700,
                              0.25000, 0.35000, 0.73000, 0.95000, 0.90000])
            C_LJ =	np.array([0.00000, 1.10000, 1.47500, 1.50000, 1.45000,
                              1.25000, 0.40000, 0.00000, -0.10000])
        else:
            print('H factor is wrong, must be either 1 or 2. Refer to ORC documentation for more info')

        # Gennaker coefficients:
        BetaG =	[0, 28, 41, 50, 60, 67, 75, 100, 115, 130, 150, 170, 180]

        C_LG =	np.array([0,0.01830, 0.73500, 0.94666, 1.08342, 1.10494, 1.09059, 0.95427,
                 0.81077, 0.60987, 0.32287, 0.10762, 0.00000])
        C_DG =	np.array([0,0.16215, 0.25184, 0.32502, 0.40897, 0.45920, 0.50225,
                          0.59839, 0.65292, 0.67086, 0.67086, 0.67086, 0.67086])

        # Cut off to zero gennaker coefficients when apparent wind angle is
        # below a threshold:
        C_LG = [0 if BetaG[i]<self.genAngle else C_LG[i] for i in range(len(BetaG))]
        C_DG = [0 if BetaG[i]<self.genAngle else C_DG[i] for i in range(len(BetaG))]

        # Interpolating functions for all sails
        C_Lm  = interp1d(BetaM, C_LM) # lift function main
        C_Lj  = interp1d(BetaJ, C_LJ) # lift function jib
        C_Dm  = interp1d(BetaM, C_DM) # drag function main
        C_Dj  = interp1d(BetaJ, C_DJ) # drag function jib
        C_Lg =	interp1d(BetaG, C_LG) # lift function gennaker
        C_Dg =	interp1d(BetaG, C_DG) # drag function gennaker

        # Compute lift and drag coefficients for datum AWA
        liftCoeff = (C_Lm(AWA)*self.Am + C_Lj(AWA)*self.Aj + C_Lg(AWA)*self.Ag) * self.F * 1/self.An
            # viscous drag
        viscDragCoeff = (C_Dm(AWA)*self.Am + C_Dj(AWA)*self.Aj + C_Dg(AWA)*self.Ag) * self.F * 1/self.An

        # Pick the correct aspect ratio
        if AWA<30: # Aspect ratio for close haul
            AR = (1.1 + (self.EHM+self.FA))**2 / self.An
        elif 30<=AWA<=180: # aspect ratio for largest runs
            AR = (1.1 + self.EHM)**2 / self.An


        indDragCoeff = (liftCoeff**2) * (1/(np.pi * AR) + 0.005) # Induced drag of sails:
        rigDragCoeff = 1.13 * (self.BMAX*self.FA + self.EHM*self.EMDC)/self.An # drag of rig
        # Total drag coefficient
        dragCoeff = viscDragCoeff + indDragCoeff + rigDragCoeff
        # Total lift coefficient:
        # liftCoeff

        # Compute coefficient of capsizing and forward-moving
        AWA = np.radians(AWA)
        capCoeff = liftCoeff*np.cos(AWA) + dragCoeff*np.sin(AWA) # capsizing coefficient
        advCoeff = liftCoeff*np.sin(AWA) - dragCoeff*np.cos(AWA) # forward movement coefficient
        return advCoeff, capCoeff
