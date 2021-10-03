# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

"""
Functions below allow for calculation of total resistance of a hull based on the
Delft model. For reference, read:
A BARE HULL RESISTANCE PREDICTION METHOD DERIVED FROM THE RESULTS
OF THE DELFT SYSTEMATIC YACHT HULL SERIES EXTENDED TO HIGHER
SPEEDS by J A Keuning and M Katgert
"""

def residuary_resistance(Boat, boatSpeed, Sea):
    a = interpolate_DSYHS_coefficients(Boat.residuaryCoefficients)
    boatSpeedFroude = boatSpeed / (Sea.gravityConstant * Boat.lengthWaterl)**0.5
    resResistance = (Boat.maxDisplacement * Sea.waterDensity * Sea.gravityConstant *
                    (a[0](boatSpeedFroude) +
                    ((
                    a[1](boatSpeedFroude) * Boat.longCentBuoy / Boat.lengthWaterl +
                    a[2](boatSpeedFroude) * Boat.prismCoeff +
                    a[3](boatSpeedFroude) * Boat.maxDisplacement**(2/3) / wetted_area_change_theta(Boat) +
                    a[4](boatSpeedFroude) * Boat.beamWaterl / Boat.lengthWaterl +
                    a[5](boatSpeedFroude) * Boat.longCentBuoy / Boat.longCentFlot +
                    a[6](boatSpeedFroude) * Boat.beamWaterl / Boat.canoeDraft +
                    a[7](boatSpeedFroude) * Boat.midshipCoeff
                    ) * Boat.maxDisplacement**(1/3) / Boat.lengthWaterl)))
    return resResistance


def frictional_resistance(Boat, boatSpeed, Sea):
    reynoldsNumber = ((boatSpeed * 0.7 * Boat.lengthWaterl) /
                      Sea.cinematicViscosity)
    fricCoefficient = 0.075 / (np.log10(reynoldsNumber) - 2)**2
    wettedAreaTheta = wetted_area_change_theta(Boat)
    fricResistance = 0.5 * Sea.waterDensity * boatSpeed**2 * wettedAreaTheta * fricCoefficient

    rugosityFactor = 0.085
    f_BS = rugosityFactor * boatSpeed - rugosityFactor
    rugResistance = f_BS * fricResistance;
    return fricResistance+rugResistance

def added_residuary_res(Boat, boatSpeed, Sea):
    """
    calcolo della resistenza residua aggiunta dovuta all'angolo di rollio
    """
    u = interpolate_DSYHS_coefficients(Boat.heelCoefficients)
    boatSpeedFroude = boatSpeed / (Sea.gravityConstant * Boat.lengthWaterl)**0.5
    deltaResidResist = (Boat.maxDisplacement * Sea.waterDensity * Sea.gravityConstant *
                       (
                       u[0](boatSpeedFroude) +
                       u[1](boatSpeedFroude) * Boat.lengthWaterl / Boat.beamWaterl +
                       u[2](boatSpeedFroude) * Boat.beamWaterl / Boat.canoeDraft +
                       u[3](boatSpeedFroude) * (Boat.beamWaterl / Boat.canoeDraft)**2 +
                       u[4](boatSpeedFroude) * Boat.longCentBuoy +
                       u[5](boatSpeedFroude) * Boat.longCentBuoy**2
                       ))
    deltaResidResistTheta = deltaResidResist * 6 * np.radians(Boat.rollAngle) **1.7
    return deltaResidResistTheta

def wetted_area_change_theta(Boat):
    """
    stima del cambiamento di area bagnata al variare dell'angolo di rollio
    """
    if Boat.rollAngle >5:
        s = interpolate_DSYHS_coefficients(Boat.deltaWettedAreaCoefficients)

        wettedAreaTheta = (Boat.wettedArea *(
                                        1 + 1/100*(s[0](Boat.rollAngle) +
                                                  s[1](Boat.rollAngle) * Boat.lengthWaterl / Boat.canoeDraft +
                                                  s[2](Boat.rollAngle) * (Boat.beamWaterl / Boat.canoeDraft)**2 +
                                                  s[3](Boat.rollAngle) * Boat.midshipCoeff)
                                        )
                           )
    else:
        wettedAreaTheta = Boat.wettedArea
    return wettedAreaTheta

############################################
# MOVE THIS TO A UTILITIES PYTHON FILE
###########################################
def interpolate_DSYHS_coefficients(coeff_table):
    x = coeff_table.columns.to_numpy()
    interpolFunctions = [interp1d(x, coeff_table.iloc[i].to_numpy(), fill_value="extrapolate") for i in range(len(coeff_table))]
    return interpolFunctions

"""
Functions below calculate stability of the hull under navigation.
Most of the code is based on basic naval stability theorems and approximations,
that can be found in any book of basics, such as Larsson's Yacht Design.
"""

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

    elif (np.radians(Boat.rollAngle) > np.arctan((h-Tc)/b)) & (np.radians(Boat.rollAngle) < np.radians(30)) :
        theta = np.radians(Boat.rollAngle)
        x = b/2 - (h-Tc)*(np.tan(theta))**(-1)
        y = b/2 + Tc*(np.tan(theta))**(-1)
        z = h
        xG = (x**2 + y**2 + x*y)/(3*(x+y))
        yG = (z*(x+y))/(3*(x+y))
        VCB = yG
        HCB = b/2 - xG

        GZ = HCB*np.cos(theta) - (VCG - VCB) * np.sin(theta)

    else:
        raise Exception('Angoli troppo alti per la stabilità di scafo')
        GZ = 0


    rightingMoment = Boat.maxDisplacement * 1000 * Sea.gravityConstant * GZ # displacement da m3 a kg

    # try:
    #     rightingMoment = rightingMoment * np.ones(shape=boatSpeed.size) # adatto a vettore per il plotting
    # except AttributeError:
    #     pass

    return rightingMoment


###############################################################################
# Master boat class below
###############################################################################

class Boat:
    def __init__(self, boatDict):
        # See the corresponding dictionary on input_data.py for information
        # on this parameters.
        self.lengthWaterl = boatDict["lengthWaterl"]
        self.beamWaterl = boatDict["beamWaterl"]
        self.canoeDraft = boatDict["canoeDraft"]
        self.midshipCoeff = boatDict["midshipCoeff"]
        self.maxDisplacement = boatDict['displacement']
        self.displacement = None
        self.prismCoeff = boatDict["prismCoeff"]
        self.longCentBuoy = boatDict["longCentBuoy"]
        self.longCentFlot = boatDict["longCentFlot"]
        self.wettedArea = boatDict["wettedArea"]
        self.leewayAngle = 0
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
