# -*- coding: utf-8 -*-
"""
funzione per calcolare la resistenza residua dello scafo basata sul modello
aggiornato di Delft. Per reference, leggi
A BARE HULL RESISTANCE PREDICTION METHOD DERIVED FROM THE RESULTS
OF THE DELFT SYSTEMATIC YACHT HULL SERIES EXTENDED TO HIGHER
SPEEDS by J A Keuning and M Katgert
"""
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

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
