# -*- coding: utf-8 -*-
"""
Functions to calculate lift, drag, righting moment and scarroccio
Created on Wed Dec 23 10:22:48 2020

@author: giova
"""
import numpy as np

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
