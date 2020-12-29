# -*- coding: utf-8 -*-
"""
Establishes cycles to go through for a static VPP.
Created on Sun Dec 27 10:36:00

@author: giova
"""
import scipy.optimize
import numpy as np
import sys


def vpp(Boat, Sails, Crew, Foil, Keel, Sea):
    AWAvector = np.arange(20, 100,1)
    TWAlist = []
    TWS = 5
    BSlist = []
    leeAngleList = []
    rollList = []
    for AWA in AWAvector:
        BSlow = 2
        BSup = 5

    try:
        BS = scipy.optimize.bisect(f,BSlow,BSup, rtol=0.001, args=(Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA))
    except ValueError:
        lowLimit = f(BS, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA)
        upLimit  = f(BS, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA)
        if np.sign(lowLimit) == -1 and np.sign(upLimit) == -1:
            raise Exception('The drag is too high. Try reducing the lower limit of the boat speed interval')sys.exit()
        elif np.sign(lowLimit) == 1 and np.sign(upLimit) == 1:
            raise Exception('The drag is too low. Try raising the upper limit of the boat speed interval')

        # Compute TWA for the polar graph:
        TWA,_ = Sails.velocity_triangle(AWA, BS, TWS)
        BSlist.append(float(BS))
        TWAlist.append(float(TWA))
        leeAngleList.append(float(Boat.leewayAngle))
        rollList.append(float(Boat.rollAngle))
    return BSlist, TWAlist, leeAngleList, rollList


def g(delta, boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA):
    Boat.leewayAngle = delta

    Boat.rollAngle = 0
    maxCrewStability = Crew.crew_stability(Boat, Sea, boatSpeed)
    cMoment = 100000
    while cMoment > maxCrewStability:
        Boat.rollAngle += 1
        TWA, AWS                        = Sails.velocity_triangle(AWA, boatSpeed, TWS)
        kLeeway, kLift, kMoment, kDrag  = Keel.keel_forces(Boat, Sea, boatSpeed)
        fMoment, fLift, fDrag, fLeeway  = Foil.foil_forces(Boat, Sea, boatSpeed)
        sMoment, sFwd, sLeeway          = Sails.sail_forces(Boat, Sea, boatSpeed, TWS, AWA)
        hMoment                         = Boat.hull_stability(Sea, boatSpeed)
        # hDrag                           = hull_resistance(boatSpeed, Sea)
        cMoment = (sMoment
        - fMoment - hMoment)
    return kLeeway - sLeeway



def f(boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA):
    deltalow = -1
    deltaup = 10
    delta = (deltalow+deltaup)/2

    try:
        Boat.leewayAngle = scipy.optimize.bisect(g,deltalow,deltaup, rtol=0.001, args=(boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA))
    except ValueError:
        lowLimit = g(deltalow, boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA)
        upLimit  = g(deltaup, boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA)
        if np.sign(lowLimit) == -1 and np.sign(upLimit) == -1:
            raise Exception('The keel does not produce enough lift. Change leeway angle '
                            + 'interval or increase keel area. Also consider raising the '
                            + 'lower limit of boat speed interval')
            sys.exit()
        elif np.sign(lowLimit) == 1 and np.sign(upLimit) == 1:
            raise Exception('The sails do not produce enough sideward force. Consider '
                            + 'checking the geometry or check xfoil results for the keel for errors.')

    sMoment, sFwd, sLeeway          = Sails.sail_forces(Boat, Sea, boatSpeed, TWS, AWA)
    kLeeway, kLift, kMoment, kDrag  = Keel.keel_forces(Boat, Sea, boatSpeed)
    fMoment, fLift, fDrag, fLeeway  = Foil.foil_forces(Boat, Sea, boatSpeed)
    hDrag                           = Boat.hull_resistance(boatSpeed, Sea)
    return sFwd - hDrag - fDrag - kDrag
