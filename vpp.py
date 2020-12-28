# -*- coding: utf-8 -*-
"""
Establishes cycles to go through for a static VPP.
Created on Sun Dec 27 10:36:00

@author: giova
"""
import scipy.optimize
import numpy as np


def vpp(Boat, Sails, Crew, Foil, Keel, Sea):
    AWAvector = np.arange(20, 100,1)
    TWAlist = []
    TWS = 6.5
    BSlist = []
    leeAngleList = []
    rollList = []
    for AWA in AWAvector:
        BSguess = 5
        BSlow = BSguess - 2
        BSup = BSguess + 2
        BS = scipy.optimize.bisect(f,BSlow,BSup, rtol=0.001, args=(Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA))

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
    deltaup = 5
    delta = (deltalow+deltaup)/2
    Boat.leewayAngle = scipy.optimize.bisect(g,deltalow,deltaup, rtol=0.001, args=(boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA))

    sMoment, sFwd, sLeeway          = Sails.sail_forces(Boat, Sea, boatSpeed, TWS, AWA)
    kLeeway, kLift, kMoment, kDrag  = Keel.keel_forces(Boat, Sea, boatSpeed)
    fMoment, fLift, fDrag, fLeeway  = Foil.foil_forces(Boat, Sea, boatSpeed)
    hDrag                           = Boat.hull_resistance(boatSpeed, Sea)
    return sFwd - hDrag - fDrag - kDrag
