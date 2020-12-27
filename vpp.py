# -*- coding: utf-8 -*-
"""
Establishes cycles to go through for a static VPP.
Created on Sun Dec 27 10:36:00

@author: giova
"""
import scipy
import numpy as np


def vpp(Boat, Sails, Crew, Foil, Keel, Sea):
    AWAvector = np.arange(20, 180,1)
    TWS = 2
    BSlist = []
    for AWA in AWAvector: 
        a = 1 # boatSpeed lower limit for bisection
        b = TWS # boatSpeed upper limit for bisection
        BS = scipy.optimize.bisect(f,a,b, xtol=2e-12, rtol=0.1, args=(Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA))
        print('AWA: ', AWA)
        print('BS: ', BS)
        print('Boat.rollAngle', Boat.rollAngle)
        print('leeway angle', Boat.leewayAngle, '\n')
        BSlist.append(float(BS))
    return BSlist, AWAvector


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
        cMoment = (sMoment + kMoment
        - fMoment - hMoment)
    return kLeeway - sLeeway



def f(boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA):
    a = 1
    b = 10
    Boat.leewayAngle = scipy.optimize.bisect(g,a,b, xtol=2e-12, rtol=0.1, args=(boatSpeed, Boat, Sea, Crew, Sails, Keel, Foil, TWS, AWA))

    sMoment, sFwd, sLeeway          = Sails.sail_forces(Boat, Sea, boatSpeed, TWS, AWA)
    kLeeway, kLift, kMoment, kDrag  = Keel.keel_forces(Boat, Sea, boatSpeed)
    fMoment, fLift, fDrag, fLeeway  = Foil.foil_forces(Boat, Sea, boatSpeed)
    hDrag                           = Boat.hull_resistance(boatSpeed, Sea)
    return sFwd - hDrag - fDrag - kDrag






