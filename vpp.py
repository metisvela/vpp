# -*- coding: utf-8 -*-
"""
Establishes cycles to go through for a static VPP.
Created on Sun Dec 27 10:36:00

@author: giova
"""
from scipy.optimize import bisect, newton
from scipy.misc import derivative as der
import numpy as np


class Vpp():
    def __init__(self,Boat, Sea, Crew, Sails, Keel, Foil, AWAvector, TWS):
        self.AWAvector = AWAvector
        self.TWS       = TWS
        self.TWAlist   = []
        self.BSlist    = []
        self.leeAngleList = []
        self.rollList  = []
        self.convergenceMethod = 'newton'
        self.Boat = Boat
        self.Sea = Sea
        self.Crew = Crew
        self.Sails = Sails
        self.Keel = Keel
        self.Foil = Foil
        return


    def vpp_static(self):

        for AWA in self.AWAvector:
            if self.convergenceMethod == 'bisection':
                BSup  = self.TWS
                BSlow = self.TWS / 2.
                boatSpeed = bisect(self.f, BSlow, BSup, rtol=0.01, args=(AWA))
            elif self.convergenceMethod == 'newton':
                boatSpeedGuess = 0.75 * self.TWS
                boatSpeed = newton(self.f, boatSpeedGuess,
                                   tol=0.01, args=(AWA,),
                                   maxiter=100)
            else:
                raise Exception(self.convergenceMethod, 'is not a valid convergence method.')

            # Compute TWA for the polar graph:
            TWA,_ = self.Sails.velocity_triangle(AWA, boatSpeed, self.TWS)
            self.BSlist.append(float(boatSpeed))
            self.TWAlist.append(float(TWA))
            self.leeAngleList.append(float(self.Boat.leewayAngle))
            self.rollList.append(float(self.Boat.rollAngle))
        return

    def f(self, boatSpeed, AWA):

        if self.convergenceMethod == 'bisection':
            deltaUp  = -2
            deltaLow = 10
            self.Boat.leewayAngle = bisect(self.g, deltaUp, deltaLow,rtol = 0.01, args=(boatSpeed, AWA))
        elif self.convergenceMethod == 'newton':
            deltaGuess = 1 
            self.Boat.leewayAngle = newton(self.g, deltaGuess,
                                           tol = 0.01, args=(boatSpeed, AWA),
                                           maxiter = 100)
        else:
            raise Exception(self.convergenceMethod, ' is not a valid convergence method.')

        sMoment, sFwd, sLeeway          = self.Sails.sail_forces(self.Boat, self.Sea, boatSpeed, self.TWS, AWA)
        kLeeway, kLift, kMoment, kDrag  = self.Keel.keel_forces(self.Boat, self.Sea, boatSpeed)
        fMoment, fLift, fDrag, fLeeway  = self.Foil.foil_forces(self.Boat, self.Sea, boatSpeed)
        self.Boat.displacement = self.Boat.maxDisplacement - fLift/10000 # in m^3
        hDrag                           = self.Boat.hull_resistance(boatSpeed, self.Sea)
        return sFwd - hDrag - fDrag - kDrag


    def g(self, delta, boatSpeed, AWA):
        self.Boat.leewayAngle = delta

        self.Boat.rollAngle = -1
        maxCrewStability = self.Crew.crew_stability(self.Boat, self.Sea, boatSpeed)
        cMoment = 100000
        while cMoment > maxCrewStability:
            self.Boat.rollAngle += 1
            TWA, AWS                        = self.Sails.velocity_triangle(AWA, boatSpeed, self.TWS)
            kLeeway, kLift, kMoment, kDrag  = self.Keel.keel_forces(self.Boat, self.Sea, boatSpeed)
            fMoment, fLift, fDrag, fLeeway  = self.Foil.foil_forces(self.Boat, self.Sea, boatSpeed)
            sMoment, sFwd, sLeeway          = self.Sails.sail_forces(self.Boat, self.Sea, boatSpeed, self.TWS, AWA)
            hMoment                         = self.Boat.hull_stability(self.Sea, boatSpeed)
            # hDrag                           = hull_resistance(boatSpeed, Sea)
            cMoment = (sMoment + kMoment - fMoment - hMoment)
        #print('delta=',delta,'for BS=',boatSpeed, 'and AWA=',AWA,'force= ',kLeeway-sLeeway)
        return kLeeway - sLeeway + fLeeway

"""
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
            raise Exception('The drag is too high. Try reducing the lower limit of the boat speed interval')
        elif np.sign(lowLimit) == 1 and np.sign(upLimit) == 1:
            raise Exception('The drag is too low. Try raising the upper limit of the boat speed interval')

        # Compute TWA for the polar graph:
        TWA,_ = Sails.velocity_triangle(AWA, BS, TWS)
        BSlist.append(float(BS))
        TWAlist.append(float(TWA))
        leeAngleList.append(float(Boat.leewayAngle))
        rollList.append(float(Boat.rollAngle))
    return BSlist, TWAlist, leeAngleList, rollList





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
"""