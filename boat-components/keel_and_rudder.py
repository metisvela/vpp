# -*- coding: utf-8 -*-
import numpy as np
import pickle
from utilities import dynamic_pressure, lift_coefficients_3D, interpolate_wing_coefficients

##############################################################################
# Master keel class starts below
##############################################################################

class Keel:

    def __init__(self,keelDict, wingProfiles):
        # See the corresponding dictionary on input_data.py for information
        # on this parameters.
        self.chord           = keelDict['chord']
        self.span            = keelDict['span']
        self.area            = self.chord * self.span
        self.aspectRatio     = self.span / self.chord
        self.profile         = keelDict['profile']
        #self.xCoord          = wingProfiles[self.profile]['x']
        #self.yCoord          = wingProfiles[self.profile]['y']
        #profilo = xfoil.model.Airfoil(self.xCoord, self.yCoord)
        #self.xf = xfoil.XFoil()
        #self.xf.airfoil = profilo
        #self.xf.n_crit = 3
        #self.xf.repanel()
        #self.xf.max_iter=50
        # Load the precomputed interpolation functions from folder
        # "foil-profiles"
        handler = open('foil-profiles/'+self.profile, 'rb')
        f_list = pickle.load(handler)
        self.liftFunction = f_list[0]
        self.dragFunction = f_list[1]
        #self.liftFunction, self.dragFunction = interpolate_wing_coefficients(self.xf)

        return

    def keel_forces(self, Boat, Sea, boatSpeed):
        AoA = Boat.leewayAngle * np.cos(np.radians(Boat.rollAngle)) # effective angle
        reynolds  = boatSpeed*self.chord / Sea.cinematicViscosity # np array storing all reynolds number
        aspectRatio = self.aspectRatio*2 # specular effect of the keel
        try:
            liftCoeff2D = self.liftFunction(AoA, reynolds.flatten())
            dragCoeff2D = self.dragFunction(AoA, reynolds.flatten())
        except AttributeError:
            liftCoeff2D = self.liftFunction(AoA, reynolds)
            dragCoeff2D = self.dragFunction(AoA, reynolds)
        liftCoeff3D, dragCoeff3D = lift_coefficients_3D(liftCoeff2D, dragCoeff2D, aspectRatio)
        dynPressure = dynamic_pressure(Sea, boatSpeed)

        lift = (dynPressure * 2*self.area * liftCoeff3D ) / 2
        drag = (dynPressure * 2*self.area * dragCoeff3D ) / 2

        # Now compute the forces in the boat frame of reference
        leewayLift = lift * np.cos(np.radians(Boat.rollAngle)) # keel force counteracting leeway
        upwdLift   = lift * np.sin(np.radians(Boat.rollAngle)) # vertical component of keel force

        # Calculate coordinates of the center of lateral resistance
        xCLR, zCLR = self.boat_CLR(Boat)
        keelMoment = lift * zCLR
        return leewayLift, upwdLift, keelMoment, drag

    def boat_CLR(self, Boat):
        """
        Calculation of the center of lateral resistance is based on the theory
        presented by Larsson in Principle of Yacht Design (2000) which is valid
        for high aspect ratio fin-keel yachts.
        """
        zCLR = 0.45*(Boat.canoeDraft + self.span)
        # xCLR = 0.25 of the chord, but a consistent coordinate system must be
        # first established. For now, let's say
        xCLR = 0
        return xCLR, zCLR
