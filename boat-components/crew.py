# -*- coding: utf-8 -*-
import numpy as np

##############################################################################
# Master keel class starts below
##############################################################################

class Crew:
    def __init__(self, crewDict):
        # See the corresponding dictionary on input_data.py for information
        # on this parameters.
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
        # crewRightMom = crewRightMom * np.ones(shape=boatSpeed.size)
        return crewRightMom
