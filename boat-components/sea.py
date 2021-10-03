# -*- coding: utf-8 -*-

##############################################################################
# Master keel class starts below
##############################################################################

class Sea:
    def __init__(self, seaDict):
        # See the corresponding dictionary on input_data.py for information
        # on this parameters.
        self.waterDensity = seaDict["waterDensity"]
        self.airDensity   = seaDict["airDensity"]
        self.cinematicViscosity = seaDict["cinematicViscosity"]
        self.gravityConstant = seaDict["gravityConstant"]
        return
