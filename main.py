# This Python file uses the following
# encoding: utf-8
"""
script per il calcolo del modello di stabilità  di tutta la barca
(scafo+equipaggio+foil)
"""
import sys
# Add packages (folders with modules in them)
sys.path.insert(0, "boat-components")

# from boatClasses import Boat, Crew, Foil, Sea, Keel, Sails

from hull import Boat
from sails import Sails
from foils import Foil
from keel_and_rudder import Keel
from sea import Sea
from crew import Crew

from interactive_plot import interactive_plot
from input_data import input_data_dictionary
import matplotlib.pyplot as plt
import numpy as np
from vpp import Vpp
from scipy.misc import derivative
import time
from postprocess import postprocess

#############################################################################
# PRE - PROCESSING

# Extract data from external dictionary to run the analysis
boatsDict, foilsDict, crewsDict, seasDict, \
    sailsDict, keelsDict, wingProfiles = input_data_dictionary()

# Init hulls
Ate = Boat(boatsDict["Ate"])#inserisco tutti i dati in qualche modo

# Init sails
Olimpics_1 = Sails(sailsDict['Set 1'])

# Init crew
Crew = Crew(crewsDict["Crew 1"])

# Init foils
noFoil = Foil(foilsDict['Foil1'], wingProfiles)

# Init rudder and keel
# TODO: add rudder functionalities
Deriva2 = Keel(keelsDict['Keel1'], wingProfiles)

# Init regatta course
Garda = Sea(seasDict["Garda"])

# Init apparent wind angles
AWAvector = np.arange(20, 140, 1)

# Init true wind speed
TWS = 2.5

# Init vpp analysis
vpp1 = Vpp(Ate, Garda, Crew, Olimpics_1, Deriva2, noFoil, AWAvector, TWS)
#vpp3 = Vpp(Ate, Garda, Crew, Olimpics_1, Deriva2, Foil3, AWAvector, TWS)
#vpp2 = Vpp(Ate, Garda, Crew, Olimpics_1, Deriva2, Foil1, AWAvector, TWS)

# Group all vpp analysis in a list
vpplist = [vpp1]#, vpp2]#, vpp3]
##############################################################################
# PROCESSING

# Old interactive plot used for educational purposes
#interactive_plot(Ate, Foil1, np.arange(1,TWS,0.1), Garda, Crew)

for vpp in vpplist:
    startTime = time.time()
    vpp.vpp_static()
    executionTime = (time.time() - startTime)
    print('time passed: ',executionTime)

##############################################################################
# POST - PROCESSING

# Choose a descriptive name for each vpp analysis
legend = ['No foil']#, 'wortmann']# 'naca64012', 'naca0021']

postprocess(vpplist, legend, TWS)
