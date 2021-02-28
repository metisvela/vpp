# This Python file uses the following
# encoding: utf-8
"""
script per il calcolo del modello di stabilità  di tutta la barca
(scafo+equipaggio+foil)
"""
from boatClasses import Boat, Crew, Foil, Sea, Keel, Sails
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

# Dizionario per lo storage dei dati dei vari scafi. Forse è più semplice salvare
# i dati in un excel e usare pandas per importarli
boatsDict, foilsDict, crewsDict, seasDict, \
    sailsDict, keelsDict, wingProfiles = input_data_dictionary()

# Inizializzo le barche
Ate = Boat(boatsDict["Ate"])#inserisco tutti i dati in qualche modo

# Inizializzo le vele
Olimpics_1 = Sails(sailsDict['Set 1'])

# Inizializzo l'equipaggio
Crew = Crew(crewsDict["Crew 1"])

# Inizializzo foil
#Foil1 = Foil(foilsDict['Foil1'], wingProfiles)
#Foil3 = Foil(foilsDict['Foil3'], wingProfiles)
noFoil = Foil(foilsDict['noFoil'], wingProfiles)

# Inizializzo deriva e timone
#Deriva1 = Keel(keelsDict['Keel1'], wingProfiles)
Deriva2 = Keel(keelsDict['Keel1'], wingProfiles)
#Timone = Keel(0.13, 1)

# Inizializzo mare
Garda = Sea(seasDict["Garda"])

AWAvector = np.arange(20, 140, 1)
TWS = 2.5


vpp1 = Vpp(Ate, Garda, Crew, Olimpics_1, Deriva2, noFoil, AWAvector, TWS)
#vpp3 = Vpp(Ate, Garda, Crew, Olimpics_1, Deriva2, Foil3, AWAvector, TWS)
#vpp2 = Vpp(Ate, Garda, Crew, Olimpics_1, Deriva2, Foil1, AWAvector, TWS)
##############################################################################
# PROCESSING

#interactive_plot(Ate, Foil1, np.arange(1,TWS,0.1), Garda, Crew)

vpplist = [vpp1]#, vpp2]#, vpp3]

for vpp in vpplist:
    startTime = time.time()
    vpp.vpp_static()
    executionTime = (time.time() - startTime)
    print('time passed: ',executionTime)


##############################################################################
# POST - PROCESSING
legend = ['No foil']#, 'wortmann']# 'naca64012', 'naca0021']
postprocess(vpplist, legend, TWS)
