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
from vpp import vpp

import time

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
Foil = Foil(foilsDict, wingProfiles)

# Inizializzo deriva e timone
Deriva1 = Keel(keelsDict['Keel1'], wingProfiles)
Deriva2 = Keel(keelsDict['Keel1'], wingProfiles)
#Timone = Keel(0.13, 1)

# Inizializzo mare
Garda = Sea(seasDict["Garda"])


startTime = time.time()



BS, TWA, delta1, roll1 = vpp(Ate, Olimpics_1, Crew, Foil, Deriva1, Garda)

executionTime = (time.time() - startTime)
print(executionTime)

plt.polar(np.radians(TWA), BS)

BS, TWA, delta2, roll2 = vpp(Ate, Olimpics_1, Crew, Foil, Deriva2, Garda)
plt.polar(np.radians(TWA), BS)
