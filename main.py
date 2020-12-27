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

# Dizionario per lo storage dei dati dei vari scafi. Forse è più semplice salvare
# i dati in un excel e usare pandas per importarli
boatsDict, foilsDict, crewsDict, seasDict, sailsDict = input_data_dictionary()

# Inizializzo le barche
Ate = Boat(boatsDict["Ate"])#inserisco tutti i dati in qualche modo

# Inizializzo le vele
Olimpics_1 = Sails(sailsDict['Set 1'])

# Inizializzo l'equipaggio
Crew = Crew(crewsDict["Crew 1"])

# Inizializzo foil
Foil = Foil(foilsDict)

# Inizializzo deriva e timone
Deriva = Keel(0.2,1.5)
Timone = Keel(0.13, 1)

# Inizializzo mare
Garda = Sea(seasDict["Garda"])

BS, AWA = vpp(Ate, Olimpics_1, Crew, Foil, Deriva, Garda)
plt.polar(np.radians(AWA), BS, '.')