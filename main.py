# This Python file uses the following
# encoding: utf-8
"""
script per il calcolo del modello di stabilità  di tutta la barca
(scafo+equipaggio+foil)
"""
from boatClasses import Boat, Crew, Foil, Sea, Keel
from interactive_plot import interactive_plot
from input_data import input_data_dictionary
import numpy as np

# Dizionario per lo storage dei dati dei vari scafi. Forse è più semplice salvare
# i dati in un excel e usare pandas per importarli
boatsDict, foilsDict, crewsDict, seasDict = input_data_dictionary()

# Inizializzo le barche
Ate = Boat(boatsDict["Ate"])#inserisco tutti i dati in qualche modo

# Inizializzo l'equipaggio
Crew = Crew(crewsDict["Crew 1"])

# Inizializzo foil
Foil = Foil(foilsDict)

# Inizializzo deriva e timone
Deriva = Keel(0.2,1.5)

# Inizializzo mare
Garda = Sea(seasDict["Garda"])

boatSpeed = np.linspace(1,6,20) # metri al secondo
rollAngle = 0 # in gradi

interactive_plot(Ate, Foil, boatSpeed, Garda, Crew)
