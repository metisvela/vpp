# -*- coding: utf-8 -*-
"""
Script per calcolare il momento raddrizzante dato da un foil a v
Created on Fri Nov 13 17:22:48 2020

@author: giova
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

def stabilita_foil(Foil, Boat, Sea, boatSpeed):
    # stabilità foil
    # Parametri di progetto
    B = Boat.foilBeam #larghezza scafo
    Hs = Boat.foilHeight # altezza buco di uscita scassa
    Tc = Boat.canoeDraft # pescaggio
    gamma1 = np.radians(Foil.gamma1) # angolo avambraccio
    gamma2 = np.radians(Foil.gamma2) # angolo braccio
    FB = Boat.freeboard # altezza bordo libero da DWL

    # costanti idrodinamiche
    pressioneDinamica = 0.5 * Sea.waterDensity * boatSpeed**2 #kg/(m2 s) pressione dinamica a 2.5 m/s
    corda = Foil.chord #cm, corda di braccio e avambraccio
    cL = Foil.cL # coefficiente di lift

    theta = np.radians(Boat.rollAngle)

    lunghezzaMassima = B/np.cos(gamma1) + (Hs - Tc)/np.sin(gamma1)

    immersione = B * np.tan(gamma1) - FB + Hs
    if immersione < 0:
        print("il foil è fuori dall'acqua! Cambia la geometria delle scasse")
        momentoTotale = 0
    else:
        spanAvambraccio = immersione / np.sin(gamma1)
        spanBraccio = immersione / np.sin(gamma2)
    
        braccioAvambraccio = lunghezzaMassima * (np.cos(gamma1))**2 + spanAvambraccio/2
        braccioBraccio = lunghezzaMassima * np.cos(gamma1) * np.sin(gamma1) + spanBraccio/2
    
        # variazione con theta delle geometrie dell'avambraccio
        deltaSpanAvambraccio = (lunghezzaMassima * np.cos(gamma1) * np.sin(theta))/(np.sin(gamma1 + theta))
        braccioAvambraccioTheta = braccioAvambraccio - deltaSpanAvambraccio/2
        spanAvambraccioTheta = spanAvambraccio + deltaSpanAvambraccio
    
    
        # variazione con theta delle geometrie del braccio
        immersioneTheta = spanAvambraccioTheta * np.sin(gamma1 + theta)
    
        spanBraccioTheta = 0.5 # immersioneTheta / (np.sin(gamma2 - theta))
        braccioBraccioTheta = lunghezzaMassima * np.cos(gamma1) * np.sin(gamma1) * spanBraccioTheta/2
    
        # calcolo forze e momento raddrizzante
        forzaAvambraccio = pressioneDinamica * corda * spanAvambraccioTheta * cL *np.cos(gamma1) * np.cos(theta)
        momentoAvambraccio = forzaAvambraccio * braccioAvambraccioTheta
    
        forzaBraccio = pressioneDinamica * corda * spanBraccioTheta * cL *np.cos(gamma2) * np.cos(theta)
        momentoBraccio = forzaBraccio * braccioBraccioTheta
    
        momentoTotale = momentoAvambraccio + momentoBraccio

    return momentoTotale
