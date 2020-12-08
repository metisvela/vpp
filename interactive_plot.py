# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 23:02:56 2020

@author: giova
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def interactive_plot(Boat, Foil, boatSpeed, Sea, Crew):
    fig, ax           = plt.subplots(1, 2, figsize=(15,8))
    resistLine,       = ax[0].plot(boatSpeed, Boat.hull_resistance(boatSpeed, Sea))
    stabFoilLine,     = ax[1].plot(boatSpeed, Foil.foil_stability(Boat, Sea, boatSpeed))
    hullRightMomLine, = ax[1].plot(boatSpeed, Boat.hull_stability(Sea, boatSpeed))
    crewRightMomLine, = ax[1].plot(boatSpeed, Crew.crew_stability(Boat, Sea, boatSpeed))
    totRightMomLine,  = ax[1].plot(boatSpeed, Foil.foil_stability(Boat, Sea, boatSpeed) + Boat.hull_stability(Sea, boatSpeed) + Crew.crew_stability(Boat, Sea, boatSpeed))
    
    


    sliderDict = {"waterline beam" : {"start" : 0.9,
                                      "end"   : 1.4,
                                      "valinit" : Boat.beamWaterl,
                                      "associated variable" : "Boat.beamWaterl"},
                  "roll angle" : {"start"   : 0,
                                  "end"     : 20,
                                  "valinit" : 0,
                                  "associated variable" : "Boat.rollAngle"},
                  "prism coeff" : {"start"  : 0.5,
                                   "end"    : 0.61,
                                   "valinit": Boat.prismCoeff,
                                   "associated variable" : "Boat.prismCoeff"},
                  "gravity"    :  {"start"  : 5,
                                   "end"    : 13,
                                   "valinit": Sea.gravityConstant,
                                   "associated variable" : "Sea.gravityConstant"},
                  "freeboard"  :  {"start"  : 0.1,
                                   "end"    : 1,
                                   "valinit": Boat.freeboard,
                                   "associated variable" : "Boat.freeboard"}}
    index = -1
    bottom = 0.1
    sliderList = []
    for slidername in sliderDict:
        index +=1
        left = 0.09
        width = 0.12
        height = 0.02
        bottom += height * 2
        start = sliderDict[slidername]["start"]
        end   = sliderDict[slidername]["end"]
        valinit = sliderDict[slidername]["valinit"]
        axes = plt.axes([left, bottom, width, height])
        sliderList.append(Slider(axes, slidername, start, end, valinit))


    fig.suptitle("Prestazioni in acqua liscia")
    ax[0].set_title("Resistenza")
    ax[0].set_xlabel("Velocità (m/s)")
    ax[0].set_ylabel("Resistenza (N)")
    ax[1].set_title("Momento raddrizzante")
    ax[1].set_xlabel("Velocità (m/s)")
    ax[1].set_ylabel("Momento (Nm)")
    plt.subplots_adjust(left = 0.3)
    plt.show()


    def update(val):

        for (variable, slider) in zip(sliderDict, sliderList):
            string = sliderDict[variable]["associated variable"]
            exec(string + "= slider.val")
        resistLine.set_ydata(Boat.hull_resistance(boatSpeed, Sea,))
        stabFoilLine.set_ydata(Foil.foil_stability(Boat, Sea, boatSpeed))
        hullRightMomLine.set_ydata(Boat.hull_stability(Sea, boatSpeed))
        crewRightMomLine.set_ydata(Crew.crew_stability(Boat, Sea, boatSpeed))
        totRightMomLine.set_ydata(Foil.foil_stability(Boat, Sea, boatSpeed) + Boat.hull_stability(Sea, boatSpeed) + Crew.crew_stability(Boat, Sea, boatSpeed))
        
        
        fig.canvas.blit(ax[0].bbox)
        fig.canvas.blit(ax[1].bbox)


    for slider in sliderList:
        slider.on_changed(update)
