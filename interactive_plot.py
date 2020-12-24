# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 23:02:56 2020

@author: giova
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def interactive_plot(Boat, Foil, boatSpeed, Sea, Crew):
    hullResist        = Boat.hull_resistance(boatSpeed, Sea)
    hullStab          = Boat.hull_stability(Sea, boatSpeed)
    crewStab          = Crew.crew_stability(Boat, Sea, boatSpeed)
    foilStab, foilLift, foilDrag, foilLeew = Foil.foil_forces(Boat, Sea, boatSpeed)
    fig, ax           = plt.subplots(1, 2, figsize=(15,8))
    resistLine,       = ax[0].plot(boatSpeed, hullResist)
    stabFoilLine,     = ax[1].plot(boatSpeed, foilStab)
    hullRightMomLine, = ax[1].plot(boatSpeed, hullStab)
    crewRightMomLine, = ax[1].plot(boatSpeed, crewStab)
    totRightMomLine,  = ax[1].plot(boatSpeed, hullStab+crewStab+foilStab.flatten())
    foilResistLine,   = ax[0].plot(boatSpeed, foilDrag)
    ax_0_legend = ['Scafo','Foil']
    ax_1_legend = ['Foil','Scafo', 'Equipaggio','Totale']



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
    ax[0].legend(ax_0_legend)
    ax[1].legend(ax_1_legend)
    plt.show()


    def update(val):

        for (variable, slider) in zip(sliderDict, sliderList):
            string = sliderDict[variable]["associated variable"]
            exec(string + "= slider.val")
        hullResist        = Boat.hull_resistance(boatSpeed, Sea)
        hullStab          = Boat.hull_stability(Sea, boatSpeed)
        crewStab          = Crew.crew_stability(Boat, Sea, boatSpeed)
        foilStab, foilLift, foilDrag, foilLeew = Foil.foil_forces(Boat, Sea, boatSpeed)
        resistLine.set_ydata(hullResist)
        stabFoilLine.set_ydata(foilStab.flatten())
        hullRightMomLine.set_ydata(hullStab)
        crewRightMomLine.set_ydata(crewStab)
        totRightMomLine.set_ydata(foilStab.flatten() + hullStab + crewStab)
        foilResistLine.set_ydata(foilDrag.flatten())

        fig.canvas.blit(ax[0].bbox)
        fig.canvas.blit(ax[1].bbox)


    for slider in sliderList:
        slider.on_changed(update)
