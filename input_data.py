"""
utilizzare questa pagina per inserire i dati che verranno utilizzati
dal programma.
"""

def input_data_dictionary():
    boatsDict = {"Ate": {"lengthWaterl" : 4.6,
                     "beamWaterl" : 1.064,
                     "canoeDraft" : 0.13,
                     "midshipCoeff" : 0.727,
                     "displacement" : 255*0.001, #m^3
                     "prismCoeff" : 0.607,
                     "longCentBuoy" : (0.2655+2.3),
                     "longCentFlot" : (0.4716+2.3),
                     "wettedArea" : 3.706,
                     "foilBeam" : 0.5,
                     "freeboard" : 0.3,
                     "foilHeight" : 0.5
                     },
             "Athena": {},
             "Arete": {},
             "Aura": {},
             "A2020": {}
             }
    foilsDict = {"gamma1" : 30,
                 "gamma2" : 40,
                 "chord"  : 0.2,
                 "cL"     : 0.8,
                 "cD"     : 0.1};
    crewsDict = {"Crew 1": {"bowmanWeight"   : 80,
                        "helmsmanWeight" : 60,
                        "bowmanHeight"   : 1.80,
                        "helmsmanHeight" : 1.65}};

    seasDict = {"Garda" : {"waterDensity"       : 1000,
                           "cinematicViscosity" : 1e-06,
                           "gravityConstant"    : 9.8,
                           "airDensity"         : 1}}

    sailsDict = {'Set 1' : {'EHM' : 7,           # Mast height from freeboard
                            'FA'  : 0.3,         # freeboard mean height
                            'BMAX': 1.2,         # max width of the yacht
                            'EMDC' : 0.12,       # Mean diameter of mast
                            'BAD' : 0.8,         # Boma height from freeboard
                            'H'   : 1,           # Cl correction factor
                            'F'   : 1,           # sail correction factor
                            'Am'  : 10.93,       # main sail area
                            'Aj'  : 4.93,        # jib area
                            'Ag'  : 17.95,       # gennaker area
                            'genAngle': 70,       # AWA at which gennaker is hoisted
                            'xMain': 2.3,         # centre of main sail, x coord
                            'zMain': 3,          # centre of main sail, z coord
                            'xJib' : 2,          # centre of jib, x coord
                            'zJib' : 1.2         # centre of jib, z coord 
                            }
                }

    return boatsDict, foilsDict, crewsDict, seasDict, sailsDict
