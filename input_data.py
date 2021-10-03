# -*- coding: utf-8 -*-
"""
Use this page to store all data that can be used by the vpp engine.
"""
import numpy as np

def input_data_dictionary():
    # Dictionary storing yacht information
    boatsDict = {"Ate"    :    {"lengthWaterl" : 4.6,          # length of the waterline [m]
                                "beamWaterl" : 1.064,          # maximum beam of waterline [m]
                                "canoeDraft" : 0.13,           # maximum draft of canoe body [m]
                                "midshipCoeff" : 0.727,        # coefficient of midship [adim]
                                "displacement" : 255*0.001,    # displacement [m^3]
                                "prismCoeff" : 0.607,          # prismatic coefficient [adim]
                                "longCentBuoy" : (0.2655+2.3), # longitudinal centre of buoyancy,
                                                               # 0 at the forward perpendicular [m]
                                "longCentFlot" : (0.4716+2.3), # longitudinal centre of flotation,
                                                               # 0 at the forward perpendicular [m]
                                "wettedArea" : 3.706,          # wetted area of canoe body [m^2]
                                "foilBeam" : 0.5,              # distance of the exit point of the foils,
                                                               # from midship [m]
                                "freeboard" : 0.3,             # freeboard height [m]
                                "foilHeight" : 0.5             # height of the exit point of the foil from LWL [m]
                               },
             "Athena"     :    {},
             "Arete"      :    {},
             "Aura"       :    {},
             "A2021"      :    {"lengthWaterl" : 4.6,          # length of the waterline [m]
                                 "beamWaterl" : 1.0299,          # maximum beam of waterline [m]
                                 "canoeDraft" : 0.1379,           # maximum draft of canoe body [m]
                                 "midshipCoeff" : 0.725,        # coefficient of midship [adim]
                                 "displacement" : 246*0.001,    # displacement [m^3]
                                 "prismCoeff" : 0.582,          # prismatic coefficient [adim]
                                 "longCentBuoy" : (0.0665+2.3), # longitudinal centre of buoyancy,
                                                                # 0 at the forward perpendicular [m]
                                 "longCentFlot" : (0.4126+2.3), # longitudinal centre of flotation,
                                                                # 0 at the forward perpendicular [m]
                                 "wettedArea" : 3.607,          # wetted area of canoe body [m^2]
                                 "foilBeam" : 0.5,              # distance of the exit point of the foils,
                                                                # from midship [m]
                                 "freeboard" : 0.3,             # freeboard height [m]
                                 "foilHeight" : 0.5             # height of the exit point of the foil from LWL [m]
                                }
             }


    # Dictionary storing foil configurations
    foilsDict = {'Foil1' :{"gamma1" : 30,
                             "gamma2" : 40,
                             "chord"  : 0.2,
                             "keying" : 2,
                             "profile" : 'wortmann'
                             },
                 'noFoil' : {'gamma1' : 40,
                             'gamma2' : 30,
                             'chord'  : 0.0001,
                             'keying' : 0,
                             'profile': 'wortmann'
                             },
                 'Foil3' :{"gamma1" : 30,
                             "gamma2" : 40,
                             "chord"  : 0.2,
                             "keying" : 2,
                             "profile" : 'naca0021'
                             }
                 }

    # Dictionary storing possible crews information
    crewsDict = {"Crew 1"  : {"bowmanWeight"   : 80,     # [kg]
                              "helmsmanWeight" : 50,     # [kg]
                              "bowmanHeight"   : 1.80,   # [m]
                              "helmsmanHeight" : 1.65}}  # [m]

    # Dictionary storing information about the regatta course
    seasDict = {"Garda" : {"waterDensity"       : 1000,  # [kg/m^3]
                           "cinematicViscosity" : 1e-06, # [m^2/s]
                           "gravityConstant"    : 9.8,   # [m/s^2]
                           "airDensity"         : 1      # [kg/m^3]
                           }
                }

    # Dictionary storing information about sail sets. The parameters are derived
    # from the ORC sail model, which uses a corrected Hazen model.
    sailsDict = {'Set 1' : {'EHM' : 7,           # Mast height from freeboard [m]
                            'FA'  : 0.3,         # freeboard mean height [m]
                            'BMAX': 1.2,         # max width of the yacht [m]
                            'EMDC' : 0.12,       # Mean diameter of mast [m]
                            'BAD' : 0.8,         # Boma height from freeboard [m]
                            'H'   : 1,           # Cl correction factor [adim]
                            'F'   : 1,           # sail correction factor [adim]
                            'Am'  : 10.93,       # main sail area [m^2]
                            'Aj'  : 4.93,        # jib area [m^2]
                            'Ag'  : 17.95,       # gennaker area [m^2]
                            'genAngle': 65,      # AWA at which gennaker is hoisted [°]
                            'xMain': 2.3,        # centre of main sail, x coord [m]
                            'zMain': 3,          # centre of main sail, z coord [m]
                            'xJib' : 2,          # centre of jib, x coord [m]
                            'zJib' : 1.2         # centre of jib, z coord [m]
                            },
                 'SetSvedese': {'EHM' : 7,           # Mast height from freeboard [m]
                                'FA'  : 0.3,         # freeboard mean height [m]
                                'BMAX': 1.2,         # max width of the yacht [m]
                                'EMDC' : 0.12,       # Mean diameter of mast [m]
                                'BAD' : 0.8,         # Boma height from freeboard [m]
                                'H'   : 1,           # Cl correction factor [adim]
                                'F'   : 1,           # sail correction factor [adim]
                                'Am'  : 15.01,       # main sail area [m^2]
                                'Aj'  : 5.7,        # jib area [m^2]
                                'Ag'  : 12,       # gennaker area [m^2]
                                'genAngle': 65,      # AWA at which gennaker is hoisted [°]
                                'xMain': 2.3,        # centre of main sail, x coord [m]
                                'zMain': 3,          # centre of main sail, z coord [m]
                                'xJib' : 2,          # centre of jib, x coord [m]
                                'zJib' : 1.2         # centre of jib, z coord [m]
                     }
                }

    # Dictionary storing information about the centreboards. Be aware that the
    # model provided in this code always assumes that the distribution of lift is
    # elliptic, and so a "quarter-chord-aligned" type of centreboard or rudder
    # is going to be designed.
    keelsDict = {'Keel1' : {'profile' : 'naca0021',              # naca profile
                            'chord'   : 0.2,                   # mean chord [m]
                            'span'    : 1.5                   # span [m]
                            },
                 'Keel2' : {'profile' : 'naca0021',
                            'chord'   : 0.2,
                            'span'    : 1.3
                            }
                 }



    # Dictionary storing the x and y coordinates of possible profiles to be used
    # in the code. Note that the frame of reference starts from x=1 at the
    # leading edge and runs over the upper part of the profile, for it to then come
    # back once passed the trailing edge through the lower part.
    # Tip: for easier implementation of profile coordinates, start from the
    # Selig file format which already follows this convention.
    wingProfiles = {'wortmann' :

                                        {'x' : np.array([1.     , 0.99893, 0.99572, 0.99039, 0.98296, 0.97347, 0.96194,
                                                            0.94844, 0.93301, 0.91573, 0.89668, 0.87592, 0.85355, 0.82967,
                                                            0.80438, 0.77779, 0.75   , 0.72114, 0.69134, 0.66072, 0.62941,
                                                            0.59755, 0.56526, 0.5327 , 0.5    , 0.4673 , 0.43474, 0.40245,
                                                            0.37059, 0.33928, 0.30866, 0.27886, 0.25   , 0.22221, 0.19562,
                                                            0.17033, 0.14645, 0.12408, 0.10332, 0.08427, 0.06699, 0.05156,
                                                            0.03806, 0.02653, 0.01704, 0.00961, 0.00428, 0.00107, 0.     ,
                                                            0.00107, 0.00428, 0.00961, 0.01704, 0.02653, 0.03806, 0.05156,
                                                            0.06699, 0.08427, 0.10332, 0.12408, 0.14645, 0.17033, 0.19562,
                                                            0.22221, 0.25   , 0.27886, 0.30866, 0.33928, 0.37059, 0.40245,
                                                            0.43474, 0.4673 , 0.5    , 0.5327 , 0.56526, 0.59755, 0.62941,
                                                            0.66072, 0.69134, 0.72114, 0.75   , 0.77779, 0.80438, 0.82967,
                                                            0.85355, 0.87592, 0.89668, 0.91573, 0.93301, 0.94844, 0.96194,
                                                            0.97347, 0.98296, 0.99039, 0.99572, 0.99893, 1.])
                                                                   ,
                                        'y': np.array([0.     ,  0.00023,  0.00086,  0.00193,  0.00341,  0.00534,
                                                         0.00766,  0.01035,  0.01342,  0.01681,  0.02053,  0.02447,
                                                         0.02864,  0.03298,  0.03747,  0.042  ,  0.04652,  0.05089,
                                                         0.05511,  0.05905,  0.06275,  0.06608,  0.06911,  0.07174,
                                                         0.07409,  0.07596,  0.0775 ,  0.07845,  0.07898,  0.07888,
                                                         0.07838,  0.0772 ,  0.07565,  0.07339,  0.07081,  0.06754,
                                                         0.06404,  0.05989,  0.05569,  0.05086,  0.04609,  0.04056,
                                                         0.03523,  0.02948,  0.02368,  0.01775,  0.01144,  0.00574,
                                                         0.     , -0.002  , -0.00436, -0.00691, -0.0097 , -0.01247,
                                                        -0.01481, -0.0173 , -0.01913, -0.0211 , -0.02246, -0.02377,
                                                        -0.02447, -0.02503, -0.025  , -0.02475, -0.02389, -0.02275,
                                                        -0.021  , -0.01895, -0.01637, -0.01357, -0.01045, -0.00731,
                                                        -0.00405, -0.00092,  0.00217,  0.00496,  0.00748,  0.00951,
                                                         0.01118,  0.01241,  0.01329,  0.01373,  0.01381,  0.01347,
                                                         0.0128 ,  0.01186,  0.01072,  0.00941,  0.00802,  0.00659,
                                                         0.00515,  0.00377,  0.00251,  0.00151,  0.0007 ,  0.00015,
                                                         0.])
                                                        },


                    'naca0021' :

                                        {'x' : np.array([1.0,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.25,0.2,0.15,
                                                                  0.1,0.075,0.05,0.025,0.0125,0.0,0.0125,0.025,
                                                                  0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,
                                                                  0.7,0.8,0.9,0.95,1.0]),

                                        'y': np.array([0.00221,0.01412,0.02534,0.04591,0.06412,0.07986,0.09265,
                                                       0.10156,0.10504,0.10397,0.1004,0.09354,0.08195,0.0735,
                                                       0.06221,0.04576,0.03315,0.0,-0.03315,-0.04576,-0.06221,-0.0735,
                                                       -0.08195,-0.09354,-0.1004,-0.10397,-0.10504,-0.10156,-0.09265,
                                                       -0.07986,-0.06412,-0.04591,-0.02534,-0.01412])

                                        },

                    'naca64012' :

                                       {'x' : np.array([1.    , 0.95  , 0.9   , 0.85  , 0.8   , 0.75  , 0.7   , 0.65  ,
                                                        0.6   , 0.55  , 0.5   , 0.45  , 0.4   , 0.35  , 0.3   , 0.25  ,
                                                        0.2   , 0.15  , 0.1   , 0.075 , 0.05  , 0.025 , 0.0125, 0.0075,
                                                        0.005 , 0.    , 0.005 , 0.0075, 0.0125, 0.025 , 0.05  , 0.075 ,
                                                        0.1   , 0.15  , 0.2   , 0.25  , 0.3   , 0.35  , 0.4   , 0.45  ,
                                                        0.5   , 0.55  , 0.6   , 0.65  , 0.7   , 0.75  , 0.8   , 0.85  ,
                                                        0.9   , 0.95  , 1.    ]),
                                        'y' : np.array([1.    , 0.95  , 0.9   , 0.85  , 0.8   , 0.75  , 0.7   , 0.65  ,
                                                       0.6   , 0.55  , 0.5   , 0.45  , 0.4   , 0.35  , 0.3   , 0.25  ,
                                                       0.2   , 0.15  , 0.1   , 0.075 , 0.05  , 0.025 , 0.0125, 0.0075,
                                                       0.005 , 0.    , 0.005 , 0.0075, 0.0125, 0.025 , 0.05  , 0.075 ,
                                                       0.1   , 0.15  , 0.2   , 0.25  , 0.3   , 0.35  , 0.4   , 0.45  ,
                                                       0.5   , 0.55  , 0.6   , 0.65  , 0.7   , 0.75  , 0.8   , 0.85  ,
                                                       0.9   , 0.95  , 1.    ])
                                        }
                    }

    return boatsDict, foilsDict, crewsDict, seasDict, sailsDict, keelsDict, wingProfiles
