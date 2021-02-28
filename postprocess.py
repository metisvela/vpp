# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 16:38:34 2021

@author: giova
"""
import numpy as np
import matplotlib.pyplot as plt


def postprocess(vpplist, legend, TWS):
    """
    Simple postprocessing of vpp results
    """
    fig = plt.Figure()
    polAxes = plt.subplot(221, projection='polar')
    linAxes = plt.subplot(222)
    ax = plt.subplot2grid((2,2), (1,0), rowspan=1, colspan=2)


    for vpp in vpplist:

        polAxes.plot(np.radians(vpp.TWAlist), vpp.BSlist)
        linAxes.plot(vpp.TWAlist, vpp.leeAngleList)
        ax.plot(vpp.TWAlist, vpp.rollList)

    polAxes.set_title('Polare')
    linAxes.set_title('Angolo di scarroccio')
    ax.set_title('Angolo di rollio')


    polAxes.legend(legend)
    linAxes.legend(legend)
    ax.legend(legend)
    title = 'Wind Speed = ' + str(TWS)
    fig.suptitle(title)
