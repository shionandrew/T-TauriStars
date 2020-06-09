# Author: Shion Andrew
# Checks to see how many of stars flagged within Stripe82 appear in published catalog

import math
import statistics
import csv
import os

# Scientific libraries
import numpy as np

# MatPlotlib
from matplotlib import pylab
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib import rc
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

# Astropy and Gaia
import astroquery
import keyring
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import stats

workingDirectory = '/Users/touatokuchi/Desktop/'
filename = 'potentialVariablesGold.csv'
filename2 = 'potentialVariablesSilverRP.csv'
filename3 = 'potentialVariablesSilverBP.csv'
filename4 = 'potentialVariablesBronze.csv'

SECONDFILE = True
THIRDFILE = True
FOURTHFILE = True


Gmags = []
Gmags1 = []
print(filename)
with open(workingDirectory + filename) as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            Gmagnitude = float(row[3])
            Gmags.append(Gmagnitude)
            Gmags1.append(Gmagnitude)

        except:
            print('error')
csvfile.close()

Gmags2 = []
if SECONDFILE == True:
    print(filename2)
    with open(workingDirectory + filename2) as csvfile:
        plots = csv.reader(csvfile, delimiter = ',')
        for row in plots:
            try:
                Gmagnitude = float(row[3])
                Gmags.append(Gmagnitude)
                Gmags2.append(Gmagnitude)

            except:
                print('error')
    csvfile.close()

Gmags3 = []
if THIRDFILE == True:
    print(filename3)
    with open(workingDirectory + filename3) as csvfile:
        plots = csv.reader(csvfile, delimiter = ',')
        for row in plots:
            try:
                Gmagnitude = float(row[3])
                Gmags.append(Gmagnitude)
                Gmags3.append(Gmagnitude)

            except:
                print('error')
    csvfile.close()

Gmags4 = []
if FOURTHFILE == True:
    print(filename4)
    with open(workingDirectory + filename4) as csvfile:
        plots = csv.reader(csvfile, delimiter = ',')
        for row in plots:
            try:
                Gmagnitude = float(row[3])
                Gmags.append(Gmagnitude)
                Gmags4.append(Gmagnitude)

            except:
                print('error')
    csvfile.close()
print(len(Gmags))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.ylabel(r'\textbf{count}')
plt.xlabel(r'\textbf{G magnitude}')
ax = plt.gca()
fig = plt.gcf()
ax.tick_params(labeltop=True, labelright=True)
 # fixed bin size
bins = np.arange(14,19.5, 0.1) # fixed bin size
GmagsSilver = Gmags2 + Gmags3

'''
plt.hist(Gmags, bins=bins, alpha=1.0, color = 'black', histtype='step')
plt.hist(Gmags1, bins=bins, alpha=1.0, color = 'gold', histtype='step')
plt.hist(GmagsSilver, bins=bins, alpha=1.0, color = 'slategrey', histtype='step')
plt.hist(Gmags4, bins=bins, alpha=1.0, color = 'brown', histtype='step')
labels= ["Gold","Silver", "Bronze"]
handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ['gold','slategrey', 'brown']]
plt.legend(handles, labels, loc='upper left')
plt.show()
'''

plt.hist(Gmags2, bins=bins, alpha=1.0, color = 'darkred', histtype='step')
plt.hist(Gmags3, bins=bins, alpha=1.0, color = 'darkblue', histtype='step')
labels= ["RP $\sigma >$ 3","BP $ \sigma >$ 3"]
handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ['darkred','darkblue']]
plt.legend(handles, labels, loc='upper left')
plt.show()
