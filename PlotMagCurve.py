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

# Astropy and Gaia
import astroquery
import keyring
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import stats

#workingDirectory = '/mnt/home/f0011229'
workingDirectory = '/Users/touatokuchi/Desktop/'


flagged = []
rand = []
count = 0

filename = raw_input("Enter filename : ")
print(workingDirectory + filename)

HEADERS = next(csv.reader(open(workingDirectory + filename)))
print (HEADERS)

SourceIdIndex = HEADERS.index('source_id')
fluxIndex =  HEADERS.index('phot_g_mean_flux')
fluxErrorIndex =  HEADERS.index('phot_g_mean_flux_error')
GmagIndex =  HEADERS.index('phot_g_mean_mag')
SourceIdIndex =  HEADERS.index('source_id')
targetSourceId = int(input("Enter target Source Id : "))

print(targetSourceId)

targetError = []
targetMag = []
magErrors = []
Gmags = []
with open(workingDirectory + filename) as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        if row[SourceIdIndex] != 'source_id':
            flux = float(row[fluxIndex])
            fluxError = float(row[fluxErrorIndex])
            magError = 1.09*fluxError/flux
            Gmag = float(row[GmagIndex])
            sourceId = int(row[SourceIdIndex])
            if sourceId == targetSourceId:
                targetError.append(magError)
                targetMag.append(Gmag)
            else:
                magErrors.append(magError)
                Gmags.append(Gmag)

csvfile.close



plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.ylabel(r'\textbf{G magnitude}')
plt.xlabel(r'\textbf{G MagError}')
ax = plt.gca()
fig = plt.gcf()
plt.gca().invert_yaxis()
ax.tick_params(labeltop=True, labelright=True)
 # fixed bin size
plt.scatter(magErrors, Gmags, color='black', s = 1)
plt.scatter(targetError, targetMag, color='red', s = 15)
plt.show()
