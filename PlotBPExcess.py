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
workingDirectory = '/Users/touatokuchi/Desktop'

dtype = [('Ra', float), ('Dec', float), ('GMagnitude', float), ('GMagError', float), ('RPMagnitude', float),('RPMagError', float),  ('BPMagnitude', float),('BPMagError', float),  ('GNumObs', int), ('RPNumObs', int), ('RBNumObs', int),('SourceId', int), ('AEN', float), ('BRexcess', float), ('w', float)]

flagged = []
rand = []
count = 0


with open(workingDirectory + '/RandomStars.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            astrometric_chi2_al = float(row[0])
            astrometric_n_good_obs_al = float(row[1])
            stat = pow(astrometric_chi2_al/(astrometric_n_good_obs_al-5),.5)
            if stat > 3:
                count += 1
            try:
                AEN = float(row[2])
            except:
                AEN = 0
            try:
                BRexcess = float(row[3])
            except:
                BRexcess = 0

            try:
                GnumObs = float(row[4])
                Gmagnitude = float(row[5])
                BRexcessboundary = 1.39 + 2.18*pow(10,-7)*math.exp(0.76*Gmagnitude)
                AENboundary = .12 + 2.66*pow(10,-6)*math.exp(0.7*Gmagnitude)
                #if AEN < AENboundary and BRexcess < BRexcessboundary:
                #    count+=1
                #    print(GnumObs)
                rand.append((0,0,Gmagnitude,0,0,0,0,0,0,0,0,0, AEN, BRexcess, stat))
            except:
                continue

        except:
            print('error')
csvfile.close


with open(workingDirectory + '/SDSS_AENBRexcess.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        if True:
            ra = float(row[0])
            dec = float(row[1])
            Gmagnitude = float(row[2])
            GMagError = float(row[3])
            RPmagnitude = float(row[4])
            RPMagError = float(row[5])
            BPmagnitude =float(row[6])
            BPMagError = float(row[7])

            GnumObs = int(row[8])
            RPnumObs = int(row[9])
            BPnumObs = int(row[10])
            sourceId = int(row[11])
            AEN = float(row[12])
            BRexcess = float(row[13])
            astrometric_chi2_al = float(row[14])
            astrometric_n_good_obs_al = float(row[15])
            stat = pow(astrometric_chi2_al/(astrometric_n_good_obs_al-5),.5)
            #if stat < 3:
            #count += 1
            BRexcessboundary = 1.39 + 2.18*pow(10,-7)*math.exp(0.76*Gmagnitude)
            AENboundary = .12 + 2.66*pow(10,-6)*math.exp(0.7*Gmagnitude)

            '''if AEN < AENboundary and BRexcess < BRexcessboundary:
                count+=1
            '''
            flagged.append((ra, dec, Gmagnitude,GMagError, RPmagnitude, RPMagError, BPmagnitude, BPMagError, GnumObs, RPnumObs, BPnumObs, sourceId, AEN, BRexcess, stat))


        else:
            print('error')
csvfile.close()


flaggedStars = np.array(flagged, dtype=dtype)
rand = np.array(rand, dtype=dtype)

print("UNDER CUT")
print(count)
print("ALL")
print(len(flaggedStars))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.ylabel(r'\textbf{count}')
plt.xlabel(r'\textbf{w}')
ax = plt.gca()
fig = plt.gcf()
ax.tick_params(labeltop=True, labelright=True)
 # fixed bin size
bins = np.arange(0,10, 0.05) # fixed bin size
plt.hist(flaggedStars['w'], bins=bins, alpha=0.5)
plt.hist(rand['w'], bins=bins, alpha=0.5)

plt.show()

'''
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.ylabel(r'\textbf{BP RP excess fator}')
plt.xlabel(r'\textbf{G magnitude}')
plt.title(r'\textbf{Random Stars}')
ax = plt.gca()
fig = plt.gcf()


G = np.arange(12, 20, 0.01)
BRboundary = []
AENboundary = []
for Gmagnitude in G:
    BRboundary.append(1.39 + 2.18*pow(10,-7)*math.exp(0.76*Gmagnitude))
    AENboundary.append(.12 + 2.66*pow(10,-6)*math.exp(0.7*Gmagnitude))
plt.scatter(G, BRboundary, s = 1, color = 'black')
plt.scatter(flaggedStars['GMagnitude'], flaggedStars['BRexcess'], s = 2, color = 'red')
plt.show()

plt.ylabel(r'\textbf{AEN}')
plt.xlabel(r'\textbf{G magnitude}')
plt.title(r'\textbf{Random Stars}')
plt.scatter(G, AENboundary, s = 1, color = 'purple')
plt.scatter(flaggedStars['GMagnitude'], flaggedStars['AEN'], s = 2, color = 'red')
plt.show()'''
