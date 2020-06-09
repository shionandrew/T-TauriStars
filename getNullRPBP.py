# Author: Shion Andrew
# Checks to see how many of stars flagged within Stripe82 appear in published catalog

import math
import statistics
import csv
import os

# Scientific libraries
import numpy as np

# Astropy and Gaia
import astroquery
import keyring
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import stats

#workingDirectory = '/mnt/home/f0011229'
workingDirectory = '/Users/touatokuchi/Desktop/MSU/HPCC/Stripe'
outputDirectory = '/Users/touatokuchi/Desktop'

dtype = [('Ra', float), ('Dec', float), ('GMagnitude', float), ('GNumObs', int), ('RPMagnitude', float), ('RPNumObs', int), ('BPMagnitude', float), ('BPNumObs', int), ('GSigmaFromBaseline', float), ('RPSigmaFromBaseline', float), ('BPSigmaFromBaseline', float), ('SourceId', int), ('AEN', float), ('BRexcess', float)]

flagged = []

count = 0
with open(workingDirectory + '/NewStripeVariablesPotentialFull.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            ra = float(row[0])
            dec = float(row[1])
            Gmagnitude = float(row[2])
            RPmagnitude = float(row[3])
            BPmagnitude = float(row[4])
            GnumObs = float(row[5])
            RPnumObs = int(row[6])
            BPnumObs = int(row[7])
            GSigmaFromBaseline = float(row[8])
            RPSigmaFromBaseline = float(row[9])
            BPSigmaFromBaseline = float(row[10])
            sourceId = int(row[11])
            falsePositive = int(row[12])
            AEN = float(row[13])
            BRexcess = float(row[14])
            AENboundary = .12 + 2.66*pow(10,-6)*math.exp(0.7*Gmagnitude)
            BRboundary = 1.39 + 2.18*pow(10,-7)*math.exp(0.76*Gmagnitude)
            count += 1
            #and GSigmaFromBaseline > 5 and RPSigmaFromBaseline > 5 and BPigmaFromBaseline > 5
            if BRexcess <  BRboundary and AEN < AENboundary:
                if falsePositive == 0:
                    if BPSigmaFromBaseline == 0 or RPSigmaFromBaseline == 0:
                        flagged.append((ra, dec, Gmagnitude, GnumObs, RPmagnitude, RPnumObs, BPmagnitude, BPnumObs, GSigmaFromBaseline, RPSigmaFromBaseline, BPSigmaFromBaseline, sourceId, AEN, BRexcess))

        except:
            continue
print(count)
csvfile.close()


flaggedStars = np.array(flagged, dtype=dtype)
print("FLAGGED STARS")
print(len(flaggedStars))


with open(outputDirectory + '/need_RP_BPinfo.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Ra', 'Dec', 'Gmagnitude', 'GNumObs', 'RPmagnitude', 'RPnumObs', 'BPmagnitude', 'BPNumObs', 'Gsigma', 'RPSigma', 'BPsigma', 'SourceID'])
    for star in range(len(flaggedStars)):
        potentialVariable = flaggedStars[star]
        writer.writerow(potentialVariable)
