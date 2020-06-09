#
# Author: Shion Andrew
#
# This program retrieves baseline curve of stars as a function of photometric observation number.
#(For a given set of stars within a given range of photometric observation numbers,
# there exists an already calculated data set is divided into magnitude bins containing the mean magnitude error and standard deviation)
# Plots surroundings of known variable stars to see if they have been flagged through this method

import math
import statistics
import csv
import os

# Scientific libraries
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

# Astropy and Gaia
from astroquery.gaia import Gaia
import astroquery
import keyring
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import stats
import sys

workingDirectory = '/mnt/home/f0011229'
outputDirectory =  '/mnt/home/f0011229/Stripe'
filename = '/Stripe/StripeNeedsNeighbors.csv'
lowerLimit = int(sys.argv[1])
amount = 500
upperLimit = lowerLimit*amount + amount
RANGE = range(lowerLimit*amount, upperLimit)

def KnownVariables():
    dtype = [('Ra', float), ('Dec', float), ('GMagnitude', float), ('RPMagnitude', float), ('BPMagnitude', float), ('GNumObs', int), ('RPNumObs', int), ('BPNumObs', int), ('GSigmaFromBaseline', float), ('RPSigmaFromBaseline', float), ('BPSigmaFromBaseline', float), ('SourceId', int), ('False positive', int)]
    flagged = []
    count = -1
    with open(workingDirectory + filename) as csvfile:

        ra = float(row[0])
        dec = float(row[1])
        Gmagnitude = float(row[2])
        RPmagnitude = float(row[3])
        BPmagnitude = float(row[4])
        GnumObs = int(row[5])
        RPnumObs = int(row[6])
        BPnumObs = int(row[7])
        Gsigma = float(row[8])
        RPsigma = float(row[9])
        BPsigma = float(row[10])
        sourceId = int(row[11])

        falsePositive = checkNeighbor(ra, dec, sourceId)
        if falsePositive == 1:
            flagged.append((ra, dec, Gmagnitude, RPmagnitude, BPmagnitude, GnumObs, RPnumObs, BPnumObs, GSigmaFromBaseline, RPSigmaFromBaseline, BPSigmaFromBaseline, sourceId, falsePositive))

    csvfile.close()

    ## create array sorted by amplitude
    potentialVariables = np.array(flagged, dtype=dtype)
    print(len(potentialVariables['SourceId']))
    with open(outputDirectory + '/CatalinaFalsePositives.csv', 'a') as csvfile:
        writer = csv.writer(csvfile)
        for star in range(len(potentialVariables)):
            potentialVariable = potentialVariables[star]
            writer.writerow(potentialVariable)


def checkNeighbor(variableRa, variableDec, sourceId):
    coord = SkyCoord(ra=variableRa, dec=variableDec, unit=(u.degree, u.degree), frame='icrs')
    radius = u.Quantity(0.00138889, u.deg)
    j = Gaia.cone_search_async(coord, radius)
    r = j.get_results()
    if len(r['ra']) > 1:
        for neighbor in range(len(r['ra'])):
            #make sure neighbor is different from target star
            if r['source_id'][neighbor] != sourceId:
                # flag star brighter than 8 magnitude within 50"
                if r['phot_g_mean_mag'][neighbor] < 8:
                    return 1
                # flag star brighter than 12 magnitude within 30"
                elif r['dist'][neighbor] < 0.00833333 and r['phot_g_mean_mag'][neighbor] < 12:
                    return 1
                # flag star brighter than 20 magnitude within 10"
                elif r['dist'][neighbor] < 0.00277778 and r['phot_g_mean_mag'][neighbor] < 20:
                    return 1
                 # flag star within 5"
                elif r['dist'][neighbor] < 0.00138889:
                    return 1

    return 0


def main():
    KnownVariables()

if __name__== "__main__":
    main()
