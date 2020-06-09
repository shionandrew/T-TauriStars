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

workingDirectory = '/mnt/home/f0011229/Stripe'
outputDirectory = '/mnt/home/f0011229/Stripe'
workingDirectory2 = '/mnt/home/f0011229/Stripe/Galaxies'
flagged = []

dtype = [('Ra', float), ('Dec', float), ('GMagnitude', float), ('RPMagnitude', float), ('BPMagnitude', float), ('GNumObs', int), ('RPNumObs', int), ('RBNumObs', int), ('GSigmaFromBaseline', float), ('RPSigmaFromBaseline', float), ('BPSigmaFromBaseline', float), ('SourceId', int), ('False positive', int)]
with open(workingDirectory + '/StripeCandidateVariablesCheckWithCatalogs.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
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
        falsePositive = float(row[12])
        #if True:
        if falsePositive == 1:
        #if falsePositive == 1:
            #if BPsigma > 3 and RPsigma > 3:
                #continue
            #elif BPsigma > 3 or RPsigma > 3:
                #continue
            #else:
            if Gsigma > 5:
                flagged.append((ra, dec, Gmagnitude, RPmagnitude, BPmagnitude, GnumObs, RPnumObs, BPnumObs, Gsigma, RPsigma, BPsigma, sourceId, falsePositive))
        #elif BPsigma > 3 and RPsigma > 3:
            #flagged.append((ra, dec, Gmagnitude, RPmagnitude, BPmagnitude, GnumObs, RPnumObs, BPnumObs, Gsigma, RPsigma, BPsigma, sourceId, falsePositive))

csvfile.close()


flaggedStars = np.array(flagged, dtype=dtype)
flaggedStars = np.sort(flaggedStars, order = 'SourceId')
print("FLAGGED STARS")
print(len(flaggedStars['SourceId']))
##############################
###### REMOVING REPEATS ######
##############################
nonRepeatingIndices = []
for star in range(len(flaggedStars['SourceId'])-1):
    currentStar = flaggedStars['SourceId'][star]
    neighborStar = flaggedStars['SourceId'][star+1]
    if currentStar != neighborStar:
        nonRepeatingIndices.append(star)
nonRepeatingIndices.append(len(flaggedStars['SourceId'])-1)
repeatsCulledStars = np.take(flaggedStars, nonRepeatingIndices)
flaggedStars = repeatsCulledStars
print(len(flaggedStars['SourceId']))
##############################
###### REMOVING GALAXIES ########
##############################
knownGalaxies = []

with open(workingDirectory2 + '/Glade_MatchedWithGaia.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[3])
            knownGalaxies.append(sourceId)
        except:
            continue
csvfile.close()
'''
with open(workingDirectory2 + '/BailerGalaxy.tsv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[2])
            knownGalaxies.append(sourceId)
        except:
            continue
csvfile.close()

with open(workingDirectory2 + '/SDSSMorphologic_MatchedWithGaia.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[2])
            knownGalaxies.append(sourceId)
        except:
            continue
csvfile.close()
'''

with open(workingDirectory2 + '/SDSS_MatchedWithGaia.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[2])
            knownGalaxies.append(sourceId)
        except:
            continue

csvfile.close()
knownGalaxies.sort()
print(len(knownGalaxies))
knownStar = 0
goodIndices = []
for star in range(len(flaggedStars['SourceId'])):
    while True:
        # reached end of galaxy list, do not increment known variables
        if knownStar >= len(knownGalaxies):
            goodIndices.append(star)
            #potentialVariables = np.append(potentialVariables, flaggedStars[star])
            break

        elif int(flaggedStars['SourceId'][star]) == knownGalaxies[knownStar]:
            knownStar += 1
            break

        # flagged star id not in galaxy
        elif int(flaggedStars['SourceId'][star]) < knownGalaxies[knownStar]:
            #potentialVariables = np.append(potentialVariables, flaggedStars[star])
            goodIndices.append(star)
            break

        # flagged star id might be further down in galaxy
        elif int(flaggedStars['SourceId'][star]) > knownGalaxies[knownStar]:
            knownStar += 1

potentialVariables = np.take(flaggedStars, goodIndices)
flaggedStars = potentialVariables
flaggedStars = np.sort(flaggedStars, order = 'SourceId')
print(len(flaggedStars['SourceId']))

##############################
###### CROSS MATCHING WITH KNOWN VARIABLES ########
##############################
knownVariables = []
with open(workingDirectory + '/StripeFull.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[3])
            knownVariables.append(sourceId)
        except:
            continue
csvfile.close()

with open(workingDirectory + '/Stripe.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[2])
            knownVariables.append(sourceId)
        except:
            continue
csvfile.close()

with open(workingDirectory + '/Atlasmatched.vsb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[21])
            knownVariables.append(sourceId)
        except:
            continue
csvfile.close()



with open(workingDirectory + '/ATLAS_ALL_MATCHED.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[0])
            knownVariables.append(sourceId)
        except:
            continue
csvfile.close()

with open('/mnt/home/f0011229/KnownVariables/CatalinaData.csv') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        try:
            sourceId = int(row[14])
            knownVariables.append(sourceId)
        except:
            continue
csvfile.close()

goodIndices = []
knownVariables.sort()
print('KNOWNVARIBLES')
print(len(knownVariables))
potentialVariables = np.array([], dtype = dtype)
knownStar = 0
falsePositives = 0
for star in range(len(flaggedStars['SourceId'])):
    while True:
        # reached end of galaxy list, do not increment known variables
        if knownStar >= len(knownVariables):
            goodIndices.append(star)
            falsePositives += 1
            #potentialVariables = np.append(potentialVariables, flaggedStars[star])
            break

        elif int(flaggedStars['SourceId'][star]) == knownVariables[knownStar]:
            knownStar += 1
            break

        # flagged star id not in galaxy
        elif int(flaggedStars['SourceId'][star]) < knownVariables[knownStar]:
            #potentialVariables = np.append(potentialVariables, flaggedStars[star])
            goodIndices.append(star)
            falsePositives += 1
            break

        # flagged star id might be further down in galaxy
        elif int(flaggedStars['SourceId'][star]) > knownVariables[knownStar]:
            knownStar += 1


print("FALSE POSITIVES")
print(falsePositives)
potentialVariables = np.take(flaggedStars, goodIndices)

with open(outputDirectory + '/StripeGold.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Ra', 'Dec', 'Gmagnitude', 'GNumObs', 'RPmagnitude', 'RPnumObs', 'BPmagnitude', 'BPNumObs', 'Gsigma', 'RPSigma', 'BPsigma', 'SourceID'])
    for star in range(len(potentialVariables)):
        potentialVariable = potentialVariables[star]
        writer.writerow(potentialVariable)
