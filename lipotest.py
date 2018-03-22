# LipoAnalysis 1.0
# This code will take the image analysis results from the Liposome Toxicity
# assay, normalize the increase in intensity in the sample image over the blank
# image utilizing the positive control (ionomycin) image,
# and organize those results into a readily plottable format,
# as well as produce a simple boxplot of the results using R, rpy2 and ggplot2.
# Troubleshooting: If using rpy2 is giving you a lot of trouble, comment out
# the portion of this code beyond 'R code starts here' and try running the
# program again.
#
# For a detailed scientific explication of the Liposome toxicity assay, please see:
# Flagmeier, P; De, S et al. Angew Chem Int. Ed. 2017. doi: 10.1002/anie.201700966
#
# Written for Python 3.6 with scipy.
# Author: Santiago Enrique Sanchez
# Author Email: ses94@cam.ac.uk
# Updated: March 22nd, 2018


import numpy as np
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import grid
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
import pandas as pd
import sklearn
from sklearn.utils.extmath import cartesian


# Initialize R instance
rprint = robjects.globalenv.get("print")
stats = importr('stats')
gridExtra = importr('gridExtra')
outliers = importr('outliers')
grdevices = importr('grDevices')
tidyr = importr('tidyr')
base = importr('base')
datasets = importr('datasets')
reshape = importr('reshape2')
#

# Input
expName = "aSyn"
expDate = "12-1-18"
workDir = 'G:\\Liposome Data ST 12 1 28\\'
sampleConcs = ["1 uM", "500 nM", "100 nM"]
timePoints = ["12 hour", "24 hour", "Day 3", "Day 5", "Day 7", "Day 11", "Day 14"]
nFOV = 9
# End Input

# These objects will be used to index the FoVs and
# individual results files read into the program
dfKeys = range(nFOV*len(timePoints)*len(sampleConcs))
sampleGrid = cartesian((timePoints, sampleConcs))
comboNum = len(sampleGrid)
breakNum = len(dfKeys)/len(timePoints)
moduloKeys = []
for p in range(len(dfKeys)):
    hold = p % 3
    moduloKeys.append(hold)


paths = []
orderTp = []
orderConc = []

# Organize all of the directory paths into a single variable, print out
# the order of the directories the program will search through.
# Timepoint is always the parent directory.
for i in timePoints:
    for j in sampleConcs:
        paths.append(workDir + '\\' + i + '\\' + j + '\\' + 'Nudged' + '\\')
        orderTp.append(i)
        orderConc.append(j)

print('this is the order: ')
print(orderTp)
print(orderConc)
print('end order')
print('this is the precise order of directories to be analysed: ')
print(paths)

# This code properly indexes each results file by concentration and timepoint,
# and aggregates all the results into a single, key-indexed dataframe.
counter = 0
beans = 0
aggRes = []
aggOnRes = []
for k in paths:
    for m in range(nFOV):
        resLoc = k
        resFile = pd.read_csv(resLoc + '\\' + 'analysis' + '\\' + 'Results' + str(m) + '.csv')
        if counter < comboNum:
            currentConc = sampleGrid[counter][1]
            currentTp = sampleGrid[counter][0]
        elif counter >= comboNum:
            cRep = counter % comboNum
            currentConc = sampleGrid[cRep][1]
            currentTp = sampleGrid[cRep][0]
        resFile['tp'] = [currentTp] * len(resFile.index)
        resFile['conc'] = [currentConc] * len(resFile.index)
        aggRes.append(resFile)
    counter = counter + 1
allRes = pd.concat(aggRes, keys=dfKeys)
allNorms = pd.DataFrame(columns=['zz', 'zzz'])
allAvg = pd.DataFrame(columns=['zz', 'zzz'])

# This code calculates the normalized intensity increase in the sample image.
for key in dfKeys:
    if key < comboNum:
        currentTp = sampleGrid[key][0]
        currentConc = sampleGrid[key][1]
    elif key >= comboNum:
        keyRep = key % comboNum
        currentTp = sampleGrid[keyRep][0]
        currentConc = sampleGrid[keyRep][1]
    tempDf = allRes.loc[key]
    # The precise 'Slice' number associated with each image will depend
    # on the order of stacking in ImageJ. The sample image should always be
    # loaded into the top slice to ensure the best results. The ordering
    # of the positive and negative controls is less important.
    sampleDf = tempDf[tempDf['Slice'] == 1]
    ionoDf = tempDf[tempDf['Slice'] == 2]
    blankDf = tempDf[tempDf['Slice'] == 3]
    ionoMean = ionoDf['Mean'].values
    sampleMean = sampleDf['Mean'].values
    blankMean = blankDf['Mean'].values
    someNorms = []
    avgInt = []
    for p in range(len(ionoDf)):
        normInt = (sampleMean[p] - blankMean[p])/(ionoMean[p] - blankMean[p])*100
        if normInt > 0 and normInt < 100:
            someNorms.append(normInt)
        else:
            continue
    avgInt.append(np.average(someNorms))
    avgIntdict = pd.DataFrame(columns=['idx'])
    tempAvg = pd.DataFrame(columns=['idx'])
    wtfTp = tempDf['tp'].values[0]
    wtfConc = tempDf['conc'].values[0]
    # Put together the final dataframe piece by piece.
    tempAvg['int'] = avgInt
    tempAvg['key'] = key*len(tempAvg.index)
    tempAvg['tp'] = [wtfTp] * len(tempAvg.index)
    tempAvg['conc'] = [wtfConc] * len(tempAvg.index)
    avgIntDf = pd.DataFrame.from_dict(avgIntdict, orient='index')
    avgIntDf['tp'] = [wtfTp] * len(avgIntDf.index)
    avgIntDf['conc'] = [wtfConc] * len(avgIntDf.index)
    allAvg = pd.concat([allAvg, tempAvg])
    allNorms = pd.concat([allNorms, tempNormsDf])
allAvg = allAvg.drop(['zz', 'zzz'], axis=1)
calcDf = pd.DataFrame(allNorms)

# Raw Normalized Intensity Output
allAvg.to_csv('LipoAnalysis' + '_' + expName + '_' + expDate + '.csv', sep=',')
print('done!')
