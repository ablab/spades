#!/usr/bin/env python2.5
# -*- coding: utf-8 -*-

import sys
import math
import random


# Define functions
def importStats(fin_stats):
    dicStats = {}
    listHeader = []

    while True:
        line = fin_stats.readline()
        if not line:
            break

        if len(dicStats) == 0:
            listHeader = line.rstrip("\n").split("\t")
            for header in listHeader:
                dicStats[header] = []
        else:
            listStats = line.rstrip("\n").split("\t")
            for i in range(len(listStats)):
                if i in [0, 1, 2, 3, 9, 10, 11]:
                    stats = int(listStats[i])
                else:
                    stats = float(listStats[i])
                dicStats[listHeader[i]].append(stats)

    return dicStats


def weightedHisto(dicStats, xMin, xMax, binWidth):
    dicHisto = {}
    listShort1Cov = dicStats["short1_cov"]
    listLgth = dicStats["lgth"]

    for x in range(xMin, xMax, binWidth):  
        dicHisto[x] = 0

    for i in range(len(listShort1Cov)):
        cov = listShort1Cov[i]
        if cov < xMin or cov >= xMax:
            continue
        for x in range(xMin, xMax+binWidth, binWidth):  
            if (cov >= x and cov < x + binWidth):
                dicHisto[x] += listLgth[i]

    return dicHisto


def smoothingHisto(dicHisto, xMin, xMax, binWidth, widthMovAve):
    dicSmoothHisto = {}
    listMovAve = []

    for x in range(xMin, xMax, binWidth):
        listMovAve.append(dicHisto[x])
        if len(listMovAve) < widthMovAve:
            continue
        dicSmoothHisto[x - binWidth * ((widthMovAve - 1) / 2)] \
            = sum(listMovAve) / float(widthMovAve)
        listMovAve.pop(0)

    return dicSmoothHisto


def printHisto(dicHisto, xMin, xMax, binWidth):
    for x in range(xMin, xMax, binWidth):
        #print str(x) + " : " + str(int(round(dicHisto[x], 0)))
        lenBar = int(round((dicHisto[x] / 20000), 0)) - 1
        print str(x) + "\t",
        for i in range(lenBar):
            print "=",
        print "\n",
    print "\n",


def setXMax(xMax, binWidth):
    return int((math.floor(xMax / binWidth)) * binWidth)
    
    
def getFirstXMax(dicStats, binWidth, thresConLen):
    listLgth = dicStats["lgth"]
    listShort1Cov = dicStats["short1_cov"]
    maxCov = 0
    subMaxCov = 0

    for i in range(len(listLgth)):
        if listLgth[i] >= thresConLen:
            if listShort1Cov[i] > maxCov:
                subMaxCov = maxCov
                maxCov = listShort1Cov[i]

    xMax = setXMax(subMaxCov, binWidth) + binWidth * 5
    return xMax


def getN50(tupleConLen):
    listSortedConLen = list(tupleConLen)
    listSortedConLen.sort()
    listSortedConLen.reverse()
    totalLen = sum(listSortedConLen)
    sumLen = 0

    for i in range(len(listSortedConLen)):
        sumLen += listSortedConLen[i]
        if sumLen >= totalLen / 2:
            return listSortedConLen[i]
    
    return -1


def setWidthByXMax(xMax):
    listWidth = [0, 0]  # [binWidth, widthMovAve]

    if xMax > 300:
        listWidth = [6, 5]
    if xMax <= 300:
        listWidth = [4, 3]
    if xMax <= 120:
        listWidth = [2, 3]
    if xMax <= 100:
        listWidth = [1, 1]

    return listWidth


def detectPeakPandS(dicHisto, xMin, xMax, binWidth, 
                    thresHeight, listPeakPandS):
    countIncrease = 0; thresIncrease = 3
    countDecrease = 0; thresDecrease = 3
    beforeHeight = -1
    flagPeakStart = False
    peakHeight = 0; peakCov = 0

    for x in range(xMax - binWidth, xMin - binWidth, -1 * binWidth):
        if beforeHeight == -1:
            beforeHeight = dicHisto[x]
            continue
        
        if not flagPeakStart:
            if dicHisto[x] >= thresHeight:
                if dicHisto[x] >= beforeHeight:
                    countIncrease += 1
                    if countIncrease >= thresIncrease:
                        countIncrease = 0
                        flagPeakStart = True
            beforeHeight = dicHisto[x]
                        
        if flagPeakStart:
            if dicHisto[x] >= peakHeight:
                peakHeight = dicHisto[x]
                peakCov = x
            else:
                countDecrease += 1
                if countDecrease >= thresDecrease:
                    for i in range(2):
                        if listPeakPandS[i] == -1:
                            tmpBias = float(binWidth) / 2
                            listPeakPandS[i] = peakCov + tmpBias
                            peakHeight = 0; peakCov = 0
                            break
                    if listPeakPandS[1] != -1:
                        return listPeakPandS
                    countDecrease = 0
                    flagPeakStart = False

    return listPeakPandS



# ---------- Main function ----------

# Import stats file
fin_stats = open(sys.argv[1], "r")
dicStats = importStats(fin_stats)

# Make weighted histogram
listPeak = []
xMin = 0
xMax = 1000
binWidth = 4
widthMovAve = 5
listPeakPandS = [-1, -1]
N50 = 0
thresHeight = 0
thresConLen = 0

while True:
    # Get N50
    if len(listPeak) == 0:
        N50 = getN50(tuple(dicStats["lgth"]))
        print "N50 : " + str(N50)
        thresConLen = N50 * 5

    # Get first xMax
    if len(listPeak) == 0:
        xMax = getFirstXMax(dicStats, binWidth, thresConLen)
        print "First xMax : " + str(xMax)
        
    # Set width and xMax
    listWidth = setWidthByXMax(xMax)
    binWidth = listWidth[0]; widthMovAve = listWidth[1]
    xMax = setXMax(xMax, binWidth)

    # Make weighted and smoothed histogram
    xMin = 0
    dicHisto = weightedHisto(dicStats, xMin, xMax, binWidth)
    dicSmoothHisto = smoothingHisto(dicHisto, xMin, xMax, 
                                    binWidth, widthMovAve)
    xMin += binWidth * ((widthMovAve - 1) / 2)
    xMax -= binWidth * ((widthMovAve - 1) / 2)

    # Get thresHeight
    if len(listPeak) == 0:
        thresHeight = dicSmoothHisto[xMax - binWidth]
        print "Thres Height : " + str(thresHeight)

    # Print histogram
    if len(listPeak) == 0:
        printHisto(dicSmoothHisto, xMin, xMax, binWidth)

    # Detect (primary and) secondary peak
    listPeakPandS = detectPeakPandS(dicSmoothHisto, xMin, xMax, binWidth, 
                                    thresHeight, listPeakPandS)

    # Record peak
    if len(listPeak) == 0:
        listPeak.append(listPeakPandS[0])
    listPeak.append(listPeakPandS[1])

    # When couldn't detect secondary peak, break
    if listPeakPandS[1] == -1:
        listPeak.pop(-1)
        print listPeak
        break

    # Prepare for next peak
    listPeakPandS[0] = listPeakPandS[1]
    listPeakPandS[1] = -1
    xMax = listPeakPandS[0]
