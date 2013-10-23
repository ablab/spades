#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import random
from xml.dom.minidom import parse, getDOMImplementation

def getInputData(finName):
    fin = open(finName)
    data = []
    line = fin.readline()
    while line:
        data.append(line.strip())
        line = fin.readline()
    fin.close()
    return data

def getAllChunkNames(data):
    names = []
    for string in data:
        for ltr in string:
            if ltr != "+" and ltr != "-":
                names.append(ltr)
    return set(names)

def generateRandomSequence(seqLen):
    seq = ""
    nucl = ["A", "C", "G", "T"]
    for x in range(seqLen):
        i = random.randint(0, 3)
        seq += nucl[i]
    return seq

def getChunkSequences(longChunkLen, shortChunkLen, chunkNames):
    seqs = {}
    for name in chunkNames:
        assert len(name) == 1
        lenToGen = longChunkLen
        if name.islower():
            lenToGen = shortChunkLen
        seqs[name] = generateRandomSequence(lenToGen)
    return seqs

def getReverseComplement(seq):
    revNucl = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}
    revSeq = ""
    for x in reversed(seq):
        revSeq += revNucl[x]
    return revSeq

def getSequenceByChunkString(chunkString, chunkSeqs):
    seq = ""
    for i in range(0, len(chunkString)-1, 2):
        ltr = chunkString[i]
        chunkName = chunkString[i+1]
        if ltr == "+":
            seq += chunkSeqs[chunkName]
        elif ltr == "-":
            seq += getReverseComplement(chunkSeqs[chunkName])
    return seq

def getAllContigData(contigNodes):
    data = []
    for contigNode in contigNodes:
        for node in contigNode.childNodes:
            if node.nodeType == node.TEXT_NODE:
                data.append(node.data)
    return data

if __name__ == "__main__":
    random.seed(239)
    if len(sys.argv) < 4:
        print("Usage: " + sys.argv[0] + " <name of input file> <name of output file> <length of chunk> (<length of chunk for lower case>)")
        sys.exit()
    #random.seed()
    chunkLen = int(sys.argv[3])
    shortChunkLen = chunkLen
    if len(sys.argv) > 4 
        shortChunkLen = int(sys.argv[4])
    dom = parse(sys.argv[1])
    contigNodes = dom.getElementsByTagName("contig")
    chunkNames = getAllChunkNames(getAllContigData(contigNodes))
    chunkSeqs = getChunkSequences(chunkLen, shortChunkLen, chunkNames)
    for node in contigNodes:
        for child in node.childNodes:
            child.data = getSequenceByChunkString(child.data, chunkSeqs)
    fout = open(sys.argv[2], "w")
    dom.writexml(fout)
    fout.close()   
