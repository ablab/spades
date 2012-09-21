#!/usr/bin/env python

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

def getChunkSequences(chunkLen, chunkNames):
    seqs = {}
    for name in chunkNames:
        seqs[name] = generateRandomSequence(chunkLen)
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
        print("Usage: " + sys.argv[0] + " <name of input file> <name of output file> <length of chunk>")
        sys.exit()
    random.seed()
    chunkLen = int(sys.argv[3])
    dom = parse(sys.argv[1])
    contigNodes = dom.getElementsByTagName("contig")
    chunkNames = getAllChunkNames(getAllContigData(contigNodes))
    chunkSeqs = getChunkSequences(chunkLen, chunkNames)
    for node in contigNodes:
        for child in node.childNodes:
            child.data = getSequenceByChunkString(child.data, chunkSeqs)
    fout = open(sys.argv[2], "w")
    dom.writexml(fout)
    fout.close()    
    
        
