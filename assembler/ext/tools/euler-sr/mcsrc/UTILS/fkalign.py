#!/usr/bin/env python

import sys
import string
import numarray


A = {'a': 0, 't': 1, 'c': 1, 'g': 1, 'n': 1}
C = {'a': 1, 't': 1, 'c': 0, 'g': 1, 'n': 1}
T = {'a': 1, 't': 0, 'c': 1, 'g': 1, 'n': 1}
G = {'a': 1, 't': 1, 'c': 1, 'g': 0, 'n': 1}
N = {'a': 1, 't': 1, 'c': 1, 'g': 1, 'n': 0}

d = {'a': A, 't': T, 'c': C, 'g': G, 'n': N}

MAX_SCORE = 99999999999;

def getAlignmentStrings(pathMat, stra, posa, strb, posb):
    r = posa
    c = posb
#    print 'starting with r: ' + str(r) + ' and ' + str(c)
    top = ''
    middle = ''
    bottom = ''
    if (posa < len(stra)):
        # a is gapped within b.
        for i in range(len(stra), r, -1):
            top = top[i-1] + top
            middle = ' ' + middle
            bottom = ' ' + bottom
    if (posb < len(strb)):
        for i in range(len(strb), c, -1):
            top = ' ' + top
            middle = ' ' + middle
            bottom = strb[i-1] + bottom
    while (r > 0):
        # check for match
        if (pathMat[r][c] == 0):
            top = stra[r-1] + top
#            print 'accessing : ' + str(r-1) + ' and ' + str(c-1)
            if (string.upper(stra[r-1]) == string.upper(strb[c-1])):
                middle = '|' + middle
            else:
                middle = 'X' + middle
            bottom = strb[c-1] + bottom
            r = r - 1
            c = c - 1
        # check for delete
        elif(pathMat[r][c] == 2):
            top = '-' + top
            middle = ' ' + middle
            bottom = strb[c-1] + bottom
            c = c  - 1
        # insert
        elif(pathMat[r][c] == 1):
            top = stra[r-1] + top
            middle = ' ' + middle
            bottom = '-' + bottom
            r = r - 1
    if (c > 0):
        #  fitting alignment ended before start of string a
        # append this.
        for i in range (c,0,-1):
            top = ' ' + top
            middle = ' ' + middle
            bottom = strb[i-1] + bottom
            
    return top, middle, bottom

def printAlignmentStrings(top, middle, bottom, width):
    while (top != ''):
        t = top[:width-1]
        m = middle[:width-1]
        b = bottom[:width-1]
        print t
        print m
        print b
        print ''
        top = top[width:]
        middle = middle[width:]
        bottom = bottom[width:]

def initKbandGrid(score, rows, cols, startCol, kband, resetScore, boundaryScore):
    # Initialize a scoring grid for k-band alignment.
    # score is the grid of dimensions rows x cols.  The dimensions are
    #    the total dimensions including any that are used for boundary
    #    conditions.
    #    
    # kband is the edit distance allowed (width = 2*kband+1),
    # resetScore is the score that resets the path on the k-band (0 for fitting alignment)
    # boundaryScore is the score of a value that should not be used in a computation.

    # initialize the first row
    for j in range(0,cols):
        score[0][j] = resetScore

    # initialize the boundarh conditions along the k-band
    startCol = startCol + 1
    leftBoundary = startCol - kband - 1
    rightBoundary = startCol + kband + 1
    for i in range(1,rows):
        if (leftBoundary < 0):
            score[i][0] = boundaryScore
        elif(leftBoundary < cols):
            score[i][leftBoundary] = boundaryScore

#        print "right boundary: ", rightBoundary, " cols: ", cols
        if (rightBoundary < cols):
            score[i][rightBoundary] = boundaryScore
            
        leftBoundary = leftBoundary + 1
        rightBoundary = rightBoundary + 1

def kbandScore(score, path, kband, startCol, seqA, seqB, rows, cols):
    for r in range(0,len(seqA)):
#        print "scoring from ", r + startCol - kband, " to ", r + startCol + kband +1
        for c in range(r + startCol - kband, r + startCol + kband +1 ):
            # if the k-bnad is off the grid, continue
            if (c >= 0 and c < len(seqB)):
                # otherwise at a valid position, compute the score
                # accordingly

                # grid indices are 1 greater than the string indices
                # because of the boundary row and columns
                gr = r+1
                gc = c+1 
                match  = score[gr-1][gc-1] + d[string.lower(seqA[r])][string.lower(seqB[c])]
                insert = score[gr-1][gc] + 1
                delete = score[gr][gc-1] + 1
                score[gr][gc] = min(match, delete, insert)
                least = min(match, insert, delete)
                if (least == match):
                    path[gr][gc] = 0
                    continue
                elif (least == insert):
                    path[gr][gc] = 1
                    continue
                else:
                    path[gr][gc] = 2
                    continue
        

def buildKmerHash(s, k, h):
    for i in range(0,len(s)-k+1):
        kmer = s[i:i+k]
        if (h.has_key(kmer)):
            h[kmer].append(i)
        else:
            h[kmer] = [i]
    
def findCommonKMer(a, b, k, p):
    # locate all k-mers in a
    akmers = {}
    bkmers = {}
#    print "k is : ", k
    buildKmerHash(a, k, akmers)
    buildKmerHash(b, k, bkmers)

    for  kmer in (akmers.keys()):
        if (bkmers.has_key(kmer)):
            for posa in (akmers[kmer]):
                for posb in (bkmers[kmer]):
                    p.append([posa, posb])
    
def delta(a, b, mat):
    return mat[a][b]

def readFasta(filename):
    file = open(filename, 'r')
    file.readline() # skip the first line
    sequences = file.readlines()
    seq = ''
    for s in sequences:
        seq = seq + s.strip()

    return seq

class kmat:
    def __init__(self,rows, cols,k):
        self.rows = rows + 1
        self.cols = cols
        self.k = k
        self.width = 2*k+3
        self.mat = []
        for i in range(0,self.rows):
            row = []
            for j in range(0,self.width):
                row.append(0)
            self.mat.append(row)

    def kget(self, i,j):
#        print 'getting ', i, ',',j, ' mapping to ', i, ',', j-i+k+1, ' size was: ', self.rows, ',', self.width, 'val: ', self.mat[i][j-i+k+1]
        return self.mat[i][j-i+k+1] # the plus 1 is for the item in the center.
            
    def kset(self,i,j,val):
        #print 'setting ', i, ',',j, ' mapping to ', i, ',', j-i+k+1, ' size was: ', self.rows, ',', self.width, ' val: ', val
        self.mat[i][j-i+k+1] = val
        
if (len(sys.argv) < 4):
    print 'usage: fkband a b kband [kmer]'
    print 'performs a fitting k-banded alignment of two sequences'
    print 'find the best fit of b inside a with at most k errors'
    sys.exit(0)


seqAName = sys.argv[1]
seqBName = sys.argv[2]

# numarray does efficient initialization of arrays, so it's ok to allocate the whole grid here 

    
kband = int(sys.argv[3])
kmer  = 10
if (len(sys.argv) == 5):
    kmer  = int(sys.argv[4])

printAlign = 0
if (len(sys.argv) == 6):
    if (sys.argv[5] == '-p'):
        printAlign = 1

printWidth = 50
seqA = readFasta(seqAName)
seqB = readFasta(seqBName)

positions = []
findCommonKMer(seqA, seqB, kmer, positions)

# initialize scoring matrix

# a in rows, b in columns.
#
rows = len(seqA)
cols = len(seqB)


MAX_SCORE = 999;

curScore = MAX_SCORE -1;
minScore = MAX_SCORE;
it = 1
# initialize values that will be computed inthe loop
srows  = rows + 1
scols  = cols + 1
global score, path, mincol
mincol = 0
score = numarray.zeros((srows, scols), numarray.Int32)  
path  = numarray.zeros((srows, scols), numarray.Int16)
#print "score1: \n", score
#print "positions: ", positions
while (curScore != minScore and len(positions) > 0):
    # pop the first item off of common
    if (curScore < minScore):
        minScore = curScore
       
    pair = positions[0]
    positions[0:] = positions[1:]

    # spurious match for sure if smaller sequence has a match after larger.
    
    if (pair[1] < pair[0]):
        continue

    startCol = pair[1] - pair[0]

    # size is rows + 1 to allow initializing the first row to 0
    # and cols+2 to allow for k-band boundaries on each side of the band
    score = numarray.zeros((srows, scols), numarray.Int32)  
    path  = numarray.zeros((srows, scols), numarray.Int16)
    initKbandGrid(score, srows, scols, startCol, kband, 0, MAX_SCORE)
#    print "score is now: \n"
#    print score
    kbandScore(score, path, kband, startCol, seqA, seqB, srows , scols)

    tempMinScore = MAX_SCORE;
    mincol = 0
    r = srows - 1
#    print "ranging from ", startCol + len(seqA) - kband, " to ", startCol + len(seqA) +kband + 1
    for c in range(startCol + len(seqA) - kband, startCol + len(seqA) +kband + 1):
        if (c > 0 and c <= len(seqB) and score[r][c] < tempMinScore ):
            tempMinScore = score[r][c]
            mincol  = c

    curScore = tempMinScore
    print 'scoring iteration: ', it, ' pair: ', pair, ' current score: ', curScore

#    print "score2:\n", score
#    print "startcol: ", startCol, " curscore ", curScore


#print 'shape of score: ', score.shape
#print "score3:\n", score
print 'got score: ', score[rows][mincol]
#@print 'got score mat: '
#print score
if (printAlign == 1):
    top, middle, bottom = getAlignmentStrings(path, seqA, len(seqA), seqB, mincol)
    printAlignmentStrings(top, middle, bottom, printWidth)

#for i in range(0,srows):
#    for j in range(0,scols):
#        sys.stdout.write(repr(score[i][j]).rjust(11) + ' ')
#    sys.stdout.write('\n')



