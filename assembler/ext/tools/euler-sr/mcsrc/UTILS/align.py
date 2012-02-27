#!/usr/bin/env python

import sys
import string

A = {'a': 0, 't': 1, 'c': 1, 'g': 1}
C = {'a': 1, 't': 1, 'c': 0, 'g': 1}
T = {'a': 1, 't': 0, 'c': 1, 'g': 1}
G = {'a': 1, 't': 1, 'c': 1, 'g': 0}

d = {'a': A, 't': T, 'c': C, 'g': G}

MAX_SCORE = 99999999999;

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


def getAlignmentStrings(pathMat, stra, strb):
    r = len(stra)
    c = len(strb)

    top = ''
    middle = ''
    bottom = ''
    while (r > 0 and c > 0):
        # check for match
        if (pathMat[r][c] == 0):
            top = stra[r-1] + top
            if (string.upper(stra[r-1]) == string.upper(strb[c-1])):
                middle = '|' + middle
            else:
                middle = 'X' + middle
            bottom = strb[c-1] + bottom
            r = r - 1
            c = c - 1
        elif(pathMat[r][c] == 2):
        # check for delete
            top = '-' + top
            middle = ' ' + middle
            bottom = strb[c-1] + bottom
            c = c  - 1
        elif(pathMat[r][c] == 1):
        # insert
            top = stra[r-1] + top
            middle = ' ' + middle
            bottom = '-' + bottom
            r = r - 1

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
        top = top[70:]
        middle = middle[70:]
        bottom = bottom[70:]
    
    
if (len(sys.argv) < 2):
    print 'usage: kband a b '
    sys.exit(0)


seqAName = sys.argv[1]
seqBName = sys.argv[2]

seqA = readFasta(seqAName)
seqB = readFasta(seqBName)

rows = len(seqA)
cols = len(seqB)

score = []
path = []
for i in range(0,rows+1):
    row = []
    rowb = []
    for j in range(0,cols +1):
        row.append(0)
        rowb.append(0)
    score.append(row)
    path.append(rowb)
    
for i in range(1, cols+1):
    score[0][i] = MAX_SCORE

for j in range(1,rows+1):
    score[j][0] = MAX_SCORE

for r in range(1,rows+1):
    for c in range(1, cols+1):  # first +1 is to count from 1 to rightcol, second is to include all k cols
        # the indices into the strings
        diag = score[r-1][c-1] + d[string.lower(seqA[r-1])][string.lower(seqB[c-1])]
        insert = score[r-1][c] + 1
        delete = score[r][c-1] + 1
        score[r][c] = min(diag, delete, insert)

        least = min(diag, insert, delete)
        if (diag == least):
            path[r][c] = 0
            continue
        elif (insert == least):
            path[r][c] = 1
            continue
        else:
            path[r][c] = 2
            continue
        
        
print 'got score: ', score[rows][cols]

top, middle, bottom = getAlignmentStrings(path, seqA, seqB)

printAlignmentStrings(top, middle, bottom, 70)
