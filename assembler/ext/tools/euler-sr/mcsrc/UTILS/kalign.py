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
        
if (len(sys.argv) < 3):
    print 'usage: kband a b k'
    sys.exit(0)


seqAName = sys.argv[1]
seqBName = sys.argv[2]


k = int(sys.argv[3])
seqA = readFasta(seqAName)
seqB = readFasta(seqBName)

#from Bio import Fasta 
#seqA = open(seqAName, 'r');
#seqB = open(seqBName, 'r'); 

#fparser = Fasta.RecordParser()
#seqAIt  = Fasta.Iterator(seqA,fparser); 
#seqBIt  = Fasta.Iterator(seqB, fparser); 


if (abs(len(seqA) - len(seqB)) > k):
    print 'dist > ', k
    sys.exit(1)

# initialize scoring matrix

# a in rows, b in columns.
#
rows = len(seqA)
cols = len(seqB)
score = kmat(rows, cols, k)
path  = kmat(rows, cols, k)

for i in range(1, k+2):
    score.kset(0,i,MAX_SCORE)
    
for r in range(1,rows+1):
    leftcol = max(r - k, 1)
    rightcol = min(r + k, cols)

    # Kset boundary conditions. These are reference on the next iteration.
    if (leftcol >= 0):
        score.kset(r, leftcol-1, MAX_SCORE)
    if (rightcol < cols):
        score.kset(r, rightcol+1, MAX_SCORE)
    minscore = 99999999
    for c in range(leftcol, rightcol+1):  # first +1 is to count from 1 to rightcol, second is to include all k cols
        # the indices into the strings
        ri = r - 1
        ci = c - 1
        #            print 'consiering ', ri, ',', ci, ' of strings lengths: ', rows, ',', cols
        #            m = score.kget(r-1,c-1) + d[string.lower(seqA[ri])][string.lower(seqB[ci])]
        #            gr = score.kget(r-1,c) + 1
        #            gc = score.kget(r,c-1) + 1
        #            if (m <= gr and m <= gc and d[string.lower(seqA[ri])][string.lower(seqB[ci])] == 1):
        #                path.kset(r,c,' \ ')
        #            if (m <= gr and m <= gc and d[string.lower(seqA[ri])][string.lower(seqB[ci])] == 0):
        #                path.kset(r,c,' M ')
        #            if (gr <= m and gr <= gc):
        #                path.kset(r,c,' |')
        #            if (gc <= m and gc <= gr):
        #                path.kset(r,c,' -')
        #
        score.kset(r,c, min(score.kget(r-1,c-1) + d[string.lower(seqA[ri])][string.lower(seqB[ci])],
                            score.kget(r-1,c) + 1,
                            score.kget(r,  c-1) + 1))

        if (score.kget(r,c) < minscore):
            minscore = score.kget(r,c)

    if (minscore >= k):
        print 'd > ', k
    #            for i in range(len(score.mat)):
    #                print score.mat[i]


    #    for i in range(len(score.mat)):
    #        print score.mat[i]
    #
    #    for i in range(len(path.mat)):
    #        print path.mat[i]
    #
print 'got score: ', score.kget(rows, cols)
