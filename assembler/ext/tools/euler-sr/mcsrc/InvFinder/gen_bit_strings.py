#!/usr/bin/env python

import sys;
import re;

if (len(sys.argv) < 2):
    print 'must enter a coords file\n'
    sys.exit(0)

coordsFile = open(sys.argv[1])
verbose = 0
if (len(sys.argv) > 2):
    opt = sys.argv[2]
    if (opt == '-v'):
        verbose = 1
    
# this line defines the pairwise comparison of two sequences
indRe = re.compile(r'(\w+) (\w+) (\-?\d+)')

line = coordsFile.readline()
refs = {}

rangeStart = []
rangeEnd   = []
rangeCount = []
rangePairwise = []
rangeP = []
binnedCoords = []
binnedCoordStart = []
binnedCoordEnd   = []
species = []
speciesBins = {}
pairwiseInversions = []

# read in the coordinates.  Each inversion is specified by a start and end
# coordinate on two separate lines.  Read in each coordinate
while (line != None and line != ""):
    s1r = 0
    s1q = 0
    indStart = 0
    s2r = 0
    s2q = 0
    indEnd = 0
    # parse the line
    indRes = indRe.search(line)
    if (indRes == None):
        print 'error parsing index line ', line
        sys.exit(0)
    else:
        # store the 
        s1r = indRes.group(1)
        s1q = indRes.group(2)
        indStart =  int(indRes.group(3))
    
    if (species.count(s1r) == 0):
        species.append(s1r)
        speciesBins[s1r] = []
    if (species.count(s1q) == 0):
        species.append(s1q)
        speciesBins[s1q] = []

    line = coordsFile.readline()

    indRes = indRe.search(line)
    if (indRes == None):
        print 'error parsing index line ', line
        sys.exit(0)
    else:
        # store the 
        s2r = indRes.group(1)
        s2q = indRes.group(2)
        indEnd =  int(indRes.group(3))
        
    # double check that the two lines correspond to the
    # same sequences
    if (s1r != s2r or s1q != s2q):
        print "paired lines do not match species ", s1r, s2r, s1q, s2q
        sys.exit(0)


    # Check to see if the inversion has been detected yet or not
    outsider = 1
    i = 0;

    
    # if the inversion was found in the reverse strand of human, the two
    # coordinates will be out of order, put them back here.
    if (indEnd < indStart):
        temp = indStart
        indStart = indEnd
        indEnd   = temp
        if (verbose):
            print 'swapped intervals '
    if (verbose):
        print 'considering ', s1r, s1q, indStart, indEnd
    while (i < len(rangeStart) and outsider == 1):
        # test 1.  Does this interval overlap with the current alignment
        if ((indStart >= rangeStart[i] and indStart <= rangeEnd[i]) or
            (indEnd >= rangeStart[i] and indEnd <= rangeEnd[i]) or
            (indStart <= rangeStart[i] and indEnd >= rangeEnd[i])):
            # this inversion overlaps with other inversions. Possibly modify the
            # interval
            if (verbose):
                print s1r, s1q, "found interval " , indStart , " " , indEnd , " overlapping " , rangeStart[i], " " , rangeEnd[i]
            # test 2.  Does this interval have the same order of magnitude of size the
            # inversion.
            lenInterval = rangeEnd[i] - rangeStart[i];
            lenInversion = indEnd - indStart;
            if ((lenInterval == lenInversion) or
                (lenInversion < lenInterval and ((lenInversion / (lenInterval*1.0)) > 0.1)) or
                (lenInversion > lenInterval and ((lenInterval  / (lenInversion*1.0))> 0.1))):
                if (verbose):
                    print ' len interval ', lenInterval, ' len inversion ', lenInversion, ' accepted at ', rangeStart[i], rangeEnd[i]
                if (indStart < rangeStart[i]):
                    rangeStart[i] = indStart
                    rangeP[i][0] = indStart
                if (indEnd > rangeEnd[i]):
                    rangeEnd[i] = indEnd
                    rangeP[i][1] = indEnd;
                outsider = 0
                rangeCount[i] = rangeCount[i] + 1
                rangePairwise[i].append([s1r, s1q])
                speciesBins[s1r].append(i)
                speciesBins[s1q].append(i)
            else:
                if (verbose):
                    print 'rejected: ', s1r, s1q, 'len interval ', lenInterval, ' len inversion ', lenInversion
                    if (lenInversion > 0):
                        print lenInterval / (1.0)*lenInversion
                    if (lenInterval > 0):
                        print lenInversion / (1.0)*lenInterval
        i = i + 1
        
    # done looking through the all intervals.  If one was not found
    # then this is an outsider, and a new interval needs to be created.
    if (outsider == 1):
        if (verbose):
            print "new interval: ", s1r, s1q, indStart, indEnd, indEnd - indStart
        rangeStart.append(indStart)
        rangeEnd.append(indEnd)
        rangeCount.append(1)
        rangeP.append([indStart, indEnd, len(rangeP)])
        rangePairwise.append([])
        rangePairwise[len(rangePairwise)-1].append([s1r, s1q])
        
    line = coordsFile.readline()

# print in a manner that phylip can understand
print len(species), len(rangeStart)

# print the coordinates of the inversion.
#import copy
#rsc = copy.copy(rangeP)
#rsc.sort()
#for i in (range(0,len(rsc))):
#    print i, rsc[i][0], rsc[i][1], rsc[i][2]

 
for i in (range(0,len(rangeStart))):
    print i, rangeStart[i] , rangeEnd[i], rangeCount[i]

#sys.exit(0)

# print each species that has the inversion
speciesBits = {}
for s in (species):
    speciesBits[s] = []
    for j in (range(0,len(rangeStart))):
        speciesBits[s].append(0)
    for bin in (speciesBins[s]):
        speciesBits[s][bin] = 1
    title = "%s" % s
    sys.stdout.write(title + " ")
    for b in (speciesBits[s]):
        sys.stdout.write("%d " % b)
    sys.stdout.write("\n")

for i in range(0,len(rangeStart)):
    print len(rangePairwise[i])*2
    for j in range(0, len(rangePairwise[i])):
        sys.stdout.write('%s ' % rangePairwise[i][j][0])
        sys.stdout.write('%s ' % rangePairwise[i][j][1])
    sys.stdout.write('\n');

