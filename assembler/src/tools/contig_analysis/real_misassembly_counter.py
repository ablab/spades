#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import itertools
import matplotlib
import shutil

########################################################################	
### CONFIG & CHECKS 
########################################################################

near_to_end = 100

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), 'libs'))
import fastaparser
from qutils import id_to_str

########################################################################

if len(sys.argv) < 4:
	print 'Real misassembly counter'
	print 'Usage: ', sys.argv[0], ' REFERENCE_COORDS MISASSEMBLIES_COORDS X_Z_DISTANCES'
	sys.exit(0)

ref_coords = open(sys.argv[1], 'r')
mis_coords = open(sys.argv[2], 'r')
x_z_dists  = open(sys.argv[3], 'w')

ref_ends = [] # [start, end] of contigs in references
for line in ref_coords:
    ref_ends.append((int(line.split()[0]), int(line.split()[1])))


## S1 E1 S2 E2
#1697963 1710156 | 12195 2 | 12194 12194 | 100.0000 | sum_contig NODE_1332_length_13042_cov_888867
#2418963 2419862 | 12142 13041 | 900 900 | 100.0000 | sum_contig NODE_1332_length_13042_cov_888867



misassemblies = 0
contig_breaks = 0
line = mis_coords.readline()
cur_contig = []
while line:
    if (line.strip() != ""):
        cur_contig.append(line)    
    else:
        mis_ends = []
        for mis in cur_contig:
            ends = mis.split('|')[1].strip().split()
            mis_ends.append(int(ends[0]))
            mis_ends.append(int(ends[1]))
        mis_ends.sort()
        
        while (len(mis_ends) >= 2 * 2):
            # find first misassembly
            first_mis_start = mis_ends[0]
            first_mis = ""
            for mis in cur_contig:
                ends = mis.split('|')[1].strip().split()
                if (first_mis_start == int(ends[0])) or (first_mis_start == int(ends[1])):
                    first_mis = mis
                    mis_ends.remove(int(ends[0]))
                    mis_ends.remove(int(ends[1]))
                    break
            
            coordX = 0
            # if asc order of S2 E2
            if (first_mis_start == first_mis.split('|')[1].strip().split()[0]):
                coordX = int(first_mis.split()[1])
            # if desc order of S2 E2       
            else:
                coordX = int(first_mis.split()[0])

            # find second misassembly
            second_mis_start = mis_ends[0]
            second_mis = ""
            for mis in cur_contig:
                ends = mis.split('|')[1].strip().split()
                if (second_mis_start == int(ends[0])) or (second_mis_start == int(ends[1])):
                    second_mis = mis
                    break
            
            coordZ = 0
            # if asc order of S2 E2
            if (second_mis_start == second_mis.split('|')[1].strip().split()[0]):
                coordZ = int(second_mis.split()[0])
            # if desc order of S2 E2       
            else:
                coordZ = int(second_mis.split()[1])

            # check lengths from coordX and coordY to Ns
            minX = sys.maxint
            minZ = sys.maxint
            for end in ref_ends:
                if (abs(end[0] - coordX) < minX):
                    minX = abs(end[0] - coordX)
                if (abs(end[1] - coordX) < minX):
                    minX = abs(end[1] - coordX)
                if (abs(end[0] - coordZ) < minZ):
                    minZ = abs(end[0] - coordZ)
                if (abs(end[1] - coordZ) < minZ):
                    minZ = abs(end[1] - coordZ)
            if (minX < near_to_end) and (minZ < near_to_end):
                contig_breaks += 1 
            else: 
                misassemblies += 1

            #test output:
            x_z_dists.write("contig=" + first_mis.split()[-1] + " ")
            x_z_dists.write("coordX=" + str(coordX) + " minX=" + str(minX) + " ")
            x_z_dists.write("coordZ=" + str(coordZ) + " minZ=" + str(minZ) + " ")
            if (minX < near_to_end) and (minZ < near_to_end):            
                x_z_dists.write(" CONTIG BREAK\n")
            else:
                x_z_dists.write(" TRUE MISASSEMBLY\n")             

        cur_contig = []

    line = mis_coords.readline()        
    
ref_coords.close()
mis_coords.close()

x_z_dists.write("Number of real misassemblies: " + str(misassemblies) + "\n")
x_z_dists.write("Number of contig breaks: " + str(contig_breaks) + "\n")
x_z_dists.close()

print "Number of real misassemblies: ", misassemblies
print "Number of contig breaks: ", contig_breaks
