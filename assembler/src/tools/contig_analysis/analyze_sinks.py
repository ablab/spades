#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys 

def parse_name(name):
      pieces = name.strip().split('_')
      return int(pieces[len(pieces) - 1])


def find_length(filename):
      inFileName = filename
      inFile = open(inFileName)
     
      length = 0
      for line in inFile:
	      parse = line.strip().split(' ')
	      if parse[0] == "Total" and parse[1] == "Region":
		      length = int(parse[3])
		      break
      print length
      return length
      
def make_snap(filename, length):
      inFileName = filename
      inFile = open(inFileName)
      nucls = [0]*(length+1)       
      contig_id = 0
      shift = 7
      for line in inFile:
	      parse = line.strip().split(' ')
	      if parse[0] == "CONTIG:":
		      name = parse[1]    
		      contig_id = parse_name(name)
	      if parse[0] == "One":
		      for i in range(int(parse[shift]),int(parse[shift + 1]) + 1):
			       nucls[i] = contig_id
			       if int(parse[shift + 3]) > int(parse[shift + 4]):     
			                  nucls[i] = nucls[i] * (-1)
			                  
	      if parse[0] == "Ambiguous" and parse[1] == "Alignment:":
		      for i in range(int(parse[2]),int(parse[3]) + 1):
				nucls[i] = contig_id
				if int(parse[5]) > int(parse[6]):     
					    nucls[i] = nucls[i] * (-1)
			                   
	      if parse[0] == "Real" and parse[1] == "Alignment":
		      for i in range(int(parse[3]),int(parse[4]) + 1):
			       nucls[i] = contig_id
			       if int(parse[6]) > int(parse[7]):     
			                   nucls[i] = nucls[i] * (-1)                      
	
	      if parse[0] == "Marking" and parse[2] == "ambiguous:":
		      for i in range(int(parse[3]),int(parse[4]) + 1):
			       nucls[i] = contig_id
			       if int(parse[6]) > int(parse[7]):     
			                   nucls[i] = nucls[i] * (-1) 
			                   
      inFile.close()
      nucls[0] = 1   
      contig_id1 = -1
      contig_id2 = -1
      result = []
      for i in range(1,length+1):
	    if (nucls[i] > 0):
	        contig_id1 = contig_id2
	        contig_id2 = nucls[i]
	        
	    if (contig_id1 > 0 and contig_id2 > 0 and not(contig_id1 == contig_id2)):
	       result.append([contig_id1, contig_id2])
	       
      return result
      
def find_first(path_id):
     inFile = open("../../../data/debruijn/ECOLI_SC_LANE_1_woHUMAN/K55/latest/paths.dat", "r")
     now = 0
     for line in inFile:
	  parse = line.strip().split(' ')
	  if (now == 1):
	      print parse[0]
	      inFile.close()
	      return parse[0]
	  if (parse[0] == "PATH" and parse[1] == str(path_id)):
	      now = 1    
     inFile.close()
      
def print_pairs_of_edges(pairs, filename):
     outFile = open("../../../sinks_with_sources.log", "w")
     for i in range(0,len(pairs)):
          inFile = open(filename, "r")
          source = '0'
          for line in inFile:
	      parse = line.strip().split(' ')
	      if (parse[0] == str(pairs[i][1])):
		 source = parse[1] 
          inFile.close()
          inFile = open("../../../sinks.log", "r")
          for line in inFile:
	      parse = line.strip().split(' ')
	      if (parse[0] == str(pairs[i][0])):
	         outFile.write(str(parse[1])+ " Must: " + str(find_first(pairs[i][1])) + " Was found: " +source + "\n")
          inFile.close()
     outFile.close()
     return 

if len(sys.argv) != 3:
        print "Print to sinks_with_sources.log next edges for sinks"
        print "Usage: " + sys.argv[0] + " <plantagora output>  scafolder_output"
        sys.exit()
       
length = find_length(sys.argv[1])       
pairs_of_contigs = make_snap(sys.argv[1], length)
print_pairs_of_edges(pairs_of_contigs, sys.argv[2])
