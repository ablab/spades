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
      
def make_snap(filename, length, file_with_genome):
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
		      #print "ID: " + str(contig_id) + " : " + parse[shift]+ " - " +parse[shift + 1]
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
      
      l = 0
      r = 0
      count = 0
      outFile = open("../../../data/debruijn/contigs.fasta", "w")
      inFile = open(file_with_genome, "r")
      inFile.readline()
      genome = []
      for line in inFile:     
	      genome.extend(line.strip())
      ind = 0
      while (ind < len(genome)):
	      if not (genome[ind] in ['A', 'C', 'G', 'T']) :
		     genome.remove(genome[ind])
 	      else:     
	             ind = ind + 1 	    
      print len(genome)
      pos = 0
      cons = 100
      for i in range(1,length+1):
	      if nucls[i] == 0 and (not nucls[i - 1] == 0):
		      l = i    
              if nucls[i] == 0 and (i + 1 > length or (not nucls[i + 1] == 0)):
		      if (l - cons > 0):
                              l = l - cons
                      else:
			      l = 1 
	              r = i      
	              if (r + cons > length):
	                      r = length
	              else:
			      r = r + cons
		      count=count + 1
		      outFile.write(">" + str(count) +  "\n")# + str(l) + " " + str(r)+"\n")
	              for j in range(l-1,r): 
		             outFile.write(genome[j])
		      outFile.write("\n")   
      print count      
      inFile.close()	      
      outFile.close()
      return nucls
                

if len(sys.argv) != 3:
        print "Saves places with gaps in genome to ../../../data/debruijn/contigs.fasta"
        print("Usage: " + sys.argv[0] + " <plantagora output> <file_with_genome>" )
        sys.exit(1)
       
length = find_length(sys.argv[1])          
gaps_before = make_snap(sys.argv[1], length, sys.argv[2])

