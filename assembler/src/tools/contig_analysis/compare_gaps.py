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
      return int(pieces[1])


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
		      print "ID: " + str(contig_id) + " : " + parse[shift]+ " - " +parse[shift + 1]
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
     
      toprint = ""
      tocount = 0
      contig_id = -1
      total_len = 0
      num = 0
      print "GAPS: "
      for i in range(1,length+1):
	      if nucls[i] == 0 and (not nucls[i - 1] == 0):
		      toprint = str(i)+" - "
		      tocount = i    
		      contig_id = nucls[i - 1]
              if nucls[i] == 0 and (i + 1 > length or (not nucls[i + 1] == 0)):
                      num = num + 1
                      print toprint+str(i)+" " + str(i-tocount+1)     
                      total_len = total_len + i - tocount + 1
      
      print "Total length: " + str(total_len)    
      print "Total number: " + str(num)
      return nucls
  
def find__gaps(nucls1, nucls2, length):  

      as_before =  False
      piece_len = 0
      total_len = 0
      for i in range(1, length + 1):    
               piece_len = piece_len + 1
               if (not nucls1[i] == nucls1[i - 1] or nucls2[i] == 1):
		       if as_before == True:
			  print " End position: "+str(i) 
			  print " Total length: " + str(piece_len)	
		       as_before = False
		       piece_len = 0
               if (nucls2[i] == 0 and (not nucls1[i] == 0) and (not as_before)):
	               print "ID: " + str(nucls1[i]) 
	               print " Start position: " + str(i)
	               as_before = True          
      
      
      

                

if len(sys.argv) != 2:
        print "Prints gaps for given plantagora output"
        print("Usage: " + sys.argv[0] + " <plantagora output>")
        sys.exit()
       
length = find_length(sys.argv[1])       
gaps_before = make_snap(sys.argv[1], length)
#print "*************************************** \nSecond:"
#gaps_after = make_snap(sys.argv[2], length)
#print "####################################### \nCounting gaps ..."
#find__gaps(gaps_before, gaps_after, length)
