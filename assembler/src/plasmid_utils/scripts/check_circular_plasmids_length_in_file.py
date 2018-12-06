
import sys
import os
import string
import re
import subprocess
import datetime
import fastaparser
from os.path import join
from genericpath import isdir, exists


if __name__ == "__main__":
	fullname = sys.argv[1]
  #  output_dir = sys.argv[2]
  #  os.mkdir(output_dir)
 #   dirs = os.listdir(main_dir)
 #   for dir in dirs:
#        fullname = join(main_dir, dir,'scaffolds.fasta')
	if exists(fullname):
            contigs = fastaparser.read_fasta(fullname)
            count = []
            #print (fullname)
            for i in range (0, 2*len(contigs)):
                 count.append(0)	
#            for contig in contigs:
             #   print contig[0]
 #               arr = contig[0].strip(';').split('_')
              #  print len (count) 
               # print arr[-1]
  #              count[int(arr[-1])] += 1
            for contig in contigs:
                arr = contig[0].strip(';').split('_')
               # comp_id = int(arr[-1])
                length = float (arr[3])
              #  print str(comp_id) + " " + str(count[comp_id])
                if  length > 500:
                    print contig[0]
		    
  #                  fastaparser.write_fasta_to_file((os.path.splitext(fullname)[0] + "_circular.fasta"), [contig])


