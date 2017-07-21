
import sys 
import os 
import string 
import re 
import subprocess 
import datetime 
import fastaparser 
from os.path import join 
from genericpath import isdir, exists

'''
Takes dir with spades result as an input, 
reads saffols.fasta, for each:
opens dir with k-mer assembly (from 129 down to 50), select largest.
Checks if it overlaps in start and end.
If yes, prints out contig name with "is circular"
'''

if __name__ == "__main__":
    file = sys.argv[1]
    contigs = fastaparser.read_fasta(file)
    count = []

    for contig in contigs:
        arr = contig[0].strip(';').split('_')
#      if float(arr[3]) > 500:
        for kval in range (200,50, -1):
#            kval = 55
            if kval >= len(contig[1]) or len(contig[1]) < 500:
                continue
            start = contig[1][:kval]
            end = contig[1][-kval:]
                   
            if start == end:
#               print (">" + contig[0][1:])
#               print (contig[1])
#                print (" k equal " + str(kval))
                print (contig[0] + " is circular")
                break   



