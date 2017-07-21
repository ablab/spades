
import sys 
import os 
import string 
import re 
import subprocess 
import datetime 
import fastaparser 
from os.path import join 
from genericpath import isdir, exists
from sets import Set
'''
Extracts suspicious linear viruses - contigs from tip to tip with coverage > cov(5) and length > l (1000)
'''

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print ("script extracts suspicious linear viruses - contigs from tip to tip with coverage > cov(5) and length > l (1000)")
        print ("Usage: " + sys.argv[0] + "<.paths> <.fastg> coverage(optional) length(optional)")
        exit(0)
    paths = sys.argv[1]
    fastg = sys.argv[2]
    tips = [Set(),Set()]
    cov_limit = 5
    l_limit = 1000
    if len(sys.argv) >3:
        cov_limit = int(sys.argv[3])
    if len(sys.argv) > 4:
        l_limit = int(sys.argv[4])
    f_fastg = open(fastg, 'r')
    for line in f_fastg:
#>EDGE_63500766_length_345_cov_10.965517:EDGE_62187737_length_61_cov_33.666667;
        if line[0] != ">":
            continue
        arr = line.split(":")
        if len(arr) == 1:
            id = int(arr[0].split("_")[1])
            
            if arr[0][-3] == "'":
                tips[0].add(id)
            else:
                tips[1].add(id)
    path_f = open(paths,'r')
    string = ""
    name = ""
    for line in path_f:
        if line.split("_")[0] == "NODE":
            if line.strip()[-1] != "'":
                if name != "":
#                    print name
#                    print string
                    arr = string.replace(";", ",").split(",")
                    
                    ends = [arr[0], arr[-1]] 
            
                    shifts = [0, 0]
                    for i in range(0, 1):
                         if ends[i][-1] == "-":
                            shifts[i] = 1
                    ids = [int(ends[0][:-1]), int(ends[1][:-1])]
#                    print ids
                    good = True
                    for i in range(0, 1):
                        good = good and ids[i] in tips[(i + shifts[i])%2]
                    if good:
                        cov = float (name.split("_")[5])
                        l = int(name.split("_")[3])
                        if cov_limit <= cov and l_limit <= l:
                            print name
                    
                name = line.strip()
            string = ""
        else:
            string+= line.strip()




