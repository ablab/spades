
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
#    main_dir = sys.argv[1]
#    dirs = os.listdir(main_dir)
    fullname = join(sys.argv[1])
    if exists(fullname):
        contigs = fastaparser.read_fasta(fullname)
        for contig in contigs:
            for kval in range (57,54, -1):
                start = contig[1][:kval]
                end = contig[1][-kval:]
                if start == end: 
                    print contig[0] + " is circular " + str(kval)     
