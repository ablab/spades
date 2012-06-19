#!/usr/bin/env python2.5
# -*- coding: utf-8 -*-

import sys
import re

listLine = ["", ""]
fin_fw = open(sys.argv[1], "r")
fin_bw = open(sys.argv[2], "r")
fout = open(sys.argv[3], "w")

while True:
    for i in range(2):
        for j in range(4):
            if i == 0: 
                line = fin_fw.readline()
            else:
                line = fin_bw.readline()
        
            if not line:
                sys.exit()

            if j >= 2:
                continue

            if i == 0 and j == 0:
                line = ">" + re.sub("^@|\n", "", line) + "/1\n"
            elif i == 1 and j == 0:
                line = ">" + re.sub("^@|\n", "", line) + "/2\n"
            
            fout.write(line)

    
