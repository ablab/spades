#!/usr/bin/python

import sys
import os
import math

def dist_from_logs(log, max_is = 1000):
    logfile = open(log, "r")

    ids = {}

    is_paired = False
    for line in logfile:
        id1 = line.split('/', 1)
        if len(id1) > 1 and id1[1][0] == '1':
            ids[id1[0]] = line
        elif len(id1) > 1 and id1[1][0] == '2':
            is_paired = True

    if not is_paired:
        print("No paired reads found.")
        return max_rl, 0.0, 0.0

    hist = [0 for i in range(0, max_is + 1)]
    max_rl = 0

    logfile = open(log, "r")
    for line in logfile:
        id2 = line.split('/', 1)
        if len(id2) > 1 and id2[1][0] == '2' and id2[0] in ids:
            line1 = ids[id2[0]]
            read1 = line1.split('\t', 5)
            read2 = line.split('\t', 5)
            pos1 = int(read1[3])
            pos2 = int(read2[3])
            len1 = len(read1[4])
            len2 = len(read2[4])

            max_rl = max(max_rl, len1, len2)
            if (read1[1] == "+" and read2[1] == "-" and pos1 < pos2) or (read1[1] == "-" and read2[1] == "+" and pos1 > pos2):
                dist = abs(pos2 - pos1)
                if dist <= max_is:
                    hist[dist] += 1

    logfile.close()

    ttl = 0
    cnt = 0
    for i in range(0, max_is + 1):
        ttl += i * hist[i]
        cnt += hist[i]

    if cnt == 0:
        print("No paired reads found.")
        return max_rl, 0.0, 0.0

    mean = float(ttl) / float(cnt)

    dev = 0.0
    for i in range(0, max_is + 1):
       	dev += (i - mean) * (i - mean) * hist[i]

    dev = math.sqrt(dev / cnt)

    return max_rl, mean, dev

