#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import math


class PairedStat:
    def __init__(self, name):
        self.name = name
        self.count = 0
        self.hist = {}

    def inc(self, distance):
        self.count += 1
        if distance in self.hist:
            self.hist[distance] += 1
        else:
            self.hist[distance] = 1

        return True


    def make_stat(self):
        if self.count == 0:
            self.mean = 0.0
            self.dev = 0.0
            return self.mean, self.dev

        ttl = 0
        for i in self.hist.iterkeys():
            ttl += i * self.hist[i]
        self.mean = float(ttl) / float(self.count)

        self.dev = 0.0
        for i in self.hist.iterkeys():
           	self.dev += (i - self.mean) * (i - self.mean) * self.hist[i]
        self.dev = math.sqrt(self.dev / self.count)

        return self.mean, self.dev

    def write_hist(self, filename):
        outf = open(filename, "w")
        for i in sorted(self.hist.keys()):
            outf.write(str(i) + " " + str(self.hist[i]) + "\n")            
        outf.close()


def stat_from_log(log, max_is= 1000000000):
    logfile = open(log, "r")

    stat = {"FR" : PairedStat("FR"), "RF" : PairedStat("RF"), "FF" : PairedStat("FF"), "AU" : PairedStat("AU"), "SP" : PairedStat("SP")}
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

    max_rl = 0
    logfile = open(log, "r")
    for line in logfile:
        id2 = line.split('/', 1)
        if len(id2) > 1 and id2[1][0] == '2':
            if id2[0] in ids:
                line1 = ids[id2[0]]
                read1 = line1.split('\t', 5)
                read2 = line.split('\t', 5)
                pos1 = int(read1[3])
                pos2 = int(read2[3])
                len1 = len(read1[4])
                len2 = len(read2[4])

                max_rl = max(max_rl, len1, len2)
                inss = 0
                if pos2 > pos1:
                    inss = max(pos2 + len2 - pos1, len1)
                else:
                    inss = max(pos1 + len1 - pos2, len2)

                if ((read1[1] == "+" and read2[1] == "-" and pos1 < pos2) or (read1[1] == "-" and read2[1] == "+" and pos1 > pos2)) and inss <= max_is:

                    stat["FR"].inc(inss)

                elif ((read1[1] == "-" and read2[1] == "+" and pos1 < pos2) or (read1[1] == "+" and read2[1] == "-" and pos1 > pos2)) and inss <= max_is:
                    stat["RF"].inc(inss)

                elif read1[1] == read2[1] and inss <= max_is:
                    stat["FF"].inc(inss)

                elif inss > max_is:
                    stat["SP"].inc(inss)
                    

                del ids[ id2[0] ]

            else:
                stat["AU"].inc(0)

    logfile.close()

    for id1 in ids.iterkeys():
        stat["AU"].inc(0)

    stat["FR"].make_stat()
    stat["RF"].make_stat()
    stat["FF"].make_stat()

    return max_rl, stat


def dist_from_log(log, max_is):
    stat = stat_from_log(log, max_is)
    return stat[0], stat[1]["FR"].mean - stat[0], stat[1]["FR"].dev

