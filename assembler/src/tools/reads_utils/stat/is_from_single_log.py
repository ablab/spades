#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import math


def round_to_zero(t):
    res = round(abs(t))
    if (t < 0.0):
        res = -res
    return res


def get_median(hist):
    S = 0.0
    for v in hist.values():
        S += v
    
    summ = S
    for k in sorted(hist.keys()):
        summ -= hist[k]
        if (summ <= S / 2):
            return k

    raise Exception    
    return 0


def get_mad(hist, median): # median absolute deviation
    hist2 = {}
    for k,v in hist.items():
        x = abs(k - round_to_zero(median))
        hist2[x] = v
    return get_median(hist2)


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
        if (len(self.hist.keys()) == 0):
            self.mean = 0.0
            self.dev = 0.0
            return self.mean, self.dev

        median = get_median(self.hist)
        mad = get_mad(self.hist, median)
        low = median - 5. * 1.4826 * mad
        high = median + 5. * 1.4826 * mad

        n = 0.0
        summ = 0.0
        sum2 = 0.0
        for k,v in self.hist.items():
            if (k < low or k > high):
                continue
            n += v
            summ += v * k  
            sum2 += v * k * k

        mean = summ / n
        delta = math.sqrt(sum2 / n - mean * mean)

        low = mean - 5 * delta
        high = mean + 5 * delta

        n = 0.0
        summ = 0.0
        sum2 = 0.0
        for k,v in self.hist.items():
            if (k < low or k > high):
               continue
            n += v
            summ += v * k
            sum2 += v * k * k

        self.mean = summ / n
        self.dev = math.sqrt(sum2 / n - self.mean * self.mean)
        
        return self.mean, self.dev


    def write_hist(self, filename):
        outf = open(filename, "w")
        for i in sorted(self.hist.keys()):
            outf.write(str(i) + " " + str(self.hist[i]) + "\n")            
        outf.close()


def stat_from_log(log, max_is= 1000000000):
    logfile = open(log, "r")

    stat = {"FR" : PairedStat("FR"), "RF" : PairedStat("RF"), "FF" : PairedStat("FF"), "AU" : PairedStat("AU"), "UU" : PairedStat("UU"),  "RL" : PairedStat("RL"), "SP" : PairedStat("SP")}
    ids = {}

    for line in logfile:
        id1 = line.split('/', 1)
        if len(id1) > 1:
            if id1[1][0] == '1':
                ids[id1[0]] = line
        else:
            id1 = line.split()
            ids[id1[0]] = line

    max_rl = 0
    logfile = open(log, "r")
    for line in logfile:
        id2 = " "
        l2 = line.split('/', 1)
        if len(l2) > 1:
            if l2[1][0] == '2':
                id2 = l2[0]
        else:
            id2 = line.split()[0]
 

        if id2 in ids:
            line1 = ids[id2]
            read1 = line1.split('\t')
            read2 = line.split('\t')
            if len(read1) < 10 or len(read2) < 10:
                print(read1, read2)
                continue

            pos1 = int(read1[3])
            pos2 = int(read2[3])
            len1 = len(read1[9])
            len2 = len(read2[9])
            str1 = int(read1[1]) & 16 == 0
            str2 = int(read2[1]) & 16 == 0
            al1 = int(read1[1]) & 4 == 0
            al2 = int(read2[1]) & 4 == 0

            max_rl = max(max_rl, len1, len2)
            inss = 0
            if pos2 > pos1:
                inss = max(pos2 + len2 - pos1, len1)
            else:
                inss = max(pos1 + len1 - pos2, len2)
             
            if not al1 and not al2:
                stat["UU"].inc(0)

            elif (al1 and not al2) or (al2 and not al1):
                stat["AU"].inc(0)

            elif ((str1 and not str2 and pos1 < pos2) or (str2 and not str1 and pos1 > pos2)) and inss <= max_is:
                stat["FR"].inc(inss)

            elif ((str2 and not str1 and pos1 < pos2) or (str1 and not str2 and pos1 > pos2)) and inss <= max_is:
                stat["RF"].inc(inss)

            elif str2 == str1 and inss <= max_is:
                stat["FF"].inc(inss)

            elif inss > max_is:
                stat["SP"].inc(inss)

            elif pos1 == pos2:
                stat["RL"].inc(inss)

            else:
                print(al1,al2,str1,str2,pos1,pos2,inss)
                    

            del ids[ id2 ]

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

