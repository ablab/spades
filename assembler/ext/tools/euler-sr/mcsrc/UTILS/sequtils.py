#!/usr/bin/env python

import sys
import getopt
import os
import string
from cStringIO import StringIO

def print_seq(outfile, seq, title, length):
    outfile.write(title)
    while (seq != ''):
        substr = seq[0:(length-1)]
        outfile.write(substr+'\n')
        seq = seq[length:]

def get_contig(infile):
    contig = ''
    file_str = StringIO()
    line = string.strip(infile.readline())
    count = 1
    while (line != '' and line[0] != '>'):
        #        contig = contig + line
        file_str.write(line)
        line = string.strip(infile.readline())
        count = count + 1
    return file_str.getvalue(), line


def extract_seq(source, start, end):
    return source[start:end]

def extract_from_fastarec(title, infile, outfile):
    repstart = -1
    repend   = -1
    line = string.strip(infile.readline())
    repeat = ''
    repNum = 0
    while (line != '' and line[0] != '>'):
        for i in range(0,len(line)):
            c = line[i]
            if (c == 'a' or c == 'c' or c == 't' or c == 'g'):
                repeat = repeat + c
            else:
                if (repeat != ''):
                    print_seq(title + ' ' + str(repNum), repeat, outfile)
                    repeat = ''
                    repNum = repNum + 1
        line = string.strip(infile.readline())
    title = line
    return title

