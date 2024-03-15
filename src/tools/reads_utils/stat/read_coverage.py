#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate coverage from raw file

import sys
import math 

def coverage(in_filename, out_filename, maxLen, bar, kmer):
    
    inFile = open(in_filename)
    outFile = open(out_filename, 'w')
    
    hist = [0] * (maxLen + 1)

    for line in inFile:
        coords = line.split()
        stpos = int(coords[0])
        rlen =  int(coords[1])
        for i in range(0, rlen - kmer + 1):
            cpos = stpos + i
            if cpos <= maxLen:
                hist[cpos] += 1

    covered = 0.0
    for i in range(0, maxLen + 1):
        if (hist[i] > 0):
            covered += 1.0

    # print("Coverage: " + str(covered/maxLen) + "\n")

    newHist = [0 for i in range(int(math.floor((maxLen + 1) / bar)) + 2)]

    for i in range(maxLen + 1):
        newHist[int(i/bar)] += hist[i]

    for i in range(int(math.floor((maxLen + 1) / bar)) + 1):
        outFile.write(str(i) + ' ' + str(int(round(newHist[i] / bar))) + '\n')

    inFile.close()
    outFile.close()

    return covered/maxLen


def read_genome(filename):
    res_seq = []

    seq = ''
    for line in open(filename):
            if line[0] == '>':
                    pass
            else:
                    seq += line.strip()

    return seq


def write_fasta(data, filename):
    outFile = open(filename, 'w')

    for seq in data:
            outFile.write('>' + seq[0] + '\n');
            i = 0
            while i < len(seq[1]):
                    outFile.write(seq[1][i:i+60] + '\n')
                    i += 60

    outFile.close()


def analyze_gaps(in_filename, out_filename, reference, out_ref, kmer):
    inFile = open(in_filename)
    outFile = open(out_filename, 'w')

    gaps_stat = {101:0, 501:0, 1001:0}
    chunks_stat = {101:0, 501:0, 1001:0}
    current = 0

    #[start,end)
    gaps = []
    chunks = []

    line = inFile.readline()
    while line:
            cov = int(line.split()[1])
            end = current
            while cov == 0 and line:
                   end += 1 
                   line = inFile.readline()
                   if line:
                          cov = int(line.split()[1])

            if end != current:
                   gaps.append((current, end));
                   length = end - current
                   for key in gaps_stat:
                          if length < key:
                                 gaps_stat[key] += 1

            current = end
            while cov > 0 and line:
                   end += 1 
                   line = inFile.readline()
                   if line:
                          cov = int(line.split()[1])

            if end != current:
                   chunks.append((current, end + kmer - 1));
                   length = end - current
                   for key in chunks_stat:
                          if length < key:
                                 chunks_stat[key] += 1

            current = end


    outFile.write("Total chunks: " + str(len(chunks)) + "\n")
    for key in chunks_stat:
            outFile.write("Chunks < " + str(key) + ": " + str(chunks_stat[key])+ "\n")

    outFile.write("Total gaps: " + str(len(gaps)) + "\n")
    for key in gaps_stat:
            outFile.write("Gaps < " + str(key) + ": " + str(gaps_stat[key])+ "\n")

    outFile.write("Gaps:\n")
    for gap in gaps:
            outFile.write(str(gap[0]) + ' ' + str(gap[1]) + "\n")
    
    outFile.close()    

#    genome = read_genome(reference, chunks)
    genome = read_genome(reference)
    ref_chunks = []
    i = 0
    for chunk in chunks:
            ref_chunks.append(("PART_" + str(i) + "_from_" + str(chunk[0]) + "_to_" + str(chunk[1]), genome[chunk[0]:chunk[1]]))
            i += 1
    write_fasta(ref_chunks, out_ref)

    return chunks, gaps

def main():

    if len(sys.argv) < 5:
        print("Usage: <coverage file> <output> <genome length> <bar width> [k = 1]");
        exit(0)

    kmer = 1
    if len(sys.argv) > 5:
        kmer = int(sys.argv[5])

    cov = coverage(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), kmer)

    print("Coverage: " + str(cov) + "\n")

if __name__ == '__main__':
    main()
