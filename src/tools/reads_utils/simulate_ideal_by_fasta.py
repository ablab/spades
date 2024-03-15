#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys

rc_map = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}

def revcomp(string):
    rc = "".join([rc_map[x] for x in string.strip()])
    return rc[::-1]


def trusted_gaps(log, gaps, last_gap, max_gap):
    logfile = open(log, "r")
    ids = {}

    trusted = set([])
    possible_gaps = set([])
    for gap in gaps:
        if gap[1] - gap[0] <= max_gap:
            possible_gaps.add(gap)

    
    is_paired = False
    for line in logfile:
        id1 = line.split('/', 1)
        if len(id1) > 1 and id1[1][0] == '1':
            ids[id1[0]] = line
        elif len(id1) > 1 and id1[1][0] == '2':
            is_paired = True

    if not is_paired:
        print("No paired reads found.")
        return trusted

    logfile = open(log, "r")
    for line in logfile:
        id2 = line.split('/', 1)
        if len(id2) > 1 and id2[1][0] == '2' and id2[0] in ids:
            line1 = ids[id2[0]]
            read1 = line1.split('\t', 5)
            read2 = line.split('\t', 5)
            pos1 = int(read1[3])
            pos2 = int(read2[3])
            
            if (read1[1] == "+" and read2[1] == "-" and pos1 < pos2) or (read1[1] == "-" and read2[1] == "+" and pos1 > pos2):
                for gap in possible_gaps:
                    if min(pos1, pos2) <= gap[0] and max(pos1, pos2) >= gap[1]:
                        trusted.add(gap)
                        possible_gaps.remove(gap)
                        break

            if ((read1[1] == "+" and read2[1] == "-" and pos1 > pos2) or (read1[1] == "-" and read2[1] == "+" and pos1 < pos2)) and min(pos1, pos2) >= last_gap[1] and max(pos1, pos2) <= last_gap[0]:
                trusted.add(last_gap)

    logfile.close()
    return trusted


def print_single(outfile, contig, rl, postfix):
    for j in range(0, len(contig) - rl + 1):
        outfile.write(">READ_" + str(j) + postfix + "\n")
        outfile.write(contig[j : j + rl] + "\n")


def print_paired(outfile, contig, rl, distance, postfix):
    for j in range(0, len(contig) - rl - distance + 1):
        outfile.write(">READ_" + str(j) + postfix + "/1\n")
        outfile.write(contig[j : j + rl] + "\n")
        outfile.write(">READ_" + str(j) + postfix + "/2\n")
        outfile.write(revcomp(contig[j + distance : j + distance + rl]) + "\n")


def simulate_single(filename, fasta, rl = 100, circular = False):
    outfile = open(filename, 'w')

    for i in range(0, len(fasta)):
        contig = fasta[i][1]
        if circular and i == len(fasta) - 1:
            contig = fasta[i][1] + fasta[0][1][:rl - 1]

        print_single(outfile, contig, min(rl, len(contig)), "_CONTIG_" + str(i))

    outfile.close()
            

def simulate_paired(filename, fasta, distance, rl = 100, circular = False):
    outfile = open(filename, 'w')

    for i in range(0, len(fasta)):
        contig = fasta[i][1]
        if circular and i == len(fasta) - 1:
            contig = fasta[i][1] + fasta[0][1][:distance + rl - 1]

        print_paired(outfile, contig, rl, distance, "_CONTIG_" + str(i))

    outfile.close()


def simulate_paired_over_gaps(filename, fasta, chunks, trusted_gaps, last_gap, ref_len, distance, rl = 100, circular = False):
    outfile = open(filename, 'w')

    for i in range(0, len(fasta)):
        if trusted_gaps[i] != (0,0) and i < len(fasta) - 1 and trusted_gaps[i][1] - trusted_gaps[i][0] <= distance - rl:
            if rl > len(fasta[i][1]) or rl > len(fasta[i + 1][1]):
                outfile.write(">READ_0" + str() + "_SGAP_" + str(i) + "/1\n")
                outfile.write(fasta[i][1][-rl:] + "\n")
                outfile.write(">READ_0" + str() + "_SGAP_" + str(i) + "/2\n")
                outfile.write(revcomp(fasta[i + 1][1][:rl]) + "\n")
                continue

            lpos = trusted_gaps[i][1] - distance - chunks[i][0]
            if lpos < 0:
                lpos = 0

            rpos = trusted_gaps[i][0] + distance
            if rpos >= chunks[i + 1][1]:
                rpos = chunks[i + 1][1]
            rpos -= chunks[i + 1][0]

            gap_contig = fasta[i][1][lpos:] + ("N" * (trusted_gaps[i][1] - trusted_gaps[i][0])) + fasta[i + 1][1][:rpos]

            print_paired(outfile, gap_contig, rl, distance, "_GAP_" + str(i))
            

    if last_gap != (0,0):
        if rl > len(fasta[0][1]) or rl > len(fasta[-1][1]):
            outfile.write(">READ_0" + str() + "_LSGAP_" + str(i) + "/1\n")
            outfile.write(fasta[-1][1][-rl:] + "\n")
            outfile.write(">READ_0" + str() + "_LSGAP_" + str(i) + "/2\n")
            outfile.write(revcomp(fasta[0][1][:rl]) + "\n")
        else:
            lpos = last_gap[1] - distance + ref_len - chunks[-1][0]
            if lpos < 0:
               lpos = 0

            rpos = last_gap[0] + distance - ref_len - chunks[0][0]
            if rpos >= chunks[0][1]:
                rpos = chunks[0][1]
            rpos -= chunks[0][0]

            gap_contig = fasta[-1][1][lpos:] + ("N" * (last_gap[1] + ref_len - last_gap[0] )) + fasta[0][1][:rpos]
            print_paired(outfile, gap_contig, rl, distance, "_LAST_GAP")

    outfile.close()


