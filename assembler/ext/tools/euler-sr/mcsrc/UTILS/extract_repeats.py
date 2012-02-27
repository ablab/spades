#!/usr/bin/env python

import sys
import string
import re
import MySQLdb
import os
from cStringIO import StringIO


import mysqlutils
import sequtils

##############
#  globals
# 

q = "'"

#   SW  perc perc perc  query                                     position in query               matching       repeat         position in  repeat
#score  div. del. ins.  sequence                                     begin      end     (left)   repeat         class/family    begin  end (left)    ID
#  308  30.7  0.0  0.0  gi|34852147|ref|NW_047596.1|Rn20_2118          162      236  (3341919) C  ORR1E          LTR/MaLR        (284)   75      1     1
#   0      1    2    3                                      4            5        6          7 8      9                10           11   12     13    14


def count_lines(filename):
    file = open(filename, 'r')
    l = file.readline()
    numLines = 0
    while (l != ''):
        numLines += 1
        l = file.readline()
    file.close()

    return numLines
    
def count_sequences(filename):
    file = open(filename, 'r')
    l = file.readline()
    numSequences = 0
    while (l != ''):
        if ((len(l) > 0) and (l[0] == '>')):
            numSequences += 1
        l = file.readline()
    file.close()

    return numSequences

def get_repmask_repeat(infile):
    line = infile.readline().strip()
    values = line.split()
#    print 'got values: ', values
    return values


def posToStrandPos(pos):
    if (pos.find('(') >= 0):
        repPosRe = re.compile(r"\((-?\d+)\)")
        repPosMatch = repPosRe.search(pos)
        if(repPosMatch== None):
            print 'error parsing position: ', pos
            sys.exit(0)
        pos = int(repPosMatch.group(1))
        strand = 'C'
    else:
        pos =  int(pos)
        strand = '+'
    return pos, strand
     
def get_ref(query_seq):
    orderRefRe = re.compile(r".*ref\|([^.]+)\.*")
    refMatch = orderRefRe.search(query_seq)
    if (refMatch != None):
        return refMatch.group(1)
    else:
        return None


def validate_repeat(repeat_str, repmask_repeat):
    # steps to validate a repeat:
    # 1.  Make sure the length of the repeat is the same as the length of the repmask repeat.
    l = len(repeat_str)
    d = int(repmask_repeat[6]) - int(repmask_repeat[5]) + 1
    if (l == d):
        return 1
    else:
        return 0

def print_seq(title, seq, outfile):
    outfile.write('>' + title + '\n')
    while (seq != ''):
        substr = seq[0:59]
        outfile.write(substr + '\n')
        seq = seq[60:]

def read_past_title(infile):
    line = infile.readline().strip()
    line2 = ''
    if (line[0] == '>'):
        line2 = infile.readline().strip()
        # always return the sequence in the first part of the tuple
        return line2, line
    else:
        # no title
        return line,  ''
            
def get_next_repeat(curline,  curindex, infile):
    repeatStarted = 0
    # advance file to first repeat.
    repeat = ''
    title = ''
    while (curline != '' and repeatStarted == 0 ):
        # look through whole, or part (curindex > 0) of first line
        for i in range(curindex, len(curline)):
            c = curline[i]
            if (c == 'a' or c == 'c' or c == 't' or c == 'g'):
                repeatStarted = 1
                break
        if (repeatStarted == 0):
            curindex = 0
            curline, title = read_past_title(infile)

    # extend for the rest of the repeat
    curindex = i
    if (repeatStarted):
        repeatEnded = 0
        while(curline != '' and repeatEnded  == 0):
            for i in range(curindex, len(curline)):
                c = curline[i]
                if (c == 'a' or c == 'c' or c == 't' or c == 'g'):
                    repeat = repeat + c
                else:
                    repeatEnded = 1
                    break

            # still looking for the end of the repeat past this line?
            if (repeatEnded == 0):
                curindex = 0
                curline, title = read_past_title(infile)

    return repeat, curline, i, title
            

def GetRepeatNameId(cursor, repmask_entry):
    rep_name = repmask_entry[9]
    select_command ='select * from reference_repeats  where name= "'+rep_name+'"'
    cursor.execute(select_command)
    res = cursor.fetchone()
#    print 'got result ', res
    if (res != None):
        return res[0]
    else:
        print 'did not find ' + rep_name
        select_command = 'select * from reference_repeats where name="NOTFOUND"'
        cursor.execute(select_command)
        res = cursor.fetchone()
        if (res == None):
            print 'Internal error with database, should have "NOTFOUND" repeat name'
            sys.exit(0)
            return None
        return res[0]
    

if (len(sys.argv) < 3):
    print 'usage: extract_repeats  genome repmask_out_file'
    sys.exit(0)


if (os.environ.has_key('SALK_DB') == False):
    print 'Must set environment variable "SALK_DB" to the database to use'
    sys.exit(1)

database = os.environ['SALK_DB']

connection = MySQLdb.connect(db=database)
cursor = connection.cursor()

total_repeats = count_lines(sys.argv[2])    
infile  = open(sys.argv[1], 'r')
repmaskFile = open(sys.argv[2], 'r')



# extract the first contig
title = infile.readline()
contig, nextTitle = get_contig(infile)
contig_ref = get_ref(title)

# get the first repeat.
repmask_repeat = get_repmask_repeat(repmaskFile)
repmask_ref    = get_ref(repmask_repeat[4])
numRepeats = 1


while (title != '' ):
    #    repeat_seq, line, index, temp_title = get_next_repeat(line, index, infile)
    while (repmask_ref == contig_ref):
        numRepeats += 1
        repeat = extract_seq(contig, int(repmask_repeat[5]), int(repmask_repeat[6]))
        valid = validate_repeat(repeat, repmask_repeat)
        name_id = GetRepeatNameId(cursor, repmask_repeat)
        locus_id = GetLocusID(cursor, repmask_ref)

        repeat_values = []
        repeat_values.append(1)  # 0, hard-wire rat species id     
        repeat_values.append(locus_id) #1           
        repeat_values.append(name_id)           # 2 name of repeat
        repeat_values.append(q+repeat+q)        # 3 repeat
        repeat_values.append(repmask_repeat[5]) # 4 genomic location of repeat start
        repeat_values.append(repmask_repeat[6]) # 5 genomic location of repeat end
        repeat_values.append(q+repmask_repeat[8]+q)   # 6 strand on genome

        pos, strand = posToStrandPos(repmask_repeat[11])
        repeat_values.append(pos)               # 7 start pos on repeat
        repeat_values.append(repmask_repeat[12]) # 8 pos of repeat
        repeat_values.append(q+strand+q)            # 9 strand of repeat
        repeat_values.append(repmask_repeat[1])  # 10 percent divergent
        repeat_values.append(repmask_repeat[2])  # 11 percent del
        repeat_values.append(repmask_repeat[3])  # 12 percent ins
        if (len(repmask_repeat) == 16): 
            repeat_values.append('1')
        else:
            repeat_values.append('0')
        
        DBStore(cursor,'genome_repeats', repeat_values )

        # move to the next repeat.
        repmask_repeat = repmaskFile.readline().strip().split()
        #print 'got repeat: ', valid, ' ',  repeat        
        repmask_ref    = get_ref(repmask_repeat[4])
        print 'repeat #' + str(numRepeats) + ' of ' + str(total_repeats) + " gs: " +str( repeat_values[4] ) + " ge: " + str(repeat_values[5])
        
    # done extracting repeats for this contig.
    # The title of the next contig is already stored in 'nextTitle', use  that
    title = nextTitle
    contig_ref = get_ref(title)

    # Get the next contig in the file.
    contig, nextTitle = get_contig(infile)
    
    numRepeats = numRepeats + 1
    if (numRepeats == 30):
        print 'bye '
        sys.exit(0)

infile.close()
repmaskFile.close()
