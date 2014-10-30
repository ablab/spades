# based on IlluminaNxSeqJunction-Split7.py, IlluminaChimera-Clean4.py and ParseFastQ.py
# all copyrights are below!

# by Scott Monsma, Copyright (c) Lucigen Corp July 2014 - based on NxSeqFOS-SplitBfa4.py
# Splits 'mates_ICC4_' files into left and right insert sequences by finding the Junction Code(s)
# usage: copy IlluminaNxSeqJunction-Split6.py and ParseFastq.py into a directory with your fastq files to process
#cd into directory with .py and .fastq
#make sure your read 1 filename contains '_R1_' and read 2 filename contains '_R2_'
#at command prompt type 'python IlluminaNxSeqJunction-Split7.py 'mates_ICC4_your-R1-filename.fastq' and hit enter
#split sequences are saved if longer than minseq
#output files are named 'R1_IJS7_mates_ICC4_your-R1-filename.fastq' and 'R2_IJS7_mates_ICC4_your-R2-filename.fastq' which are the trimmed mate pairs, and
#'unsplit_IJS7_yourfilename.fastq' which contains interleaved reads where no junction was found.

#IlluminaChimera-Clean4 by Scott Monsma, Lucigen Corp Copyright (C) July 2014
# usage: copy IlluminaChimera-Clean4.py and ParseFastq.py into a directory with your fastq file to process
#cd into directory with .py and .fastq
#at command prompt type 'python IlluminaChimera-Clean4.py yourfilename.fastq' and hit enter
#four new files will be created, 'mates_ICC4_your-R1-filename.fastq' and 'mates_ICC4_your-R2-filename.fastq' containing the
#true mate pairs with matching chimera codes, and 'non-mates_ICC4_your-R1-filename.fastq' and 'non-mates_ICC4_your-R2-filename.fastq'
#containing the chimera read pairs and unidentified read pairs


import os
import time
import string
import support

try:
    import regex
except ImportError:
    support.error("Can't process Lucigen NxMate reads! Python module regex is not installed!")

minseq = 25  #minimum length sequence to keep after trimming


class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""

    def __init__(self, filePath, headerSymbols=['@', '+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        #  if filePath.endswith('.gz'):
        #     self._file = gzip.open(filePath)
        # else:
        self._file = open(filePath, 'rU')  #filePath, 'rU') test with explicit filename
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self

    def next(self):  # for both Python2 and Python3
        return self.__next__()

    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1  ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else:
                elemList.append(None)

        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4, \
            "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]), \
            "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._hdSyms[0], self._currentLineNumber)
        assert elemList[2].startswith(self._hdSyms[1]), \
            "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._hdSyms[1], self._currentLineNumber)
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
        self._currentLineNumber)

        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)


def chimera_clean(infilename1, infilename2, dst, log, silent=True):
    starttime = time.clock()

    #open four  outfiles
    outfilename1 = os.path.join(dst, 'mates_ICC4_' + os.path.basename(infilename1))
    outfile1 = open(outfilename1, 'w')

    slagfilename1 = os.path.join(dst, 'non-mates_ICC4_' + os.path.basename(infilename1))
    slagfile1 = open(slagfilename1, 'w')

    outfilename2 = os.path.join(dst, 'mates_ICC4_' + os.path.basename(infilename2))
    outfile2 = open(outfilename2, 'w')

    slagfilename2 = os.path.join(dst, 'non-mates_ICC4_' + os.path.basename(infilename2))
    slagfile2 = open(slagfilename2, 'w')

    #set up regular expression patterns for chimera codes- for illumin use the  reverse complements of right codes
    csslist1 = ['(TGGACTCCACTGTG){e<=1}', '(ACTTCGCCACTGTG){e<=1}', '(TGAGTCCCACTGTG){e<=1}', '(TGACTGCCACTGTG){e<=1}',
                '(TCAGGTCCACTGTG){e<=1}', '(ATGTCACCACTGTG){e<=1}', '(GTATGACCACTGTG){e<=1}', '(GTCTACCCACTGTG){e<=1}',
                '(GTTGGACCACTGTG){e<=1}', '(CGATTCCCACTGTG){e<=1}', '(GGTTACCCACTGTG){e<=1}', '(TCACCTCCACTGTG){e<=1}']

    csslist2 = ['(TCCAGACCAATGTG){e<=1}', '(ACATCACCAATGTG){e<=1}', '(TCACGACCAATGTG){e<=1}', '(TAGCACCCAATGTG){e<=1}',
                '(AACCTCCCAATGTG){e<=1}', '(ACAACTCCAATGTG){e<=1}', '(GTCTAACCAATGTG){e<=1}', '(TACACGCCAATGTG){e<=1}',
                '(GAGAACCCAATGTG){e<=1}', '(GAGATTCCAATGTG){e<=1}', '(GACCTACCAATGTG){e<=1}', '(AGACTCCCAATGTG){e<=1}']

    csscounter = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  #12 slots for the correct code combinations

    #PARSE both files in tuples of 4 lines
    parserR1 = ParseFastQ(infilename1)
    parserR2 = ParseFastQ(infilename2)

    #init counters
    readcounter = 0
    matecounter = 0  #for pairs where both trimmed reads are equal to or longer than minseq
    TOTALmatecounter = 0
    slagcounter = 0

    #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
    for recR1 in parserR1:
        readcounter += 1
        recR2 = parserR2.__next__()

        #check if rec.seqStr contains match to chimera pattern
        for cssindex, css1 in enumerate(csslist1):
            m = regex.search(css1, recR1[1])
            css2 = csslist2[cssindex]
            n = regex.search(css2, recR2[1])

            if m and n:  # a true mate pair! write out to mates files
                TOTALmatecounter += 1
                #NOTE  TAKE THIS OPPORTUNITY TO RECORD CSS CODE AND TRUNCATE READS
                #need to trim additional 9+4 nts from end of match to remove css, Bst, barcode (9) and CGAT (4) linker
                csscounter[cssindex] += 1  #increment the appropriate css counter
                R1matches = m.span()
                mend = R1matches[1]
                mend = mend + 13
                mySeq = recR1[1]
                myR1 = mySeq[mend:]  #trim the left end off of Read1
                myQual1 = recR1[3]
                myR1Qual = myQual1[mend:]  #trim the left end off of Read1 quality string

                R2matches = n.span()
                nend = R2matches[1]
                nend = nend + 13
                mySeq2 = recR2[1]
                myR2 = mySeq2[nend:]  #trim the left end off of Read2
                myQual2 = recR2[3]
                myR2Qual = myQual2[nend:]  #trim the left end off of Read2 quality string

                if (len(myR1) >= minseq) and (len(myR2) >= minseq):  #and if one or other is too short, toss both
                    matecounter += 1
                    outfile1.write(recR1[0])
                    outfile1.write('\n')
                    outfile1.write(myR1)
                    outfile1.write('\n')
                    outfile1.write(recR1[2])
                    outfile1.write('\n')
                    outfile1.write(myR1Qual)
                    outfile1.write('\n')

                    outfile2.write(recR2[0])
                    outfile2.write('\n')
                    outfile2.write(myR2)
                    outfile2.write('\n')
                    outfile2.write(recR2[2])
                    outfile2.write('\n')
                    outfile2.write(myR2Qual)
                    outfile2.write('\n')
                break  # found it, go on to next rec
            else:  # no chimera code in R1 so can't be a mate pair; write out to slag files if this is the last one
                if cssindex == 11:
                    slagcounter += 1
                    slagfile1.write(recR1[0])
                    slagfile1.write('\n')
                    slagfile1.write(recR1[1])
                    slagfile1.write('\n')
                    slagfile1.write(recR1[2])
                    slagfile1.write('\n')
                    slagfile1.write(recR1[3])
                    slagfile1.write('\n')

                    slagfile2.write(recR2[0])
                    slagfile2.write('\n')
                    slagfile2.write(recR2[1])
                    slagfile2.write('\n')
                    slagfile2.write(recR2[2])
                    slagfile2.write('\n')
                    slagfile2.write(recR2[3])
                    slagfile2.write('\n')
    #all done
    outfile1.close()
    slagfile1.close()
    outfile2.close()
    slagfile2.close()
    if TOTALmatecounter+slagcounter != readcounter:
        support.error("lucigen_nxmate.py, chimera_clean: error in the script somewhere! unequal read counts!", log)
    if not silent:
        #print some stats
        percentmates = 100. * matecounter / readcounter
        percentslag = 100. * slagcounter / readcounter
        log.info("==== chimera_clean info: %d reads processed, %d true mate reads (%.2f %%) "
                 "and %d non-mates/chimeras (%.2f %%)."
                 % (readcounter, matecounter, percentmates, slagcounter, percentslag))
        shortmates = TOTALmatecounter - matecounter
        log.info("==== chimera_clean info: %d mates too short to keep after trimming" % shortmates)
        elapsedtime = time.clock() - starttime
        log.info("==== chimera_clean info: processed %d reads in %.1f seconds" % (readcounter, elapsedtime))
        log.info("==== chimera_clean info: " + str(csscounter))
    return outfilename1, outfilename2


def nx_seq_junstion(infilename1, infilename2, dst, log, silent=True):
    starttime = time.clock()

    #open three  outfiles
    splitfilenameleft = os.path.join(dst, 'R1_IJS7_' + os.path.basename(infilename1))
    splitfile1 = open(splitfilenameleft, 'w')

    splitfilenameright = os.path.join(dst, 'R2_IJS7_' + os.path.basename(infilename2))
    splitfile2 = open(splitfilenameright, 'w')

    unsplitfilename = os.path.join(dst, 'unsplit_IJS7_' + os.path.basename(infilename1).replace('_R1_', '_R1R2_'))
    unsplitfile = open(unsplitfilename, 'w')

    #jctstr = '(GGTTCATCGTCAGGCCTGACGATGAACC){e<=4}' # JS7 24/28 required results in ~92% detected in ion torrent
    # from NextClip: --adaptor_sequence GTTCATCGTCAGG -e --strict_match 22,11 --relaxed_match 20,10 eg strict 22/26 = 4 errors, relaxed 20/26 = 6 errors
    jctstr = '(GTTCATCGTCAGGCCTGACGATGAAC){e<=4}'  # try 22/26 to match NextClip strict (e<=6 for relaxed)

    #PARSE both files in tuples of 4 lines, don't have to explicitly open them
    parserR1 = ParseFastQ(infilename1)
    parserR2 = ParseFastQ(infilename2)

    #init counters
    readcounter = 0
    jctcounter = 0
    splitcounter = 0
    bothjctcounter = 0
    R1jctcounter = 0
    R2jctcounter = 0
    R1R2jctcounter = 0

    #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
    for recR1 in parserR1:
        readcounter += 1
        recR2 = parserR2.__next__()

        m = regex.search(jctstr, recR1[1])
        n = regex.search(jctstr, recR2[1])
        if m and n:  #found jctstr in both reads; need to save left part of R1 and LEFT part of R2
            bothjctcounter += 1
            matches = m.span()
            start = matches[0]
            mySeq = recR1[1]  #get the left part of Read1
            myLeft = mySeq[:start]
            myQual = recR1[3]
            myLeftQual = myQual[:start]  #get left part of Read1 quality string

            nmatches = n.span()
            nstart = nmatches[0]
            nend = nmatches[1]
            mySeq2 = recR2[1]

            myRight2 = mySeq2[:nstart]  #get left part of Read2
            myQual2 = recR2[3]
            myRightQual2 = myQual2[:nstart]  #get left part of Read2 quality string

            #only write out as split if both pieces are big enough
            if (len(myLeft) > minseq) and (len(myRight2) > minseq):
                splitcounter += 1
                R1R2jctcounter += 1
                splitfile1.write(recR1[0])
                splitfile1.write('\n')
                #write left sequence
                splitfile1.write(myLeft)
                splitfile1.write('\n')
                splitfile1.write(recR1[2])
                splitfile1.write('\n')
                #write quality
                splitfile1.write(myLeftQual)
                splitfile1.write('\n')

                splitfile2.write(recR2[0])
                splitfile2.write('\n')
                #write right sequence
                splitfile2.write(myRight2)
                splitfile2.write('\n')
                splitfile2.write(recR2[2])
                splitfile2.write('\n')
                #write quality
                splitfile2.write(myRightQual2)
                splitfile2.write('\n')
        elif n:  #junction only in R2, so save entire R1 and LEFT part of R2 IFF R2 long enough
            nmatches = n.span()
            nstart = nmatches[0]
            nend = nmatches[1]
            mySeq2 = recR2[1]
            myRight2 = mySeq2[:nstart]
            myQual2 = recR2[3]
            myRightQual2 = myQual2[:nstart]
            if (len(myRight2) > minseq):
                splitcounter += 1

                splitfile2.write(recR2[0])
                splitfile2.write('\n')
                #write right sequence
                splitfile2.write(myRight2)
                splitfile2.write('\n')
                splitfile2.write(recR2[2])
                splitfile2.write('\n')
                #write quality
                splitfile2.write(myRightQual2)
                splitfile2.write('\n')

                splitfile1.write(recR1[0])
                splitfile1.write('\n')
                #write left sequence
                splitfile1.write(recR1[1])
                splitfile1.write('\n')
                splitfile1.write(recR1[2])
                splitfile1.write('\n')
                #write quality
                splitfile1.write(recR1[3])
                splitfile1.write('\n')

                jctcounter += 1
                R2jctcounter += 1
        elif m:  #junction only in R1, save left part of R1 and entire R2, IFF R1 is long enough
            matches = m.span()
            start = matches[0]
            mySeq = recR1[1]
            myLeft = mySeq[:start]
            myQual = recR1[3]
            myLeftQual = myQual[:start]
            if (len(myLeft) > minseq):
                splitcounter += 1
                splitfile1.write(recR1[0])
                splitfile1.write('\n')
                #write left sequence
                splitfile1.write(myLeft)
                splitfile1.write('\n')
                splitfile1.write(recR1[2])
                splitfile1.write('\n')
                #write quality
                splitfile1.write(myLeftQual)
                splitfile1.write('\n')

                splitfile2.write(recR2[0])
                splitfile2.write('\n')
                #write left sequence
                splitfile2.write(recR2[1])
                splitfile2.write('\n')
                splitfile2.write(recR2[2])
                splitfile2.write('\n')
                #write quality
                splitfile2.write(recR2[3])
                splitfile2.write('\n')

                jctcounter += 1
                R1jctcounter += 1
        else:  #no junctions, save for frag use, as is 'unsplit'; note this file will be interleaved R1 R2 R1 R2...
            unsplitfile.write(recR1[0])
            unsplitfile.write('\n')
            #write left sequence
            unsplitfile.write(recR1[1])
            unsplitfile.write('\n')
            unsplitfile.write(recR1[2])
            unsplitfile.write('\n')
            #write quality
            unsplitfile.write(recR1[3])
            unsplitfile.write('\n')

            unsplitfile.write(recR2[0])
            unsplitfile.write('\n')
            #write left sequence
            unsplitfile.write(recR2[1])
            unsplitfile.write('\n')
            unsplitfile.write(recR2[2])
            unsplitfile.write('\n')
            #write quality
            unsplitfile.write(recR2[3])
            unsplitfile.write('\n')
    #all done
    splitfile1.close()
    splitfile2.close()
    unsplitfile.close()

    if not silent:
        #print some stats
        percentsplit = 100 * splitcounter / readcounter
        percentR1R2 = 100 * R1R2jctcounter / splitcounter
        percentR1 = 100 * R1jctcounter / splitcounter
        percentR2 = 100 * R2jctcounter / splitcounter
        log.info("==== nx_seq_junstion info: %d reads processed" % (readcounter))
        log.info("==== nx_seq_junstion info: %d total split pairs (%.2f %% of processed reads))"
                 % (splitcounter, percentsplit))
        log.info("==== nx_seq_junstion info: %d junctions in both R1 and R2 (%.2f %% of split junctions))"
                 % (R1R2jctcounter, percentR1R2))
        log.info("==== nx_seq_junstion info: %d split junctions are in Read1 (%.2f %% of split junctions))"
                 % (R1jctcounter, percentR1))
        log.info("==== nx_seq_junstion info: %d split junctions are in Read2 (%.2f %% of split junctions))"
                 % (R2jctcounter, percentR2))
        elapsedtime = time.clock() - starttime
        log.info("==== nx_seq_junstion info: processed %d reads in %.1f seconds" % (readcounter, elapsedtime))
    return splitfilenameleft, splitfilenameright, unsplitfilename


def process_reads(left_reads_fpath, right_reads_fpath, dst, log):
    log.info("== Processing Lucigen NxMate reads (" + left_reads_fpath + " and " +
             os.path.basename(right_reads_fpath) + " (results are in " + dst + " directory)")
    cleaned_filename1, cleaned_filename2 = chimera_clean(left_reads_fpath, right_reads_fpath, dst, log, silent=True)
    split_filename1, split_filename2, unsplit_filename = nx_seq_junstion(cleaned_filename1, cleaned_filename2, dst, log, silent=True)
    return split_filename1, split_filename2, unsplit_filename
