############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

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
import support
import gzip
import itertools
import sys
from site import addsitedir
import spades_init
import options_storage

try:
    import regex
except ImportError:
    support.error("Can't process Lucigen NxMate reads! Python module regex is not installed!")

addsitedir(spades_init.ext_python_modules_home)
if sys.version.startswith('2.'):
    from joblib2 import Parallel, delayed
elif sys.version.startswith('3.'):
    from joblib3 import Parallel, delayed

    
# CONSTANTS
READS_PER_THREAD = 25000
READS_PER_BATCH = READS_PER_THREAD * options_storage.threads  # e.g. 100000 for 4 threads
minseq = 25  # minimum length sequence to keep after trimming


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
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
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

    def close(self):
        if self._file:
            self._file.close()


def write_to_files(file_handlers, record_lists):
    for file_handler, record_list in zip(file_handlers, record_lists):
        for record in record_list:
            for line in record:
                file_handler.write(line + '\n')


def split_into_chunks(l, n):
    avg = len(l) / float(n)
    out = []
    last = 0.0
    while last < len(l):
        out.append(l[int(last):int(last + avg)])
        last += avg
    return out


class CleanStats(object):
    def __init__(self):
        self.readcounter = 0
        self.matecounter = 0  # for pairs where both trimmed reads are equal to or longer than minseq
        self.TOTALmatecounter = 0
        self.slagcounter = 0
        self.csscounter = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 12 slots for the correct code combinations

    def __add__(self, other):
        self.readcounter += other.readcounter
        self.matecounter += other.matecounter
        self.TOTALmatecounter += other.TOTALmatecounter
        self.slagcounter += other.slagcounter
        self.csscounter = [x + y for x, y in zip(self.csscounter, other.csscounter)]
        return self


def chimera_clean_process_batch(reads, csslist1, csslist2):
    stats = CleanStats()
    processed_out1 = []
    processed_out2 = []
    processed_slag1 = []
    processed_slag2 = []

    #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
    for recR1, recR2 in reads:
        stats.readcounter += 1

        #check if rec.seqStr contains match to chimera pattern
        for cssindex, css1 in enumerate(csslist1):
            m = regex.search(css1, recR1[1])
            css2 = csslist2[cssindex]
            n = regex.search(css2, recR2[1])

            if m and n:  # a true mate pair! write out to mates files
                stats.TOTALmatecounter += 1
                #NOTE  TAKE THIS OPPORTUNITY TO RECORD CSS CODE AND TRUNCATE READS
                #need to trim additional 9+4 nts from end of match to remove css, Bst, barcode (9) and CGAT (4) linker
                stats.csscounter[cssindex] += 1  #increment the appropriate css counter
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
                    stats.matecounter += 1
                    processed_out1.append([recR1[0], myR1, recR1[2], myR1Qual])
                    processed_out2.append([recR2[0], myR2, recR2[2], myR2Qual])
                break  # found it, go on to next rec
            else:  # no chimera code in R1 so can't be a mate pair; write out to slag files if this is the last one
                if cssindex == 11:
                    stats.slagcounter += 1
                    processed_slag1.append([recR1[0], recR1[1], recR1[2], recR1[3]])
                    processed_slag2.append([recR2[0], recR2[1], recR2[2], recR2[3]])
    return [processed_out1, processed_out2, processed_slag1, processed_slag2], stats


def chimera_clean(infilename1, infilename2, dst, log, silent=True):
    starttime = time.time()

    basename1 = os.path.basename(infilename1)
    if os.path.splitext(basename1)[1] == '.gz':
        basename1 = os.path.splitext(basename1)[0]
    basename2 = os.path.basename(infilename2)
    if os.path.splitext(basename2)[1] == '.gz':
        basename2 = os.path.splitext(basename2)[0]
    #open four outfiles
    outfilename1 = os.path.join(dst, 'mates_ICC4_' + basename1)
    outfile1 = open(outfilename1, 'w')

    slagfilename1 = os.path.join(dst, 'non-mates_ICC4_' + basename1)
    slagfile1 = open(slagfilename1, 'w')

    outfilename2 = os.path.join(dst, 'mates_ICC4_' + basename2)
    outfile2 = open(outfilename2, 'w')

    slagfilename2 = os.path.join(dst, 'non-mates_ICC4_' + basename2)
    slagfile2 = open(slagfilename2, 'w')

    #set up regular expression patterns for chimera codes- for illumin use the  reverse complements of right codes
    csslist1 = ['(TGGACTCCACTGTG){e<=1}', '(ACTTCGCCACTGTG){e<=1}', '(TGAGTCCCACTGTG){e<=1}', '(TGACTGCCACTGTG){e<=1}',
                '(TCAGGTCCACTGTG){e<=1}', '(ATGTCACCACTGTG){e<=1}', '(GTATGACCACTGTG){e<=1}', '(GTCTACCCACTGTG){e<=1}',
                '(GTTGGACCACTGTG){e<=1}', '(CGATTCCCACTGTG){e<=1}', '(GGTTACCCACTGTG){e<=1}', '(TCACCTCCACTGTG){e<=1}']

    csslist2 = ['(TCCAGACCAATGTG){e<=1}', '(ACATCACCAATGTG){e<=1}', '(TCACGACCAATGTG){e<=1}', '(TAGCACCCAATGTG){e<=1}',
                '(AACCTCCCAATGTG){e<=1}', '(ACAACTCCAATGTG){e<=1}', '(GTCTAACCAATGTG){e<=1}', '(TACACGCCAATGTG){e<=1}',
                '(GAGAACCCAATGTG){e<=1}', '(GAGATTCCAATGTG){e<=1}', '(GACCTACCAATGTG){e<=1}', '(AGACTCCCAATGTG){e<=1}']

    #PARSE both files in tuples of 4 lines
    parserR1 = ParseFastQ(infilename1)
    parserR2 = ParseFastQ(infilename2)

    all_stats = CleanStats()
    n_jobs = options_storage.threads
    while True:
        # prepare input
        reads1 = list(itertools.islice(parserR1, READS_PER_BATCH))
        reads2 = list(itertools.islice(parserR2, READS_PER_BATCH))
        if len(reads1) != len(reads2):
            support.error("lucigen_nxmate.py, chimera_clean: "
                          "number of left reads (%d) is not equal to number of right reads (%d)!"
                          % (len(reads1), len(reads2)), log)
        if not reads1:
            break
        chunks = split_into_chunks(list(zip(reads1, reads2)), n_jobs)
        # processing
        outputs = Parallel(n_jobs=n_jobs)(delayed(chimera_clean_process_batch)(reads, csslist1, csslist2)
                                          for reads in chunks)
        results, stats = [x[0] for x in outputs], [x[1] for x in outputs]
        # writing results
        for result, stat in zip(results, stats):
            write_to_files([outfile1, outfile2, slagfile1, slagfile2], result)
            all_stats += stat
        if not silent:
            log.info("==== chimera_clean progress: reads processed: %d, time elapsed: %s"
                     % (all_stats.readcounter, time.strftime('%H:%M:%S', time.gmtime(time.time() - starttime))))
    parserR1.close()
    parserR2.close()

    outfile1.close()
    slagfile1.close()
    outfile2.close()
    slagfile2.close()

    if all_stats.TOTALmatecounter + all_stats.slagcounter != all_stats.readcounter:
        support.error("lucigen_nxmate.py, chimera_clean: error in the script somewhere! Unequal read counts!", log)
    if all_stats.readcounter == 0:
        support.error("lucigen_nxmate.py, chimera_clean: error in input data! Number of processed reads is 0!", log)
    if not silent:
        #print some stats
        percentmates = 100. * all_stats.matecounter / all_stats.readcounter
        percentslag = 100. * all_stats.slagcounter / all_stats.readcounter
        log.info("==== chimera_clean info: processing finished!")
        log.info("==== chimera_clean info: %d reads processed, %d true mate reads (%.2f %%) "
                 "and %d non-mates/chimeras (%.2f %%)."
                 % (all_stats.readcounter, all_stats.matecounter, percentmates, all_stats.slagcounter, percentslag))
        shortmates = all_stats.TOTALmatecounter - all_stats.matecounter
        log.info("==== chimera_clean info: %d mates too short to keep after trimming" % shortmates)
        elapsedtime = time.strftime('%H:%M:%S', time.gmtime(time.time() - starttime))
        log.info("==== chimera_clean info: time elapsed: %s" % (elapsedtime))
        log.info("==== chimera_clean info: " + str(all_stats.csscounter))
    return outfilename1, outfilename2


class JunctionStats(object):
    def __init__(self):
        self.readcounter = 0
        self.jctcounter = 0
        self.splitcounter = 0
        self.bothjctcounter = 0
        self.R1jctcounter = 0
        self.R2jctcounter = 0
        self.R1R2jctcounter = 0

    def __add__(self, other):
        self.readcounter += other.readcounter
        self.jctcounter += other.jctcounter
        self.splitcounter += other.splitcounter
        self.bothjctcounter += other.bothjctcounter
        self.R1jctcounter += other.R1jctcounter
        self.R2jctcounter += other.R2jctcounter
        self.R1R2jctcounter += other.R1R2jctcounter
        return self


def nx_seq_junction_process_batch(reads, jctstr):
    stats = JunctionStats()
    processed_split1 = []
    processed_split2 = []
    processed_unsplit = []

    for recR1, recR2 in reads:
        stats.readcounter += 1

        m = regex.search(jctstr, recR1[1])
        n = regex.search(jctstr, recR2[1])
        if m and n:  #found jctstr in both reads; need to save left part of R1 and LEFT part of R2
            stats.bothjctcounter += 1
            matches = m.span()
            start = matches[0]
            mySeq = recR1[1]  #get the left part of Read1
            myLeft = mySeq[:start]
            myQual = recR1[3]
            myLeftQual = myQual[:start]  #get left part of Read1 quality string

            nmatches = n.span()
            nstart = nmatches[0]
            mySeq2 = recR2[1]
            myRight2 = mySeq2[:nstart]  #get left part of Read2
            myQual2 = recR2[3]
            myRightQual2 = myQual2[:nstart]  #get left part of Read2 quality string

            #only write out as split if both pieces are big enough
            if (len(myLeft) > minseq) and (len(myRight2) > minseq):
                stats.splitcounter += 1
                stats.R1R2jctcounter += 1
                processed_split1.append([recR1[0], myLeft, recR1[2], myLeftQual])
                processed_split2.append([recR2[0], myRight2, recR2[2], myRightQual2])
        elif n:  #junction only in R2, so save entire R1 and LEFT part of R2 IFF R2 long enough
            nmatches = n.span()
            nstart = nmatches[0]
            mySeq2 = recR2[1]
            myRight2 = mySeq2[:nstart]
            myQual2 = recR2[3]
            myRightQual2 = myQual2[:nstart]
            if (len(myRight2) > minseq):
                stats.splitcounter += 1
                processed_split2.append([recR2[0], myRight2, recR2[2], myRightQual2])
                processed_split1.append([recR1[0], recR1[1], recR1[2], recR1[3]])
                stats.jctcounter += 1
                stats.R2jctcounter += 1
        elif m:  #junction only in R1, save left part of R1 and entire R2, IFF R1 is long enough
            matches = m.span()
            start = matches[0]
            mySeq = recR1[1]
            myLeft = mySeq[:start]
            myQual = recR1[3]
            myLeftQual = myQual[:start]
            if (len(myLeft) > minseq):
                stats.splitcounter += 1
                processed_split1.append([recR1[0], myLeft, recR1[2], myLeftQual])
                processed_split2.append([recR2[0], recR2[1], recR2[2], recR2[3]])
                stats.jctcounter += 1
                stats.R1jctcounter += 1
        else:  #no junctions, save for frag use, as is 'unsplit'; note this file will be interleaved R1 R2 R1 R2...
            processed_unsplit.append([recR1[0], recR1[1], recR1[2], recR1[3]])
            processed_unsplit.append([recR2[0], recR2[1], recR2[2], recR2[3]])
    return [processed_split1, processed_split2, processed_unsplit], stats


def nx_seq_junction(infilename1, infilename2, dst, log, silent=True):
    starttime = time.time()

    basename1 = os.path.basename(infilename1)
    if os.path.splitext(basename1)[1] == '.gz':
        basename1 = os.path.splitext(basename1)[0]
    basename2 = os.path.basename(infilename2)
    if os.path.splitext(basename2)[1] == '.gz':
        basename2 = os.path.splitext(basename2)[0]
    #open three outfiles
    splitfilenameleft = os.path.join(dst, 'R1_IJS7_' + basename1)
    splitfile1 = open(splitfilenameleft, 'w')

    splitfilenameright = os.path.join(dst, 'R2_IJS7_' + basename2)
    splitfile2 = open(splitfilenameright, 'w')

    unsplitfilename = os.path.join(dst, 'unsplit_IJS7_' + basename1.replace('_R1_', '_R1R2_'))
    unsplitfile = open(unsplitfilename, 'w')

    #jctstr = '(GGTTCATCGTCAGGCCTGACGATGAACC){e<=4}' # JS7 24/28 required results in ~92% detected in ion torrent
    # from NextClip: --adaptor_sequence GTTCATCGTCAGG -e --strict_match 22,11 --relaxed_match 20,10 eg strict 22/26 = 4 errors, relaxed 20/26 = 6 errors
    jctstr = '(GTTCATCGTCAGGCCTGACGATGAAC){e<=4}'  # try 22/26 to match NextClip strict (e<=6 for relaxed)

    #PARSE both files in tuples of 4 lines
    parserR1 = ParseFastQ(infilename1)
    parserR2 = ParseFastQ(infilename2)

    all_stats = JunctionStats()
    n_jobs = options_storage.threads
    while True:
        # prepare input
        reads1 = list(itertools.islice(parserR1, READS_PER_BATCH))
        reads2 = list(itertools.islice(parserR2, READS_PER_BATCH))
        if len(reads1) != len(reads2):
            support.error("lucigen_nxmate.py, nx_seq_junction: "
                          "number of left reads (%d) is not equal to number of right reads (%d)!"
                          % (len(reads1), len(reads2)), log)
        if not reads1:
            break
        chunks = split_into_chunks(list(zip(reads1, reads2)), n_jobs)
        # processing
        outputs = Parallel(n_jobs=n_jobs)(delayed(nx_seq_junction_process_batch)(reads, jctstr)
                                          for reads in chunks)
        results, stats = [x[0] for x in outputs], [x[1] for x in outputs]
        # writing results
        for result, stat in zip(results, stats):
            write_to_files([splitfile1, splitfile2, unsplitfile], result)
            all_stats += stat
        if not silent:
            log.info("==== nx_seq_junction progress: reads processed: %d, time elapsed: %s"
                     % (all_stats.readcounter, time.strftime('%H:%M:%S', time.gmtime(time.time() - starttime))))
    parserR1.close()
    parserR2.close()

    splitfile1.close()
    splitfile2.close()
    unsplitfile.close()

    if all_stats.readcounter == 0:
        support.error("lucigen_nxmate.py, nx_seq_junction: error in input data! Number of processed reads is 0!", log)
    if all_stats.splitcounter == 0:
        support.error("lucigen_nxmate.py, nx_seq_junction: error in input data! Number of split pairs is 0!", log)
    if not silent:
        #print some stats
        percentsplit = 100 * all_stats.splitcounter / all_stats.readcounter
        percentR1R2 = 100 * all_stats.R1R2jctcounter / all_stats.splitcounter
        percentR1 = 100 * all_stats.R1jctcounter / all_stats.splitcounter
        percentR2 = 100 * all_stats.R2jctcounter / all_stats.splitcounter
        log.info("==== nx_seq_junction info: processing finished!")
        log.info("==== nx_seq_junction info: %d reads processed" % (all_stats.readcounter))
        log.info("==== nx_seq_junction info: %d total split pairs (%.2f %% of processed reads))"
                 % (all_stats.splitcounter, percentsplit))
        log.info("==== nx_seq_junction info: %d junctions in both R1 and R2 (%.2f %% of split junctions))"
                 % (all_stats.R1R2jctcounter, percentR1R2))
        log.info("==== nx_seq_junction info: %d split junctions are in Read1 (%.2f %% of split junctions))"
                 % (all_stats.R1jctcounter, percentR1))
        log.info("==== nx_seq_junction info: %d split junctions are in Read2 (%.2f %% of split junctions))"
                 % (all_stats.R2jctcounter, percentR2))
        elapsedtime = time.strftime('%H:%M:%S', time.gmtime(time.time() - starttime))
        log.info("==== nx_seq_junction info: time elapsed: %s" % (elapsedtime))
    parserR1.close()
    parserR2.close()
    return splitfilenameleft, splitfilenameright, unsplitfilename


def process_reads(left_reads_fpath, right_reads_fpath, dst, log):
    log.info("== Processing Lucigen NxMate reads (" + left_reads_fpath + " and " +
             os.path.basename(right_reads_fpath) + " (results are in " + dst + " directory)")
    cleaned_filename1, cleaned_filename2 = chimera_clean(left_reads_fpath, right_reads_fpath, dst, log, silent=False)
    split_filename1, split_filename2, unsplit_filename = nx_seq_junction(cleaned_filename1, cleaned_filename2, dst, log, silent=False)
    return split_filename1, split_filename2, unsplit_filename
