############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import moleculo_filter_contigs
import break_by_coverage
import SeqIO
import sys
import generate_quality
import sam_parser

pattern = "TACGCTTGCAT"
rc_pattern = "ATGCAAGCGTA"

def SplitAndFilter(contigs, coverage_breaker, length_filter, n_breaker, pattern_breaker, pattern_filter):
    result = []
    for contig in contigs:
        if (pattern_filter.Filter(contig)):
            for subcontig in coverage_breaker.Break(contig):
                for subsubcontig in pattern_breaker.Break(subcontig):
                    for subsubsubcontig in n_breaker.Break(subsubcontig):
                        if length_filter.Filter(subsubsubcontig):
                            result.append(subsubsubcontig)
    return result


def OutputResults(output_file, format, result):
    output = open(output_file + "." + format, "w")
    for contig in result:
        SeqIO.write(contig, output, format)
    output.close()


def moleculo_postprocessing(contigs_file, output_file, sam_files, log):
    log.info("===== Starting postprocessing based on read alignment")
    log.info("Processing scaffolds from " + contigs_file)
    log.info("Using read alignments to break and filter scaffolds")
    contigs = list(SeqIO.parse(open(contigs_file, "rU"), "fasta"))
    sam = sam_parser.SamChain([sam_parser.Samfile(sam_file) for sam_file in sam_files])
    generate_quality.GenerateQuality(contigs, sam)
    pattern_filter = moleculo_filter_contigs.PatternContigFilter(contigs, sam, pattern, rc_pattern)
    length_filter = moleculo_filter_contigs.ContigLengthFilter(1500)
    coverage_breaker = break_by_coverage.ContigBreaker(contigs, sam, 100, 50)
    pattern_breaker = break_by_coverage.PatternBreaker(pattern, rc_pattern, 150)
    n_breaker = break_by_coverage.NBreaker(3)
    result = SplitAndFilter(contigs, coverage_breaker, length_filter, n_breaker, pattern_breaker, pattern_filter)
    OutputResults(output_file, "fasta", result)
    OutputResults(output_file, "fastq", result)
    log.info("===== Postprocessing finished. Results can be found in " + output_file + ".fastq")

