

Fragment assembly using EULER-SR

December 14, 2008 Contents 1 Installing EULER-SR. The installation of EULER-SR has been made more straightforward from theprevious versions. You probably will only need to set one environment variable. The first step is to unpack euler-sr.tgz, probably already done since you arereading this.

> tar zxvf euler-sr.tgz. 1.1 Required software. Software: You must have g++ and gmake installed on your system.Hardware: EULER-SR is only tested on an x86 64 linux system. It should

run on the intel macs, and has limited support for running on a powerbook.

1.2 Environment variables. You need to set EUSRC to the full path to where euler-sr was unzipped, andMACHTYPE to one of x86 64, i686, or powerpc:

If you are running bash: > tar zxvf euler-sr.tgz > cd euler-sr > export EUSRC=`pwd` > export MACHTYPE=x86\_64

If you are running tcsh: > tar zxvf euler-sr.tgz > cd euler-sr > setenv EUSRC `pwd` > setenv MACHTYPE x86\_64

1

1.3 Build euler-sr. If all environment variables are set, this is a straightforward make. >cd euler-sr >gmake

1.4 Running programs. When ran without arguments (or when there is an error on the command line),all programs produce a brief help file.

1.5 Sample data. Provided with the assembly is a sample paired end reads file (reads.fasta), arule file for this (readtitles.rules), and a sample file that contains reads from

two haplotypes (reads.variants.fasta).To assemble the reads file:

${EUSRC}/assembly/Assemble.pl reads.fasta 25 readtitle.rules. Instructions to extract the indel variants from the haplotypes file are givenbelow.

2 CLEANING UP DATA The assembly quality is strongly dependent on how clean the data is. Qualityvalues may be used to trim/remove low quality reads. If you are assembling Illumina reads, particularly those produced by the first generation sequencer,some extra programs are provided to filter reads that have errors particular to the system (discussed below).

2.1 Translating data. EULER-SR takes as input a single fasta file with ALL reads, an optional matefile that specifies mate-pairing (which may be generated), and a rule file that

describes how to map a fasta title to a clone type.Utilities are provided to translate from: .sff (the 454 binary output file). .fastq (Sanger format) .fastq (Illumina quality values). 2.1.1 FASTQ -? Fasta (Sanger) To translate fastq to fasta, use the quality trimmer, given reads.fastq

${EUSRC}/assembly/x86\_64/qualityTrimmer -fasta reads.fastq -outFasta reads.fasta

Options may be given to quality trimmer to tune the trimming, and aredescribed in the output help.

2

2.1.2 FASTQ -? Fasta (Illumina) Illumina uses a different scaling of quality values than the normal phred scores.Although qualityTrimmer attempts to guess the type (Sanger vs. Illumina),

this may be forced on the command line: ${EUSRC}/assembly/x86\_64/qualityTrimmer -fasta reads.fastq -outFasta reads.fasta -type illumina

2.1.3 sff -? fasta The 454 output files by default are in a binary format that contains both thefasta, quality, and flow values. Newbler can take advantage of the flow values,

but euler-sr does not. To translate from sff to fasta, use sff2fasta

${EUSRC}/assembly/x86\_64/sff2fasta reads.sff reads.fasta

Some diagnostic information will be printed regarding the parsing. If theparser is able to read the binary format, the diagnostic information will be

sensical (below). If 454 has updated the sff format, and sff2fasta is out of date,the following sample diagnostic output will look like junk:

magic number: 779314790 version: 1 index offset: 535227192 index length: 6438410 number of reads: 321891 header length: 440 key length: 4 flowsPerRead: 400 format code: 1 flow chars: TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG

2.2 Quality filtering. If you have quality values, you can trim low quality sequences:

${EUSRC}/assembly/${MACHTYPE}/qualityTrimmer -fasta reads.fasta -qual reads.fasta.qual -outFasta reads.trimmed The sensitivity of the trimmer may be tuned with -span S -minQual Q, wherereads are trimmed so that the entire read is of min quality Q for each consecutive

S nucleotides.

2.3 Filtering bad illumina reads: Illumina reads, particularly ones from the first generation machine, have par-ticular output values. An example bad read is:

3

>read1 ACGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

These may be filtered using filterIlluminaReads

${EUSRC}/assembly/${MACHTYPE}/filterIlluminaReads reads.fasta reads.filt.fasta 2.4 Preparing paired-end 454 reads. 454 Paired-ends are prepared by linking ends of a clone with a linker, shearingthe clones, pulling out the linked sequences (hopefully linked), then sequencing:

READ1-linker-READ2.To split these, you will need to know the sequence of the linker, here called linker.fasta

${EUSRC}/assembly/${MACHTYPE}/splitLinkedClones reads.454.fasta linker.fasta reads.454.fasta.split

To print reads that are not mated to a separate file, use: ${EUSRC}/assembly/${MACHTYPE}/splitLinkedClones reads.454.fasta linker.fasta reads.454.fasta.split -singletons reads.454.fasta.singletons

2.5 Cleaning vector and masked sequences from SANGER

reads.

The vector and masked sequences in Sanger reads can cause problems for errorcorrection and assembly. These can be filtered by piping them through two

utilities cat reads.fasta ${EUSRC}/assembly\_utils/RemoveMaskedSequence.pl | ${EUSRC}/assembly\_utils/RemoveVectorSequence.pl > reads.trimmed

3 Running assemblies. 3.1 Preparing your input file. If you have multiple input files, concatenate them into one file: >cat file1.fasta file2.fasta file3.fasta > reads.fasta

3.2 Running a default assembly. You can run a default assembly using:

${EUSRC}/assembly/Assemble.pl reads.fa 25

4

3.3 Running a paired-end assembly. You will have to generate a "rule file" that describes the mates you are using.This is a file with a regular expression in quotes, followed by two required

keywords, CloneLength and CloneVar.An example line is:

"([^/]*)/([12])" CloneLength=200 CloneVar=50

In general, the format of a pair of reads should be:

CLONE\_NAME.DIR1 CLONE\_NAME.DIR2

where CLONE NAME is the unique name of each clone, and DIR1 and DIR2name what side from the clone the reads are from.

The regular expressions follow the extended regular expression format, de-scribed at:

http://opengroup.org/onlinepubs/007908775/xbd/re.html#tag\_007\_004

Each regular expression needs two sub matches, one for the clone name, andthe other for the direction. No more than two reads should have the same clone

name.The CloneLength keyword is the expected gap from the end of the first read to the beginning of the second read, or the length of the DNA fragment fromwhich both reads are sequenced - the sum of the two read lengths.

The CloneVar keyword is expected window of variation of clone lengths. Thegap separating two reads in an assembly is estimated, and if it is greater than CloneLength + CloneVar, or less than CloneLength - CloneVar, the mate-pairis considered invalid (it may be a chimeric clone), and it is not used.

Examples of mate-files are:454 clones:

>000556\_0922\_1963.a ACTGGCGAGAGCCCAGACGT... >000556\_0922\_1963.b GCCCGAGACCGGACTGGGAT...

The rule line for this is:

"([0-9]+\_[0-9]+\_[0-9]+)\.([ab])" CloneLength=2500 CloneVar=500 Type=6

Illumina clones:

>SLXA-EAS1\_89:3:1:715:750/1 GTCTTGAAAGCTATGATGTCAAGATTAATTTAATC >SLXA-EAS1\_89:3:1:715:750/2 GTGTATTGCTCAATCTTCGAACGGGGGGAGGATTG

5

The rule line for this is: "([^/]*)/([12])" CloneLength=200 CloneVar=50 Type=1 # Illumina mate-pair

3.4 Nonstandard assemblies. 3.4.1 Detecting variation. To detect variants, not all the steps to run an assembly need to be performed.Since variant detection is still being developed, you have to manually run most

assembly steps, however a final post-processing tool exists that prints likelyvariants.

To detect variants using a file reads.fasta: 1. Fix errors:

${EUSRC}/assembly/FixErrors.pl reads.fasta 25

2. Build the de bruijn graph: ${EUSRC}/assembly/assemblesec.pl reads.fasta.fixed -vertexSize 25

3. Remove some errors from the graph mkdir simple; ${EUSRC}/assembly/${MACHTYPE}/simplifyGraph reads.fasta.fixed simple/reads.fasta -minEdgeLength 75 -removeLowCoverage 5 3

4. Run post-processing to detect variants. ${EUSRC}/assembly/${MACHTYPE}/printVariants simple/reads.fasta reads.variants.fasta

This will produce a fasta file, reads.variants.fasta. You can use the fastatitles to collect the variants. The format is:

>SHARED\_EDGE INDEX (length, read support) VARIANT\_INDEX (length,

support) [... VARIANT\_INDEX2 (length2, support2)] SHARED\_EDGE2(length,support)

This outputs a variant and the flanking sequences so that the variant may beeasily aligned to a reference. The numbers in the FASTA titles are: EDGE INDEX(LENGTH,SUPPORT), where EDGE INDEX is the index of the edge in the ".edge" file (the same as acontig, but edges exist for both the forward and reverse orientation), the length of the edge (the middle lengths represent the lengths of the variants), and finallythe support, the number of reads used to build the edge in the assembly. Very low support indicates that the variant is probably erroneous. In the examplebelow, there is an indel of length 6 in the assembly. The alignment below is the bl2seq alignment of the two variants.

6

>0(1855, 2957) 21(49, 44) 8(2274, 3692) AGAAAGGGATTTTAGTTTGTAATATCGCAGCAAGTCGATTGATTTTACCGTCTCCCAATGCATTCAAAGATAGTATTGTAAAAATCTCAGTTGGTGAAGAATATGATCAACACGCGTTTATCCATCAGTTAAAGGAAAATGGCTATCGAAAAGTTACTCAAGTACAAACTCAGGGCGAATTTAGTCTTCGAGGAGATATTATTTTTGAAATATCCCAGTTAGAACCTTGTCGAATTGAGTTTTTTGGTGATGAAATTGATGGTATCAGGTCATTTGAAGTAGAAACACAATTATCGAAAGAAAATAAGACAGAACTCACTATCTTTCCAGCTAGTGATATGCTTTTGAGAGAAAAGGATTATCAACGAGGACAGTCAGCTTTAGAAAAACAAATTTCAA >0(1855, 2957) 22(55, 39) 8(2274, 3692) AGAAAGGGATTTTAGTTTGTAATATCGCAGCAAGTCGATTGATTTTACCGTCTCCCAATGCATTCAAAGATAGTATTGTAAAAATCTCAGTTGGTGAAGAATATGATCAACACGCGTTTATCCATCAGTTAAAGGAAAATGGCTATCGAAAAGTTACTCAAGTACAAACTCAGGGCGAATTTAGTCTTCGAGGAGATATTTTAGATATTTTTGAAATATCCCAGTTAGAACCTTGTCGAATTGAGTTTTTTGGTGATGAAATTGATGGTATCAGGTCATTTGAAGTAGAAACACAATTATCGAAAGAAAATAAGACAGAACTCACTATCTTTCCAGCTAGTGATATGCTTTTGAGAGAAAAGGATTATCAACGAGGACAGTCAGCTTTAGAAAAACAAATTTCAA

Query 1 AGAAAGGGATTTTAGTTTGTAATATCGCAGCAAGTCGATTGATTTTACCGTCTCCCAATG 60

|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| Sbjct 1 AGAAAGGGATTTTAGTTTGTAATATCGCAGCAAGTCGATTGATTTTACCGTCTCCCAATG 60

Query 61 CATTCAAAGATAGTATTGTAAAAATCTCAGTTGGTGAAGAATATGATCAACACGCGTTTA 120

|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| Sbjct 61 CATTCAAAGATAGTATTGTAAAAATCTCAGTTGGTGAAGAATATGATCAACACGCGTTTA 120

Query 121 TCCATCAGTTAAAGGAAAATGGCTATCGAAAAGTTACTCAAGTACAAACTCAGGGCGAAT 180

|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| Sbjct 121 TCCATCAGTTAAAGGAAAATGGCTATCGAAAAGTTACTCAAGTACAAACTCAGGGCGAAT 180

Query 181 TTAGTCTTCGAGGAGATAT------TATTTTTGAAATATCCCAGTTAGAACCTTGTCGAA 234

||||||||||||||||||| ||||||||||||||||||||||||||||||||||| Sbjct 181 TTAGTCTTCGAGGAGATATTTTAGATATTTTTGAAATATCCCAGTTAGAACCTTGTCGAA 240

Query 235 TTGAGTTTTTTGGTGATGAAATTGATGGTATCAGGTCATTTGAAGTAGAAACACAATTAT 294

|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| Sbjct 241 TTGAGTTTTTTGGTGATGAAATTGATGGTATCAGGTCATTTGAAGTAGAAACACAATTAT 300

Query 295 CGAAAGAAAATAAGACAGAACTCACTATCTTTCCAGCTAGTGATATGCTTTTGAGAGAAA 354

|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| Sbjct 301 CGAAAGAAAATAAGACAGAACTCACTATCTTTCCAGCTAGTGATATGCTTTTGAGAGAAA 360

Query 355 AGGATTATCAACGAGGACAGTCAGCTTTAGAAAAACAAATTTCAA 399

||||||||||||||||||||||||||||||||||||||||||||| Sbjct 361 AGGATTATCAACGAGGACAGTCAGCTTTAGAAAAACAAATTTCAA 405

3.4.2 Fast pooling of reads. If the error rate of the reads is low, and all that is required is a clustering ofreads that have shared k-mers, simply running assemblesec is sufficient.

${EUSRC}/assembly/assemblesec.pl reads.fasta -vertexSize 25

You can increase specificity by increasing the vertex size, although runningwith vertex sizes greater than 32 is VERY slow. You will have to parse the file "reads.fasta.intv" file to determine the readsthat overlap. The format of this file is:

7

EDGE E Length L Multiplicity M INTV 1 0 50 135 INTV 8 0 50 140 INTV a b c d ... EDGE ... INTV ...

The edge E is an edge in the assembly, corresponding to either the forwardor reverse contig (edges for both are created). Length is the length of the edge,

and M is the number of reads included in the assembly of the edge.each INTV corresponds to a portion of a read that maps to the edge. In INTV a b c d, "a" is the index of the read. Odd numbered indices are thereverse complement of a read included in the data. So read 0 is a read in your input file, and read 1 is the reverse complement of 0, 2 is a read, 3 is 2's reversecomplement, and so on. "b" is the starting position in the read mapped to the edge, "c" is the length of the portion of the read mapped to the edge, and "d"is the position along the edge where the read maps.

When "c" is the length of the read, the full read maps to the edge, when itis less, the read either maps to two edges, or part of the read does not map to the assembly, usually when there are errors in the read.

4 Examining output and refining assemblies. Some assemblies may be improved by iteratively re-simplifying the graph, orexamined to determine if the data are insufficient for producing a valid assembly.

4.1 Detecting fragmented assemblies. If you have mate-pairs, the assembly will be in the 'matetransformed' directory,otherwise the 'transformed' directory. For examples here the 'transformed' directory is used.Print all 'dead-end' edges:

${EUSRC}/assembly/${MACHTYPE}/printGraphSummary transformed/reads.fasta -sources

This will print edge indices followed by their lengths. If there is a very highnumber of dead-end edges, the coverage is probably low causing a fragmented

assembly. Some sequencing platforms such as Illumina have a lot of "coveragedeserts", such as in GC-rich regions, where the coverage is very low, causing fragmentation of the assembly.

4.2 Joining fragmented assemblies. If there are a lot of dead end edges, you can try and fix them by joining themwhen the overlap is less than the vertex size, but greater than some minimal

8

length, say 10. ${EUSRC}/assembly/${MACHTYPE}/joinsas transformed/reads.fasta 10 transformed/reads.j ${EUSRC}/assembly/${MACHTYPE}/printContigs transformed/reads.j cp transformed/reads.j.contig .

This will look for unambiguous joins of edges that are at least 10 exactlymatching nucleotides.

5 Release notes. 5.1 Known issues. When using mate-pairs to resolve repeats, the sequence of the repeat will bea consensus, so there will likely be many mismatches in repeat regions. The

overall accuracy of the layout of the alignment is usually correct though.Support for running on a rocks cluster has been removed. Since EULER-SR is 100X faster than the first release, running on a cluster is less important. Still,a new release will use multiple cores for error correction.

5.2 Notes. This release contains most new functionality described in a recently acceptedpaper: De novo fragment assembly with short mate-paired reads: does the read

length matter? The significant improvements over the previous release are:- Mate-paired assembly is used. - Multiple clone libraries are allowed. - Repeats are resolved using longer reads by default (previously this was justdone as part of a "beta" assembly step). - Speed. While there are some areas for improvement, the speed of the method is much faster than before.If speed is still an issue, more improvements will be made, however they will take lower priority to new functionality.- Easier installation. Now only one (possibly two) environment variables need to be set to install and run EULER-SR. Before, dynamic linking was usedto save disk space with executables, however since most data sets are in the 10's of gigabytes now, a few megs of executables will probably fit anywhere.- Easier running. Previously different scripts were used to assemble illumina and 454 (or any non illumina) reads. This is now merged into one script, al-though pre-processing is done by hand at an earlier step. This is necessary for heterogeneous data sets.

9
