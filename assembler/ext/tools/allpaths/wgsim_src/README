Introduction
============

Wgsim is a small tool for simulating sequence reads from a reference genome.
It is able to simulate diploid genomes with SNPs and insertion/deletion (INDEL)
polymorphisms, and simulate reads with uniform substitution sequencing errors.
It does not generate INDEL sequencing errors, but this can be partly
compensated by simulating INDEL polymorphisms.

Wgsim outputs the simulated polymorphisms, and writes the true read coordinates
as well as the number of polymorphisms and sequencing errors in read names.
One can evaluate the accuracy of a mapper or a SNP caller with wgsim_eval.pl
that comes with the package.


Compilation
===========

gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm


History
=======

Wgsim was modified from MAQ's read simulator by dropping dependencies to other
source codes in the MAQ package and incorporating patches from Colin Hercus
which allow to simulate INDELs longer than 1bp. Wgsim was originally released
in the SAMtools software package. I forked it out in 2011 as a standalone
project. A few improvements were also added in this course.


Evaluation
==========

Simulation and evaluation
-------------------------

The command line for simulation:

  wgsim -Nxxx -1yyy -d0 -S11 -e0 -rzzz hs37m.fa yyy-zzz.fq /dev/null

where yyy is the read length, zzz is the error rate and $xxx * $yyy = 10000000.
By default, 15% of polymorphisms are INDELs and their lengths are drawn from a
geometric distribution with density 0.7*0.3^{l-1}.

The command line for evaluation:

  wgsim_eval.pl unique aln.sam | wgsim_eval.pl alneval -g 20

The '-g' option may be changed with mappers.


System
------

GCC: 4.1.2
CPU: AMD Opteron 8350 @ 2.0GHz
Mem: 128GB


Results
-------

==================================================================================================================
                          100bp              200bp              500bp              1000bp            10000bp
                   ------------------  -----------------  -----------------  -----------------  -----------------
 Program  Metrics     2%    5%   10%     2%    5%   10%     2%    5%   10%     2%    5%   10%     2%    5%   10%
------------------------------------------------------------------------------------------------------------------
            CPU      249   198   136    325   262   163    332   243   232    320   235   215    235   197   189
 BWA-SW     Q20%    85.1  63.6  21.4   93.7  88.9  53.5   96.4  95.7  89.2   96.6  96.2  95.1   97.7  98.3  97.7
            err%    0.01  0.06  0.20   0.00  0.01  0.14   0.00  0.01  0.01   0.00  0.00  0.01   0.00  0.00  0.00
            one%    94.6  77.4  35.7   97.5  95.1  67.6   98.6  98.5  93.4   99.0  98.9  98.3   99.7  99.8  99.7
------------------------------------------------------------------------------------------------------------------
            CPU                                            302   484  1060    330   352   607    381   480   919
 AGILE      Q20%                                          98.6  98.4  98.4   98.4  98.4  98.6   98.2  98.6  99.3
            err%                                          0.66  0.69  2.31   0.34  0.40  0.70   0.10  0.00  0.20
            one%                                           100  99.4     0    100   100   100    100   100   100
==================================================================================================================

1) AGILE throws "Floating point exception" halfway for 100/200bp reads.  The
   default output is supposed to be PSL, but actually has an additional "score"
   column. AGILE is reportedly faster than BWA-SW for 1000bp reads. It is
   slower here possibly because of suboptimal command line options.

2) Gassst uses over 27GB memory in 20 minutes. The memory then quickly
   increases to over 40GB. It gets killed.

3) Lastz complains: "FAILURE: bad fasta character in hs37m.fa ...".

4) Pash only gives 'unique mapping'. Its unique mapping is better than BWA-SW's
   Q1 mapiping. It is very slow, though, possibly because of suboptimal
   options.

