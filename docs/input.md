# SPAdes basic input

SPAdes takes as input paired-end reads, mate-pairs and single (unpaired) reads in FASTA and FASTQ (can be gzipped).
For IonTorrent data SPAdes also supports unpaired reads in unmapped BAM format (like the one produced by Torrent Server).
However, in order to run read error correction, reads should be in FASTQ or BAM format.
Sanger, Oxford Nanopore and PacBio CLR reads can be provided in both formats since SPAdes does not run error correction for these types of data.

To run SPAdes you need at least one library of the following types:

-   Illumina paired-end/high-quality mate-pairs/unpaired reads
-   IonTorrent paired-end/high-quality mate-pairs/unpaired reads
-   PacBio CCS reads

Illumina and IonTorrent libraries should not be assembled together. All other types of input data are compatible.
SPAdes should not be used if only PacBio CLR, Oxford Nanopore, Sanger reads or additional contigs are available.

SPAdes supports mate-pair only assembly. However, we recommend to use only high-quality mate-pair libraries in this case
(e.g. that do not have a paired-end part). We tested the mate-pair-only pipeline using Illumina Nextera mate-pairs. See more [here](running.md#specifying-multiple-libraries).

Notes:

-   It is strongly recommended to provide multiple paired-end and mate-pair libraries according to their insert size (from smallest to longest).
-   It is not recommended to run SPAdes on PacBio reads with low coverage (less than 5).
-   We suggest not to run SPAdes on PacBio reads for large genomes.

## Paired read libraries

By using command line interface, you can specify (1) paired-end, (2) standard mate-pair and (3) high-quality mate-pair libraries.
You can provide up to nine libraries of each type. Libraries can be used in any combination, but recommend not to assemble low-quality mate-pairs alone.
If you wish to use more libraries, you can use [YAML data set file](running.md#specifying-multiple-libraries-with-yaml-data-set-file).

By default, SPAdes assumes that paired-end and high-quality mate-pair reads have forward-reverse (fr) orientation and usual mate-pairs have reverse-forward
(rf) orientation. However, different orientations can be indicated for any library by using SPAdes options.

To distinguish reads in pairs we refer to them as left and right reads. For forward-reverse orientation, the forward reads correspond to the left reads and the reverse reads, to the right. Similarly, in reverse-forward orientation left and right reads correspond to reverse and forward reads, respectively, etc.

Each paired-end or mate-pair library can be stored in several files or several pairs of files. Paired reads can be organized in two different ways:

-   In file pairs. In this case left and right reads are placed in different files and must go in the same order.
I.e. for every left read at line X in the first file the corresponding right read from the pair must be at line X in the second file.
-   In interleaved files. In this case, the reads are interlaced, so that each right read goes after the corresponding paired left read.

For example, Illumina produces paired-end reads in two files: `R1.fastq` and `R2.fastq`. If you choose to store reads in file pairs make sure that for every read from `R1.fastq` the corresponding paired read from `R2.fastq` is placed in the respective paired file on the same line number. If you choose to use interleaved files, every read from `R1.fastq` should be followed by the corresponding paired read from `R2.fastq`.

If adapter and/or quality trimming software has been used prior to assembly, files with the orphan reads can be provided as "single read files" for the corresponding read-pair library.

If you have merged some of the reads from your paired-end (not mate-pair or high-quality mate-pair) library (using tools s.a. [BBMerge](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) or [STORM](https://bitbucket.org/yaoornl/align_test/overview)), you should provide the file with resulting reads as a "merged read file" for the corresponding library.
Note that non-empty files with the remaining unmerged left/right reads (separate or interlaced) **must** be provided for the same library (for SPAdes to correctly detect the original read length).

In an unlikely case some of the reads from your mate-pair (or high-quality mate-pair) library are "merged", you should provide the resulting reads as a SEPARATE single-read library.

See [examples](running.md#examples).

## Unpaired (single-read) libraries

By using the command line interface, you can specify up to nine different single-read libraries. To input more libraries, you can use [YAML data set file](running.md#specifying-multiple-libraries-with-yaml-data-set-file).

Single libraries are assumed to have high quality and reasonable coverage. For example, you can provide PacBio CCS reads as a single-read library.

Note, that you should not specify PacBio CLR, Sanger reads or additional contigs as single-read libraries, each of them has a separate [option](running.md#input-data).

See [examples](running.md#examples).