# Stand-alone tools released within SPAdes package


## k-mer counter

To provide input data to SPAdes k-mer counting tool `spades-kmercounter ` you may just specify files in [SPAdes-supported formats](running.md#spades-input) without any flags (after all options) or provide dataset description file in [YAML format](running.md#specifying-multiple-libraries-with-yaml-data-set-file).

Output: `<output_dir>/final_kmers` - unordered set of kmers in binary format. Kmers from both forward and reverse-complementary reads are taken into account.

Output format: All kmers are written sequentially without any separators. Each k-mer takes the same number of bits. One k-mer of length K takes 2*K bits. Kmers are aligned by 64 bits. For example, one kmer with length=21 takes 8 bytes, with length=33 takes 16 bytes, and with length=55 takes 16 bytes. Each nucleotide is coded with 2 bits: 00 - A, 01 - C, 10 - G, 11 - T.
                                                 
Example:

        For k-mer: AGCTCT
        Memory: 6 bits * 2 = 12, 64 bits (8 bytes)
        Let’s describe bytes:
        data[0] = AGCT -> 11 01 10 00 -> 0xd8
        data[1] = CT00 -> 00 00 11 01 -> 0x0d
        data[2] = 0000 -> 00 00 00 00 -> 0x00
        data[3] = 0000 -> 00 00 00 00 -> 0x00
        data[4] = 0000 -> 00 00 00 00 -> 0x00
        data[5] = 0000 -> 00 00 00 00 -> 0x00
        data[6] = 0000 -> 00 00 00 00 -> 0x00
        data[7] = 0000 -> 00 00 00 00 -> 0x00

Synopsis: `spades-kmercount [OPTION...] <input files>`

The options are:

`-d, --dataset file <file name> `
    dataset description (in YAML format), input files ignored

`-k, --kmer <int> `
    k-mer length (default: 21)

`-t, --threads <int> `
    number of threads to use (default: number of CPUs)

`-w, --workdir <dir name> `
    working directory to use (default: current directory)

`-b, --bufsize <int> `
    sorting buffer size in bytes, per thread (default 536870912)

`-h, --help `
    print help message


## k-mer coverage read filter

`spades-read-filter` is a tool for filtering reads with median kmer coverage less than threshold.

To provide input data to SPAdes k-mer read filter tool `spades-read-filter ` you should provide dataset description file in [YAML format](running.md#specifying-multiple-libraries-with-yaml-data-set-file).

Synopsis: `spades-read-filter [OPTION...] -d <yaml>`

The options are:

`-d, --dataset file <file name> `
    dataset description (in YAML format)

`-k, --kmer <int> `
    k-mer length (default: 21)

`-t, --threads <int> `
    number of threads to use (default: number of CPUs)

`-o, --outdir <dir> `
    output directory to use (default: current directory)

`-c, --cov <value> `
    median kmer count threshold (read pairs, s.t. kmer count median for BOTH reads LESS OR EQUAL to this value will be ignored)

`-h, --help `
    print help message


## k-mer cardinality estimating

`spades-kmer-estimating` is a tool for estimating the approximate number of unique k-mers in the provided reads. Kmers from reverse-complementary reads aren"t taken into account for k-mer cardinality estimating.

To provide input data to SPAdes k-mer cardinality estimating tool `spades-kmer-estimating ` you should provide dataset description file in [YAML format](running.md#specifying-multiple-libraries-with-yaml-data-set-file).

Synopsis: `spades-kmer-estimating [OPTION...] -d <yaml>`

The options are:

`-d, --dataset file <file name> `
    dataset description (in YAML format)

`-k, --kmer <int> `
    k-mer length (default: 21)

`-t, --threads <int> `
    number of threads to use (default: number of CPUs)

`-h, --help `
    print help message


## Graph construction
Graph construction tool `spades-gbuilder ` has two mandatory options: dataset description file in [YAML format](running.md#specifying-multiple-libraries-with-yaml-data-set-file) and an output file name.

Synopsis: `spades-gbuilder <dataset description (in YAML)> <output filename> [-k <value>] [-t <value>] [-tmpdir <dir>] [-b <value>] [-unitigs|-fastg|-gfa|-spades]`

Additional options are:

`-k <int> `
    k-mer length used for construction (must be odd)

`-t <int> `
    number of threads

`-tmp-dir <dir_name>  `
    scratch directory to use

`-b <int> `
    sorting buffer size (per thread, in bytes)

`-unitigs `
    k-mer length used for construction (must be odd)

`-fastg `
    output graph in FASTG format

`-gfa `
    output graph in GFA1 format

`-spades `
    output graph in SPAdes internal format

## Graph simplification

Graph simplification tool `spades-gsimplifier` has four mandatory options: the first one is an input graph in GFA format, or a prefix of the SPAdes internal graph pack format created by setting [checkpoint options](running.md#pipeline-options). The second one is the prefix of the output simplified graph. The last two are the k-mer size and the read length.

Synopsis: `spades-gsimplifier <graph. In GFA (ending with .gfa) or prefix to SPAdes graph pack> <output prefix> [--gfa] [--spades-gp] [--use-cov-ratios] -k <value> --read-length <value> [OPTION...]`

Additional options are:

`--gfa`
    produce GFA output (default: true)

`--spades-gp` 
    produce output graph pack in the SPAdes internal format (default: false)

`--use-cov-ratios` 
    enable simplification procedures based on unitig coverage ratios (default: false)

`-k <value>`
    k-mer length to use

`--read-length <value>`
    read length to use

`-c, --coverage <coverage>`
    estimated average (k+1-mer) bin coverage (default: 0.) or 'auto' (works only with '-d/--dead-ends' provided)

`-t, --threads <value>`
    number of threads to use (default: max_threads / 2)

`-p, --profile <file>`
    file with edge coverage profiles across multiple samples

`-s, --stop-codons <file>`
    file with stop codon positions

`-d, --dead-ends <file>`
    while processing a subgraph -- file listing edges which are dead-ends in the
    original graph


## Long read to graph alignment


### hybridSPAdes aligner
A tool `spades-gmapper ` gives the opportunity to extract long read alignments generated with hybridSPAdes pipeline options. It has three mandatory options: dataset description file in [YAML format](running.md#specifying-multiple-libraries-with-yaml-data-set-file), graph file in GFA format and an output file name.

Synopsis: `spades-gmapper <dataset description (in YAML)> <graph (in GFA)> <output filename> [-k <value>] [-t <value>] [-tmpdir <dir>]`

Additional options are:

`-k <int> `
    k-mer length that was used for graph construction

`-t <int> `
    number of threads

`-tmpdir <dir_name>  `
    scratch directory to use

While `spades-gmapper` is a solution for those who work on hybridSPAdes assembly and want to get exactly its intermediate results, [SPAligner](standalone.md#spaligner) is an end-product application for sequence-to-graph alignment with tunable parameters and output types.


### SPAligner
A tool for fast and accurate alignment of nucleotide sequences to assembly graphs. It takes file with sequences (in fasta/fastq format) and assembly in GFA format and outputs long read to graph alignment in various formats (such as tsv, fasta and [GPA](https://github.com/ocxtal/gpa "GPA-format spec")).

Synopsis: `spaligner src/projects/spaligner_config.yaml -d <value> -s <value> -g <value> -k <value> [-t <value>] [-o <value>]`

Parameters are:

`-d <type> `
    long reads type: nanopore, pacbio

`-s <filename> `
    file with sequences (in fasta/fastq)

`-g <filename> `
    file with graph (in GFA)

`-k <int> `
    k-mer length that was used for graph construction

`-t <int> `
    number of threads (default: 8)

`-o, --outdir <dir> `
    output directory to use (default: spaligner_result/)

For more information on parameters and options please refer to the main SPAligner manual (assembler/src/projects/spaligner/README.md).

Also if you want to align protein sequences please refer to our [pre-release version](https://github.com/ablab/spades/releases/tag/spaligner-paper).

Note that in order you use SPAligner one needs either to use pre-built binaries or compile SPAdes from sources using the additional `-DSPADES_ENABLE_PROJECTS=spaligner` option.
