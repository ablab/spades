# parasail applications

Though the parasail library is the focus of this work, we do provide one example application that may be of immediate use.

## parasail_aligner

The aligner tool will take a FASTA- or FASTQ-formatted database file as input and align all of the database sequences against a set of query sequences found in a second FASTA- or FASTQ-formatted database file.  Alternatively, if only one file is supplied, all of the sequences in the file will be compared against themselves.  If parasail was linked against libz, the input files can also be gzip compressed.

The parasail_aligner can also receive an input file piped from another program in order to use the aligner as part of a larger workflow. If the database file is specified with `-q`, and input is detected on stdin, then this input is assumed to be the query file. If neither the database file or the query file is specified, and input is detected on stdin, then this input is assumed to be the single database file.

### Command-Line Interface

```bash
usage: parasail_aligner [-a funcname] [-c cutoff] [-x] [-e gap_extend] [-o gap_open] [-m matrix] [-t threads] [-d] [-M match] [-X mismatch] [-k band size (for nw_banded)] [-l AOL] [-s SIM] [-i OS] [-v] [-V] -f file [-q query_file] [-g output_file] [-O output_format {EMBOSS,SAM,SAMH,SSW}] [-b batch_size] [-r memory_budget]

Defaults:
     funcname: sw_stats_striped_16
       cutoff: 7, must be >= 1, exact match length cutoff
           -x: if present, don't use suffix array filter
   gap_extend: 1, must be >= 0
     gap_open: 10, must be >= 0
       matrix: blosum62
           -d: if present, assume DNA alphabet
        match: 1, must be >= 0
     mismatch: 0, must be >= 0
      threads: Warning: ignored; OpenMP was not supported by your compiler
          AOL: 80, must be 0 <= AOL <= 100, percent alignment length
          SIM: 40, must be 0 <= SIM <= 100, percent exact matches
           OS: 30, must be 0 <= OS <= 100, percent optimal score
                                           over self score
           -v: verbose output, report input parameters and timing
           -V: verbose memory output, report memory use
         file: no default, must be in FASTA format
   query_file: no default, must be in FASTA format
  output_file: parasail.csv
output_format: no deafult, must be one of {EMBOSS,SAM,SAMH,SSW}
   batch_size: 0 (calculate based on memory budget),
               how many alignments before writing output
memory_budget: 2GB or half available from system query (X.XXX GB)
```

#### Using the Enhanced Suffix Array Filter

One feature of this tool is its ability to filter out sequence pairs based on an exact-match cutoff.  Using the cutoff paramter (`-c`), the filter will keep only those pairs of sequences which contain an exact match of length greater than or equal to the cutoff.  The assumption is that any pair of sequences which are highly similar should also contain an exact-matching k-mer of length at least c, our cutoff.  This is similar to the seed and extend model of sequence alignment found in other tools, however, our filter allows for arbitrarily long exact-matches (aka seeds) and once a match is found the entire alignment is performed rather than extending the seed.  The filter is turned on by default but can be disabled using the -x command-line parameter.

#### Alignment Function

Please review the function naming conventions of the top-level parasail README.

The alignment routine to use defaults to one of the stats Smith-Waterman routines, but any of the parasail routines (including the profile-based routines and the global banded routine) may be selected using the appropriate command-line parameter (`-a`).  Follow the naming conventions for parasail functions in order to select the desired alignment routine.

#### DNA Mode

The aligner assumes amino acid sequences; use `-d` to indicate DNA sequences as well as `-M` and `-X` to indicate the match and mismatch scores, respectively.

#### Substitution Matrix Selection

Please review the substitution matrix section of the top-level parasail README.

Using the `-m` parameter, you can select from a variety of built-in matrices, e.g., blosum50, pam100. You can also specify a matrix filename if the file is of the appropriate format (see README). For DNA alignments using a simple match/mismatch scoring criteria, use `-M` and `-X` to create the simple substitution matrix (see above).

#### Threading

If your compiler supports OpenMP and parasail's configure script is able to detect the support, then the parasail_aligner will be compiled with OpenMP parallel for loops. You can control the number of threads using the `-t` parameter.

If you are using threading in conjunction with batch mode, threading occurs within a batch of alignments.

#### Memory Concerns and Batch Mode

For very large inputs or when using traceback-capable routines, the parasail_aligner can use a significant amount of memory. Instead of storing all alignment results until the end, the aligner will complete a batch of alignments, write the results, and free their memory before attempting the next batch of alignments.  There are two modes for using batches, either explicitly specifying how many alignments per batch `-b` or specifying a memory budget `-r`.

The `-b` parameter indicates how many alignments to perform before writing results to output and freeing their memory. The default is 0 indicating the memory budget will be used instead.

The `-r` parameter indicates how much memory can be used. By default, it will query the system for the amount of physical memory and set the limit to half of the physical memory.

The larger the batch size, the better the runtime performance. This is a tuning parameter to balance between memory requirements and performance. Ideally, you will not need to specify either batch size or memory budget; the default settings are sufficient for most cases.

### Output

#### Traceback Output

If a trace-capable alignment function is used, e.g., `sw_trace`, then the output will be to stdout by default.  You can optionally redirect to a file using the `-g` option.  You must select an output format from one of SAM, SAMH, EMBOSS, or SSW.  This is what the outputs look like using the single target and single query sequence shown below, using the default values for the aligner and the `sw_trace` function.

Target:

```
>AF0017_1 COG1250 # Protein_GI_number: 11497638 # Func_class: I Lipid transport and metabolism  # Function: 3-hydroxyacyl-CoA dehydrogenase # Organism: Archaeoglobus fulgidus
MMVLEIRNVAVIGAGSMGHAIAEVVAIHGFNVKLMDVSEDQLKRAMEKIEEGLRKSYERGYISEDPEKVLKRIEATADLIEVAKDADLVIEAIPEIFDLKKKVFSEIEQYCPDHTIFATNTSSLSITKLAEATKRPEKFIGMHFFNPPKILKLLEIVWGEKTSEETIRIVEDFARKIDRIIIHVRKDVPGFIVNRIFVTMSNEASWAVEMGEGTIEEIDSAVKYRLGLPMGLFELHDVLGGGSVDVSYHVLEYYRQTLGESYRPSPLFERLFKAGHYGKKTGKGFYDWSEGKTNEVPLRAGANFDLLRLVAPAVNEAAWLIEKGVASAEEIDLAVLHGLNYPRGLLRMADDFGIDSIVKKLNELYEKYNGEERYKVNPVLQKMVEEGKLGRTTGEGFYKYGD
```

Query:

```
>AF0017_2 COG1024 # Protein_GI_number: 11497638 # Func_class: I Lipid transport and metabolism  # Function: Enoyl-CoA hydratase/carnithine racemase # Organism: Archaeoglobus fulgidus
GNYEFVKVEKEGKVGVLKLNRPRRANALNPTFLKEVEDALDLLERDEEVRAIVIAGEGKNFCAGADIAMFASGRPEMVTEFSQLGHKVFRKIEMLSKPVIAAIHGAAVGGGFELAMACDLRVMSERAFLGLPELNLGIIPGWGGTQRLAYYVGVSKLKEVIMLKRNIKPEEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPLAVKYLKKVIALGTMPALETGNLAESEAGAVIALTDDVAEGIQAFNYRRKPNFRGR
```

##### SAMH:

```
@HD	VN:1.4	SO:queryname
@SQ	SN:AF0017_1	LN:402
AF0017_2	0	AF0017_1	81	255	169S1=1X2=3X2=1X1=2X1=1X2I2X1=3X1=3X1=5X1=56S	*	GNYEFVKVEKEGKVGVLKLNRPRRANALNPTFLKEVEDALDLLERDEEVRAIVIAGEGKNFCAGADIAMFASGRPEMVTEFSQLGHKVFRKIEMLSKPVIAAIHGAAVGGGFELAMACDLRVMSERAFLGLPELNLGIIPGWGGTQRLAYYVGVSKLKEVIMLKRNIKPEEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPLAVKYLKKVIALGTMPALETGNLAESEAGAVIALTDDVAEGIQAFNYRRKPNFRGR	*	AS:i:37	NM:i:23	
```

##### SAM

```
AF0017_2	0	AF0017_1	81	255	169S1=1X2=3X2=1X1=2X1=1X2I2X1=3X1=3X1=5X1=56S	*	GNYEFVKVEKEGKVGVLKLNRPRRANALNPTFLKEVEDALDLLERDEEVRAIVIAGEGKNFCAGADIAMFASGRPEMVTEFSQLGHKVFRKIEMLSKPVIAAIHGAAVGGGFELAMACDLRVMSERAFLGLPELNLGIIPGWGGTQRLAYYVGVSKLKEVIMLKRNIKPEEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPLAVKYLKKVIALGTMPALETGNLAESEAGAVIALTDDVAEGIQAFNYRRKPNFRGR	*	AS:i:37	NM:i:23	
```

##### EMBOSS

```
AF0017_1            81 EVAKDADLVIEAIPE--IFDLKKKVFSEIEQYCP     112
                       |.||:..||.|..|:  .:|...|:..|:.:..|
AF0017_2           170 EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPP     203

Length: 34
Identity:        11/34 (32.4%)
Similarity:      17/34 (50.0%)
Gaps:             2/34 ( 5.9%)
Score: 37
```

##### SSW

```
target_name: AF0017_1
query_name: AF0017_2
optimal_alignment_score: 37	strand: +	target_begin: 81	target_end: 112	query_begin: 170	query_end: 203

Target:         81 EVAKDADLVIEAIPE--IFDLKKKVFSEIEQYCP     112
                   |*||***||*|**|*  **|***|***|*****|
Query:         170 EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPP     203
```


#### CSV Output

If a non-trace function is used, the output is a comma-separated values (CSV) file.  However, the exact output will depend on how the tool is used.  At minimum, there will be seven values per line in the file.

index1, index2, length1, length2, score, end_query, end_ref[, matches, similarities, length]

If a query file was used, then index1 is the query index and index2 is the database index.

If a query file was not used, then both index1 and index2 refer to the same single FASTA/FASTQ file that was supplied.

If a statistics-calculating function is used, for example 'sw_stats_striped_16', then the number of exact matches, similarities, and alignment length are also computed and returned.

### Generating a Homology Graph

The parasail_aligner already can take a FASTA- or FASTQ-formatted set of sequences and all of the sequences in the file will be compared against themselves.  If the 'edge' parameter (`-E`) or 'graph' parameter (`-G`) in combination with any of the statistics-calculating parasail routines is selected, this changes the output calculation.  The reason statistics must be calculated is that the output depends on them.  This application is used in a metagenomics workflow, creating a homology graph as output which is later processed by a community detection application.  The 'edges' in the graph consist of any highly similar pair of sequences such that their alignment meets certain criteria.  An 'edge' is only output if it meets the following criteria.

 * AOL = percent alignment length -- the alignment must cover at least XX percent of the longer sequence.
 * SIM = percent exact matches -- the alignment must contain at least XX percent exact character matches.
 * OS = percent optimal score -- the calculated score must be XX percent of the longer sequence's self score.

The 'edge' (`-E`) output is comma-separated values (CSV).

index1, index2, length/max_lengh, matches/length, score/self_score

The 'graph' (`-G`) output is a METIS graph file.
