# plasmidVerify: plasmid contig verification tool

plasmidVerify classifies contigs (output of metaplasmidSPAdes or other assemblers) as plasmidic or non-plasmidic, based on gene content. 


plasmidVerify predicts genes in the contig using Prodigal in the metagenomic mode, runs hmmsearch on the predicted proteins, 
and classifies the contig as plasmidic or chromosomal by applying the Naive Bayes classifier (NBC). 
For the set of predicted HMMs, plasmidVerify uses trained NBC to classify this set to be plasmidic or chromosomal. 


## Installation

### Requirements

plasmidVerify is a Python script, thus, installation is not required. However, it has following dependencies:

* Python 2.7+,
* Prodigal (https://github.com/hyattpd/Prodigal, available via conda),
* hmmsearch (from the hmmer package, http://hmmer.org/download.html),
* Plasmid- and chromosome-specific HMMs database (you can download it here:<LINK>  . Also you can use recent release of Pfam-A database, ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/)

In order to work properly, plasmidVerify require Prodigal and hmmsearch in your PATH environment variable.


### Optional BLAST verification

You can verify your output by BLAST to check if you found novel plasmid. In this case you need to have blastn in your $PATH and provide path to local copy of nr database. 

## Usage 

    ./plasmidverify.py 
            -f Input fasta file
            -o output_directory 
            --hmm HMM    Path to Pfam-A HMM database

            Optional arguments:
            -h, --help  Show the help message and exit
            -b          Run BLAST on input contigs
            --db DB      Path to BLAST db


Output file input_file_result_table.csv contains comma-separated table with predicted HMMs and verification result.
