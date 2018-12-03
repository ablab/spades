

# plasmidVerify: plasmid contig verification tool

plasmidVerify classifies contigs (output of metaplasmidSPAdes or other assemblers) as plasmidic or non-plasmidic, based on gene content. 


plasmidVerify predicts genes in the contig using Prodigal in the metagenomic mode, runs hmmsearch on the predicted proteins, 
and classifies the contig as plasmidic or chromosomal by applying the Naive Bayes classifier (NBC). 
For the set of predicted HMMs, plasmidVerify uses trained NBC to classify this set to be plasmidic or chromosomal. 


## Installation

### Requirements

plasmidVerify is a Python script, thus, installation is not required. However, it has following dependencies:

* Python 2.7+,
* Biopython,
* Prodigal,
* hmmsearch (from the hmmer package),
* Pfam-A database (check most recent version - http://pfam.xfam.org, FTP tab)

### Optional BLAST verification

You can verify your output by BLAST to check if you found novel plasmid. In this case you need to have blastn in your $PATH and provide path to local copy of nr database. 

## Usage 

    ./plasmidverify.py -f input_file.fasta 
            -o output_directory 
            -hmm path_to_plasmid_HMMs_folder 
            -b -db path_to_BLAST_database


Output file input_file_result_table.csv contain comma-separated table with predicted HMMs and verification result.

