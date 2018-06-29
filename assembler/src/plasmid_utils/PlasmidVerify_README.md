

# PlasmidVerify: plasmid contig verification tool

PlasmidVerify classifies contigs (output of metaplasmidSPAdes or other assemblers) as plasmidic or non-plasmidic. 

Contig is classified as plasmidic (chromosomal) if it contains (does not contain) plasmid-specific genes. Since plasmids harbor a large variety of genes, plasmidVerify uses a plasmid-specific profile-HMM database to detect remote similarities to known plasmid-specific genes.


## Installation

### Requirements

PlasmidVerify is a Python script, thus, installation is not required. However, it has following dependencies:

* Python 2.7+,
* Biopython,
* Prodigal,
* hmmsearch (from the hmmer package),
* Plasmid-specific HMMs (you can download it [here](/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/hmms.tar.gz )).

### Optional BLAST verification

You can verify your output by BLAST to check if you found novel plasmid. In this case you need to have blastn in your $PATH and provide path to local copy of nr database. 

## Usage 

    ./plasmidverify.py -f input_file.fasta 
            -o output_directory 
            -hmm path_to_plasmid_HMMs_folder 
            -b -db path_to_BLAST_database


Output file input_file_result_table.tsv contain tab-separated table with predicted plasmid-specific HMMs and verification result.

