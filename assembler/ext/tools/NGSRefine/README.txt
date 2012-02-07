#########################################################################
########################## NGS-Refine README ############################
#########################################################################

NGS-Refine is a positional reassembly tool that corrects errors in NGS assemblies.

To run NGS-Refine:
    1) prep.pl to pre-process the data (aligns reads to contigs using BWA and 
       prepaires the input-dir for NGS-Refine).
    2) NGSRefine.jar to correct errors in all contigs in input dir using positional 
       reassembly.

#########################################################################
################################# PREP ##################################
#########################################################################

    SYNOPSIS 
        prep.pl <reads1.fq> <reads2.fq> <contigs.fa> <out-dir> <min-len> <num-threads>

    INPUT 
        - reads1/reads2     paired-end reads from NGS (FASTQ/A)
        - contigs.fa        contigs from assembly (FASTA). Headers should be either 
                            Vlevet-style, or have 1st word be a unique ID (INT).
	
    OUTPUT
        Directory with what is needed for NGS-Refine to run:
        1) all reads aligned to all contigs (BWA must be pre-installed).
        2) individual contig files (FASTA).
        3) individual alignment file (SAM) per contig.

    NOTES
        The 'bwa' executable (v0.6.*) should be placed in PATH prior to running prep.pl

#########################################################################
############################## NGS-REFINE ###############################
#########################################################################

    SYNOPSIS
        java -Xmx16g -jar NGSRefine.jar -A PATH [OPTIONS...]

    INPUT
        PATH    path to the input-dir from prep.pl
	
    OPTIONS
        java -jar NGSRefine.jar for more information.

    OUTPUT
        Directory (NGSRefine_out_YYYMMDD) containing the refined contigs.fa file, and, 
        for each contig, an ID.improved.fa (where ID is the unique contig ID from the 
        original contigs file). Some contigs will be unchanged, where NGS-Refine simply 
        agreed with the assembled contig, or did not have enough information to refine. 
        Others will be refined. 

    NOTES:
        When using the -g option (for evaluation of the refinement when a reference exists), the 
        'blat' executable and 'blat_wrapper.pl' script must be placed in user's PATH.

