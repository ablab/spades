README_velvetsc.TXT

VELVET-SC SOURCE 
March 11, 2011
Hamidreza Chitsaz

> SUMMARY
	* A/ Credits
	* B/ Velvet-SC documentation
	* C/ Usage

----------------------------------------------------------------------------------
A/ Credits

Hamidreza Chitsaz, UCSD made changes to Velvet version 0.7.62 to transform it into Velvet-SC version 0.7.62. Velvet-SC source code is distributed under the same license as the original Velvet source code. 

Reference

Hamidreza Chitsaz, Joyclyn L. Yee-Greenbaum, Glenn Tesler, Mary-Jane Lombardo, Christopher L. Dupont, Jonathan H. Badger, Mark Novotny, Douglas B. Rusch, Louise J. Fraser, Niall A. Gormley, Ole Schulz-Trieglaff, Geoffrey P. Smith, Dirk J. Evers, Pavel A. Pevzner, Roger S. Lasken.
De Novo Assembly of Bacterial Genomes from Single Cells.
Submitted. 2011. 	 

----------------------------------------------------------------------------------
B/ Velvet-SC documentation

Our changes to Velvet, that constitute Velvet-SC, are specified by
/* Added by Hamidreza Chitsaz, UCSD, Date */
/* Begin */
...
/* End */
in the source code. In particular, file run2.c contains the main algorithm of Velvet-SC; there are three constants defined in run2.c: VELVETSC_STARTING_COV, VELVETSC_INC_COV, and VELVETSC_MAX_LENGTH that specify respectively the starting average coverage, increment of average coverage in each pruning iteration, and the maximum length of contigs subject to removal if their average coverage is below the threshold. Please refer to the above reference for details.

----------------------------------------------------------------------------------
C/ Usage

Velvet-SC does not have any specific command line options/parameters. However, the meaning of -cov_cutoff in Velvet is different from that in Velvet-SC. Velvet-SC iteratively eliminates those contigs that are shorter than VELVETSC_MAX_LENGTH and whose average coverage is less than covCutoff, starting with covCutoff = VELVETSC_STARTING_COV, incrementing covCutoff by VELVETSC_INC_COV until covCutoff is more than the value specified by '-cov_cutoff' command line parameter. If '-cov_cutoff auto' is used, the value is computed automatically by the Velvet's coverage estimation algorithm. 
NOTE: for non-uniform coverage assembly, such as single cell genome sequencing, '-cov_cutoff auto' is particularly discouraged. A higher value in '-cov_cutoff' gives larger contigs, less total assembled nucleotides, and more misassemblies. Manual tuning of '-cov_cutoff' in the range 0.1-0.5 times the estimated coverage is expected to yield the best result for MDA amplified Illumina GAIIx reads.


