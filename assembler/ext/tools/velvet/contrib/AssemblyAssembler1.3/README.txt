NAME
====

AssemblyAssembler1.3.py

VERSION
=======

Version 1.3

LICENCE
=======

Copyright 2010 - Jacob Crawford - Cornell University.
	
jc598@cornell.edu

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
    
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
        
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.


INTRODUCTION
============

The AssemblyAssembler (AA) is designed to automate a directed series of 
assemblies using the Velvet assembler (Daniel Zerbino, EBI UK).  The 
assembly routine was initially developed in an attempt to improve de novo 
transcriptome assemblies using Velvet (prior to the release of Oases 
(Daniel Zerbino, EBI UK)).  The AA conducts Velvet assemblies using 
default parameter values across a user specified range of kmer values.  
A coverage cutoff can be included if desired (see description of -c). 
It then determines which kmer values produced the best assembly (largest 
max contig) and conducts additional assemblies in that region of the 
kmer-parameter space.  Lastly, the AA takes contigs from all previous 
assemblies and uses them as input for a final assembly.  The final contig 
set is stored in Finalcontigs.fa and all Velvet operations recorded in 
GrandVelvetLog.txt, both stored in the directory FinalDir.  


DISCLAIMER
==========
This script comes with NO GUARANTEE. Misassemblies are possible in Velvet
assemblies and are possible using this script.  It is up to the user to 
determine the quality of the assembly.  

PREREQUISITES
=============

Python >= 2.5
Python modules sys and os

COMMAND LINE
============
	
	AssemblyAssembler.py -s INT -e INT -v /path/to/velvetdir -f 
	'velveth input line'
  
  	-s	[Required] The starting (lowest) hash value (no default). 
		Enter integer value (no quotes).  It is not recommended 
		that you attempt small values (< k=15) for this parameter.
  	
  	-e	[Required] The end (highest) hash value (no default). 
		Enter integer value (no quotes).  The difference between 
		the highest and lowest hash values must be at least 16.
  		I found it useful to set this parameter to (length of
		read-1).
  	
  	-v	[Required] The full path to the Velvet directory containing 
		velveth and velvetg binaries (executable files). 
		Example: /programs/velvet_0.7.61 (no quotes !! ).
  			
  	-f	[Required] The velveth parameters you would enter for a 
		standard velveth run.  This entry MUST be INSIDE QUOTES.
 		You MUST enter 1) the input file type (e.g. -fasta or 
		-fastq), 2) the read type (e.g. -short or =shortPaired) and
 		3) the full path to the sequence file(s). Example: ‘-fastq 
		-short /myfiles/Illumina/s_8_1_sequence.txt’).

	-c	[Optional] The coverage cutoff value to be used if desired.  
		It was noted by Dieter Best that including -cov_cutoff when 
		assembling at higher kmer values provides a better
 		environment for assembly of highly expressed genes, but
 		lowly expressed genes best assembled with no -cov_cutoff at 
		lower kmer values.  Thus, the -cov_cutoff must be
 		accompanied by a minimum threshold kmer value (-t) above 
		which the -cov_cutoff will be included in the velvetg runs.
 		You can enter either an integer value or auto if you'd like
 		Velvet to estimate a coverage cutoff, although the auto
 		version is not recommended for transcriptome assemblies.
	
	-t	[Optional] Minimum kmer value to include -cov_cutoff if
 		desired.  See description of -c.   
  	
	-i	[Optional] Insert length for paired-reads.  This parameter
 		is required to activate the use of paired reads during
 		assembly.  You must also have flagged the reads as paired
 		in your Velveth call line (-f)

	-m	[Optional] Expected kmer coverage.  If you prefer, you can
 		enter 'auto' (no quotes) to have Velvet estimate it from
 		the data.  This parameter is required to activate the use
 		of paired reads during assembly.  You must also have
 		flagged the reads as paired in your Velveth call line (-f)

	-r 	[Optional] enter either 'y' or 'n' (no quotes) to indicate
 		whether you would like short reads to be included in the
 		final summary assembly (either as paired or unpaired as
 		indicated in Velveth call -f)
		

EXAMPLES
========

Assemble a lane of Illumina single-end reads by searching across the range
k = 19 to k = 35:

	AssemblyAssembler1.2.py -s 19 -e 35 -v /programs/velvet_0.7.61 -f 
	'-short -fastq s_8_1_sequence.txt'

Assemble two lanes of Illumina paired-end reads by searching across the
range k = 19 to k = 59:

	AssemblyAssembler1.2.py -s 19 -e 35 -v /programs/velvet_0.7.61 -f 
	‘-fastq -shortPaired s_7_shuffled.txt -shortPaired2 
	s_8_shuffled.txt’ -i 250 -m auto

Assemble two lanes of Illumina paired-end reads by searching across the range
k = 19 to k = 59 and include a coverage cutoff of 8 at kmer values above 31:

	AssemblyAssembler1.2.py -s 19 -e 35 -v /programs/velvet_0.7.61 -f 
	‘-fastq -shortPaired s_7_shuffled.txt -shortPaired2
 	s_8_shuffled.txt’ -c 8 -t 31 -i 250 -m auto

Assemble one paired-end Illumina lane by searching across the range k = 21
to k = 59 with coverage cutoff of 8 initiated at kmer values above 45 and 
raw short reads included in the final summary assemblies. 

	AssemblyAssembler1.2.py -s 21 -e 59 -v /home/fs01/jc598/data/
	velvet_0.7.58 -f "-fastq -shortPaired /tmp/jc598/
	JC_paired_8both_clean.txt" -c 8 -t 45 -i 250 -m auto -r y

BUGS
=======

* None that I am aware of.

CONTACT
=======

Jacob Crawford <jc598@cornell.edu>

CHANGE LOG
=======

Version 1.2
- added check to see make sure user is using Python version 2.5 or above- added feature to allow -cov_cotoff to be implemented in specific kmer range (Per Dieter Best's recommendation)- added switch to trigger use of paired-end reads- added the ability to include raw short reads in the final summary assemblies. Joe Fass found this feature may be useful for whole genome assemblies. 

Version 1.3
- Included default value ('n') for the -r parameter to make it truly an optional parameter.  Thanks to Anar Khan for pointing this out. 
 
TROUBLESHOOTING
=======

If you get a 'Permission Denied' error, change the permissions associated with 
the AssemblyAssembler1.2.py file by entering the directory containing the file
and typing…

chmod 777 AssemblyAssembler1.2.py 

