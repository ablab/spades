/* shuffle-fast - interleave two fasta or fastq files 
*
* Version 1.0 9/16/2009 by  Eric L. Cabot, Genome Center of Wisconsin 
*
*   Usage:  shuffleSeqs-fast [file1] [file2]
*
*   Example: ./shuffleSeqs-fast s_1_1.fas s_1_2.fas > bigger.fas
*
*   
*   Notes: 
*     1. file1 and file2 are either fasta or fastq 
*        formatted. The two need not be of the same type.
*
*     2. Output goes to stdout.
*   
*     3. Reading gzipped 
*
*    This program is based on H. Li's kseq_test.c, modified 
*    to mimic the behavior of Velvet's shuffleSequences-fasta.pl
*	   and d shuffleSequences-fastq.pl. 
*
*    Please see kseq.h for  Li's copyright notice.
**********************************************************************
*
*	Installation
*
*   1. Place this file and kseq.h in the same directory*
*
*   2. Make sure that zlib.h is available on your system
*      (there is a copy of zlib shipped with Velvet)
*
*   3. Typical complilation on linux/unix:
*    
*       gcc -g -O2 shuffleSeqs-fasta.c -lz -o shuffleSeqs-fast
*
**********************************************************************/

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

KSEQ_INIT(gzFile,gzread)


int main (int argc, char *argv[]) 
{
	gzFile  fpA,fpB;
  kseq_t  *seqA,*seqB;
  int     lA,lB;

	if (argc < 3 ) {
		fprintf(stderr,"Usage: %s <in.seq> <in.seq>\n", argv[0]);
		return 1;
	}


  fpA = gzopen(argv[1],"r");
	fpB = gzopen(argv[2],"r");

  if ( fpA == NULL || fpB == NULL ) {
		fprintf(stderr, "Could not open both input files\n");
    return 2;
	}

  seqA = kseq_init(fpA);
  seqB = kseq_init(fpB);
  
	while ((lA = kseq_read(seqA)) >= 0) {
		if (lB = kseq_read(seqB)) {
			if (seqA->qual.l) 
				printf("@%s\n%s\n+\n%s\n", seqA->name.s, seqA->seq.s, seqA->qual.s);
			else
				printf(">%s\n%s\n", seqA->name.s, seqA->seq.s);
						
			if (seqB->qual.l)
				printf("@%s\n%s\n+\n%s\n", seqB->name.s, seqB->seq.s, seqB->qual.s);
			else
				printf(">%s\n%s\n", seqB->name.s, seqB->seq.s);
			
		}   /*  if   */
	}    /* while */

	kseq_destroy(seqA);
	kseq_destroy(seqB);
	
	gzclose(fpA);
	gzclose(fpB);
	return 0;
}
