/***************************************************************************
 * Title:          FilterLowComplexityReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {

	std::string readsInName, readsOutName, readsDiscardName;
	std::ofstream readsOut, discardOut;

	if (argc < 2) {
		std::cout << "usage: filterLCReads readsIn readsOut [composition=0.8]" << std::endl;
		std::cout << "        Remove any reads if one or two nucleotides makes up more than" << std::endl
							<< "        'composition' percent of the read." << std::endl;
	}

	double comp = 0.8;
	readsInName = argv[1];
	readsOutName = readsInName + ".hc";
	readsDiscardName = readsInName + ".lc";
	if (argc > 2) {
		comp = atof(argv[2]);
	}

	DNASequenceList reads;
	ReadDNASequences(readsInName, reads);
	openck(readsOutName, readsOut, std::ios::out);
	openck(readsDiscardName, discardOut, std::ios::out);
	ssize_t actg[4];
	ssize_t r;
	ssize_t numFiltered = 0;
	ssize_t p;
	ssize_t lc;
	for (r = 0; r < reads.size(); r++) {
		actg[0] = actg[1] = actg[2] = actg[3] = 0;
		for (p = 0; p < reads[r].length; p++) {
			actg[(unsigned char) unmasked_nuc_index[reads[r].seq[p]]]++;
		}

		lc = 0;
		ssize_t n, n2;
		/*		for (n = 0; n < 4 and !lc; n++ ){
			std::cout << " " << actg[n] ;
		}

		std::cout << std::endl;
		*/		
		for (n = 0; n < 4 and !lc; n++ ){

			if (double(actg[n])/reads[r].length > comp) {
				lc= 1;
				break;
			}
		}
		for (n = 0; n < 3 and !lc; n++) {
			for (n2 = n+1; n2 < 4; n2++) {
				if ((double(actg[n] + actg[n2]))/reads[r].length > comp) {
					lc = 1;
					break;
				}
			}
		}

		if (!lc) {
			reads[r].PrintSeq(readsOut);
			readsOut << std::endl;
		}
		else {
			reads[r].PrintSeq(discardOut);
			discardOut << std::endl;
			numFiltered++;
		}			
	}
	std::cout << "removed: " << numFiltered << " reads." << std::endl;
	return 0;
}

