#include "AlternativeEdge.h"
#include "ParseTitle.h"
#include "utils.h"
#include <assert.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include "SeqReader.h"

void AppendAlternativeEdges(AlternativeEdgeVect &src, 
														ssize_t offset,
														AlternativeEdgeVect &dest) {
	if (src.size() == 0)
		return;
	ssize_t d = dest.size();
	dest.resize(src.size() + dest.size());
	ssize_t s;
	for (s = 0; s < src.size(); s++, d++) {
		dest[d].offset = offset;
		dest[d].seq = src[s].seq;
	}
	ssize_t i;
	for (i = 0; i < dest.size(); i++) {
		assert(dest[i].offset != -1);
	}
	
}

void AppendAlternativeEdge(AlternativeEdgeVect &src,
													 SimpleSequence &seq,
													 ssize_t offset) {
	ssize_t p;
	p = src.size();
	src.resize(p+1);
	src[p].seq.seq = new unsigned char[seq.length];
	memcpy(src[p].seq.seq, seq.seq, seq.length);
	src[p].seq.length = seq.length;
	src[p].offset = offset;
	ssize_t i;
	for (i = 0; i < src.size(); i++) {
		assert(src[i].offset != -1);
	}
}

void WriteAltEdges(AlternativeEdgeVect &vect,
									 ssize_t edgeIndex, ofstream &seqOut) {
	
	ssize_t e;
	stringstream titleStrm;
	for (e = 0; e < vect.size(); e++ ){
		titleStrm.str("");
		titleStrm << e << " DestEdge=" << edgeIndex << " Offset=" 
							<< vect[e].offset;
		vect[e].seq.PrintSeq(seqOut, titleStrm.str());
	}
}


ssize_t ReadAltEdge(ifstream &seqIn,
								DNASequence &seq,
								ssize_t &destEdge,
								ssize_t &offset) {
	ssize_t status = 0;
	if (SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
		if (ParseKeyword(seq.namestr, "DestEdge", destEdge) and
				ParseKeyword(seq.namestr, "Offset", offset)) {
			status = 1;
		}
	}	
	return status;
}		

