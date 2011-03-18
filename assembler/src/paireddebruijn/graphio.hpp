#ifndef IOPROCEDURES_HPP_
#define IOPROCEDURES_HPP_
#include "common.hpp"
#include "pairedGraph.hpp"
using namespace paired_assembler;


void codeRead(char *read, char *code);

void outputLongEdges(longEdgesMap &longEdges);

inline bool nextReadPair(char * &read1, char * &read2) {
	return (scanf("%s %s", read1, read2)== 2);
}

ll extractMer(char *read, int shift, int length);

string decompress(ll a, int l);

#endif /* IOPROCEDURES_HPP_ */
