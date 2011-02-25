#include "simpleGraph.hpp"
#include "constructHashTable.hpp"

char *tmpFile = "data/reads.out";

void clusterReadPaires(PairedGraph *g, char *tmpFile) {

}

void threadReads(PairedGraph *g, char *inputFile) {

}

void constructGraph(PairedGraph *g, char *inputFile) {
	readsToPairs(inputFile, tmpFile);
	clusterReadPaires(g, tmpFile);
	threadReads(g, inputFile);
}

