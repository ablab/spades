#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"

void constructGraph(PairedGraph *g, char *inputFile);

void clusterReadPaires(PairedGraph *g, char *tmpFile);

void threadReads(PairedGraph *g, char *inputFile);

#endif /*GRAPHCONSTRUCTION_H_*/
