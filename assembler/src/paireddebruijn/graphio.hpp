#ifndef IOPROCEDURES_HPP_
#define IOPROCEDURES_HPP_
#include "common.hpp"
#include "pairedGraph.hpp"
using namespace paired_assembler;


void outputLongEdges(longEdgesMap &longEdges);
void outputLongEdgesThroughGenome(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount);
void codeRead(char *read, char *code);

inline bool nextReadPair(char * &read1, char * &read2) {
	return (scanf("%s %s", read1, read2)== 2);
}

ll extractMer(char *read, int shift, int length);

string decompress(ll a, int l);


class DataPrinter {
	FILE *f_;
public:
	DataPrinter(char *fileName);
	void outputInt(int a);
	void outputEdge(Edge *edge);
	void outputSequence(Sequence *sequence);
	void outputLongEdgesMap(longEdgesMap &map);
	void outputIntArray(int *array, int length);
	void outputIntArray(int *array, int length, int width);
	void close();
};

class DataReader {
	FILE *f_;
public:
	DataReader(char *fileName);
	void readInt(int &a);
	void readEdge(Edge * &edge);
	void readSequence(Sequence * &sequence);
	void readLongEdgesMap(longEdgesMap &map);
	void readIntArray(int *array, int length);
	void readIntArray(int *array, int length, int width);
	void close();
};

void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges, int &VertexCount, int EdgeId);

void load(char *fileName, PairedGraph &g, longEdgesMap &longEdges, int &VertexCount, int EdgeId);

#endif /* IOPROCEDURES_HPP_ */
