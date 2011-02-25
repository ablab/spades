#include <vector>
#include "hashTable.h"
using namespace std;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

namespace paired_assembler {

typedef int Kmer;

class Sequence {
	char *_nucleotides;
	short _length;
public:
	Sequence(char *nucleotides, short length) : _nucleotides(nucleotides), _length(length)
	{}
	char operator[](const int &index);
};

class Vertex;

struct Arc {
	short _coverage;
	Vertex* _head;
};

class Vertex {
	int _coverage;
	Sequence *_upper;
	Sequence *_lower;
	int _neighbours_count;
	Arc* _neighbours;
	short _delta_d;
	Vertex *real_vertex;
public:
	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);

	int coverage() {return _coverage;};

	int neighbours_count() {return _neighbours_count;};

	int addEdge(Vertex *neighbour, short coverage);

	void getEdges(Arc *begin, Arc *end);

	Kmer getKmer(int position);
};

class Graph {
	HashTable map;

	int merge(Vertex *u, Vertex* v);

	int split(Vertex *u, Vertex* v, short position);
};

class GraphIterator {
public:
	GraphIterator(Graph *graph);

	Vertex *nextVertex();

	bool hasNext();
};

}
#endif /* CONDENSED_GRAPH_H_ */
