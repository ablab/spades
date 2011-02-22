#include <vector>

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

class Arc {
	int _coverage;
	Vertex* _head;
public:
	Arc(int coverage, Vertex* head) : _coverage(coverage), _head(head)
	{}

	int coverage() {return _coverage;};
	Vertex* head() {return _head;};
};

class Vertex {
	int _coverage;
	Sequence *_upper;
	Sequence *_lower;
	int _neighbours_count;
	Vertex* _neighbours;
	short _delta_d;
	Vertex *real_vertex;
public:
	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);

	int coverage() {return _coverage;};

	int neighbours_count() {return _neighbours_count;};

	int addEdge(Vertex *neighbour);

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

};
#endif /* CONDENSED_GRAPH_H_ */
