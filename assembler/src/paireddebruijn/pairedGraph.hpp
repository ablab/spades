#include "vector"
#include "sequence.hpp"
//#include "hashTable.h"
using namespace std;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

namespace paired_assembler {

//typedef int Kmer;

//class Sequence {
//	char *_nucleotides;
//	short _length;
//public:
//	Sequence(char *nucleotides, short length) : _nucleotides(nucleotides), _length(length)
//	{}
//	char operator[](const int &index);
//};

class Vertex;

//struct Arc {
//	short _coverage;
//	Vertex* _head;
//};
class VertexPrototype {
	public:
	VertexPrototype(Sequence *lower_, int start_) {
		lower = lower_;
		start = start_;
		used = false;
	}
	Sequence *lower;
	int start;
//	int finish;
	bool used;
};

class Vertex {
	//	int _coverage;
	Sequence *upper;
	Sequence *lower;
	vector<Vertex*> neighbours;
	//	Arc* _neighbours;
	//	int _neighbours_count;
	//	short _delta_d;
	//	Vertex *real_vertex;
public:
	//	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);
	Vertex(Sequence *up, Sequence *low) {
		upper = up;
		lower = low;
	}

	void glue(Sequence *up, Sequence *low, int glueDepth);
	//
	//	int coverage() {return _coverage;};
	//
	//	int neighbours_count() {return _neighbours_count;};
	//
	//	int addEdge(Vertex *neighbour, short coverage);
	//
	//	vector<Vertex*> getEdges();
	//
	//	Kmer getKmer(int position);

};

class Graph {
	//	HashTable map;
	//
	//	int merge(Vertex *u, Vertex* v);
	//
	//	int split(Vertex *u, Vertex* v, short position);
public:
	vector<Vertex *> vertices;
	int addVertex(Sequence *upper, Sequence *lower) {
		Vertex *newVertex = new Vertex(upper, lower);
		vertices.push_back(newVertex);
		return vertices.size()-1;
	}
};

//class GraphIterator {
//public:
//	GraphIterator(Graph *graph);
//
//	Vertex *nextVertex();
//
//	bool hasNext();
//};

}
#endif /* CONDENSED_GRAPH_H_ */
