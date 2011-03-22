#include "vector"
#include "sequence.hpp"
#include "common.hpp"
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
		VertexId = start_;
		used = false;
	}
	Sequence *lower;
	int VertexId;
	bool used;
};

class Edge {
	//	int _coverage;
public:
	Sequence *upper;
	Sequence *lower;
	int length;
	int FromVertex;
	int ToVertex;
	int EdgeId;
	//	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);
	void ExpandRight(Edge newRigth) {
		ToVertex = newRigth.ToVertex;
		if (newRigth.length > 0) {
			length = length + newRigth.length;
			string toOut = newRigth.upper->str();
			//		cerr <<endl << "str:: "<< upper->str() <<" "<< toOut <<" "<< k <<endl;
			assert(k-1 < toOut.length());
			upper = new Sequence(
					upper->str() + newRigth.upper->Subseq(k - 1).str());

			toOut = newRigth.lower->str();
			//		cerr <<endl << l-1 <<" "<< toOut.length()<<" "<< toOut<< endl;
			assert(l-1 < toOut.length());
			//		cerr <<endl << "strL:: "<< lower->str() <<" "<< toOut <<" "<< l <<endl;
			lower = new Sequence(
					lower->str() + newRigth.lower->Subseq(l - 1).str());
			//		cerr<< endl << "expanded" << endl;
		}
	}
	void ExpandLeft(Edge newLeft) {
		FromVertex = newLeft.FromVertex;
		if (newLeft.length > 0) {
			length = length + newLeft.length;
			upper = new Sequence(
					newLeft.upper->Subseq(0, newLeft.length).str()
							+ upper->str());
			lower = new Sequence(
					newLeft.lower->Subseq(0, newLeft.length).str()
							+ lower->str());
		}
	}
	Edge(Sequence *up, Sequence *low, int from, int to, int len, int id) {
		upper = up;
		lower = low;
		FromVertex = from;
		ToVertex = to;
		length = len;
		EdgeId = id;
	}
};

inline int edgeRealId(int id, longEdgesMap &longEdges){
	int res = id;
	while (longEdges[res]->EdgeId !=res){
		res = longEdges[res]->EdgeId;
	}
	return res;
}

class PairedGraph {
public:
	//0 - in-degrees
	//1 -out-degrees
	int degrees[MAX_VERT_NUMBER][2];//, outD[MAX_VERT_NUMBER][2];
	int outputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
	int inputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
//	int firstDiff[Msv]
};
//
//class Vertex {
//	//	int _coverage;
//	Sequence *upper;
//	Sequence *lower;
//	vector<Vertex*> neighbours;
//	//	Arc* _neighbours;
//	//	int _neighbours_count;
//	//	short _delta_d;
//	//	Vertex *real_vertex;
//public:
//	//	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);
//	Vertex(Sequence *up, Sequence *low) {
//		upper = up;
//		lower = low;
//	}
//
//	void glue(Sequence *up, Sequence *low, int glueDepth);
//	//
//	//	int coverage() {return _coverage;};
//	//
//	//	int neighbours_count() {return _neighbours_count;};
//	//
//	//	int addEdge(Vertex *neighbour, short coverage);
//	//
//	//	vector<Vertex*> getEdges();
//	//
//	//	Kmer getKmer(int position);
//
//};
//
//class Graph {
//	//	HashTable map;
//	//
//	//	int merge(Vertex *u, Vertex* v);
//	//
//	//	int split(Vertex *u, Vertex* v, short position);
//public:
//	vector<Vertex *> vertices;
//	int addVertex(Sequence *upper, Sequence *lower) {
//		Vertex *newVertex = new Vertex(upper, lower);
//		vertices.push_back(newVertex);
//		return vertices.size()-1;
//	}
//};
//
////class GraphIterator {
////public:
////	GraphIterator(Graph *graph);
////
////	Vertex *nextVertex();
////
////	bool hasNext();
////};

}
#endif /* CONDENSED_GRAPH_H_ */
