#include "vector"
#include "sequence.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"
//#include "hashTable.h"
using namespace std;

#ifndef PAIREDGRAPH_H_
#define PAIREDGRAPH_H_

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
	VertexPrototype(Sequence *lower_, int start_, int coverage_=1) {
		lower = lower_;
		VertexId = start_;
		used = false;
		coverage = coverage_;
	}
	Sequence *lower;
	int VertexId;
	bool used;
	int coverage;
};

class EdgePrototype {
public:
	EdgePrototype(Sequence *lower_, int start_, int coverage_=1) {
		lower = lower_;
		VertexId = start_;
		used = false;
		coverage = coverage_;
	}
	Sequence *lower;
	int VertexId;
	bool used;
	int coverage;
};
/*
 * length- including one vertex.
 * upper and lower including both in and out vertex
 *
 *
 *
 */
class Edge {
	//	int _coverage;
public:
	Sequence *upper;
	Sequence *lower;
	int length;
	int FromVertex;
	int ToVertex;
	int EdgeId;
	int coverage;
	Edge(Edge &e){
		length = e.length;
		FromVertex = e.FromVertex;
		ToVertex = e.ToVertex;
		EdgeId=e.EdgeId;
		coverage= e.coverage;
		upper = new Sequence(const_cast<Sequence&> (*e.upper));
		lower = new Sequence(const_cast<Sequence&> (*e.lower));
	}
	//	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);
	void ExpandRight(Edge &newRigth) {
		ToVertex = newRigth.ToVertex;
		if (newRigth.length > 0) {
			length = length + newRigth.length;
			string toOut = newRigth.upper->str();
			assert(k-1 < toOut.length());
			upper = new Sequence(
					upper->str() + newRigth.upper->Subseq(k - 1).str());

			toOut = newRigth.lower->str();
			assert(l-1 < toOut.length());
			lower = new Sequence(
					lower->str() + newRigth.lower->Subseq(l - 1).str());

		}
	}
	void shortenEdge(int toCut, int direction){
		if (toCut < length){
			if (direction == OUT_EDGE) {
				upper = new Sequence(upper->Subseq(toCut, length));
				lower = new Sequence(lower->Subseq(toCut, length));
			} else {
				upper = new Sequence(upper->Subseq(0, upper->size() - toCut));
				lower = new Sequence(lower->Subseq(0, lower->size() - toCut));
			}
			length -= toCut;
		}
	}
	void ExpandLeft(Edge &newLeft) {
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
	Edge(Sequence *up, Sequence *low, int from, int to, int len, int id, int cov = 1) {
		upper = up;
		lower = low;
		FromVertex = from;
		ToVertex = to;
		length = len;
		EdgeId = id;
		coverage = cov;
	}

	~Edge() {
//		cerr << "destructing" << upper->str() << endl;
		if (upper != lower) {
			delete upper;
			delete lower;
		} else
			delete upper;
	}
};

inline int edgeRealId(int id, longEdgesMap &longEdges) {
	int res = id;
	while (longEdges[res]->EdgeId != res) {
		res = longEdges[res]->EdgeId;
	}
	return res;
}

class PairedGraph {
public:
	//0 - in-degrees
	//1 -out-degrees
	int degrees[MAX_VERT_NUMBER][2];//, outD[MAX_VERT_NUMBER][2];
	int edgeIds[MAX_VERT_NUMBER][MAX_DEGREE][2];
	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
	longEdgesMap longEdges;
	verticesMap verts;
	int VertexCount;
	int EdgeId;
	PairedGraph(){
		cerr<<"VAH Paired created"<<endl;
	}
	void RebuildVertexMap(void);

};

int storeVertex(gvis::GraphPrinter<int> &g, PairedGraph &graph, ll newKmer, Sequence* newSeq);
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq);
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, int VertNum);
void resetVertexCount(PairedGraph &graph);


}
#endif /* PAIREDGRAPH_H_ */
