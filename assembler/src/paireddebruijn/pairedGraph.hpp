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
	VertexPrototype(Sequence *lower_, int start_, int coverage_ = 1) {
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
	EdgePrototype(Sequence *lower_, int start_, int coverage_ = 1) {
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
	Edge(Edge &e) {
		length = e.length;
		FromVertex = e.FromVertex;
		ToVertex = e.ToVertex;
		EdgeId = e.EdgeId;
		coverage = e.coverage;
		upper = new Sequence(const_cast<Sequence&> (*e.upper));
		lower = new Sequence(const_cast<Sequence&> (*e.lower));
	}
	//	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);
	void ExpandRight(Edge &newRight) {
		ToVertex = newRight.ToVertex;
		if (newRight.length > 0) {
			length = length + newRight.length;
			string toOut = newRight.upper->str();
			assert(k-1 < toOut.length());
			upper = new Sequence(
					upper->str() + newRight.upper->Subseq(k - 1).str());

			toOut = newRight.lower->str();
			assert(l-1 < toOut.length());
			lower = new Sequence(
					lower->str() + newRight.lower->Subseq(l - 1).str());

		}
	}
	void shortenEdge(int toCut, int direction) {
		if (toCut < length) {
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
	Edge(Sequence *up, Sequence *low, int from, int to, int len, int id,
			int cov = 1) {
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

template<typename tVertex>
class IVertexIterator {
	virtual bool hasNext() = 0;
	virtual tVertex next() = 0;
};

template<typename tVertex, typename tEdge>
class IPairedGraph {
public:
	/*Better way to do it is to implement begin() and end() functions which would return pointers to start and
	 * end of lists of neighbours but it is not possible because of the order of dimentions in edgeIds and
	 * because those are edge ids in array instead of pointers
	 */
	virtual int rightDegree(tVertex vertex) = 0;
	virtual int leftDegree(tVertex vertex) = 0;
	int degree(tVertex vertex, int direction) {
		if (direction == RIGHT)//RIGHT = 1
			return rightDegree(vertex);
		else if (direction == LEFT)
			return leftDegree(vertex);//LEFT = -1
	}
	virtual tEdge rightNeighbour(tVertex vertex, int number) = 0;
	virtual tEdge leftNeighbour(tVertex vertex, int number) = 0;
	tEdge neighbour(tVertex vertex, int number, int direction) {
		if (direction == RIGHT)//RIGHT = 1
			return getRightNeighbour(vertex);
		else if (direction == LEFT)
			return getRightNeighbour(vertex);//LEFT = -1
	}

	virtual IVertexIterator<tVertex> *vertexIterator() = 0;

	//In order to add edge to graph one should create this edge first!
	virtual void addEdge(tEdge newEdge) = 0;
	virtual void removeEdge(tEdge edge) = 0;

	virtual void addVertex(tVertex vertex) = 0;
	//add adjecent edges should be removed as well
	virtual void removeVertex(tVertex vertex) = 0;
	virtual tEdge merge(tEdge edge1, tEdge edge2) = 0;
	virtual pair<tEdge, tEdge> splitEdge(tEdge edge, int position) = 0;

	//glue edges, there start and end vertices
	virtual tEdge glueEdges(tEdge edge1, tEdge edge2) = 0;
	//seperate edges  adjecent to the vertex
	virtual pair<tEdge, tEdge> unGlueEdgesLeft(tVertex vertex) = 0;
};

class PairedGraphData {
public:
	//0 - in-degrees
	//1 -out-degrees
	int degrees[MAX_VERT_NUMBER][2];//, outD[MAX_VERT_NUMBER][2];
	int edgeIds[MAX_VERT_NUMBER][MAX_DEGREE][2];
	//	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
	longEdgesMap longEdges;verticesMap verts;
	int VertexCount;
	int EdgeId;
	PairedGraphData() {
		cerr << "VAH Paired created" << endl;
	}
	//	void RebuildVertexMap(void);
};

class VertexIterator: public IVertexIterator<int> {
	friend class PairedGraphData;
private:
	int currentVertex_;
	PairedGraphData *graph_;
public:
	VertexIterator(PairedGraphData *graph) {
		currentVertex_ = 0;
	}

	virtual bool hasNext() {
		while (currentVertex_ < graph_->VertexCount
				&& graph_->degrees[currentVertex_][0]
						+ graph_->degrees[currentVertex_][1]
						== 0) {
			currentVertex_++;
		}
		return currentVertex_ < graph_->VertexCount;
	}

	virtual int next() {
		assert(!hasNext());
		return currentVertex_;
	}
};

class PairedGraph: public PairedGraphData, public IPairedGraph<int, Edge *> {
public:
	virtual int rightDegree(int vertex) {
		return degrees[vertex][1];
	}
	virtual int leftDegree(int vertex) {
		return degrees[vertex][0];
	}
	virtual Edge *rightNeighbour(int vertex, int number) {
		assert(number < degrees[vertex][1]);
		return longEdges[edgeRealId(edgeIds[vertex][number][1], longEdges)];
	}
	virtual Edge *leftNeighbour(int vertex, int number) {
		assert(number >= degrees[vertex][0]);
		return longEdges[edgeRealId(edgeIds[vertex][number][0], longEdges)];
	}

	//This is very bad method!!!!
	virtual VertexIterator *vertexIterator() {
//		return new VertexIterator(this);
	}

	virtual void addEdge(Edge *newEdge) {
		longEdges.insert(make_pair(newEdge->EdgeId, newEdge));
		edgeIds[newEdge->FromVertex][degrees[newEdge->FromVertex][1]][1]
				= newEdge->EdgeId;
		degrees[newEdge->FromVertex][1]++;
		edgeIds[newEdge->FromVertex][degrees[newEdge->FromVertex][0]][0]
				= newEdge->EdgeId;
		degrees[newEdge->FromVertex][0]++;
	}

	virtual void removeEdge(Edge *edge) {
		assert(false);
	}

	virtual void addVertex(int vertex) {
		assert(false);
	}

	virtual void removeVertex(int vertex) {
		assert(false);
	}

	virtual Edge * merge(Edge *edge1, Edge *edge2) {
		assert(false);
	}

	virtual pair<Edge *, Edge *> splitEdge(Edge *edge, int position) {
		assert(false);
	}

	//glue edges, there start and end vertices
	virtual Edge *glueEdges(Edge *edge1, Edge *edge2) {
		assert(false);
	}
	//seperate edges  adjecent to the vertex
	virtual pair<Edge *, Edge *> unGlueEdgesLeft(int vertex) {
		assert(false);
	}

	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
	void RebuildVertexMap(void);
};


int storeVertex(gvis::GraphPrinter<int> &g, PairedGraph &graph, ll newKmer,
		Sequence* newSeq);
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq);
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, int VertNum);
void resetVertexCount(PairedGraph &graph);

}
#endif /* PAIREDGRAPH_H_ */
