#include "vector"
#include "sequence.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"
#include "logging.hpp"
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
	ll upper;
	Sequence *lower;
	int VertexId;
	bool used;
	int coverage;
	VertexPrototype(Sequence *lower_, int id, int coverage_ = 1);
	VertexPrototype(ll upper_, Sequence *lower_, int id, int coverage_ = 1);
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
		else if (direction == LEFT)//LEFT = -1
			return leftDegree(vertex);
	}
	virtual tEdge rightNeighbour(tVertex vertex, int number) = 0;
	virtual tEdge leftNeighbour(tVertex vertex, int number) = 0;
	tEdge neighbour(tVertex vertex, int number, int direction) {
		if (direction == RIGHT)//RIGHT = 1
			return getRightNeighbour(vertex);
		else if (direction == LEFT)//LEFT = -1
			return getRightNeighbour(vertex);
	}

	virtual IVertexIterator<tVertex> *vertexIterator() = 0;

	//In order to add edge to graph one should create this edge first!
	virtual tEdge addEdge(tEdge newEdge) = 0;
	virtual void removeEdge(tEdge edge) = 0;

	//create ne vertex, adds it to graph and return
	virtual tVertex addVertex(tVertex) = 0;
	//add adjecent edges should be removed as well
	virtual void removeVertex(tVertex vertex) = 0;
	virtual tEdge merge(tEdge edge1, tEdge edge2, int direction = RIGHT) = 0;
	virtual pair<tEdge, tEdge> splitEdge(tEdge edge, int position) = 0;

	virtual tVertex glueVertices(tVertex vertex1, tVertex vertex2) = 0;
	//glue edges, there start and end vertices
	virtual tEdge glueEdges(tEdge edge1, tEdge edge2) = 0;
	//seperate edges  adjecent to the vertex
	virtual void unGlueEdges(tVertex vertex) = 0;
	virtual void unGlueEdgesLeft(tVertex vertex) = 0;
	virtual void unGlueEdgesRight(tVertex vertex) = 0;
};

class PairedGraphData {
public:
	//0 - in-degrees
	//1 -out-degrees
	int degrees[MAX_VERT_NUMBER][2];//, outD[MAX_VERT_NUMBER][2];
	int edgeIds[MAX_VERT_NUMBER][MAX_DEGREE][2];
	//	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
	vector<VertexPrototype *> vertexList_;
	longEdgesMap longEdges;
	verticesMap verts;
	int VertexCount;
	int EdgeId;
	PairedGraphData() {
		cerr << "VAH Paired created" << endl;
	}
	//	void RebuildVertexMap(void);
};

class VertexIterator: public IVertexIterator<VertexPrototype *> {
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
						+ graph_->degrees[currentVertex_][1] == 0) {
			currentVertex_++;
		}
		return currentVertex_ < graph_->VertexCount;
	}

	virtual VertexPrototype *next() {
		if(!hasNext())
			assert(false);
		VertexPrototype *result = graph_->vertexList_[currentVertex_];
		currentVertex_++;
		return result;
	}
};

class PairedGraph: public PairedGraphData, public IPairedGraph<VertexPrototype *, Edge *> {
private:
	//direction is 1 or -1. index is 0 or 1.
	inline int directionToIndex(int direction) {
		assert(direction != 0);
		return (direction + 1) >> 1;
	}

	void removeEdgeVertexAdjacency(VertexPrototype *vertex, Edge *edge, int direction) {
		int index = directionToIndex(direction);
		int current = 0;
		while (edgeIds[vertex->VertexId][current][index] != edge->EdgeId) {
			current++;
		}
		while (current + 1 < degrees[vertex->VertexId][index]) {
			edgeIds[vertex->VertexId][current][index]
					= edgeIds[vertex->VertexId][current + 1][index];
			current++;
		}
	}

	void addEdgeVertexAdjacency(VertexPrototype *vertex, Edge *edge, int direction) {
		int index = directionToIndex(direction);
		edgeIds[vertex->VertexId][degrees[vertex->VertexId][index]][index] = edge->EdgeId;
		degrees[vertex->VertexId][index]++;
	}
public:
	virtual int rightDegree(VertexPrototype *vertex) {
		return degrees[vertex->VertexId][1];
	}

	virtual int leftDegree(VertexPrototype *vertex) {
		return degrees[vertex->VertexId][0];
	}

	virtual Edge *rightNeighbour(VertexPrototype *vertex, int number) {
		assert(number < degrees[vertex->VertexId][1]);
		return longEdges[edgeRealId(edgeIds[vertex->VertexId][number][1], longEdges)];
	}
	virtual Edge *leftNeighbour(VertexPrototype *vertex, int number) {
		assert(number >= degrees[vertex->VertexId][0]);
		return longEdges[edgeRealId(edgeIds[vertex->VertexId][number][0], longEdges)];
	}

	//This is very bad method!!!!
	virtual IVertexIterator<VertexPrototype *> *vertexIterator() {
		return new VertexIterator(this);
	}

	virtual Edge *addEdge(Edge *newEdge) {
		newEdge->EdgeId = EdgeId;
		EdgeId++;
		longEdges.insert(make_pair(newEdge->EdgeId, newEdge));
		edgeIds[newEdge->FromVertex][degrees[newEdge->FromVertex][1]][1]
				= newEdge->EdgeId;
		degrees[newEdge->FromVertex][1]++;
		edgeIds[newEdge->FromVertex][degrees[newEdge->FromVertex][0]][0]
				= newEdge->EdgeId;
		degrees[newEdge->FromVertex][0]++;
		return newEdge;
	}

	virtual void removeEdge(Edge *edge) {
		if (edge = longEdges[edge->EdgeId]) {
			removeEdgeVertexAdjacency(vertexList_[edge->FromVertex], edge, 1);
			removeEdgeVertexAdjacency(vertexList_[edge->ToVertex], edge, -1);
		}
		delete edge;
	}

	virtual VertexPrototype *addVertex(VertexPrototype *vertex) {
		int vertexIndex = VertexCount;
		VertexCount++;
		degrees[VertexCount][0] = 0;
		degrees[VertexCount][1] = 0;
		vertex->VertexId = vertexIndex;
		vertexList_.push_back(vertex);
		verts[vertex->upper].push_back(vertex);
		return vertex;
	}

	virtual void removeVertex(VertexPrototype *vertex) {
		for (int index = 0; index <= 1; index++) {
			while (degrees[VertexCount][0] > 0) {
				removeEdge(longEdges[edgeIds[vertex->VertexId][0][index]]);
			}
		}
	}

	virtual Edge *merge(Edge *edge1, Edge *edge2, int direction = RIGHT) {
		Sequence *upper = new Sequence(
				edge1->upper->Subseq(0, edge1->length) + *(edge2->upper));
		Sequence *lower = new Sequence(
				edge1->lower->Subseq(0, edge1->length) + *(edge2->lower));
		Edge *edge = new Edge(upper, lower, edge1->FromVertex, edge2->ToVertex,
				edge1->length + edge2->length, 0, 0);
		addEdge(edge);
		removeEdge(edge1);
		removeEdge(edge2);
	}

	virtual pair<Edge *, Edge *> splitEdge(Edge *edge, int position) {
		assert(position > 0 && position < edge->length);
		Sequence *vertexUpper = new Sequence(edge->upper->Subseq(position, position + k - 1));
		Sequence *vertexLower = new Sequence(edge->lower->Subseq(position, position + k - 1));
		VertexPrototype *newVertex = addVertex(new VertexPrototype(vertexLower, 0));
		Sequence *upper1 = new Sequence(
				edge->upper->Subseq(0, position + k - 1));
		Sequence *upper2 = new Sequence(
				edge->upper->Subseq(position, edge->length + k - 1));
		Sequence *lower1 = new Sequence(
				edge->lower->Subseq(0, position + k - 1));
		Sequence *lower2 = new Sequence(
				edge->lower->Subseq(position, edge->length + k - 1));
		Edge *edge1 = new Edge(upper1, lower1, edge->FromVertex, newVertex->VertexId,
				position, 0, 0);
		Edge *edge2 = new Edge(upper2, lower2, newVertex->VertexId, edge->ToVertex,
				edge->length - position, 0, 0);
		removeEdge(edge);
		addEdge(edge1);
		addEdge(edge2);
		return make_pair(edge1, edge2);
	}

	virtual VertexPrototype *glueVertices(VertexPrototype *vertex1, VertexPrototype *vertex2) {
		if (vertex1 != vertex2) {
			return vertex1;
			for (int direction = -1; direction <= 1; direction += 2) {
				int index = directionToIndex(direction);
				for (int i = 0; i < degrees[vertex2->VertexId][index]; i++) {
					addEdgeVertexAdjacency(vertex1,
							longEdges[edgeIds[vertex2->VertexId][i][index]], direction);
				}
			}
			degrees[vertex2->VertexId][0] = 0;
			degrees[vertex2->VertexId][1] = 0;
			removeVertex(vertex2);
		}
		return vertex1;
	}

	//glue edges, there start and end vertices
	virtual Edge *glueEdges(Edge *edge1, Edge *edge2) {
		int fromVertex = edge2->FromVertex;
		int toVertex = edge2->ToVertex;
		removeEdge(edge2);
		glueVertices(vertexList_[edge1->FromVertex], vertexList_[fromVertex]);
		glueVertices(vertexList_[edge1->ToVertex], vertexList_[toVertex]);
		return edge1;
	}

	virtual void unGlueEdges(VertexPrototype *vertex) {
		if (degrees[vertex->VertexId][0] == 1) {
			unGlueEdgesLeft(vertex);
		} else if (degrees[vertex->VertexId][1] == 1) {
			unGlueEdgesRight(vertex);
		} else {
			assert(false);
		}
	}

	virtual void unGlueEdgesLeft(VertexPrototype *vertex) {
		assert(degrees[vertex->VertexId][0] == 1);
		Edge *leftEdge = longEdges[edgeIds[vertex->VertexId][0][0]];
		for (int i = 0; i < degrees[vertex->VertexId][1]; i++) {
			Edge *rightEdge = longEdges[edgeIds[vertex->VertexId][i][1]];
			rightEdge->ExpandLeft(*leftEdge);
			addEdgeVertexAdjacency(vertexList_[rightEdge->FromVertex], rightEdge, RIGHT);
			removeEdgeVertexAdjacency(vertex, rightEdge, RIGHT);
		}
		removeVertex(vertex);
	}

	virtual void unGlueEdgesRight(VertexPrototype *vertex) {
		assert(degrees[vertex->VertexId][0] == 1);
		Edge *rightEdge = longEdges[edgeIds[vertex->VertexId][0][1]];
		for (int i = 0; i < degrees[vertex->VertexId][0]; i++) {
			Edge *leftEdge = longEdges[edgeIds[vertex->VertexId][i][0]];
			leftEdge->ExpandLeft(*rightEdge);
			addEdgeVertexAdjacency(vertexList_[leftEdge->ToVertex], rightEdge, LEFT);
			removeEdgeVertexAdjacency(vertex, leftEdge, LEFT);
		}
		removeVertex(vertex);
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
