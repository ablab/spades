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
	cerr<<"realId for "<<id<<endl;
	while (longEdges[res]->EdgeId != res) {
		res = longEdges[res]->EdgeId;
		cerr<<"possible "<<res<<endl;
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
			return rightNeighbour(vertex, number);
		else if (direction == LEFT)
			return leftNeighbour(vertex, number);//LEFT = -1
	}

	virtual IVertexIterator<tVertex> *vertexIterator() = 0;

	//In order to add edge to graph one should create this edge first!
	virtual void addEdge(tEdge newEdge) = 0;
	virtual void removeEdge(tEdge edge) = 0;

	//create ne vertex, adds it to graph and return
	virtual tVertex addVertex() = 0;
	//add adjecent edges should be removed as well
	virtual void removeVertex(tVertex vertex) = 0;
	virtual tEdge concat(tEdge edge1, tEdge edge2, int direction = RIGHT) = 0;
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
	longEdgesMap longEdges;
	verticesMap verts;
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
		cerr<<"hasNext "<<currentVertex_<<endl;
		if (graph_==NULL) cerr<<"Iterator has not graph"<<endl;
		while (currentVertex_ < graph_->VertexCount){
			if (graph_->degrees[currentVertex_][0]+graph_->degrees[currentVertex_][1] > 0) break;
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
private:
	//direction is 1 or -1. index is 0 or 1.
	inline int directionToIndex(int direction) {
		assert(direction != 0);
		return (direction + 1) >> 1;
	}

	void removeEdgeVertexAdjacency(int vertex, Edge *edge, int direction) {
		cerr<<"removeEdgeVertexAjacency vert "<<vertex<<" edge "<<edge->EdgeId<<" dir "<<direction<<endl;
		int index = directionToIndex(direction);
		int current = 0;
		while (edgeIds[vertex][current][index] != edge->EdgeId) {
			cerr<<"index "<<current<<" edge "<<edgeIds[vertex][current][index]<<endl;
			current++;
			assert(current<degrees[vertex][index]);
		}
		while (current + 1 < degrees[vertex][index]) {
			edgeIds[vertex][current][index]
					= edgeIds[vertex][current + 1][index];
			current++;
		}
		degrees[vertex][index]--;
		cerr<<"new degree "<<degrees[vertex][index]<<endl;
	}

	void addEdgeVertexAdjacency(int vertex, Edge *edge, int direction) {
		int index = directionToIndex(direction);
		cerr<<"add vertex adjacency"<<endl;
		edgeIds[vertex][degrees[vertex][index]][index] = edge->EdgeId;
		degrees[vertex][index]++;
		cerr<<"new degree "<<degrees[vertex][index]<<"for dir "<<direction<<endl;
	}
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
		return new VertexIterator(this);
	}

	virtual void addEdge(Edge *newEdge) {

		newEdge->EdgeId = EdgeId;
		cerr<<"addEdge: new id "<<newEdge->EdgeId<< "("<<newEdge->FromVertex<<"->"<<newEdge->ToVertex<<")"<<endl;
		EdgeId++;
		longEdges.insert(make_pair(newEdge->EdgeId, newEdge));
		edgeIds[newEdge->FromVertex][degrees[newEdge->FromVertex][1]][1]
				= newEdge->EdgeId;
		degrees[newEdge->FromVertex][1]++;
		cerr<<"Vert "<<newEdge->FromVertex<< " out degree "<<degrees[newEdge->FromVertex][1]<<endl;

		edgeIds[newEdge->ToVertex][degrees[newEdge->ToVertex][0]][0]
				= newEdge->EdgeId;
		degrees[newEdge->ToVertex][0]++;
		cerr<<"Vert "<<newEdge->ToVertex<< " in  degree "<<degrees[newEdge->ToVertex][0]<<endl;
	}

	virtual void removeEdge(Edge *edge) {
		cerr<<"removeEdge "<<edge->EdgeId<<endl;
		int edgeId = edge->EdgeId;
		if (edge == longEdges[edge->EdgeId]) {
			cerr<<"real removeEdge "<<edge->EdgeId<<endl;
			removeEdgeVertexAdjacency(edge->FromVertex, edge, 1);
			removeEdgeVertexAdjacency(edge->ToVertex, edge, -1);

		}
		cerr<<"delete ";
		delete edge;
		longEdges.erase(edgeId);
		cerr<<" ok "<<endl;
	}

	virtual int addVertex() {
		degrees[VertexCount][0] = 0;
		degrees[VertexCount][1] = 0;
		VertexCount++;
		return VertexCount-1;
	}

	virtual void removeVertex(int vertex) {
		for (int index = 0; index <= 1; index++) {
			while (degrees[VertexCount][0] > 0) {
				removeEdge(longEdges[edgeIds[vertex][0][index]]);
			}
		}
	}

	virtual Edge *concat(Edge *edge1, Edge *edge2, int direction = RIGHT) {
		Sequence *upper = new Sequence(
				edge1->upper->Subseq(0, edge1->length) + *(edge2->upper));
		Sequence *lower = new Sequence(
				edge1->lower->Subseq(0, edge1->length) + *(edge2->lower));
		Edge *edge = new Edge(upper, lower, edge1->FromVertex, edge2->ToVertex,
				edge1->length + edge2->length, 0, 0);
		addEdge(edge);
		removeEdge(edge1);
		removeEdge(edge2);
		return edge;
	}

	virtual pair<Edge *, Edge *> splitEdge(Edge *edge, int position) {
		assert(position > 0 && position < edge->length);
		int newVertex = addVertex();
		Sequence *upper1 = new Sequence(
				edge->upper->Subseq(0, position + k - 1));
		Sequence *upper2 = new Sequence(
				edge->upper->Subseq(position, edge->length + k - 1));
		Sequence *lower1 = new Sequence(
				edge->lower->Subseq(0, position + k - 1));
		Sequence *lower2 = new Sequence(
				edge->lower->Subseq(position, edge->length + k - 1));
		Edge *edge1 = new Edge(upper1, lower1, edge->FromVertex, newVertex,
				position, 0, 0);
		Edge *edge2 = new Edge(upper2, lower2, newVertex, edge->ToVertex,
				edge->length - position, 0, 0);
		removeEdge(edge);
		addEdge(edge1);
		addEdge(edge2);
		return make_pair(edge1, edge2);
	}

	virtual int glueVertices(int vertex1, int vertex2) {
		cerr<<"Glue vertices "<<vertex1<<" "<<vertex2<<endl;

		if (vertex1 != vertex2) {
//			return vertex1;
			for (int direction = -1; direction <= 1; direction += 2) {
				int index = directionToIndex(direction);
				for (int i = 0; i < degrees[vertex2][index]; i++) {
					addEdgeVertexAdjacency(vertex1,
							longEdges[edgeIds[vertex2][i][index]], direction);
					if (direction == 1) longEdges[edgeIds[vertex2][i][index]]->FromVertex = vertex1;
					if (direction == -1) longEdges[edgeIds[vertex2][i][index]]->ToVertex = vertex1;
				}
			}
			degrees[vertex2][0] = 0;
			degrees[vertex2][1] = 0;
			removeVertex(vertex2);
		}
		return vertex1;
	}

	//glue edges, there start and end vertices
	virtual Edge *glueEdges(Edge *edge1, Edge *edge2) {
		cerr<<"Glue edge "<<edge1->EdgeId<<" "<<edge2->EdgeId<<endl;
		int fromVertex = edge2->FromVertex;
		int toVertex = edge2->ToVertex;
		removeEdge(edge2);
		glueVertices(edge1->FromVertex, fromVertex);
		glueVertices(edge1->ToVertex, toVertex);
		return edge1;
	}

	virtual void unGlueEdges(int vertex) {
		if (degrees[vertex][0] == 1) {
			unGlueEdgesLeft(vertex);
		} else if (degrees[vertex][1] == 1) {
			unGlueEdgesRight(vertex);
		} else {
			assert(false);
		}
	}

	virtual void unGlueEdgesLeft(int vertex) {
		assert(degrees[vertex][0] == 1);
		Edge *leftEdge = longEdges[edgeIds[vertex][0][0]];
		for (int i = 0; i < degrees[vertex][1]; i++) {
			Edge *rightEdge = longEdges[edgeIds[vertex][i][1]];
			rightEdge->ExpandLeft(*leftEdge);
			addEdgeVertexAdjacency(rightEdge->FromVertex, rightEdge, RIGHT);
			removeEdgeVertexAdjacency(vertex, rightEdge, RIGHT);
		}
		removeVertex(vertex);
	}

	virtual void unGlueEdgesRight(int vertex) {
		assert(degrees[vertex][0] == 1);
		Edge *rightEdge = longEdges[edgeIds[vertex][0][1]];
		for (int i = 0; i < degrees[vertex][0]; i++) {
			Edge *leftEdge = longEdges[edgeIds[vertex][i][0]];
			leftEdge->ExpandLeft(*rightEdge);
			addEdgeVertexAdjacency(leftEdge->ToVertex, rightEdge, LEFT);
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
