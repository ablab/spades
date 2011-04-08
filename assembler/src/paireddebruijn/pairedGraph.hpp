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



class Vertex;

struct constructingEdge {
	string upper;
	string lower;
	int coverage;
	short deltaShift;
};

class VertexPrototype {
public:
	ll upper;
	Sequence *lower;
	int VertexId;
	bool used;
	int coverage;
	int position;
	int deltaShift;
	VertexPrototype(Sequence *lower_, int id, int coverage_ = 1, int position = 0);
	VertexPrototype(ll upper_, Sequence *lower_, int id, int coverage_ = 1, int deltaShift_ = 0, int position_ = 0);
	~VertexPrototype() {
		delete lower;
	}
};

class EdgePrototype {
public:
	EdgePrototype(Sequence *lower_, int start_, int coverage_ = 1) {
		lower = lower_;
		VertexId = start_;
		used = false;
		looped = 0;
		coverage = coverage_;
	}
	Sequence *lower;
	int VertexId;
	bool used;
	int coverage;
	char looped;

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
	//length of upper part of edge, including one of adjacent vertex.
	int length;
	int FromVertex;
	int ToVertex;
	int EdgeId;
	int coverage;
	int deltaShift;
	Edge(Edge &e) {
		length = e.length;
		FromVertex = e.FromVertex;
		ToVertex = e.ToVertex;
		EdgeId = e.EdgeId;
		coverage = e.coverage;
		upper = new Sequence(const_cast<Sequence&> (*e.upper));
		lower = new Sequence(const_cast<Sequence&> (*e.lower));
		deltaShift = e.deltaShift;
	}
	//	Vertex(int coverage, int length, Sequence *kmer, Sequence *pair, bool direction, int delta_d);
	void ExpandRight(Edge &newRight) {
		ToVertex = newRight.ToVertex;
		if (newRight.length > 0) {
			length = length + newRight.length;
			string toOut = newRight.upper->str();
			assert(k-1 < (int)toOut.length());
			upper = new Sequence(
					upper->str() + newRight.upper->Subseq(k - 1).str());

			toOut = newRight.lower->str();
			assert(l-1 < (int)toOut.length());
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

	int computeInsertLength() {
		assert(0);
		return 0;
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
			int cov = 1, int dShift = 0) {
		upper = up;
		lower = low;
		FromVertex = from;
		ToVertex = to;
		length = len;
		EdgeId = id;
		coverage = cov;
		deltaShift = dShift;
	}
	Edge(constructingEdge protoEdge, int from, int to, int id,
				 int dShift = 0) {
		upper = new Sequence(protoEdge.upper);
		lower = new Sequence(protoEdge.lower);
		FromVertex = from;
		ToVertex = to;
		length = protoEdge.upper.length() - (k - 1);
		EdgeId = id;
		coverage = protoEdge.coverage;
		if (length > insertLength + k) {
			//TODO: recompute dShift

		}
		deltaShift = dShift;
	}

	void clearData() {
		delete upper;
		delete lower;
		upper = NULL;
		lower = NULL;
	}

	~Edge() {
		//		cerr << "destructing" << upper->str() << endl;
		if (upper != NULL) {
			if (upper != lower) {
				delete upper;
				delete lower;
			} else
				delete upper;
		}
	}
};

inline int edgeRealId(int id, longEdgesMap &longEdges) {
	int res = id;
	cerr << "realId for " << id << endl;
	while (longEdges[res]->EdgeId != res) {
		res = longEdges[res]->EdgeId;
		cerr << "possible " << res << endl;
	}
	return res;
}

template<typename tVertex>
class IVertexIterator {
	virtual bool hasNext() = 0;
	virtual tVertex next() = 0;
};

//TODO Replace all Sequence * with just Sequence
template<typename tVertex, typename tEdge, typename tJVertexIterator,
		typename tVertexIterator, typename tEdgeIterator>
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

	virtual tEdge neighbourEdge(tVertex vertex, int number, int direction) = 0;

	tEdge rightEdge(tVertex vertex, int number) {
		return this->neighbourEdge(vertex, number, RIGHT);
	}

	tEdge leftEdge(tVertex vertex, int number) {
		return this->neighbourEdge(vertex, number, LEFT);
	}

	virtual tEdgeIterator beginEdge(tVertex vertex, int direction) = 0;
	virtual tEdgeIterator endEdge(tVertex vertex, int direction) = 0;

	virtual tEdgeIterator beginRightEdge(tVertex vertex) {
		return this->beginEdge(vertex, RIGHT);
	}

	virtual tEdgeIterator endRightEdge(tVertex vertex) {
		return this->endEdge(vertex, RIGHT);
	}

	virtual tEdgeIterator beginLeftEdge(tVertex vertex) {
		return this->beginEdge(vertex, LEFT);
	}

	virtual tEdgeIterator endLeftEdge(tVertex vertex) {
		return this->endEdge(vertex, LEFT);
	}

	//	virtual tJVertexIterator jVertexIterator() = 0;

	virtual tVertexIterator beginVertex() = 0;
	virtual tVertexIterator endVertex() = 0;

	//In order to add edge to graph one should create this edge first!
	virtual tEdge addEdge(tEdge newEdge) = 0;
	virtual void removeEdge(tEdge edge) = 0;

	virtual tVertex leftEnd(tEdge edge) = 0;
	virtual tVertex rightEnd(tEdge edge) = 0;
	virtual tVertex end(tEdge edge, int direction) = 0;

	//create ne vertex, adds it to graph and return
	virtual tVertex addVertex(tVertex) = 0;
	//add adjecent edges should be removed as well
	virtual void removeVertex(tVertex vertex) = 0;
	virtual tEdge concat(tEdge edge1, tEdge edge2) = 0;
	virtual pair<tEdge, tEdge> splitEdge(tEdge edge, int position,
			int direction) = 0;

	virtual tVertex glueVertices(tVertex vertex1, tVertex vertex2) = 0;
	//glue edges, there start and end vertices
	virtual tEdge glueEdges(tEdge edge1, tEdge edge2) = 0;
	//seperate edges  adjecent to the vertex
	virtual bool unGlueEdgesIfPossible(tVertex vertex) = 0;
	virtual bool unGlueEdgesLeft(tVertex vertex) = 0;
	virtual bool unGlueEdgesRight(tVertex vertex) = 0;

	void unGlueEdges(tVertex vertex, int direction) {
		if (direction == RIGHT)
			unGlueEdgesRight(vertex);
		else if (direction == LEFT)
			unGlueEdgesLeft(vertex);
		else
			assert(false);
	}

	//This method has bad parameters but still nothing to replace them
	virtual tVertex findVertex(ll kmer, Sequence *s) = 0;
	virtual vector<tVertex> findVertices(ll kmer) = 0;
};

class PairedGraphData {
	//TODO write destructor which would delete all VertexPrototypes * and all Edge *
public:
	//0 - in-degrees
	//1 -out-degrees
	bool isUpToDate;
	int degrees[MAX_VERT_NUMBER][2];//, outD[MAX_VERT_NUMBER][2];
	int edgeIds[MAX_VERT_NUMBER][MAX_DEGREE][2];
	//	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
	//	vector<VertexPrototype *> vertexList_;
	longEdgesMap longEdges;verticesMap verts;
	int VertexCount;
	int EdgeId;
	PairedGraphData() {
		isUpToDate = false;
		VertexCount = 0;
		EdgeId = 0;
		cerr<<"Paired created"<<endl;
	}
	//	void RebuildVertexMap(void);
};

class JVertexIterator {
	friend class PairedGraphData;
private:
	int currentVertex_;
	PairedGraphData *graph_;
public:
	JVertexIterator(PairedGraphData *graph);
	JVertexIterator(const JVertexIterator &iterator);

	/*
	 * commented while merge. Possible we need it.
	 *
	 virtual bool hasNext() {
	 cerr<<"hasNext "<<currentVertex_<<endl;
	 if (graph_==NULL) cerr<<"Iterator has not graph"<<endl;
	 while (currentVertex_ < graph_->VertexCount){
	 if (graph_->degrees[currentVertex_][0]+graph_->degrees[currentVertex_][1] > 0) break;
	 currentVertex_++;
	 }
	 return currentVertex_ < graph_->VertexCount;
	 }
	 */
	bool hasNext();

	int next();
};

class VertexIterator {
	friend class PairedGraphData;
private:
	PairedGraphData *graph_;
	int currentVertex_;
	void findNext();
public:
	VertexIterator(PairedGraphData *graph, int current = 0);

	VertexIterator(const VertexIterator &iterator);

	VertexIterator &operator++();

	int &operator*();

	bool operator==(const VertexIterator &other);

	bool operator!=(const VertexIterator &other);
};

class EdgeIterator {
	friend class PairedGraphData;
private:
	PairedGraphData *graph_;
	int vertex_;
	int index_;
	int current_;
public:
	EdgeIterator(PairedGraphData *graph, int vertex, int index, int current);

	EdgeIterator(const EdgeIterator &iterator);

	EdgeIterator &operator++();

	Edge *&operator*();

	bool operator==(const EdgeIterator &other);

	bool operator!=(const EdgeIterator &other);
};

class PairedGraph: public PairedGraphData, public IPairedGraph<int, Edge *,
		JVertexIterator, VertexIterator, EdgeIterator> {
private:
	//TODO add vertexNumber and edgeNumber fields and getters for them
	/**Method takes direction as RIGHT or LEFT, which are +1 and -1 and returns corresponding index
	 * for arrays in PairedGraphData which is 0 for LEFT and 1 for RIGHT
	 * @param direction LEFT or RIGHT
	 */
	inline int directionToIndex(int direction);

	//TODO replace these methods with transferEdgeAdjesency!!
	/**
	 * Method deletes @edge from the list of edges adjecent to @vertex
	 * @param vertex Vertex to delete adjacent edge from
	 * @param edge Edge to delete
	 * @direction Direction from which edge should be deleted. @direction can be LEFT or RIGHT.
	 */
	void removeEdgeVertexAdjacency(int vertex, Edge *edge, int direction);
	//	void removeEdgeVertexAdjacency(int vertex, Edge *edge, int direction);

	/**
	 * Method adds @edge to the list of edges adjecent to @vertex. This method works in O(number of neighbours)
	 * time.
	 * @param vertex Vertex to add new adjacent edge to
	 * @param edge Edge to add
	 * @direction Direction from which edge is added to vertex. @direction can be LEFT or RIGHT.
	 */
	void addEdgeVertexAdjacency(int vertex, Edge *edge, int direction);
	//	void addEdgeVertexAdjacency(int vertex, Edge *edge, int direction);
public:

	//TODO C-style outcoming and incoming edge iterators should be added!!!

	/**
	 *Method returns number of outgoing edges for @vertex
	 */
	virtual int rightDegree(int vertex);
	//	virtual int rightDegree(int vertex) {
	//		return degrees[vertex][1];
	//	}

	virtual int degree(int vertexId, int direction) {
		int index = directionToIndex(direction);
		return degrees[vertexId][index];

	}

	/**
	 *Method returns number of incoming edges for @vertex
	 */
	virtual int leftDegree(int vertex);
	//	virtual int leftDegree(int vertex) {
	//		return degrees[vertex][0];
	//	}

	/**
	 * Method returns incoming edge for given vertex with a given number in the given direction
	 */
	virtual Edge *neighbourEdge(int vertex, int number, int direction);

	//TODO Replace with C-stype iterator
	//TODO Make it return iterator instead of iterator *.
	//This is java style iterator. Better use newly created c style iterator VertexIterator
	//	virtual JVertexIterator jVertexIterator();

	virtual VertexIterator beginVertex();
	virtual VertexIterator endVertex();

	virtual EdgeIterator beginEdge(int vertex, int direction);
	virtual EdgeIterator endEdge(int vertex, int direction);


	/**
	 * Method adds edge to graph and updates all data stored in graph correspondingly.
	 *@param newEdge edge with any id value.
	 *@return the same Edge. After the edge is added to graph it is assigned with new id value
	 */
	virtual Edge *addEdge(Edge *newEdge);

	/*
	 * Method removes edge from graph, deletes @edge object and all its contents including Sequences stored in it.
	 * @param edge edge to be deleted
	 */
	virtual void removeEdge(Edge *edge);

	/**
	 * Method returns start vertex of edge
	 * @param edge to find start of
	 */
	virtual int leftEnd(Edge *edge);

	/**
	 * Method returns end vertex of edge
	 * @param edge to find end of
	 */
	virtual int rightEnd(Edge *edge);

	/**
	 * Method returns start or finish of the edge depending on direction
	 */
	virtual int end(Edge *edge, int direction);

	/**
	 * Method adds vertex to graph and updates all data stored in graph correspondingly. @newVertex is supposed
	 * to have upper sequence defined otherwise graph data may become inconsistent.
	 *@param newVertex vertex with any id value.
	 *@return the same vertex. After the vertex is added to graph it is assigned with new id value
	 */
	virtual int addVertex(int newVertex);
	int addVertex();

	/*
	 * Method removes vertex from graph, deletes @vertex object and all its contents including Sequences stored
	 * in it.
	 * @param vertex vertex to be deleted
	 */
	virtual void removeVertex(int vertex);
	//	virtual void removeVertex(int vertex);

	/**
	 * Method concatenate two edges into new one and deletes old edges.
	 * @param edge1 left edge to be concatenated
	 * @param edge2 right edge to be concatenated
	 * @return concatenated edge
	 */
	virtual Edge *concat(Edge *edge1, Edge *edge2);

	/**
	 * Method splits edge in two at given position
	 * @param edge edge to be split
	 * @param position position to split edge. Can not be less or equal to 0 or larger or equal than length of
	 * @return Two edges created
	 */
	virtual pair<Edge *, Edge *> splitEdge(Edge *edge, int position,
			int direction = RIGHT);

	/**
	 * Method transfers all connections of @vertex2 to @vertex1 and removes @vertex2
	 * @param vertex1 vertex to be glued to
	 * @param vertex2 vertex to be removed
	 * @return vertex which is @vertex1 glued to @vertex2. In this implementation it is @vertex1
	 */
	virtual int glueVertices(int vertex1, int vertex2);
	//	virtual int glueVertices(int vertex1, int vertex2);

	/**
	 * Method glues edges, there start and end points. Currently implemented as gluing of start and end vertices
	 * and then deletion of @edges
	 * @param edge1 first edge to be glued
	 * @param edge2 second edge to be glued
	 * @return resulting edge
	 */
	virtual Edge *glueEdges(Edge *edge1, Edge *edge2);

	/**
	 * Method checks if given vertex has incoming or outcoming edge equal to 1 and in this case removes this
	 * vertex from graph and creates new edges which are concatenations of incoming and outcoming edges of
	 * the vertex.
	 *@vertex vertex edges adjecent to which should be processed
	 *@return true if either incoming or outcoming degree of the vertex is equal to 1 and false otherwise
	 */
	virtual bool unGlueEdgesIfPossible(int vertex);

	/**
	 * Method checks if given vertex has incoming degree equal to 1 and in this case removes this
	 * vertex from graph and creates new edges which are concatenations of the incoming edge with all
	 * outcoming edges of the vertex.
	 *@vertex vertex edges adjecent to which should be processed
	 *@return true if either incoming degree of the vertex is equal to 1 and false otherwise
	 */
	virtual bool unGlueEdgesLeft(int vertex);

	/**
	 * Method checks if given vertex has outcoming degree equal to 1 and in this case removes this
	 * vertex from graph and creates new edges which are concatenations of the all incoming edges with the
	 * outcoming edge of the vertex.
	 *@vertex vertex edges adjecent to which should be processed
	 *@return true if either outcoming degree of the vertex is equal to 1 and false otherwise
	 */
	virtual bool unGlueEdgesRight(int vertex);

	/**
	 * This method finds vertex in graph with specific parameters.
	 * @param kmer upper sequence
	 * @param s lower sequence
	 * @return vertex with given parameters if found and -1 otherwise.
	 */
	virtual int findVertex(ll kmer, Sequence *s);

	/**
	 * This method finds all vertices in graph with given k-mer as upper sequence.
	 * @param kmer upper sequence
	 * @return vector of vertices found
	 */
	virtual vector<int> findVertices(ll kmer);

	void recreateVerticesInfo(int vertCount, longEdgesMap &longEdges);
	void removeLowCoveredEdges(longEdgesMap &longEdges, int CoverageThreshold = 1);
	void RebuildVertexMap(void);
};

int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq);
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, bool ensureNew);
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, int VertNum);
void resetVertexCount(PairedGraph &graph);

}
#endif /* PAIREDGRAPH_H_ */
