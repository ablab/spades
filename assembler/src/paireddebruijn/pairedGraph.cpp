#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphio.hpp"
//
//int isDeleted(Node* A) {
//	return A->upperSize == 0;
//}
//int deleteNode(Node* A) {
//	return A->upperSize = 0;
//}
LOGGER("p.pairedGraph");

namespace paired_assembler {

VertexPrototype::VertexPrototype(Sequence *lower_, int id, int coverage_) {
	upper = NULL;
	TRACE("Upper sequence is not defined");
	lower = lower_;
	VertexId = id;
	used = false;
	coverage = coverage_;
}

VertexPrototype::VertexPrototype(ll upper_, Sequence *lower_, int id,
		int coverage_) {
	upper = upper_;
	lower = lower_;
	VertexId = id;
	used = false;
	coverage = coverage_;
}

void PairedGraph::recreateVerticesInfo(int vertCount, longEdgesMap &longEdges) {
	forn(i, vertCount) {
		forn(j, 2)
			degrees[i][j] = 0;
	}
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			edgeIds[it->second->FromVertex][degrees[it->second->FromVertex][1]++][OUT_EDGE]
					= it->first;
			edgeIds[it->second->ToVertex][degrees[it->second->ToVertex][0]++][IN_EDGE]
					= it->first;
		}
	}
}

//todo: Complete this
void PairedGraph::RebuildVertexMap(void) {
	verts.clear();
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			ll kmer = 0;
			forn(i,k-1)
				kmer = (kmer << 2) | (it->second->upper->operator[](i));
			Sequence * seq = new Sequence(it->second->lower->Subseq(0, l - 1));
			storeVertex(*this, kmer, seq, it->second->FromVertex);
			kmer = 0;
			forn(i,k-1)
				kmer = (kmer << 2) | (it->second->upper->operator[](
						i + it->second->length));
			seq = new Sequence(
					it->second->lower->Subseq(it->second->length,
							it->second->length + l - 1));
			storeVertex(*this, kmer, seq, it->second->ToVertex);
		}
	}
}

verticesMap::iterator addKmerToMap(verticesMap &verts, ll kmer) {
	verticesMap::iterator position = verts.find(kmer);
	if (position == verts.end()) {
		vector<VertexPrototype *> prototypes;
		return verts.insert(make_pair(kmer, prototypes)).first;
	} else {
		return position;
	}
}

/*First argument of result is id of the vertex. Second argument is true if new entry
 * was created and false otherwise
 *
 */
pair<int, bool> addVertexToMap(PairedGraph &graph, ll newKmer, Sequence* newSeq) {
	verticesMap::iterator position = addKmerToMap(graph.verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		Sequence *otherSequence = (*it)->lower;
		if (newSeq->similar(*otherSequence, minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, graph.VertexCount));
	graph.VertexCount++;
	return make_pair(graph.VertexCount - 1, true);
}
pair<int, bool> addVertexToMap(PairedGraph &graph, ll newKmer,
		Sequence* newSeq, int VertNum) {
	verticesMap::iterator position = addKmerToMap(graph.verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		Sequence *otherSequence = (*it)->lower;
		if (newSeq->similar(*otherSequence, minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, VertNum));
	return make_pair(VertNum, true);
}

int storeVertex(gvis::GraphPrinter<int> &g, PairedGraph &graph, ll newKmer,
		Sequence* newSeq) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq);
	if (addResult.second) {
		g.addVertex(addResult.first, decompress(newKmer, k - 1));
	}
	return addResult.first;
}
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq);
	return addResult.first;
}
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, int VertNum) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq, VertNum);
	return addResult.first;
}

void resetVertexCount(PairedGraph &graph) {
	graph.VertexCount = 0;
}

inline int PairedGraph::directionToIndex(int direction) {
	assert(direction != 0);
	return (direction + 1) >> 1;
}

void PairedGraph::removeEdgeVertexAdjacency(VertexPrototype *vertex,
		Edge *edge, int direction) {
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
void PairedGraph::removeEdgeVertexAdjacency(int vertex, Edge *edge, int direction) {
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


void PairedGraph::addEdgeVertexAdjacency(VertexPrototype *vertex, Edge *edge,
		int direction) {
	int index = directionToIndex(direction);
	edgeIds[vertex->VertexId][degrees[vertex->VertexId][index]][index]
			= edge->EdgeId;
	degrees[vertex->VertexId][index]++;
	if (direction == RIGHT)
		edge->FromVertex = vertex->VertexId;
	else
		edge->ToVertex = vertex->VertexId;
}

void PairedGraph::addEdgeVertexAdjacency(int vertex, Edge *edge, int direction) {
	int index = directionToIndex(direction);
	cerr<<"add vertex adjacency"<<endl;
	edgeIds[vertex][degrees[vertex][index]][index] = edge->EdgeId;
	degrees[vertex][index]++;
	cerr<<"new degree "<<degrees[vertex][index]<<"for dir "<<direction<<endl;
}

int PairedGraph::rightDegree(VertexPrototype *vertex) {
	return degrees[vertex->VertexId][1];
}

int PairedGraph::leftDegree(VertexPrototype *vertex) {
	return degrees[vertex->VertexId][0];
}

Edge *PairedGraph::rightEdge(VertexPrototype *vertex, int number) {
	assert(number < degrees[vertex->VertexId][1]);
	return longEdges[edgeRealId(edgeIds[vertex->VertexId][number][1], longEdges)];
}
Edge *PairedGraph::leftEdge(VertexPrototype *vertex, int number) {
	assert(number < degrees[vertex->VertexId][0]);
	return longEdges[edgeRealId(edgeIds[vertex->VertexId][number][0], longEdges)];
}

//This is very bad method!!!!
VertexIterator *PairedGraph::vertexIterator() {
	return new VertexIterator(this);
}

Edge *PairedGraph::addEdge(Edge *newEdge) {

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
	return newEdge;
}


/*
void PairedGraph::removeEdge(Edge *edge) {
	if (edge = longEdges[edge->EdgeId]) {
		removeEdgeVertexAdjacency(vertexList_[edge->FromVertex], edge, 1);
		removeEdgeVertexAdjacency(vertexList_[edge->ToVertex], edge, -1);
	}
	delete edge;
}
*/

void PairedGraph::removeEdge(Edge *edge) {
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


VertexPrototype *PairedGraph::addVertex(VertexPrototype *vertex) {
	int vertexIndex = VertexCount;
	degrees[vertexIndex][0] = 0;
	degrees[vertexIndex][1] = 0;
	vertex->VertexId = vertexIndex;
	vertexList_.push_back(vertex);
	verts[vertex->upper].push_back(vertex);
	VertexCount++;
	return vertex;
}
int PairedGraph::addVertex() {
	degrees[VertexCount][0] = 0;
	degrees[VertexCount][1] = 0;
	VertexCount++;
	return VertexCount-1;
}


void PairedGraph::removeVertex(VertexPrototype *vertex) {
	for (int index = 0; index <= 1; index++) {
		while (degrees[vertex->VertexId][index] > 0) {
			removeEdge(longEdges[edgeIds[vertex->VertexId][index][index]]);
		}
	}
	vertexList_[vertex->VertexId] = NULL;
	delete vertex;
}

void PairedGraph::removeVertex(int vertex) {
	for (int index = 0; index <= 1; index++) {
		while (degrees[vertex][index] > 0) {
			removeEdge(longEdges[edgeIds[vertex][index][index]]);
		}
	}
}

Edge *PairedGraph::concat(Edge *edge1, Edge *edge2) {
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

pair<Edge *, Edge *> PairedGraph::splitEdge(Edge *edge, int position, int direction) {
	assert(position > 0 && position < edge->length);
//todo: Possible it can be simplified!

	if (direction == RIGHT){
		Sequence *vertexUpper = new Sequence(
				edge->upper->Subseq(position, position + k - 1));
		Sequence *vertexLower = new Sequence(
				edge->lower->Subseq(position, position + l - 1));
		VertexPrototype *newVertex = addVertex(new VertexPrototype(vertexLower, 0));
		Sequence *upper1 = new Sequence(edge->upper->Subseq(0, position + k - 1));
		Sequence *upper2 = new Sequence(
				edge->upper->Subseq(position, edge->length + k - 1));
		Sequence *lower1 = new Sequence(edge->lower->Subseq(0, position + l - 1));
		Sequence *lower2 = new Sequence(
				edge->lower->Subseq(position, edge->length + l - 1));
		Edge *edge1 = new Edge(upper1, lower1, edge->FromVertex,
				newVertex->VertexId, position, 0, edge->coverage);
		Edge *edge2 = new Edge(upper2, lower2, newVertex->VertexId, edge->ToVertex,
				edge->length - position, 0, edge->coverage);
		removeEdge(edge);
		addEdge(edge1);
		addEdge(edge2);
		return make_pair(edge1, edge2);
	}
	else if (direction == LEFT){
		Sequence *vertexUpper = new Sequence(
				edge->upper->Subseq(edge->upper->size() - k + 1));
		Sequence *vertexLower = new Sequence(
				edge->lower->Subseq(edge->lower->size() - l + 1));
		VertexPrototype *newVertex = addVertex(new VertexPrototype(vertexLower, 0));
		Sequence *upper1 = new Sequence(edge->upper->Subseq(edge->upper->size() - position - k + 1));
		Sequence *upper2 = new Sequence(
				edge->upper->Subseq(0, edge->upper->size() - position));
		Sequence *lower1 = new Sequence(edge->lower->Subseq(edge->lower->size() - position - l + 1));
		Sequence *lower2 = new Sequence(
				edge->lower->Subseq(0, edge->lower->size() - position));
		Edge *edge1 = new Edge(upper1, lower1, newVertex->VertexId, edge->ToVertex,
				 position, 0, edge->coverage);
		Edge *edge2 = new Edge(upper2, lower2, edge->FromVertex, newVertex->VertexId,
				edge->length - position, 0, edge->coverage);
		removeEdge(edge);
		addEdge(edge1);
		addEdge(edge2);
		return make_pair(edge1, edge2);

	}
	else {ERROR("Incorrect direction in split edge!"); assert(0); };


}

VertexPrototype *PairedGraph::glueVertices(VertexPrototype *vertex1,
		VertexPrototype *vertex2) {
	if (vertex1 != vertex2) {
		for (int direction = LEFT; direction <= RIGHT; direction += 2) {
			int index = directionToIndex(direction);
			for (int i = 0; i < degrees[vertex2->VertexId][index]; i++) {
				addEdgeVertexAdjacency(vertex1,
						longEdges[edgeIds[vertex2->VertexId][i][index]],
						direction);
			}
		}
		degrees[vertex2->VertexId][0] = 0;
		degrees[vertex2->VertexId][1] = 0;
		removeVertex(vertex2);
	}
	return vertex1;
}

int PairedGraph::glueVertices(int vertex1, int vertex2) {
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
/*Edge *PairedGraph::glueEdges(Edge *edge1, Edge *edge2) {
	int fromVertex = edge2->FromVertex;
	int toVertex = edge2->ToVertex;
	removeEdge(edge2);
	glueVertices(vertexList_[edge1->FromVertex], vertexList_[fromVertex]);
	glueVertices(vertexList_[edge1->ToVertex], vertexList_[toVertex]);
	return edge1;
}*/

Edge *PairedGraph::glueEdges(Edge *edge1, Edge *edge2) {
	cerr<<"Glue edge "<<edge1->EdgeId<<" "<<edge2->EdgeId<<endl;
	edge1->coverage += edge2->coverage;
	int fromVertex = edge2->FromVertex;
	int toVertex = edge2->ToVertex;
	removeEdge(edge2);
	glueVertices(edge1->FromVertex, fromVertex);
	glueVertices(edge1->ToVertex, toVertex);
	return edge1;
}



bool PairedGraph::unGlueEdgesIfPossible(VertexPrototype *vertex) {
	if (degrees[vertex->VertexId][0] == 1) {
		unGlueEdgesLeft(vertex);
	} else if (degrees[vertex->VertexId][1] == 1) {
		unGlueEdgesRight(vertex);
	} else {
		assert(false);
	}
}

bool PairedGraph::unGlueEdgesLeft(VertexPrototype *vertex) {
	if (degrees[vertex->VertexId][0] != 1) {
		return false;
	}
	Edge *leftEdge = longEdges[edgeIds[vertex->VertexId][0][0]];
	for (int i = 0; i < degrees[vertex->VertexId][1]; i++) {
		Edge *rightEdge = longEdges[edgeIds[vertex->VertexId][i][1]];
		rightEdge->ExpandLeft(*leftEdge);
		addEdgeVertexAdjacency(vertexList_[rightEdge->FromVertex], rightEdge,
				RIGHT);
		removeEdgeVertexAdjacency(vertex, rightEdge, RIGHT);
	}
	removeVertex(vertex);
	return true;
}

bool PairedGraph::unGlueEdgesRight(VertexPrototype *vertex) {
	if (degrees[vertex->VertexId][0] != 1)
		return false;
	Edge *rightEdge = longEdges[edgeIds[vertex->VertexId][0][1]];
	for (int i = 0; i < degrees[vertex->VertexId][0]; i++) {
		Edge *leftEdge = longEdges[edgeIds[vertex->VertexId][i][0]];
		leftEdge->ExpandLeft(*rightEdge);
		addEdgeVertexAdjacency(vertexList_[leftEdge->ToVertex], rightEdge, LEFT);
		removeEdgeVertexAdjacency(vertex, leftEdge, LEFT);
	}
	removeVertex(vertex);
	return true;
}

VertexPrototype *PairedGraph::findVertex(ll kmer, Sequence *s) {
	verticesMap::iterator it = verts.find(kmer);
	if (it == verts.end())
		return NULL;
	for (int i = 0; i < it->second.size(); i++) {
		if (*s == *(it->second[i]->lower)) {
			return it->second[i];
		}
	}
	return NULL;
}

VertexIterator::VertexIterator(PairedGraphData *graph) {
	graph_ = graph;
	currentVertex_ = 0;
}

bool VertexIterator::hasNext() {
	while (currentVertex_ < graph_->VertexCount
			&& graph_->degrees[currentVertex_][0]
					+ graph_->degrees[currentVertex_][1] == 0) {
		currentVertex_++;
	}
	return currentVertex_ < graph_->VertexCount;
}

VertexPrototype *VertexIterator::next() {
	if (!hasNext())
		assert(false);
	VertexPrototype *result = graph_->vertexList_[currentVertex_];
	currentVertex_++;
	return result;
}

}
