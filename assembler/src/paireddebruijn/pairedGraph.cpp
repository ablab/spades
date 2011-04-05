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

VertexPrototype::VertexPrototype(Sequence *lower_, int id, int coverage_,int position_) {
	upper = NULL;
	TRACE("Upper sequence is not defined");
	lower = lower_;
	VertexId = id;
	used = false;
	coverage = coverage_;

	position = (lower->size()- k + 1) /2 ;
}

VertexPrototype::VertexPrototype(ll upper_, Sequence *lower_, int id,
		int coverage_, int deltaShift_, int position_) {
	upper = upper_;
	lower = lower_;
	VertexId = id;
	used = false;
	coverage = coverage_;
	deltaShift = deltaShift_;
	position = position_;
}

void PairedGraph::recreateVerticesInfo(int vertCount, longEdgesMap &longEdges) {
	INFO("recreateVerticesInfo");
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
void PairedGraph::removeLowCoveredEdges(longEdgesMap &longEdges, int CoverageThreshold){
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
//			if ((degrees[it->second->FromVertex][1]>1)&&(degrees[it->second->ToVertex][0]>1))
			if (it->second->coverage <= CoverageThreshold) longEdges.erase(it--);
		}
	}
}


//todo: Complete this
void PairedGraph::RebuildVertexMap(void) {
	INFO("RebuildVertexMap");
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
	isUpToDate = true;
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
pair<int, bool> addVertexToMap(PairedGraph &graph, ll newKmer, Sequence* newSeq, bool ensureNew) {
	verticesMap::iterator position = addKmerToMap(graph.verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	if(!ensureNew)
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		Sequence *otherSequence = (*it)->lower;
		if (newSeq->similar(*otherSequence, minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, graph.VertexCount));
	cerr<<"added vertex "<<graph.VertexCount<<" kmer "<<newKmer<<" seq "<<newSeq->str()<<endl;
	graph.VertexCount++;
	return make_pair(graph.VertexCount - 1, true);
}
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
	cerr<<"added vertex "<<graph.VertexCount<<" kmer "<<newKmer<<" seq "<<newSeq->str()<<endl;
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
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, bool ensureNew) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq, ensureNew);
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

void PairedGraph::removeEdgeVertexAdjacency(int vertex, Edge *edge,
		int direction) {
	int index = directionToIndex(direction);
	int current = 0;
	while (edgeIds[vertex][current][index] != edge->EdgeId) {
		current++;
	}
	while (current + 1 < degrees[vertex][index]) {
		edgeIds[vertex][current][index] = edgeIds[vertex][current + 1][index];
		current++;
	}
	degrees[vertex][index]--;
}

//void PairedGraph::removeEdgeVertexAdjacency(int vertex, Edge *edge,
//		int direction) {
//	removeEdgeVertexAdjacency(vertexList_[vertex], edge, direction);
//	//	cerr<<"removeEdgeVertexAjacency vert "<<vertex<<" edge "<<edge->EdgeId<<" dir "<<direction<<endl;
//	//	int index = directionToIndex(direction);
//	//	int current = 0;
//	//	while (edgeIds[vertex][current][index] != edge->EdgeId) {
//	//		cerr<<"index "<<current<<" edge "<<edgeIds[vertex][current][index]<<endl;
//	//		current++;
//	//		assert(current<degrees[vertex][index]);
//	//	}
//	//	while (current + 1 < degrees[vertex][index]) {
//	//		edgeIds[vertex][current][index]
//	//				= edgeIds[vertex][current + 1][index];
//	//		current++;
//	//	}
//	//	degrees[vertex][index]--;
//	//	cerr<<"new degree "<<degrees[vertex][index]<<endl;
//}

//void PairedGraph::removeEdgeVertexAdjacency(int vertex, Edge *edge, int direction) {
//
//	cerr<<"removeEdgeVertexAdjacency for vert "<<vertex<<endl;
//	cerr<<"check vertexList "<<vertexList_[vertex]->VertexId<<endl;
//	removeEdgeVertexAdjacency(vertexList_[vertex], edge, direction);
////	cerr<<"removeEdgeVertexAjacency vert "<<vertex<<" edge "<<edge->EdgeId<<" dir "<<direction<<endl;
////	int index = directionToIndex(direction);
////	int current = 0;
////	while (edgeIds[vertex][current][index] != edge->EdgeId) {
////		cerr<<"index "<<current<<" edge "<<edgeIds[vertex][current][index]<<endl;
////		current++;
////		assert(current<degrees[vertex][index]);
////	}
////	while (current + 1 < degrees[vertex][index]) {
////		edgeIds[vertex][current][index]
////				= edgeIds[vertex][current + 1][index];
////		current++;
////	}
////	degrees[vertex][index]--;
////	cerr<<"new degree "<<degrees[vertex][index]<<endl;
//}

/*
 void PairedGraph::addEdgeVertexAdjacency(VertexPrototype *vertex, Edge *edge,
 int direction) {
 int index = directionToIndex(direction);
 edgeIds[vertex][degrees[vertex][index]][index]
 = edge->EdgeId;
 degrees[vertex][index]++;
 if (direction == RIGHT)
 edge->FromVertex = vertex;
 else
 edge->ToVertex = vertex;
 }
 */
void PairedGraph::addEdgeVertexAdjacency(int vertex, Edge *edge, int direction) {
	int index = directionToIndex(direction);
	cerr << "add vertex adjacency" << endl;
	edgeIds[vertex][degrees[vertex][index]][index] = edge->EdgeId;
	degrees[vertex][index]++;
	cerr << "new degree " << degrees[vertex][index] << "for dir " << direction
			<< endl;
	if (direction == RIGHT)
		edge->FromVertex = vertex;
	else
		edge->ToVertex = vertex;
}

int PairedGraph::rightDegree(int vertex) {
	return degrees[vertex][1];
}

int PairedGraph::leftDegree(int vertex) {
	return degrees[vertex][0];
}

Edge *PairedGraph::neighbourEdge(int vertex, int number, int direction) {
	int index = directionToIndex(direction);
	assert(number < degrees[vertex][index]);
	return longEdges[edgeRealId(edgeIds[vertex][number][index], longEdges)];
}

//This is very bad method!!!!
//JVertexIterator PairedGraph::jVertexIterator() {
//	return JVertexIterator(this);
//}

VertexIterator PairedGraph::beginVertex() {
	return VertexIterator(this, 0);
}

VertexIterator PairedGraph::endVertex() {
	return VertexIterator(this, -1);
}

EdgeIterator PairedGraph::beginEdge(int vertex, int direction) {
	int index = directionToIndex(direction);
	return EdgeIterator(this, vertex, index, 0);
}

EdgeIterator PairedGraph::endEdge(int vertex, int direction) {
	int index = directionToIndex(direction);
	return EdgeIterator(this, vertex, index, -1);
}

Edge *PairedGraph::addEdge(Edge *newEdge) {
	newEdge->EdgeId = EdgeId;
	cerr << "addEdge: new id " << newEdge->EdgeId << "(" << newEdge->FromVertex
			<< "->" << newEdge->ToVertex << ")" << endl;
	EdgeId++;
	longEdges.insert(make_pair(newEdge->EdgeId, newEdge));
	addEdgeVertexAdjacency(newEdge->FromVertex, newEdge, RIGHT);
	addEdgeVertexAdjacency(newEdge->ToVertex, newEdge, LEFT);
	//	edgeIds[newEdge->FromVertex][degrees[newEdge->FromVertex][1]][1]
	//			= newEdge->EdgeId;
	//	degrees[newEdge->FromVertex][1]++;
	//	cerr<<"Vert "<<newEdge->FromVertex<< " out degree "<<degrees[newEdge->FromVertex][1]<<endl;
	//
	//	edgeIds[newEdge->ToVertex][degrees[newEdge->ToVertex][0]][0]
	//			= newEdge->EdgeId;
	//	degrees[newEdge->ToVertex][0]++;
	//	cerr<<"Vert "<<newEdge->ToVertex<< " in  degree "<<degrees[newEdge->ToVertex][0]<<endl;
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
	cerr << "removeEdge " << edge->EdgeId << endl;
	if (edge == longEdges[edge->EdgeId]) {
		cerr << "real removeEdge " << edge->EdgeId << endl;
		removeEdgeVertexAdjacency(edge->FromVertex, edge, 1);
		removeEdgeVertexAdjacency(edge->ToVertex, edge, -1);

	}
	cerr << "delete ";
	longEdges.erase(edge->EdgeId);
	delete edge;
	//	edge->clearData();
	//	edge->EdgeId = -1;
	cerr << " ok " << endl;
}

int PairedGraph::leftEnd(Edge *edge) {
	return edge->FromVertex;
}

int PairedGraph::rightEnd(Edge *edge) {
	return edge->ToVertex;
}

int PairedGraph::end(Edge *edge, int direction) {
	if (direction == RIGHT)
		return rightEnd(edge);
	else if (direction == LEFT)
		return leftEnd(edge);
	else
		assert(false);
}

int PairedGraph::addVertex(int vertex) {
	isUpToDate = false;
	int vertexIndex = VertexCount;
	degrees[vertexIndex][0] = 0;
	degrees[vertexIndex][1] = 0;
	VertexCount++;
	return vertexIndex;
}

int PairedGraph::addVertex() {
	return addVertex(0);
}

void PairedGraph::removeVertex(int vertex) {
	isUpToDate = false;
	for (int index = 0; index <= 1; index++) {
		while (degrees[vertex][index] > 0) {
			removeEdge(longEdges[edgeIds[vertex][0][index]]);
		}
	}
}

//void PairedGraph::removeVertex(int vertex) {
//	removeVertex(vertexList_[vertex]);
//	//	for (int index = 0; index <= 1; index++) {
//	//		while (degrees[vertex][index] > 0) {
//	//			removeEdge(longEdges[edgeIds[vertex][index][index]]);
//	//		}
//	//	}
//}

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

Edge *subEdge(Edge *edge, int startPosition, int endPosition, int leftVertexId,
		int rightVertexId) {
	Sequence *upper = new Sequence(
			edge->upper->Subseq(startPosition, endPosition + k - 1));
	Sequence *lower = new Sequence(
			edge->lower->Subseq(startPosition, endPosition + l - 1));
	return new Edge(upper, lower, leftVertexId, rightVertexId,
			endPosition - startPosition, 0, edge->coverage);
}

pair<Edge *, Edge *> PairedGraph::splitEdge(Edge *edge, int position,
		int direction) {
	assert(position > 0 && position < edge->length);
	assert(direction == RIGHT || direction == LEFT);
	int realPosition = 0;
	if (direction == RIGHT)
		realPosition = position;
	else {
		realPosition = edge->length - position;
	}
	//	Sequence *vertexUpper = new Sequence(
	//			edge->upper->Subseq(realPosition, realPosition + k - 1));
	//	Sequence *vertexLower = new Sequence(
	//			edge->lower->Subseq(position, realPosition + l - 1));
	int newVertex = addVertex();
	//	Sequence *upper1 = new Sequence(
	//			edge->upper->Subseq(0, realPosition + k - 1));
	//	Sequence *upper2 = new Sequence(
	//			edge->upper->Subseq(realPosition, edge->upper->size()));
	//
	//	Sequence *lower1 = new Sequence(edge->lower->Subseq(0, position + l - 1));
	//	Sequence *lower2 = new Sequence(
	//			edge->lower->Subseq(position, edge->lower->size()));
	Edge *edge1 = subEdge(edge, 0, realPosition, edge->FromVertex, newVertex);
	Edge *edge2 = subEdge(edge, realPosition, edge->length, newVertex,
			edge->ToVertex);
	removeEdge(edge);
	addEdge(edge1);
	addEdge(edge2);
	if (direction == RIGHT) {
		return make_pair(edge1, edge2);
	} else if (direction == LEFT) {
		return make_pair(edge2, edge1);
	}
}

int PairedGraph::glueVertices(int vertex1, int vertex2) {
	if (vertex1 != vertex2) {
		for (int direction = LEFT; direction <= RIGHT; direction += 2) {
			int index = directionToIndex(direction);
			for (int i = 0; i < degrees[vertex2][index]; i++) {
				addEdgeVertexAdjacency(vertex1,
						longEdges[edgeIds[vertex2][i][index]], direction);
			}
		}
		degrees[vertex2][0] = 0;
		degrees[vertex2][1] = 0;
		removeVertex(vertex2);
	}
	return vertex1;
}
//
//int PairedGraph::glueVertices(int vertex1, int vertex2) {
//	return glueVertices(vertexList_[vertex1], vertexList_[vertex2])->VertexId;
//	//	cerr<<"Glue vertices "<<vertex1<<" "<<vertex2<<endl;
//	//
//	//	if (vertex1 != vertex2) {
//	////			return vertex1;
//	//		for (int direction = -1; direction <= 1; direction += 2) {
//	//			int index = directionToIndex(direction);
//	//			for (int i = 0; i < degrees[vertex2][index]; i++) {
//	//				addEdgeVertexAdjacency(vertex1,
//	//						longEdges[edgeIds[vertex2][i][index]], direction);
//	//				if (direction == 1) longEdges[edgeIds[vertex2][i][index]]->FromVertex = vertex1;
//	//				if (direction == -1) longEdges[edgeIds[vertex2][i][index]]->ToVertex = vertex1;
//	//			}
//	//		}
//	//		degrees[vertex2][0] = 0;
//	//		degrees[vertex2][1] = 0;
//	//		removeVertex(vertex2);
//	//	}
//	//	return vertex1;
//}

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
	cerr << "Glue edge " << edge1->EdgeId << " " << edge2->EdgeId << endl;
	edge1->coverage += edge2->coverage;
	int fromVertex = edge2->FromVertex;
	int toVertex = edge2->ToVertex;
	removeEdge(edge2);
	glueVertices(edge1->FromVertex, fromVertex);
	glueVertices(edge1->ToVertex, toVertex);
	return edge1;
}

bool PairedGraph::unGlueEdgesIfPossible(int vertex) {
	if (degrees[vertex][0] == 1) {
		unGlueEdgesLeft(vertex);
		return true;
	} else if (degrees[vertex][1] == 1) {
		unGlueEdgesRight(vertex);
		return true;
	} else {
		return false;
	}
}

bool PairedGraph::unGlueEdgesLeft(int vertex) {
	if (degrees[vertex][0] != 1) {
		return false;
	}
	Edge *leftEdge = longEdges[edgeIds[vertex][0][0]];
	for (int i = 0; i < degrees[vertex][1]; i++) {
		Edge *rightEdge = longEdges[edgeIds[vertex][i][1]];
		rightEdge->ExpandLeft(*leftEdge);
		addEdgeVertexAdjacency(rightEdge->FromVertex, rightEdge, RIGHT);
		removeEdgeVertexAdjacency(vertex, rightEdge, RIGHT);
	}
	removeVertex(vertex);
	return true;
}

bool PairedGraph::unGlueEdgesRight(int vertex) {
	if (degrees[vertex][0] != 1)
		return false;
	Edge *rightEdge = longEdges[edgeIds[vertex][0][1]];
	for (int i = 0; i < degrees[vertex][0]; i++) {
		Edge *leftEdge = longEdges[edgeIds[vertex][i][0]];
		leftEdge->ExpandLeft(*rightEdge);
		addEdgeVertexAdjacency(leftEdge->ToVertex, rightEdge, LEFT);
		removeEdgeVertexAdjacency(vertex, leftEdge, LEFT);
	}
	removeVertex(vertex);
	return true;
}

int PairedGraph::findVertex(ll kmer, Sequence *s) {
	if (!isUpToDate)
		RebuildVertexMap();
	verticesMap::iterator it = verts.find(kmer);
	if (it == verts.end())
		return -1;
	for (int i = 0; i < it->second.size(); i++) {
		if (*s == *(it->second[i]->lower)) {
			return it->second[i]->VertexId;
		}
	}
	return -1;
}

vector<int> PairedGraph::findVertices(ll kmer) {
	if (!isUpToDate)
		RebuildVertexMap();
	verticesMap::iterator it = verts.find(kmer);
	vector<int> result;
	if (it != verts.end())
		for (int i = 0; i < it->second.size(); i++) {
			result.push_back(it->second[i]->VertexId);
		}
	return result;
}

JVertexIterator::JVertexIterator(PairedGraphData *graph) {
	graph_ = graph;
	currentVertex_ = 0;
}

JVertexIterator::JVertexIterator(const JVertexIterator &iterator) {
	graph_ = iterator.graph_;
	currentVertex_ = iterator.currentVertex_;
}

VertexIterator::VertexIterator(PairedGraphData *graph, int current) {
	graph_ = graph;
	currentVertex_ = current;
	findNext();
}

VertexIterator::VertexIterator(const VertexIterator &iterator) {
	graph_ = iterator.graph_;
	currentVertex_ = iterator.currentVertex_;
	findNext();
}

void VertexIterator::findNext() {
	while (currentVertex_ < graph_->VertexCount
			&& graph_->degrees[currentVertex_][0]
					+ graph_->degrees[currentVertex_][1] == 0) {
		currentVertex_++;
	}
	if (currentVertex_ >= graph_->VertexCount) {
		currentVertex_ = -1;
	}
}

VertexIterator &VertexIterator::operator++() {
	if (currentVertex_ == -1)
		return *this;
	currentVertex_++;
	findNext();
	return *this;
}

int &VertexIterator::operator*() {
	if (currentVertex_ == -1) {
		INFO("End of vertices reached");
		assert(false);
	}
	return currentVertex_;
}

bool VertexIterator::operator==(const VertexIterator &other) {
	return this->graph_ == other.graph_ && this->currentVertex_
			== other.currentVertex_;
}

bool VertexIterator::operator!=(const VertexIterator &other) {
	return this->graph_ != other.graph_ || this->currentVertex_
			!= other.currentVertex_;
}

bool JVertexIterator::hasNext() {
	while (currentVertex_ < graph_->VertexCount
			&& graph_->degrees[currentVertex_][0]
					+ graph_->degrees[currentVertex_][1] == 0) {
		currentVertex_++;
	}
	return currentVertex_ < graph_->VertexCount;
}

int JVertexIterator::next() {
	if (!hasNext())
		assert(false);
	currentVertex_++;
	return currentVertex_ - 1;
}

EdgeIterator::EdgeIterator(PairedGraphData *graph, int vertex, int index,
		int current) {
	graph_ = graph;
	vertex_ = vertex;
	index_ = index;
	current_ = current;
}

EdgeIterator::EdgeIterator(const EdgeIterator &iterator) {
	graph_ = iterator.graph_;
	vertex_ = iterator.vertex_;
	index_ = iterator.index_;
	current_ = iterator.current_;
}

EdgeIterator &EdgeIterator::operator++() {
	if (current_ == -1)
		return *this;
	current_++;
	if (current_ >= graph_->degrees[vertex_][index_])
		current_ = -1;
	return *this;
}

Edge *&EdgeIterator::operator*() {
	if (current_ == -1) {
		INFO("edges ended!");
		assert(false);
	}
	return graph_->longEdges[graph_->edgeIds[vertex_][current_][index_]];
}

bool EdgeIterator::operator==(const EdgeIterator &other) {
	if (current_ == -1)
		return other.current_ == -1;
	return graph_ == other.graph_ && vertex_ == other.vertex_ && index_
			== other.index_ && current_ == other.current_;
}

bool EdgeIterator::operator!=(const EdgeIterator &other) {
	return !this->operator ==(other);
}

}
