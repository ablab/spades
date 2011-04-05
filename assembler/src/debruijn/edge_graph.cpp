#include "edge_graph.hpp"
#include "logging.hpp"

namespace edge_graph {

Sequence EdgeGraph::vertexNucls(const Vertex *v) const {
	if (v->outgoing_edges_.size() > 0) {
		return v->outgoing_edges_[0]->nucls().Subseq(0, k_);
	} else if (v->complement_->outgoing_edges_.size() > 0) {
		return !vertexNucls(v->complement_);
	}
	assert(false);
	//	return new Sequence("");
}

bool EdgeGraph::CheckIfNoIncoming(Vertex* v) const {
	return v->begin() == v->end();
}

bool EdgeGraph::CanBeDeleted(Vertex* v) const {
	return CheckIfNoIncoming(v) && CheckIfNoIncoming(v->complement());
}

Edge* EdgeGraph::AddSingleEdge(Vertex* v1, Vertex* v2, const Sequence& s) {
	Edge *newEdge = new Edge(s, v2);
	v1->AddOutgoingEdge(newEdge);
	return newEdge;
}

//void EdgeGraph::DeleteSingleEdge(const Edge* edge) {
//	Vertex *v = edgeStart(edge);
//	v->RemoveOutgoingEdge(edge);
//}

void EdgeGraph::OutgoingEdges(const Vertex* v, Vertex::EdgeIterator &begin,
		Vertex::EdgeIterator &end) const {
	begin = v->begin();
	end = v->end();
}

Vertex* EdgeGraph::AddVertex() {
	Vertex* v1 = new Vertex();
	Vertex* v2 = new Vertex();
	v1->set_complement(v2);
	v2->set_complement(v1);
	vertices_.insert(v1);
	vertices_.insert(v2);
	action_handler_->HandleAdd(v1);
	return v1;
}

void EdgeGraph::DeleteVertex(Vertex* v) {
	assert(IsDeadEnd(v) && IsDeadStart(v));
	Vertex* complement = v->complement();
	vertices_.erase(v);
	vertices_.erase(complement);
	action_handler_->HandleDelete(v);
	delete v;
	delete complement;
}

//TODO Method needs a tiny bit of refactoring
void EdgeGraph::ForceDeleteVertex(Vertex* v) {
	vector<Edge *> toDelete;
	Vertex::EdgeIterator begin, end;
	OutgoingEdges(v, begin, end);
	toDelete.insert(toDelete.end(), begin, end);
	Vertex* complement = v->complement();
	OutgoingEdges(complement, begin, end);
	toDelete.insert(toDelete.end(), begin, end);
	for (vector<Edge *>::iterator it = toDelete.begin(); it != toDelete.end(); ++it) {
		DeleteEdge(*it);
	}
	DeleteVertex(v);
}

Edge* EdgeGraph::AddEdge(Vertex* v1, Vertex* v2, const Sequence &nucls) {
	assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
	assert(nucls.size() >= k_ + 1);
	Edge *result = AddSingleEdge(v1, v2, nucls);
	AddSingleEdge(v2->complement(), v1->complement(), !nucls);
	action_handler_->HandleAdd(result);
	return result;
}

void EdgeGraph::DeleteEdge(Edge* edge) {
	const Edge *rcEdge = ComplementEdge(edge);
	Vertex *rcStart = ComplementVertex(edge->end());
	Vertex *start = ComplementVertex(rcEdge->end());
	rcStart->RemoveOutgoingEdge(rcEdge);
	start->RemoveOutgoingEdge(edge);
	action_handler_->HandleDelete(edge);
	delete edge;
	delete rcEdge;
}

bool EdgeGraph::AreLinkable(Vertex* v1, Vertex* v2, const Sequence &nucls) const {
	return vertexNucls(v1) == nucls.Subseq(0, k_) && vertexNucls(
			v2->complement()) == (!nucls).Subseq(0, k_);
}

Edge* EdgeGraph::OutgoingEdge(const Vertex* v, char nucl) const {
	for (Vertex::EdgeIterator iter = v->begin(); iter != v->end(); ++iter) {
		char lastNucl = (*iter)->nucls()[k_];
		if (lastNucl == nucl) {
			return *iter;
		}
	}
	return NULL;
}

Edge *EdgeGraph::ComplementEdge(const Edge* edge) const {
	Sequence s = !(edge->nucls());
	char nucl = s[k_];
	Vertex *v = edge->end();
	Edge *result = OutgoingEdge(v->complement(), nucl);
	assert(result != NULL);
	return result;
}

Vertex *EdgeGraph::edgeStart(const Edge *edge) const {
	return ComplementEdge(edge)->end()->complement();
}

Vertex *EdgeGraph::edgeEnd(const Edge *edge) const {
	return edge->end();
}

}
