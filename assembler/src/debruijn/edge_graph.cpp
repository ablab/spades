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

Edge* EdgeGraph::AddSingleEdge(Vertex* v1, Vertex* v2, const Sequence& s) {
	Edge *newEdge = new Edge(s, v2);
	v1->AddOutgoingEdge(newEdge);
	return newEdge;
}

//void EdgeGraph::DeleteSingleEdge(const Edge* edge) {
//	Vertex *v = edgeStart(edge);
//	v->RemoveOutgoingEdge(edge);
//}

void EdgeGraph::DeleteAllOutgoing(Vertex *v) {
	vector<Edge *> out = v->outgoing_edges_;
	for (vector<Edge *>::iterator it = out.begin(); it != out.end(); ++it) {
		DeleteEdge(*it);
	}
}

void EdgeGraph::OutgoingEdges(const Vertex* v, Vertex::EdgeIterator &begin,
		Vertex::EdgeIterator &end) const {
	begin = v->begin();
	end = v->end();
}

const vector<Edge *> EdgeGraph::OutgoingEdges(const Vertex* v) const {
	return v->outgoing_edges_;
}

const vector<Edge *> EdgeGraph::IncomingEdges(const Vertex* v) const {
	vector<Edge *> result;
	Vertex *rcv = ComplementVertex(v);
	for (Vertex::EdgeIterator it = rcv->begin(); it != rcv->end(); ++it) {
		result.push_back(*it);
	}
	return result;
}

Vertex* EdgeGraph::AddVertex() {
	Vertex* v1 = new Vertex();
	Vertex* v2 = new Vertex();
	v1->set_complement(v2);
	v2->set_complement(v1);
	vertices_.insert(v1);
	vertices_.insert(v2);
	for (vector<GraphActionHandler<EdgeGraph> *>::iterator it =
			action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
		(*it)->HandleAdd(v1);
	}
	return v1;
}

void EdgeGraph::DeleteVertex(Vertex* v) {
	assert(IsDeadEnd(v) && IsDeadStart(v));
	assert(v != NULL);
	for (vector<GraphActionHandler<EdgeGraph> *>::iterator it =
			action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
		(*it)->HandleDelete(v);
	}
	Vertex* complement = v->complement();
	vertices_.erase(v);
	delete v;
	vertices_.erase(complement);
	delete complement;
}

void EdgeGraph::ForceDeleteVertex(Vertex* v) {
	DeleteAllOutgoing(v);
	DeleteAllOutgoing(v->complement());
	DeleteVertex(v);
}

Edge* EdgeGraph::AddEdge(Vertex* v1, Vertex* v2, const Sequence &nucls) {
	assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
	assert(nucls.size() >= k_ + 1);
	Edge *result = AddSingleEdge(v1, v2, nucls);
	if (nucls != !nucls)
		AddSingleEdge(v2->complement(), v1->complement(), !nucls);
	for (vector<GraphActionHandler<EdgeGraph> *>::iterator it =
			action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
		(*it)->HandleAdd(result);
	}
	return result;
}

void EdgeGraph::DeleteEdge(Edge* edge) {
	const Edge *rcEdge = ComplementEdge(edge);
	Vertex *rcStart = ComplementVertex(edge->end());
	Vertex *start = ComplementVertex(rcEdge->end());
	start->RemoveOutgoingEdge(edge);
	rcStart->RemoveOutgoingEdge(rcEdge);
	for (vector<GraphActionHandler<EdgeGraph> *>::iterator it =
			action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
		(*it)->HandleDelete(edge);
	}
	delete edge;
	if (edge != rcEdge)
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

bool EdgeGraph::CanCompressVertex(Vertex *v) const {
	return v->OutgoingEdgeCount() == 1 && v->complement()->OutgoingEdgeCount()
			== 1;
}

Edge *EdgeGraph::CompressVertex(Vertex *v) {
	assert(v->OutgoingEdgeCount() == 1 && v->complement()->OutgoingEdgeCount() == 1);
	Edge *edge1 = GetUniqueIncomingEdge(v);
	Edge *edge2 = GetUniqueOutgoingEdge(v);
	Sequence nucls = edge1->nucls() + edge2->nucls().Subseq(k_);
	Vertex *v1 = edgeStart(edge1);
	Vertex *v2 = edgeEnd(edge2);
	DeleteEdge(edge1);
	DeleteEdge(edge2);
	DeleteVertex(v);
	return AddEdge(v1, v2, nucls);
}

Edge *EdgeGraph::CompressPath(const vector<Vertex *> path) {
	assert(!path.empty());
	SequenceBuilder sb;
	assert(CheckUniqueIncomingEdge(path[0]));
	sb.append(GetUniqueIncomingEdge(path[0])->nucls());
	Vertex *v1 = edgeStart(GetUniqueIncomingEdge(path[0]));
	Vertex *v2 = edgeEnd(GetUniqueOutgoingEdge(path[path.size() - 1]));
	for (vector<Vertex *>::const_iterator it = path.begin(); it != path.end(); ++it) {
		sb.append(GetUniqueOutgoingEdge(*it)->nucls().Subseq(k_));
		ForceDeleteVertex(*it);
	}
	return AddEdge(v1, v2, sb.BuildSequence());
}

void EdgeGraph::CompressAllVertices() {
	SmartVertexIterator<EdgeGraph> end = this->SmartVertexEnd();
	for (SmartVertexIterator<EdgeGraph> it = this->SmartVertexBegin(); it
			!= end; ++it) {
		Vertex *v = *it;
		if(CheckUniqueOutgiongEdge(v) && CheckUniqueIncomingEdge(v)) {
			while(CheckUniqueOutgiongEdge(v))
				v = edgeEnd(GetUniqueOutgoingEdge(v));
			vector<Vertex *> compressList;
			v = ComplementVertex(v);
			while(CheckUniqueOutgiongEdge(v)) {
				compressList.push_back(v);
				v = edgeEnd(GetUniqueOutgoingEdge(v));
			}
			CompressPath(compressList);
		}
	}
}

}
