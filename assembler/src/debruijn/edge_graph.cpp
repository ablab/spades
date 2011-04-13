#include "edge_graph.hpp"
#include "logging.hpp"

namespace edge_graph {

Sequence EdgeGraph::VertexNucls(VertexId v) const {
	if (v->outgoing_edges_.size() > 0) {
		return v->outgoing_edges_[0]->nucls().Subseq(0, k_);
	} else if (v->complement_->outgoing_edges_.size() > 0) {
		return !VertexNucls(v->complement_);
	}
	assert(false);
	//	return new Sequence("");
}

EdgeId EdgeGraph::AddSingleEdge(VertexId v1, VertexId v2, const Sequence& s) {
	EdgeId newEdge = new Edge(s, v2);
	v1->AddOutgoingEdge(newEdge);
	return newEdge;
}

//void EdgeGraph::DeleteSingleEdge(const Edge* edge) {
//	Vertex *v = edgeStart(edge);
//	v->RemoveOutgoingEdge(edge);
//}

void EdgeGraph::DeleteAllOutgoing(Vertex *v) {
	vector<EdgeId> out = v->outgoing_edges_;
	for (vector<EdgeId>::iterator it = out.begin(); it != out.end(); ++it) {
		DeleteEdge(*it);
	}
}

void EdgeGraph::OutgoingEdges(VertexId v, EdgeIterator& begin,
		EdgeIterator& end) const {
	begin = v->begin();
	end = v->end();
}

const vector<EdgeId> EdgeGraph::OutgoingEdges(VertexId v) const {
	return v->outgoing_edges_;
}

const vector<EdgeId> EdgeGraph::IncomingEdges(VertexId v) const {
	vector<EdgeId> result;
	VertexId rcv = Complement(v);
	for (EdgeIterator it = rcv->begin(); it != rcv->end(); ++it) {
		result.push_back(*it);
	}
	return result;
}

VertexId EdgeGraph::AddVertex() {
	VertexId v1 = new Vertex();
	VertexId v2 = new Vertex();
	v1->set_complement(v2);
	v2->set_complement(v1);
	vertices_.insert(v1);
	vertices_.insert(v2);
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		(*it)->HandleAdd(v1);
	}
	return v1;
}

void EdgeGraph::DeleteVertex(VertexId v) {
	assert(IsDeadEnd(v) && IsDeadStart(v));
	assert(v != NULL);
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		(*it)->HandleDelete(v);
	}
	VertexId complement = v->complement();
	vertices_.erase(v);
	delete v;
	vertices_.erase(complement);
	delete complement;
}

void EdgeGraph::ForceDeleteVertex(VertexId v) {
	DeleteAllOutgoing(v);
	DeleteAllOutgoing(v->complement());
	DeleteVertex(v);
}

EdgeId EdgeGraph::AddEdge(VertexId v1, VertexId v2, const Sequence &nucls) {
	assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
	assert(nucls.size() >= k_ + 1);
	assert(OutgoingEdge(v1, nucls[k_]) == NULL);
	EdgeId result = AddSingleEdge(v1, v2, nucls);
	if (nucls != !nucls)
		AddSingleEdge(v2->complement(), v1->complement(), !nucls);
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		(*it)->HandleAdd(result);
	}
	return result;
}

void EdgeGraph::DeleteEdge(EdgeId edge) {
	EdgeId rcEdge = Complement(edge);
	VertexId rcStart = Complement(edge->end());
	VertexId start = Complement(rcEdge->end());
	start->RemoveOutgoingEdge(edge);
	rcStart->RemoveOutgoingEdge(rcEdge);
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		(*it)->HandleDelete(edge);
	}
	delete edge;
	if (edge != rcEdge)
		delete rcEdge;
}

bool EdgeGraph::AreLinkable(VertexId v1, VertexId v2, const Sequence &nucls) const {
	return VertexNucls(v1) == nucls.Subseq(0, k_) && VertexNucls(
			v2->complement()) == (!nucls).Subseq(0, k_);
}

EdgeId EdgeGraph::OutgoingEdge(VertexId v, char nucl) const {
	for (EdgeIterator iter = v->begin(); iter != v->end(); ++iter) {
		char lastNucl = (*iter)->nucls()[k_];
		if (lastNucl == nucl) {
			return *iter;
		}
	}
	return NULL;
}

EdgeId EdgeGraph::Complement(EdgeId edge) const {
	Sequence s = !(edge->nucls());
	char nucl = s[k_];
	VertexId v = edge->end();
	EdgeId result = OutgoingEdge(v->complement(), nucl);
	assert(result != NULL);
	return result;
}

VertexId EdgeGraph::EdgeStart(EdgeId edge) const {
	return Complement(edge)->end()->complement();
}

VertexId EdgeGraph::EdgeEnd(EdgeId edge) const {
	return edge->end();
}

bool EdgeGraph::CanCompressVertex(VertexId v) const {
	return v->OutgoingEdgeCount() == 1 && v->complement()->OutgoingEdgeCount()
			== 1;
}

EdgeId EdgeGraph::CompressVertex(VertexId v) {
	assert(v->OutgoingEdgeCount() == 1 && v->complement()->OutgoingEdgeCount() == 1);
	EdgeId edge1 = GetUniqueIncomingEdge(v);
	EdgeId edge2 = GetUniqueOutgoingEdge(v);
	Sequence nucls = edge1->nucls() + edge2->nucls().Subseq(k_);
	VertexId v1 = EdgeStart(edge1);
	VertexId v2 = EdgeEnd(edge2);
	DeleteEdge(edge1);
	DeleteEdge(edge2);
	DeleteVertex(v);
	return AddEdge(v1, v2, nucls);
}

EdgeId EdgeGraph::CompressPath(const vector<VertexId>& path) {
	assert(!path.empty());
	SequenceBuilder sb;
	assert(CheckUniqueIncomingEdge(path[0]));
	sb.append(GetUniqueIncomingEdge(path[0])->nucls());
	VertexId v1 = EdgeStart(GetUniqueIncomingEdge(path[0]));
	VertexId v2 = EdgeEnd(GetUniqueOutgoingEdge(path[path.size() - 1]));
	for (vector<Vertex *>::const_iterator it = path.begin(); it != path.end(); ++it) {
		sb.append(GetUniqueOutgoingEdge(*it)->nucls().Subseq(k_));
		ForceDeleteVertex(*it);
	}
	return AddEdge(v1, v2, sb.BuildSequence());
}

bool EdgeGraph::GoUniqueWay(VertexId &v) {
	VertexId u = EdgeEnd(GetUniqueOutgoingEdge(v));
	if (!CheckUniqueOutgiongEdge(u) || !CheckUniqueIncomingEdge(u))
		return false;
	v = u;
	return true;
}

void EdgeGraph::CompressAllVertices() {
	SmartVertexIterator<EdgeGraph> end = SmartVertexEnd();
	for (SmartVertexIterator<EdgeGraph> it = SmartVertexBegin(); it != end; ++it) {
		VertexId v = *it;
		if (CheckUniqueOutgiongEdge(v) && CheckUniqueIncomingEdge(v)) {
			while (GoUniqueWay(v))
				;
			//				v = EdgeEnd(GetUniqueOutgoingEdge(v));
			vector<VertexId> compressList;
			v = Complement(v);
			do
				compressList.push_back(v);
			while (GoUniqueWay(v));
			CompressPath(compressList);
		}
	}
}

void SimpleGraphVisualizer::Visualize(const EdgeGraph& g) {
	VisHandler h(g, gp_);
	de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
	gp_.output();
}

void ComplementGraphVisualizer::Visualize(const EdgeGraph& g) {
	ComplementVisHandler h(g, gp_);
	de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
	gp_.output();
}

void WriteToFile(const string& file_name, const string& graph_name,
		const EdgeGraph& g) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<VertexId> gp("simulated_data_graph", filestr);
	ComplementGraphVisualizer gv(gp);
	gv.Visualize(g);
	filestr.close();

}

}
