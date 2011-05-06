#include "edge_graph.hpp"
#include "logging.hpp"
#include "visualization_utils.hpp"

namespace edge_graph {

Sequence EdgeGraph::VertexNucls(VertexId v) const {
	if (v->outgoing_edges_.size() > 0) {
		return v->outgoing_edges_[0]->nucls().Subseq(0, k_);
	} else if (v->complement_->outgoing_edges_.size() > 0) {
		return !VertexNucls(v->complement_);
	}
	assert(false);
}

EdgeId EdgeGraph::AddSingleEdge(VertexId v1, VertexId v2, const Sequence& s,
		size_t coverage) {
	EdgeId newEdge = new Edge(s, v2, coverage);
	v1->AddOutgoingEdge(newEdge);
	return newEdge;
}

void EdgeGraph::DeleteAllOutgoing(Vertex *v) {
	vector<EdgeId> out = v->outgoing_edges_;
	for (vector<EdgeId>::iterator it = out.begin(); it != out.end(); ++it) {
		DeleteEdge(*it);
	}
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

void EdgeGraph::FireAddVertex(VertexId v) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplyAdd(*it, v);
	}
}

void EdgeGraph::FireAddEdge(EdgeId edge) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplyAdd(*it, edge);
	}
}

void EdgeGraph::FireDeleteVertex(VertexId v) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplyDelete(*it, v);
	}
}

void EdgeGraph::FireDeleteEdge(EdgeId edge) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplyDelete(*it, edge);
	}
}

void EdgeGraph::FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplyMerge(*it, oldEdges, newEdge);
	}
}

void EdgeGraph::FireGlue(EdgeId edge1, EdgeId edge2) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplyGlue(*it, edge1, edge2);
	}
}

void EdgeGraph::FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_->ApplySplit(*it, edge, newEdge1, newEdge2);
	}
}

VertexId EdgeGraph::HiddenAddVertex() {
	VertexId v1 = new Vertex();
	VertexId v2 = new Vertex();
	v1->SetComplement(v2);
	v2->SetComplement(v1);
	vertices_.insert(v1);
	vertices_.insert(v2);
	return v1;
}

VertexId EdgeGraph::AddVertex() {
	VertexId result = HiddenAddVertex();
	FireAddVertex(result);
	return result;
}

void EdgeGraph::DeleteVertex(VertexId v) {
	assert(IsDeadEnd(v) && IsDeadStart(v));
	assert(v != NULL);
	FireDeleteVertex(v);
	VertexId complement = v->Complement();
	vertices_.erase(v);
	delete v;
	vertices_.erase(complement);
	delete complement;
}

void EdgeGraph::ForceDeleteVertex(VertexId v) {
	DeleteAllOutgoing(v);
	DeleteAllOutgoing(v->Complement());
	DeleteVertex(v);
}

EdgeId EdgeGraph::HiddenAddEdge(VertexId v1, VertexId v2,
		const Sequence &nucls, size_t coverage) {
	assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
	assert(nucls.size() >= k_ + 1);
	//	assert(OutgoingEdge(v1, nucls[k_]) == NULL);
	EdgeId result = AddSingleEdge(v1, v2, nucls, coverage);
	EdgeId rcEdge = result;
	if (nucls != !nucls) {
		rcEdge = AddSingleEdge(v2->Complement(), v1->Complement(), !nucls,
				coverage);
	}
	result->SetComplement(rcEdge);
	rcEdge->SetComplement(result);
	return result;
}

EdgeId EdgeGraph::AddEdge(VertexId v1, VertexId v2, const Sequence &nucls,
		size_t coverage) {
	EdgeId result = HiddenAddEdge(v1, v2, nucls, coverage);
	FireAddEdge(result);
	return result;
}

void EdgeGraph::DeleteEdge(EdgeId edge) {
	FireDeleteEdge(edge);
	EdgeId rcEdge = Complement(edge);
	VertexId rcStart = Complement(edge->end());
	VertexId start = Complement(rcEdge->end());
	start->RemoveOutgoingEdge(edge);
	rcStart->RemoveOutgoingEdge(rcEdge);
	delete edge;
	if (edge != rcEdge) {
		delete rcEdge;
	}
}

bool EdgeGraph::AreLinkable(VertexId v1, VertexId v2, const Sequence &nucls) const {
	return VertexNucls(v1) == nucls.Subseq(0, k_) && VertexNucls(
			v2->Complement()) == (!nucls).Subseq(0, k_);
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
	return edge->Complement();
}

VertexId EdgeGraph::EdgeStart(EdgeId edge) const {
	return Complement(edge)->end()->Complement();
}

VertexId EdgeGraph::EdgeEnd(EdgeId edge) const {
	return edge->end();
}

bool EdgeGraph::CanCompressVertex(VertexId v) const {
	return v->OutgoingEdgeCount() == 1 && v->Complement()->OutgoingEdgeCount()
			== 1;
}

void EdgeGraph::CompressVertex(VertexId v) {
	//assert(CanCompressVertex(v));
	if (CanCompressVertex(v)) {
		Merge(GetUniqueIncomingEdge(v), GetUniqueOutgoingEdge(v));
	}
}

void EdgeGraph::Merge(EdgeId edge1, EdgeId edge2) {
	assert(EdgeEnd(edge1) == EdgeStart(edge2));
	vector<EdgeId> toCompress;
	toCompress.push_back(edge1);
	toCompress.push_back(edge2);
	MergePath(toCompress);
}

EdgeId EdgeGraph::MergePath(const vector<EdgeId>& path) {
	assert(!path.empty());
	SequenceBuilder sb;
	//	sb.append(GetUniqueIncomingEdge(path[0])->nucls());
	VertexId v1 = EdgeStart(path[0]);
	VertexId v2 = EdgeEnd(path[path.size() - 1]);
	for (vector<EdgeId>::const_iterator it = path.begin(); it != path.end(); ++it) {
		sb.append(EdgeNucls(*it));
	}
	EdgeId newEdge = HiddenAddEdge(v1, v2, sb.BuildSequence());
	FireMerge(path, newEdge);
	DeleteEdge(path[0]);
	for (size_t i = 0; i + 1 < path.size(); i++) {
		VertexId v = EdgeEnd(path[i]);
		DeleteEdge(path[i + 1]);
		DeleteVertex(v);
	}
	FireAddEdge(newEdge);
	return newEdge;
}

bool EdgeGraph::GoUniqueWay(EdgeId &e) {
	VertexId u = EdgeEnd(e);
	if (!CheckUniqueOutgiongEdge(u) || !CheckUniqueIncomingEdge(u)) {
		return false;
	}
	e = GetUniqueOutgoingEdge(u);
	return true;
}

void EdgeGraph::CompressAllVertices() {
	SmartVertexIterator<EdgeGraph> end = SmartVertexEnd();
	for (SmartVertexIterator<EdgeGraph> it = SmartVertexBegin(); it != end; ++it) {
		VertexId v = *it;
		if (CheckUniqueOutgiongEdge(v) && CheckUniqueIncomingEdge(v)) {
			EdgeId e = GetUniqueOutgoingEdge(v);
			while (GoUniqueWay(e)) {
			}
			vector<EdgeId> mergeList;
			e = Complement(e);
			do {
				mergeList.push_back(e);
			} while (GoUniqueWay(e));
			MergePath(mergeList);
		}
	}
}

pair<EdgeId, EdgeId> EdgeGraph::SplitEdge(EdgeId edge, size_t position) {
	assert(position >= 1 && position < length(edge));
	assert(edge != Complement(edge));
	Sequence s1 = EdgeNucls(edge).Subseq(0, position + k_);
	Sequence s2 = EdgeNucls(edge).Subseq(position);
	Sequence newSequence = s1 + s2.Subseq(k_);
	VertexId splitVertex = HiddenAddVertex();
	EdgeId newEdge1 = HiddenAddEdge(this->EdgeStart(edge), splitVertex, s1);
	EdgeId newEdge2 = HiddenAddEdge(splitVertex, this->EdgeEnd(edge), s2);
	FireSplit(edge, newEdge1, newEdge2);
	FireAddVertex(splitVertex);
	FireAddEdge(newEdge1);
	FireAddEdge(newEdge2);
	DeleteEdge(edge);
	return make_pair(newEdge1, newEdge2);
}

void EdgeGraph::GlueEdges(EdgeId edge1, EdgeId edge2) {
	FireDeleteEdge(edge2);
	FireGlue(edge1, edge2);
	FireAddEdge(edge2);
	VertexId start = EdgeStart(edge1);
	VertexId end = EdgeEnd(edge1);
	DeleteEdge(edge1);
	if (IsDeadStart(start) && IsDeadEnd(start)) {
		DeleteVertex(start);
	}
	if (IsDeadStart(end) && IsDeadEnd(end)) {
		DeleteVertex(end);
	}
}

}
