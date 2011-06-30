#include "debruijn_graph.hpp"
#include "logging.hpp"
#include "visualization_utils.hpp"

namespace debruijn_graph {

Sequence OldDeBruijnGraph::VertexNucls(OldDeBruijnGraph::VertexId v) const {
	if (v->outgoing_edges_.size() > 0) {
		return v->outgoing_edges_[0]->nucls().Subseq(0, k_);
	} else if (v->conjugate_->outgoing_edges_.size() > 0) {
		return !VertexNucls(v->conjugate_);
	}
	assert(false);
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::AddSingleEdge(OldDeBruijnGraph::VertexId v1, OldDeBruijnGraph::VertexId v2,
		const Sequence& s, size_t coverage) {
	OldDeBruijnGraph::EdgeId newEdge = new Edge(s, v2, coverage);
	v1->AddOutgoingEdge(newEdge);
	return newEdge;
}

void OldDeBruijnGraph::DeleteAllOutgoing(Vertex *v) {
	vector<OldDeBruijnGraph::EdgeId> out = v->outgoing_edges_;
	for (vector<OldDeBruijnGraph::EdgeId>::iterator it = out.begin(); it != out.end(); ++it) {
		DeleteEdge(*it);
	}
}

const vector<OldDeBruijnGraph::EdgeId> OldDeBruijnGraph::OutgoingEdges(OldDeBruijnGraph::VertexId v) const {
	return v->outgoing_edges_;
}

const vector<OldDeBruijnGraph::EdgeId> OldDeBruijnGraph::IncomingEdges(OldDeBruijnGraph::VertexId v) const {
	vector<OldDeBruijnGraph::EdgeId> result;
	OldDeBruijnGraph::VertexId rcv = conjugate(v);
	vector<OldDeBruijnGraph::EdgeId> edges = rcv->OutgoingEdges();
	for (EdgeIterator it = edges.begin(); it != edges.end(); ++it) {
		result.push_back(conjugate(*it));
	}
	return result;
}

const vector<OldDeBruijnGraph::EdgeId> OldDeBruijnGraph::IncidentEdges(OldDeBruijnGraph::VertexId v) const {
	vector<OldDeBruijnGraph::EdgeId> result;
	DEBUG("Incident for vert: "<< v);
	for (EdgeIterator it = v->begin(); it != v->end(); ++it) {
		DEBUG("out:"<< *it);
		result.push_back(*it);
	}
	OldDeBruijnGraph::VertexId rcv = conjugate(v);

	for (EdgeIterator it = rcv->begin(); it != rcv->end(); ++it) {
		int fl = 1;
		for (int j = 0, sz = result.size(); j < sz; j++) {
			if (result[j] == *it) {
				fl = 0;
				break;
			}
		}
		if (fl) {
			DEBUG("in:"<< *it);
			result.push_back(conjugate(*it));
		}
	}
	return result;
}

const vector<OldDeBruijnGraph::EdgeId> OldDeBruijnGraph::NeighbouringEdges(OldDeBruijnGraph::EdgeId e) const {
	OldDeBruijnGraph::VertexId v_out = EdgeEnd(e);
	OldDeBruijnGraph::VertexId v_in = EdgeStart(e);
	vector<OldDeBruijnGraph::EdgeId> result = OldDeBruijnGraph::IncidentEdges(v_in);
	vector<OldDeBruijnGraph::EdgeId> out_res = OldDeBruijnGraph::IncidentEdges(v_out);
	// these vectors are small, and linear time is less than log in this case.
	for (vector<OldDeBruijnGraph::EdgeId>::iterator it = out_res.begin(); it != out_res.end(); ++it) {
		int fl = 1;
		for (int j = 0, sz = result.size(); j < sz; j++)
			if (result[j] == *it) {
				fl = 0;
				break;
			}

		if (fl)
			result.push_back(*it);
	}
	DEBUG(result.size());
	return result;
}

void OldDeBruijnGraph::FireAddVertex(OldDeBruijnGraph::VertexId v) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplyAdd(*it, v);
	}
}

void OldDeBruijnGraph::FireAddEdge(OldDeBruijnGraph::EdgeId edge) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplyAdd(*it, edge);
	}
}

void OldDeBruijnGraph::FireDeleteVertex(OldDeBruijnGraph::VertexId v) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplyDelete(*it, v);
	}
}

void OldDeBruijnGraph::FireDeleteEdge(OldDeBruijnGraph::EdgeId edge) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplyDelete(*it, edge);
	}
}

void OldDeBruijnGraph::FireMerge(vector<OldDeBruijnGraph::EdgeId> oldEdges, OldDeBruijnGraph::EdgeId newEdge) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplyMerge(*it, oldEdges, newEdge);
	}
}

void OldDeBruijnGraph::FireGlue(OldDeBruijnGraph::EdgeId new_edge, OldDeBruijnGraph::EdgeId edge1, OldDeBruijnGraph::EdgeId edge2) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplyGlue(*it, new_edge, edge1, edge2);
	}
}

void OldDeBruijnGraph::FireSplit(OldDeBruijnGraph::EdgeId edge, OldDeBruijnGraph::EdgeId newEdge1, OldDeBruijnGraph::EdgeId newEdge2) {
	for (vector<ActionHandler*>::iterator it = action_handler_list_.begin(); it
			!= action_handler_list_.end(); ++it) {
		applier_.ApplySplit(*it, edge, newEdge1, newEdge2);
	}
}

OldDeBruijnGraph::VertexId OldDeBruijnGraph::HiddenAddVertex() {
	OldDeBruijnGraph::VertexId v1 = new Vertex();
	OldDeBruijnGraph::VertexId v2 = new Vertex();
	v1->Setconjugate(v2);
	v2->Setconjugate(v1);
	vertices_.insert(v1);
	vertices_.insert(v2);
	return v1;
}

OldDeBruijnGraph::VertexId OldDeBruijnGraph::AddVertex() {
	OldDeBruijnGraph::VertexId result = HiddenAddVertex();
	FireAddVertex(result);
	return result;
}

void OldDeBruijnGraph::DeleteVertex(OldDeBruijnGraph::VertexId v) {
	TRACE("OldDeBruijnGraph::DeleteVertex "<<v);
	assert(IsDeadEnd(v) && IsDeadStart(v));
	TRACE("OldDeBruijnGraph::DeleteVertex "<<v);
	assert(v != NULL);
	TRACE("OldDeBruijnGraph::DeleteVertex "<<v);
	FireDeleteVertex(v);
	TRACE("OldDeBruijnGraph::DeleteVertex "<<v);
	OldDeBruijnGraph::VertexId conjugate = v->conjugate();
	TRACE("OldDeBruijnGraph::DeleteVertex "<<v);
	vertices_.erase(v);
	TRACE("OldDeBruijnGraph::DeleteVertex "<<v);
	delete v;
	TRACE("OldDeBruijnGraph::DeleteVertex ");
	vertices_.erase(conjugate);
	TRACE("OldDeBruijnGraph::DeleteVertex ");
	delete conjugate;
	TRACE("OldDeBruijnGraph::DeleteVertex Ok");
}

void OldDeBruijnGraph::ForceDeleteVertex(OldDeBruijnGraph::VertexId v) {
	DeleteAllOutgoing(v);
	DeleteAllOutgoing(v->conjugate());
	DeleteVertex(v);
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::HiddenAddEdge(OldDeBruijnGraph::VertexId v1, OldDeBruijnGraph::VertexId v2,
		const Sequence &nucls, size_t coverage) {
	assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
	assert(nucls.size() >= k_ + 1);
	//	assert(OutgoingEdge(v1, nucls[k_]) == NULL);
	OldDeBruijnGraph::EdgeId result = AddSingleEdge(v1, v2, nucls, coverage);
	OldDeBruijnGraph::EdgeId rcEdge = result;
	if (nucls != !nucls) {
		rcEdge = AddSingleEdge(v2->conjugate(), v1->conjugate(), !nucls,
				coverage);
	}
	result->set_conjugate(rcEdge);
	rcEdge->set_conjugate(result);
	return result;
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::AddEdge(OldDeBruijnGraph::VertexId v1, OldDeBruijnGraph::VertexId v2, const Sequence &nucls,
		size_t coverage) {
	OldDeBruijnGraph::EdgeId result = HiddenAddEdge(v1, v2, nucls, coverage);
	FireAddEdge(result);
	return result;
}

void OldDeBruijnGraph::DeleteEdge(OldDeBruijnGraph::EdgeId edge) {
	FireDeleteEdge(edge);
	OldDeBruijnGraph::EdgeId rcEdge = conjugate(edge);
	OldDeBruijnGraph::VertexId rcStart = conjugate(edge->end());
	OldDeBruijnGraph::VertexId start = conjugate(rcEdge->end());
	start->RemoveOutgoingEdge(edge);
	rcStart->RemoveOutgoingEdge(rcEdge);
	if (edge != rcEdge) {
		delete rcEdge;
	}
	delete edge;
}

bool OldDeBruijnGraph::AreLinkable(OldDeBruijnGraph::VertexId v1, OldDeBruijnGraph::VertexId v2, const Sequence &nucls) const {
	return VertexNucls(v1) == nucls.Subseq(0, k_) && VertexNucls(
			v2->conjugate()) == (!nucls).Subseq(0, k_);
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::OutgoingEdge(OldDeBruijnGraph::VertexId v, char nucl) const {
	vector<OldDeBruijnGraph::EdgeId> edges = v->OutgoingEdges();
	for (EdgeIterator iter = edges.begin(); iter != edges.end(); ++iter) {
		char lastNucl = (*iter)->nucls()[k_];
		if (lastNucl == nucl) {
			return *iter;
		}
	}
	return NULL;
}

OldDeBruijnGraph::VertexId OldDeBruijnGraph::conjugate(OldDeBruijnGraph::VertexId v) const {
	return v->conjugate();
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::conjugate(OldDeBruijnGraph::EdgeId edge) const {
	return edge->conjugate();
}

OldDeBruijnGraph::VertexId OldDeBruijnGraph::EdgeStart(OldDeBruijnGraph::EdgeId edge) const {
	return conjugate(edge)->end()->conjugate();
}

OldDeBruijnGraph::VertexId OldDeBruijnGraph::EdgeEnd(OldDeBruijnGraph::EdgeId edge) const {
	return edge->end();
}

bool OldDeBruijnGraph::CanCompressVertex(OldDeBruijnGraph::VertexId v) const {
	return v->OutgoingEdgeCount() == 1 && v->conjugate()->OutgoingEdgeCount()
			== 1;
}

void OldDeBruijnGraph::CompressVertex(OldDeBruijnGraph::VertexId v) {
	//assert(CanCompressVertex(v));
	if (CanCompressVertex(v)) {
		Merge(GetUniqueIncomingEdge(v), GetUniqueOutgoingEdge(v));
	}
}

void OldDeBruijnGraph::Merge(OldDeBruijnGraph::EdgeId edge1, OldDeBruijnGraph::EdgeId edge2) {
	assert(EdgeEnd(edge1) == EdgeStart(edge2));
	vector<OldDeBruijnGraph::EdgeId> toCompress;
	toCompress.push_back(edge1);
	toCompress.push_back(edge2);
	MergePath(toCompress);
}

vector<OldDeBruijnGraph::EdgeId> OldDeBruijnGraph::CorrectMergePath(const vector<OldDeBruijnGraph::EdgeId>& path) {
	vector<OldDeBruijnGraph::EdgeId> result;
	for (size_t i = 0; i < path.size(); i++) {
		if (path[i] == conjugate(path[i])) {
			if (i < path.size() - 1 - i) {
				for (size_t j = 0; j < path.size(); j++)
					result.push_back(conjugate(path[path.size() - 1 - j]));
				i = path.size() - 1 - i;
			} else {
				result = path;
			}
			size_t size = 2 * i + 1;
			for (size_t j = result.size(); j < size; j++) {
				result.push_back(conjugate(result[size - 1 - j]));
			}
			return result;
		}
	}
	return path;
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::MergePath(const vector<OldDeBruijnGraph::EdgeId>& path) {
	vector<OldDeBruijnGraph::EdgeId> correctedPath = CorrectMergePath(path);
	assert(!correctedPath.empty());
	OldDeBruijnGraph::EdgeId newEdge = AddMergedEdge(correctedPath);
	FireMerge(correctedPath, newEdge);
	DeletePath(correctedPath);
	FireAddEdge(newEdge);
	return newEdge;
}

OldDeBruijnGraph::EdgeId OldDeBruijnGraph::AddMergedEdge(const vector<OldDeBruijnGraph::EdgeId> &path) {
	OldDeBruijnGraph::VertexId v1 = EdgeStart(path[0]);
	OldDeBruijnGraph::VertexId v2 = EdgeEnd(path[path.size() - 1]);
	SequenceBuilder sb;
	sb.append(EdgeNucls(path[0]).Subseq(0, k_));
	for (vector<OldDeBruijnGraph::EdgeId>::const_iterator it = path.begin(); it != path.end(); ++it) {
		sb.append(EdgeNucls(*it).Subseq(k_));
//		sb.append(EdgeNucls(*it));
	}
	return HiddenAddEdge(v1, v2, sb.BuildSequence());
}

void OldDeBruijnGraph::DeletePath(const vector<OldDeBruijnGraph::EdgeId> &path) {
	set<OldDeBruijnGraph::EdgeId> edgesToDelete;
	set<OldDeBruijnGraph::VertexId> verticesToDelete;
	edgesToDelete.insert(path[0]);
	for (size_t i = 0; i + 1 < path.size(); i++) {
		OldDeBruijnGraph::EdgeId e = path[i + 1];
		if (edgesToDelete.find(conjugate(e)) == edgesToDelete.end())
			edgesToDelete.insert(e);
		OldDeBruijnGraph::VertexId v = EdgeStart(e);
		if (verticesToDelete.find(conjugate(v)) == verticesToDelete.end())
			verticesToDelete.insert(v);
	}
	for (auto it = edgesToDelete.begin(); it != edgesToDelete.end(); ++it)
		DeleteEdge(*it);
	for (auto it = verticesToDelete.begin(); it != verticesToDelete.end(); ++it)
		DeleteVertex(*it);
}

pair<OldDeBruijnGraph::EdgeId, OldDeBruijnGraph::EdgeId> OldDeBruijnGraph::SplitEdge(OldDeBruijnGraph::EdgeId edge, size_t position) {
	assert(position >= 1 && position < length(edge));
	assert(edge != conjugate(edge));
	Sequence s1 = EdgeNucls(edge).Subseq(0, position + k_);
	Sequence s2 = EdgeNucls(edge).Subseq(position);
	assert(s1 + s2.Subseq(k_) == EdgeNucls(edge));
//	Sequence newSequence = s1 + s2.Subseq(k_);
	OldDeBruijnGraph::VertexId splitVertex = HiddenAddVertex();
	OldDeBruijnGraph::EdgeId newEdge1 = HiddenAddEdge(this->EdgeStart(edge), splitVertex, s1);
	OldDeBruijnGraph::EdgeId newEdge2 = HiddenAddEdge(splitVertex, this->EdgeEnd(edge), s2);
	FireSplit(edge, newEdge1, newEdge2);
	FireAddVertex(splitVertex);
	FireAddEdge(newEdge1);
	FireAddEdge(newEdge2);
	DeleteEdge(edge);
	return make_pair(newEdge1, newEdge2);
}

void OldDeBruijnGraph::GlueEdges(OldDeBruijnGraph::EdgeId edge1, OldDeBruijnGraph::EdgeId edge2) {
	OldDeBruijnGraph::EdgeId newEdge = HiddenAddEdge(EdgeStart(edge2), EdgeEnd(edge2), EdgeNucls(edge2));
	FireGlue(newEdge, edge1, edge2);
	FireDeleteEdge(edge1);
	FireDeleteEdge(edge2);
	FireAddEdge(newEdge);
	OldDeBruijnGraph::VertexId start = EdgeStart(edge1);
	OldDeBruijnGraph::VertexId end = EdgeEnd(edge1);
	DeleteEdge(edge1);
	if (IsDeadStart(start) && IsDeadEnd(start)) {
		DeleteVertex(start);
	}
	if (IsDeadStart(end) && IsDeadEnd(end)) {
		DeleteVertex(end);
	}
}

}
