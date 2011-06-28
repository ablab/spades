#ifndef ABSTRACT_CONJUGATE_GRAPH_HPP_
#define ABSTRACT_CONJUGATE_GRAPH_HPP_

#include <vector>
#include <set>
#include <cstring>
#include "seq.hpp"
#include "sequence.hpp"
#include "logging.hpp"
#include "nucl.hpp"
#include "strobe_read.hpp"
#include "omni_utils.hpp"
#include "observable_graph.hpp"
#include "coverage.hpp"

namespace omnigraph {

template<typename VertexData, typename EdgeData, class DataMaster>
class AbstractConjugateGraph;

template<typename VertexData, typename EdgeData, class DataMaster>
class PairedEdge;

template<typename VertexData, typename EdgeData, class DataMaster>
class PairedVertex {
private:
	typedef PairedVertex<VertexData, EdgeData, DataMaster>* VertexId;
	typedef PairedEdge<VertexData, EdgeData, DataMaster>* EdgeId;

	friend class AbstractConjugateGraph<VertexData, EdgeData, DataMaster>;

	vector<EdgeId> outgoing_edges_;

	VertexId conjugate_;

	VertexData data_;

	void set_conjugate(VertexId conjugate) {
		conjugate_ = conjugate;
	}

	size_t OutgoingEdgeCount() const {
		return outgoing_edges_.size();
	}

	const vector<EdgeId> OutgoingEdges() const {
		return outgoing_edges_;
	}

	const vector<EdgeId> OutgoingEdgesTo(VertexId v) const {
		vector<EdgeId> result;
		for (auto it = outgoing_edges_.begin(); it != outgoing_edges_.end(); ++it) {
			if ((*it)->end() == v) {
				result.push_back(*it);
			}
		}
		return result;
	}

	size_t IncomingEdgeCount() const {
		return conjugate_->OutgoingEdgeCount();
	}

	const vector<EdgeId> IncomingEdges() const {
		vector<EdgeId> result = conjugate_->OutgoingEdges();
		for (size_t i = 0; i < result.size(); i++) {
			result[i] = result[i]->conjugate();
		}
		return result;
	}

	PairedVertex(VertexData data) :
		data_(data) {
	}

	VertexData &data() {
		return data_;
	}

	void set_data(VertexData data) {
		data_ = data;
	}

	bool IsDeadend() {
		return outgoing_edges_.size() == 0;
	}

	void AddOutgoingEdge(EdgeId e) {
		outgoing_edges_.push_back(e);
	}

	bool RemoveOutgoingEdge(const EdgeId e) {
		auto it = outgoing_edges_.begin();
		while (it != outgoing_edges_.end() && *it != e) {
			++it;
		}
		if (it == outgoing_edges_.end()) {
			return false;
		}
		outgoing_edges_.erase(it);
		return true;
	}

	VertexId conjugate() const {
		return conjugate_;
	}

	~PairedVertex() {
		TRACE("PairedVertex destructor");
		assert(outgoing_edges_.size() == 0);
		TRACE("PairedVertex destructor ok");
	}
};

template<typename VertexData, typename EdgeData, class DataMaster>
class PairedEdge {
private:
	typedef PairedVertex<VertexData, EdgeData, DataMaster>* VertexId;
	typedef PairedEdge<VertexData, EdgeData, DataMaster>* EdgeId;

	friend class AbstractConjugateGraph<VertexData, EdgeData, DataMaster> ;
	friend class PairedVertex<VertexData, EdgeData, DataMaster> ;
	VertexId end_;

	EdgeData data_;

	EdgeId conjugate_;

	PairedEdge(VertexId end, const EdgeData &data) :
		end_(end), data_(data) {
	}

	EdgeData &data() {
		return data_;
	}

	void set_data(EdgeData &data) {
		data_ = data;
	}

	VertexId end() const {
		return end_;
	}

	EdgeId conjugate() {
		return conjugate_;
	}

	void set_conjugate(EdgeId conjugate) {
		conjugate_ = conjugate;
	}

	~PairedEdge() {
	}
};

template<typename VertexData, typename EdgeData, class DataMaster>
class AbstractConjugateGraph: public ObservableGraph<PairedVertex<VertexData,
		EdgeData, DataMaster>*, PairedEdge<VertexData, EdgeData, DataMaster>*> {
public:
	typedef PairedVertex<VertexData, EdgeData, DataMaster>* VertexId;
	typedef set<VertexId> Vertices;
	typedef typename Vertices::const_iterator VertexIterator;
	typedef PairedEdge<VertexData, EdgeData, DataMaster>* EdgeId;
	typedef vector<EdgeId> Edges;
	typedef typename Edges::const_iterator EdgeIterator;
	typedef ObservableGraph<VertexId, EdgeId> super;
	typedef SmartVertexIterator<ObservableGraph<VertexId, EdgeId> >
			SmartVertexItarator;
	typedef SmartEdgeIterator<ObservableGraph<VertexId, EdgeId>>
			SmartEdgeItarator;

	set<VertexId> vertices_;

	DataMaster master_;

	VertexId HiddenAddVertex(const VertexData &data1, const VertexData &data2) {
		VertexId v1 =
				new PairedVertex<VertexData, EdgeData, DataMaster> (data1);
		VertexId v2 =
				new PairedVertex<VertexData, EdgeData, DataMaster> (data2);
		v1->set_conjugate(v2);
		v2->set_conjugate(v1);
		vertices_.insert(v1);
		vertices_.insert(v2);
		return v1;
	}

	VertexId HiddenAddVertex(const VertexData &data) {
		return HiddenAddVertex(data, master_.conjugate(data));
	}

	EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
		EdgeId result = AddSingleEdge(v1, v2, data);
		if (master_.isSelfConjugate(data)) {
			result->set_conjugate(result);
			return result;
		}
		EdgeId rcEdge = AddSingleEdge(v2->conjugate(), v1->conjugate(),
				master_.conjugate(data));
		result->set_conjugate(rcEdge);
		rcEdge->set_conjugate(result);
		return result;
	}

	EdgeId AddSingleEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		EdgeId newEdge = new PairedEdge<VertexData, EdgeData, DataMaster> (v2,
				data);
		v1->AddOutgoingEdge(newEdge);
		return newEdge;
	}

	void DeleteAllOutgoing(VertexId v) {
		vector<EdgeId> out = v->outgoing_edges_;
		TRACE("DeleteAllOutgoing "<<out.size());
		for (auto it = out.begin(); it != out.end(); ++it) {
			TRACE("DeleteOutgoing "<<*it);
			DeleteEdge(*it);
			TRACE("DeleteOutgoing ok");
		}
	}

public:

	VertexIterator begin() const {
		return vertices_.begin();
	}

	VertexIterator end() const {
		return vertices_.end();
	}

	size_t size() {
		return vertices_.size();
	}

	AbstractConjugateGraph(DataMaster master) :
				ObservableGraph<VertexId, EdgeId> (
						new PairedHandlerApplier<AbstractConjugateGraph<
								VertexData, EdgeData, DataMaster> > (*this)),
				master_(master) {
	}

	virtual ~AbstractConjugateGraph() {
		while (!vertices_.empty()) {
			ForceDeleteVertex(*vertices_.begin());
		}
	}

	vector<EdgeId> OutgoingEdges(VertexId v) const {
		return v->OutgoingEdges();
	}

	const vector<EdgeId> IncomingEdges(VertexId v) const {
		return v->IncomingEdges();
	}

	size_t OutgoingEdgeCount(VertexId v) const {
		return v->OutgoingEdgeCount();
	}

	size_t IncomingEdgeCount(VertexId v) const {
		return v->conjugate()->OutgoingEdgeCount();
	}

	bool CheckUniqueOutgoingEdge(VertexId v) const {
		return v->OutgoingEdgeCount() == 1;
	}

	EdgeId GetUniqueOutgoingEdge(VertexId v) const {
		assert(CheckUniqueOutgoingEdge(v));
		return (v->OutgoingEdges())[0];
	}

	bool CheckUniqueIncomingEdge(VertexId v) const {
		return CheckUniqueOutgoingEdge(v->conjugate());
	}

	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		return conjugate(GetUniqueOutgoingEdge(v->conjugate()));
	}

	vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) {
		return v->OutgoingEdgesTo(u);
	}

	const EdgeData& data(EdgeId edge) const {
		return edge->data();
	}

	const VertexData& data(VertexId v) const {
		return v->data();
	}

	std::string str(const EdgeId edge) const {
		return master_.str(data(edge));
	}

	std::string str(const VertexId v) const {
		return master_.str(data(v));
	}

	size_t length(const EdgeId edge) const {
		return master_.length(data(edge));
	}

	size_t length(const VertexId v) const {
		return master_.length(data(v));
	}

	VertexId AddVertex(const VertexData& data) {
		VertexId result = HiddenAddVertex(data);
		FireAddVertex(result);
		return result;
	}

	VertexId AddVertex(const VertexData& data1, const VertexData& data2) {
		VertexId result = HiddenAddVertex(data1, data2);
		FireAddVertex(result);
		return result;
	}

	void DeleteVertex(VertexId v) {
		TRACE("ab_conj DeleteVertex "<<v);
		assert(IsDeadEnd(v) && IsDeadStart(v));
		TRACE("ab_conj DeleteVertex "<<v);
		assert(v != NULL);
		TRACE("ab_conj DeleteVertex "<<v);
		FireDeleteVertex(v);
		TRACE("ab_conj DeleteVertex "<<v);
		VertexId conjugate = v->conjugate();
		TRACE("ab_conj DeleteVertex "<<v<<" and conj "<<conjugate);
		vertices_.erase(v);
		TRACE("ab_conj delete "<<v);
		delete v;
		TRACE("ab_conj erase "<<conjugate);
		vertices_.erase(conjugate);
		TRACE("ab_conj delete "<<conjugate);
		delete conjugate;
		TRACE("ab_conj delete FINISHED");
	}

	void ForceDeleteVertex(VertexId v) {
		TRACE("ForceDeleteVertex "<<v);
		DeleteAllOutgoing(v);
		TRACE("DeleteAllOutgoing OK");
		DeleteAllOutgoing(v->conjugate());
		TRACE("DeleteAllOutgoing  conj OK");
		DeleteVertex(v);
		TRACE("ForceDeleteVertex OK");
	}

	EdgeId AddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		EdgeId result = HiddenAddEdge(v1, v2, data);
		FireAddEdge(result);
		return result;
	}

	bool HasEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		for (auto it = v1->outgoing_edges_.begin(); it
				!= v1->outgoing_edges_.end(); ++it) {
			if (((*it)->end() == v2) && (master_.equals((*it)->data(), data))) {
				return true;
			}
		}
		return false;
	}

	EdgeId GetEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		for (auto it = v1->outgoing_edges_.begin(); it
				!= v1->outgoing_edges_.end(); ++it) {
			if (((*it)->end() == v2) && (master_.equals((*it)->data(), data))) {
				return *it;
			}
		}
		return NULL;
	}

	void DeleteEdge(EdgeId edge) {
		FireDeleteEdge(edge);
		EdgeId rcEdge = conjugate(edge);
		VertexId rcStart = conjugate(edge->end());
		VertexId start = conjugate(rcEdge->end());
		start->RemoveOutgoingEdge(edge);
		rcStart->RemoveOutgoingEdge(rcEdge);
		if (edge != rcEdge) {
			delete rcEdge;
		}
		delete edge;
	}

	bool IsDeadEnd(VertexId v) const {
		return v->IsDeadend();
	}

	bool IsDeadStart(VertexId v) const {
		return IsDeadEnd(v->conjugate());
	}

	VertexId EdgeStart(EdgeId edge) const {
		return edge->conjugate()->end()->conjugate();
	}

	VertexId EdgeEnd(EdgeId edge) const {
		return edge->end();
	}

	VertexId conjugate(VertexId v) const {
		return v->conjugate();
	}

	EdgeId conjugate(EdgeId edge) const {
		return edge->conjugate();
	}

	bool CanCompressVertex(VertexId v) const {
		return v->OutgoingEdgeCount() == 1
				&& v->conjugate()->OutgoingEdgeCount() == 1;
	}

	void CompressVertex(VertexId v) {
		//assert(CanCompressVertex(v));
		if (CanCompressVertex(v)) {
			vector<EdgeId> toMerge;
			toMerge.push_back(GetUniqueIncomingEdge(v));
			toMerge.push_back(GetUniqueOutgoingEdge(v));
			MergePath(toMerge);
		}
	}

	vector<EdgeId> CorrectMergePath(const vector<EdgeId>& path) {
		vector<EdgeId> result;
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

	EdgeId MergePath(const vector<EdgeId>& path) {
		assert(!path.empty());
		vector<EdgeId> correctedPath = CorrectMergePath(path);
		SequenceBuilder sb;
		VertexId v1 = EdgeStart(correctedPath[0]);
		VertexId v2 = EdgeEnd(correctedPath[correctedPath.size() - 1]);
		vector<EdgeData*> toMerge;
		for (auto it = correctedPath.begin(); it != correctedPath.end(); ++it) {
			toMerge.push_back(&((*it)->data()));
		}
		EdgeId newEdge = HiddenAddEdge(v1, v2, master_.MergeData(toMerge));
		FireMerge(correctedPath, newEdge);
		DeletePath(correctedPath);
		FireAddEdge(newEdge);
		return newEdge;
	}

	void DeletePath(const vector<EdgeId> &path) {
		set<EdgeId> edgesToDelete;
		set<VertexId> verticesToDelete;
		edgesToDelete.insert(path[0]);
		for (size_t i = 0; i + 1 < path.size(); i++) {
			EdgeId e = path[i + 1];
			if (edgesToDelete.find(conjugate(e)) == edgesToDelete.end())
				edgesToDelete.insert(e);
			VertexId v = EdgeStart(e);
			if (verticesToDelete.find(conjugate(v)) == verticesToDelete.end())
				verticesToDelete.insert(v);
		}
		for (auto it = edgesToDelete.begin(); it != edgesToDelete.end(); ++it)
			DeleteEdge(*it);
		for (auto it = verticesToDelete.begin(); it != verticesToDelete.end(); ++it)
			DeleteVertex(*it);
	}

	pair<EdgeId, EdgeId> SplitEdge(EdgeId edge, size_t position) {
		assert(edge != conjugate(edge));
		pair<VertexData, pair<EdgeData, EdgeData>> newData = master_.SplitData(
				edge->data(), position);
		VertexId splitVertex = HiddenAddVertex(newData.first);
		EdgeId newEdge1 = HiddenAddEdge(this->EdgeStart(edge), splitVertex,
				newData.second.first);
		EdgeId newEdge2 = HiddenAddEdge(splitVertex, this->EdgeEnd(edge),
				newData.second.second);
		FireSplit(edge, newEdge1, newEdge2);
		FireAddVertex(splitVertex);
		FireAddEdge(newEdge1);
		FireAddEdge(newEdge2);
		DeleteEdge(edge);
		return make_pair(newEdge1, newEdge2);
	}

	void GlueEdges(EdgeId edge1, EdgeId edge2) {
		FireDeleteEdge(edge2);
		FireGlue(edge1, edge2);
		edge2->set_data(master_.GlueData(edge1->data(), edge2->data()));
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

private:
	DECL_LOGGER("AbstractConjugateGraph")
};

}
#endif /* ABSTRUCT_CONJUGATE_GRAPH_HPP_ */
