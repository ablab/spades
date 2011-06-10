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

namespace omnigraph {

template<typename VertexData, typename EdgeData, class DataMaster>
class AbstractConjugateGraph {
public:
	class Vertex;
	typedef Vertex* VertexId;
	typedef set<VertexId> Vertices;
	typedef typename Vertices::const_iterator VertexIterator;
	class Edge;
	typedef Edge* EdgeId;
	typedef vector<EdgeId> Edges;
	typedef typename Edges::const_iterator EdgeIterator;
	typedef ActionHandler<VertexId, EdgeId> Handler;

	class Vertex {
	private:
		friend class AbstractConjugateGraph<VertexData, EdgeData, DataMaster> ;

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
				result[i] = result[i].conjugate();
			}
			return result;
		}

		Vertex(VertexData data) :
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

		~Vertex() {
			assert(outgoing_edges_.size() == 0);
		}
	};

	class Edge {
	private:
		friend class AbstractConjugateGraph<VertexData, EdgeData, DataMaster> ;
		VertexId end_;

		EdgeData data_;

		EdgeId conjugate_;

		Edge(Vertex* end, const EdgeData &data) :
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

		~Edge() {
		}
	};

	const PairedHandlerApplier<AbstractConjugateGraph<VertexData, EdgeData,
			DataMaster> > applier_;

	vector<Handler*> action_handler_list_;

	set<Vertex*> vertices_;

	DataMaster master_;

	VertexId HiddenAddVertex(const VertexData &data1, const VertexData &data2) {
		VertexId v1 = new Vertex(data1);
		VertexId v2 = new Vertex(data2);
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
		EdgeId newEdge = new Edge(v2, data);
		v1->AddOutgoingEdge(newEdge);
		return newEdge;
	}

	void DeleteAllOutgoing(VertexId v) {
		vector<EdgeId> out = v->outgoing_edges_;
		for (auto it = out.begin(); it != out.end(); ++it) {
			DeleteEdge(*it);
		}
	}

	void FireAddVertex(VertexId v) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyAdd(*it, v);
		}
	}

	void FireAddEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyAdd(*it, edge);
		}
	}

	void FireDeleteVertex(VertexId v) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyDelete(*it, v);
		}
	}

	void FireDeleteEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyDelete(*it, edge);
		}
	}

	void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyMerge(*it, oldEdges, newEdge);
		}
	}

	void FireGlue(EdgeId edge1, EdgeId edge2) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyGlue(*it, edge1, edge2);
		}
	}

	void FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplySplit(*it, edge, newEdge1, newEdge2);
		}
	}

public:

	VertexIterator begin() const {
		return vertices_.begin();
	}

	VertexIterator end() const {
		return vertices_.end();
	}

	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<AbstractConjugateGraph, Comparator> SmartVertexBegin(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<AbstractConjugateGraph, Comparator> (*this,
				true, comparator);
	}

	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<AbstractConjugateGraph, Comparator> SmartVertexEnd(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<AbstractConjugateGraph, Comparator> (*this,
				false, comparator);
	}

	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<AbstractConjugateGraph, Comparator> SmartEdgeBegin(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<AbstractConjugateGraph, Comparator> (*this,
				true, comparator);
	}

	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<AbstractConjugateGraph, Comparator> SmartEdgeEnd(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<AbstractConjugateGraph, Comparator> (*this,
				false, comparator);
	}

	size_t size() {
		return vertices_.size();
	}

	AbstractConjugateGraph(DataMaster master) :
		applier_(*this), master_(master) {
	}

	~AbstractConjugateGraph() {
		while (!vertices_.empty()) {
			ForceDeleteVertex(*vertices_.begin());
		}
	}

	void AddActionHandler(Handler* action_handler) {
		TRACE("Action handler added");
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(Handler* action_handler) {
		TRACE("Trying to remove action handler");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				action_handler_list_.erase(it);
				TRACE("Action handler removed");
				return true;
			}
		}
		return false;
	}

	const vector<EdgeId> OutgoingEdges(VertexId v) const {
		return v->OutgoingEdges();
	}

	const vector<EdgeId> IncomingEdges(VertexId v) const {
		return v->IncomingEdges();
	}

	size_t OutgoingEdgeCount(VertexId v) const {
		return v->OutgoingEdgeCount();
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
		assert(IsDeadEnd(v) && IsDeadStart(v));
		assert(v != NULL);
		FireDeleteVertex(v);
		VertexId conjugate = v->conjugate();
		vertices_.erase(v);
		delete v;
		vertices_.erase(conjugate);
		delete conjugate;
	}

	void ForceDeleteVertex(VertexId v) {
		DeleteAllOutgoing(v);
		DeleteAllOutgoing(v->conjugate());
		DeleteVertex(v);
	}

	EdgeId AddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		EdgeId result = HiddenAddEdge(v1, v2, data);
		FireAddEdge(result);
		return result;
	}

	bool HasEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		for (auto it = v1->outgoing_edges_.begin(); it != v1->outgoing_edges_.end(); ++it) {
			if (((*it)->end() == v2) && (master_.equals((*it)->data(), data))) {
				return true;
			}
		}
		return false;
	}

	EdgeId GetEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		for (auto it = v1->outgoing_edges_.begin(); it != v1->outgoing_edges_.end(); ++it) {
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
			Merge(GetUniqueIncomingEdge(v), GetUniqueOutgoingEdge(v));
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
