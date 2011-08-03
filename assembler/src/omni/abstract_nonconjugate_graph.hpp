#ifndef ABSTRACT_GRAPH_HPP_
#define ABSTRACT_GRAPH_HPP_

#include <vector>
#include <set>
#include <cstring>
#include "seq.hpp"
#include "sequence.hpp"
#include "logging.hpp"
#include "nucl.hpp"
//#include "strobe_read.hpp"
#include "common/io/paired_read.hpp"
#include "omni_utils.hpp"
#include "observable_graph.hpp"

namespace omnigraph {

template<typename VertexData, typename EdgeData, class DataMaster>
class AbstractNonconjugateGraph;

template<typename VertexData, typename EdgeData, class DataMaster>
class SingleEdge;

template<typename VertexData, typename EdgeData, class DataMaster>
class SingleVertex {
private:
	typedef SingleVertex<VertexData, EdgeData, DataMaster>* VertexId;
	typedef SingleEdge<VertexData, EdgeData, DataMaster>* EdgeId;

	friend class AbstractNonconjugateGraph<VertexData, EdgeData, DataMaster> ;

	vector<EdgeId> outgoing_edges_;

	vector<EdgeId> incoming_edges_;

	VertexData data_;

	size_t OutgoingEdgeCount() const {
		return outgoing_edges_.size();
	}

	const vector<EdgeId> OutgoingEdges() const {
		return outgoing_edges_;
	}

	size_t IncomingEdgeCount() const {
		return incoming_edges_.size();
	}

	const vector<EdgeId> IncomingEdges() const {
		return incoming_edges_;
	}

	SingleVertex(VertexData data) :
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

	bool IsDeadstart() {
		return incoming_edges_.size() == 0;
	}

	void AddIncomingEdge(EdgeId e) {
		incoming_edges_.push_back(e);
	}

	bool RemoveIncomingEdge(const EdgeId e) {
		auto it = incoming_edges_.begin();
		while (it != incoming_edges_.end() && *it != e) {
			++it;
		}
		if (it == incoming_edges_.end()) {
			return false;
		}
		incoming_edges_.erase(it);
		return true;
	}

	~SingleVertex() {
		assert(outgoing_edges_.size() == 0);
	}
};

template<typename VertexData, typename EdgeData, class DataMaster>
class SingleEdge {
private:
	typedef SingleVertex<VertexData, EdgeData, DataMaster>* VertexId;
	typedef SingleEdge<VertexData, EdgeData, DataMaster>* EdgeId;

	friend class AbstractNonconjugateGraph<VertexData, EdgeData, DataMaster> ;

	VertexId start_;
	VertexId end_;

	EdgeData data_;

	SingleEdge(VertexId start, VertexId end, const EdgeData &data) :
		start_(start), end_(end), data_(data) {
	}

	EdgeData &data() {
		return data_;
	}

	void set_data(EdgeData &data) {
		data_ = data;
	}

	VertexId start() const {
		return start_;
	}

	VertexId end() const {
		return end_;
	}

	~SingleEdge() {
	}
};

template<typename VertexData, typename EdgeData, class DataMaster>
class AbstractNonconjugateGraph: public ObservableGraph<SingleVertex<
		VertexData, EdgeData, DataMaster>*, SingleEdge<VertexData, EdgeData,
		DataMaster>*> {
public:
	typedef SingleVertex<VertexData, EdgeData, DataMaster>* VertexId;
	typedef set<VertexId> Vertices;
	typedef typename Vertices::const_iterator VertexIterator;
	typedef SingleEdge<VertexData, EdgeData, DataMaster>* EdgeId;
	typedef vector<EdgeId> Edges;
	typedef typename Edges::const_iterator EdgeIterator;

	set<VertexId> vertices_;

	DataMaster master_;

	VertexId HiddenAddVertex(const VertexData &data1) {
		VertexId v = new SingleVertex<VertexData, EdgeData, DataMaster> (data1);
		vertices_.insert(v);
		return v;
	}

	EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		assert(vertices_.find(v1) != vertices_.end() && vertices_.find(v2) != vertices_.end());
		EdgeId newEdge = new SingleEdge<VertexData, EdgeData, DataMaster> (v1,
				v2, data);
		v1->AddOutgoingEdge(newEdge);
		v2->AddIncomingEdge(newEdge);
		return newEdge;
	}

	void DeleteAllOutgoing(VertexId v) {
		vector < EdgeId > out = v->outgoing_edges_;
		for (auto it = out.begin(); it != out.end(); ++it) {
			DeleteEdge(*it);
		}
	}

	void DeleteAllIncoming(VertexId v) {
		vector < EdgeId > out = v->incoming_edges_;
		for (auto it = out.begin(); it != out.end(); ++it) {
			DeleteEdge(*it);
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

	AbstractNonconjugateGraph(DataMaster master) :
				ObservableGraph<VertexId, EdgeId> (
						new SimpleHandlerApplier<AbstractNonconjugateGraph<
								VertexData, EdgeData, DataMaster> > ()),
				master_(master) {
	}

	virtual ~AbstractNonconjugateGraph() {
		while (!vertices_.empty()) {
			ForceDeleteVertex(*vertices_.begin());
		}
	}

	size_t length(const EdgeId edge) const {
		return master_.length(data(edge));
	}

	vector<EdgeId> OutgoingEdges(VertexId v) const {
		return v->OutgoingEdges();
	}

	vector<EdgeId> IncomingEdges(VertexId v) const {
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

	size_t IncomingEdgeCount(VertexId v) const {
		return v->IncomingEdgeCount();
	}

	bool CheckUniqueIncomingEdge(VertexId v) const {
		return v->IncomingEdgeCount() == 1;
	}

	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		assert(CheckUniqueIncomingEdge(v));
		return (v->IncomingEdges())[0];
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

	void DeleteVertex(VertexId v) {
		assert(IsDeadEnd(v) && IsDeadStart(v));
		assert(v != NULL);
		FireDeleteVertex(v);
		vertices_.erase(v);
		delete v;
	}

	void ForceDeleteVertex(VertexId v) {
		DeleteAllOutgoing(v);
		DeleteAllIncoming(v);
		DeleteVertex(v);
	}

	virtual EdgeId AddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		EdgeId result = HiddenAddEdge(v1, v2, data);
		FireAddEdge(result);
		return result;
	}

	void DeleteEdge(EdgeId edge) {
		FireDeleteEdge(edge);
		VertexId start = edge->start();
		VertexId end = edge->end();
		start->RemoveOutgoingEdge(edge);
		end->RemoveIncomingEdge(edge);
		delete edge;
	}

	bool IsDeadEnd(VertexId v) const {
		return v->IsDeadend();
	}

	bool IsDeadStart(VertexId v) const {
		return v->IsDeadstart();
	}

	VertexId EdgeStart(EdgeId edge) const {
		return edge->start();
	}

	VertexId EdgeEnd(EdgeId edge) const {
		return edge->end();
	}

	bool CanCompressVertex(VertexId v) const {
		return v->OutgoingEdgeCount() == 1 && v->IncomingEdgeCount() == 1;
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

	EdgeId MergePath(const vector<EdgeId>& path) {
		assert(!path.empty());
		SequenceBuilder sb;
		VertexId v1 = EdgeStart(path[0]);
		VertexId v2 = EdgeEnd(path[path.size() - 1]);
		vector<EdgeData*> toMerge;
		for (auto it = path.begin(); it != path.end(); ++it) {
			toMerge.push_back(&((*it)->data()));
		}
		EdgeId newEdge = HiddenAddEdge(v1, v2, master_.MergeData(toMerge));
		FireMerge(path, newEdge);
		DeletePath(path);
		FireAddEdge(newEdge);
		return newEdge;
	}

	void DeletePath(const vector<EdgeId> &path) {
		set < EdgeId > edgesToDelete;
		set < VertexId > verticesToDelete;
		edgesToDelete.insert(path[0]);
		for (size_t i = 1; i < path.size(); i++) {
			edgesToDelete.insert(path[i]);
			verticesToDelete.insert(EdgeStart(path[i]));
		}
		for (auto it = edgesToDelete.begin(); it != edgesToDelete.end(); ++it)
			DeleteEdge(*it);
		for (auto it = verticesToDelete.begin(); it != verticesToDelete.end(); ++it)
			DeleteVertex(*it);
	}

	pair<EdgeId, EdgeId> SplitEdge(EdgeId edge, size_t position) {
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
		EdgeId newEdge = HiddenAddEdge(EdgeStart(edge2), EdgeEnd(edge2), master_.GlueData(edge1->data(), edge2->data()));
		FireGlue(newEdge, edge1, edge2);
		FireDeleteEdge(edge1);
		FireDeleteEdge(edge2);
		FireAddEdge(newEdge);
		VertexId start = EdgeStart(edge1);
		VertexId end = EdgeEnd(edge1);
		DeleteEdge(edge1);
		if (IsDeadStart(start) && IsDeadEnd(start)) {
			DeleteVertex(start);
		}
		if (IsDeadStart(end) && IsDeadEnd(end)) {
			DeleteVertex(end);
		}
/*


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
*/
	}

private:
	DECL_LOGGER("AbstractNonconjugateGraph")
};

}
#endif /* ABSTRACT_GRAPH_HPP_ */
