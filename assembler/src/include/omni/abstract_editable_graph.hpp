#pragma once

#include <vector>
#include <set>
#include <cstring>
#include "sequence/seq.hpp"
#include "sequence/sequence.hpp"
#include "logging.hpp"
#include "sequence/nucl.hpp"
#include "omni_utils.hpp"
#include "observable_graph.hpp"

namespace omnigraph {

template<typename VertexIdT, typename EdgeIdT, class DataMasterT,
		typename VertexIt>
class AbstractEditableGraph: public ObservableGraph<VertexIdT, EdgeIdT,
		VertexIt> {
	typedef ObservableGraph<VertexIdT, EdgeIdT, VertexIt> base;
	//todo maybe rename template params themselves???
public:
	typedef VertexIdT VertexId;
	typedef EdgeIdT EdgeId;
	typedef DataMasterT DataMaster;
	typedef typename DataMaster::VertexData VertexData;
	typedef typename DataMaster::EdgeData EdgeData;
	typedef VertexIt VertexIterator;
private:
	//todo think of necessity to pull these typedefs through hierarchy
	DataMaster master_;

	virtual VertexId HiddenAddVertex(const VertexData &data) = 0;

	virtual void HiddenDeleteVertex(VertexId v) = 0;

	virtual EdgeId
	HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data) = 0;

	virtual void HiddenDeleteEdge(EdgeId edge) = 0;

	virtual vector<EdgeId> CorrectMergePath(const vector<EdgeId>& path) = 0;

	virtual vector<EdgeId> EdgesToDelete(const vector<EdgeId> &path) = 0;

	virtual vector<VertexId> VerticesToDelete(const vector<EdgeId> &path) = 0;

	void DeleteAllOutgoing(VertexId v) {
		TRACE("DeleteAllOutgoing " << OutgoingEdgeCount(v));
		while (OutgoingEdgeCount(v) > 0) {
			EdgeId edge = OutgoingEdges(v)[0];
			TRACE("DeleteOutgoing " << edge);
			DeleteEdge(edge);
			TRACE("DeleteOutgoing ok");
		}
	}

	void DeleteAllIncoming(VertexId v) {
		TRACE("DeleteAllIncoming " << IncomingEdgeCount(v));
		while (IncomingEdgeCount(v) > 0) {
			EdgeId edge = IncomingEdges(v)[0];
			TRACE("DeleteIncoming " << edge);
			DeleteEdge(edge);
			TRACE("DeleteIncoming ok");
		}
		TRACE("DeleteAllIncoming ok");
	}

	void FireDeletePath(const vector<EdgeId> &edgesToDelete, const vector<
			VertexId> &verticesToDelete) {
		for (auto it = edgesToDelete.begin(); it != edgesToDelete.end(); ++it)
			FireDeleteEdge(*it);
		for (auto it = verticesToDelete.begin(); it != verticesToDelete.end(); ++it)
			FireDeleteVertex(*it);
	}

	void HiddenDeletePath(const vector<EdgeId> &edgesToDelete, const vector<
			VertexId> &verticesToDelete) {
		for (auto it = edgesToDelete.begin(); it != edgesToDelete.end(); ++it)
			HiddenDeleteEdge(*it);
		for (auto it = verticesToDelete.begin(); it != verticesToDelete.end(); ++it)
			HiddenDeleteVertex(*it);
	}

public:
	typedef typename base::SmartVertexIt SmartVertexIt;
	typedef typename base::SmartEdgeIt SmartEdgeIt;

	AbstractEditableGraph(HandlerApplier<VertexId, EdgeId>* applier,
			const DataMaster& master) :
		base(applier), master_(master) {
	}

	virtual ~AbstractEditableGraph() {
		TRACE("~AbstractEditableGraph");
		//		doesn't work this way because call to virtual function is needed
		//		for (auto it = this->SmartVertexBegin(); !it.IsEnd(); ++it) {
		//			ForceDeleteVertex(*it);
		//		}
	}

	const DataMaster& master() const {
		return master_;
	}

	virtual const EdgeData& data(EdgeId edge) const = 0;
	virtual const VertexData& data(VertexId v) const = 0;

	virtual const vector<EdgeId> OutgoingEdges(VertexId v) const = 0;

	virtual const vector<EdgeId> IncomingEdges(VertexId v) const = 0;

	virtual size_t OutgoingEdgeCount(VertexId v) const = 0;

	virtual size_t IncomingEdgeCount(VertexId v) const = 0;

	virtual vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) const = 0;

	virtual VertexId EdgeStart(EdgeId edge) const = 0;

	virtual VertexId EdgeEnd(EdgeId edge) const = 0;

	virtual bool RelatedVertices(VertexId v1, VertexId v2) const = 0;

	bool CheckUniqueOutgoingEdge(VertexId v) const {
		return OutgoingEdgeCount(v) == 1;
	}

	EdgeId GetUniqueOutgoingEdge(VertexId v) const {
		VERIFY(CheckUniqueOutgoingEdge(v));
		return OutgoingEdges(v)[0];
	}

	bool CheckUniqueIncomingEdge(VertexId v) const {
		return IncomingEdgeCount(v) == 1;
	}

	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		VERIFY(CheckUniqueIncomingEdge(v));
		return IncomingEdges(v)[0];
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
		TRACE("Adding vertex");
		VertexId v = HiddenAddVertex(data);
		FireAddVertex(v);
		TRACE("Vertex " << PrintVertex(v) << " added");
		return v;
	}

	void DeleteVertex(VertexId v) {
		VERIFY(IsDeadEnd(v) && IsDeadStart(v));
		VERIFY(v != NULL);
		TRACE("Deleting vertex " << PrintVertex(v));
		FireDeleteVertex(v);
		HiddenDeleteVertex(v);
		TRACE("Vertex " << v << " deleted");
	}

	void ForceDeleteVertex(VertexId v) {
		TRACE("Forcing deletion of vertex " << PrintVertex(v));
		DeleteAllOutgoing(v);
		DeleteAllIncoming(v);
		DeleteVertex(v);
		TRACE("Vertex " << v << " force-deleted");
	}

	EdgeId AddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		TRACE("Adding edge connecting " << v1 << " and " << v2)
		EdgeId e = HiddenAddEdge(v1, v2, data);
		FireAddEdge( e);
		TRACE("Added edge " << PrintEdge(e) << " connecting " << v1 << " and "
				<< v2);
		return e;
	}

	//todo delete if not used
	bool HasEdge(VertexId v1, VertexId v2, const EdgeData &data) {
		vector<EdgeId> out = OutgoingEdges(v1);
		for (auto it = out.begin(); it != out.end(); ++it) {
			if ((EdgeEnd(*it) == v2) && (master_.equals(data(*it), data))) {
				return true;
			}
		}
		return false;
	}

	//todo delete if not used
	EdgeId GetEdge(VertexId v1, VertexId v2, const EdgeData &edgeData) {
		vector<EdgeId> out = OutgoingEdges(v1);
		for (auto it = out.begin(); it != out.end(); ++it) {
			if ((EdgeEnd(*it) == v2) && (master_.equals(data(*it), edgeData))) {
				return *it;
			}
		}
		return NULL;
	}

	void DeleteEdge(EdgeId e) {
		TRACE("Deleting edge " << PrintEdge(e));
		FireDeleteEdge(e);
		HiddenDeleteEdge(e);
		TRACE("Edge " << e << " deleted");
	}

	bool IsDeadEnd(VertexId v) const {
		return OutgoingEdgeCount(v) == 0;
	}

	bool IsDeadStart(VertexId v) const {
		return IncomingEdgeCount(v) == 0;
	}

	virtual bool AdditionalCompressCondition(VertexId v) const {
		return true;
	}

	bool CanCompressVertex(VertexId v) const {
		return OutgoingEdgeCount(v) == 1 && IncomingEdgeCount(v) == 1 /*one-in one-out*/
		&& GetUniqueOutgoingEdge(v) != GetUniqueIncomingEdge(v) /*not loop*/;
	}

	void CompressVertex(VertexId v) {
		TRACE("Trying to compress vertex " << v);
		//VERIFY(CanCompressVertex(v));
		if (CanCompressVertex(v) && AdditionalCompressCondition(v)) {
			TRACE("Compressing vertex " << v);
			vector<EdgeId> edges_to_merge;
			edges_to_merge.push_back(GetUniqueIncomingEdge(v));
			edges_to_merge.push_back(GetUniqueOutgoingEdge(v));
			MergePath(edges_to_merge);
			TRACE("Vertex compressed");
		} else {
			TRACE("Vertex " << v << " can't be compressed");
		}
	}

	//todo remove after debug
	virtual std::string PrintDetailedPath(const vector<EdgeId>& path) const {
		VERIFY(false);
		return "";
	}

	//todo remove after debug
	std::string PrintConjugatePath(const vector<EdgeId>& path) const {
		vector<EdgeId> conjugate_path;
		for (int i = path.size() - 1; i >= 0; --i) {
			conjugate_path.push_back(path[i]);
		}
		return PrintDetailedPath(conjugate_path);
	}

	//todo remove after debug
	virtual std::string PrintDetailedVertexInfo(VertexId v) const {
		VERIFY(false);
		return "";
	}

	//todo remove after debug
	virtual std::string PrintEdges(const vector<EdgeId>& path) const {
		VERIFY(false);
		return "";
	}

	//todo remove after debug
	template<class T>
	std::string SimplePrint(const vector<T>& v) const {
		stringstream ss;
		for (auto it = v.begin(); it != v.end(); ++it) {
			ss << *it << ", ";
		}
		return ss.str();
	}

	//todo remove after debug
	virtual std::string PrintEdge(EdgeId e) const {
		return "";
	}

	//todo remove after debug
	virtual std::string PrintVertex(VertexId v) const {
		return "";
	}

	//todo remove after debug
	std::string PrintVertices(const vector<VertexId>& path) const {
		stringstream ss;
		ss << "Vertices ";
		for (auto it = path.begin(); it != path.end(); ++it) {
			ss << *it << ", ";
		}
		return ss.str();
	}

	EdgeId MergePath(const vector<EdgeId>& path) {
		VERIFY(!path.empty());
		for (size_t i = 0; i < path.size(); i++)
			for (size_t j = i + 1; j < path.size(); j++) {
				VERIFY(path[i] != path[j]);
			}
		if (path.size() == 1) {
			TRACE("Path of single edge " << PrintEdge(*(path.begin()))
					<< ". Nothing to merge.");
		}
		TRACE("Merging path " << PrintEdges(path));

		//		cerr << "Merging " << PrintDetailedPath(pObservableGraph<VertexIdT, EdgeIdT, VertexIt>ath) << endl;
		//		cerr << "Conjugate " << PrintConjugatePath(path) << endl;
		vector<EdgeId> corrected_path = CorrectMergePath(path);
		TRACE("Corrected path " << PrintEdges(corrected_path));
		VertexId v1 = EdgeStart(corrected_path[0]);
		VertexId v2 = EdgeEnd(corrected_path[corrected_path.size() - 1]);
		vector<const EdgeData*> to_merge;
		for (auto it = corrected_path.begin(); it != corrected_path.end(); ++it) {
			to_merge.push_back(&(data(*it)));
		}
		EdgeId new_edge = HiddenAddEdge(v1, v2, master_.MergeData(to_merge));
		FireMerge(corrected_path, new_edge);

		//		cerr << "Corrected " << PrintDetailedPath(corrected_path) << endl;
		//		cerr << "Corrected conjugate " << PrintConjugatePath(corrected_path) << endl;
		vector<EdgeId> edges_to_delete = EdgesToDelete(corrected_path);
		//		cerr << "To delete " << PrintEdges(edges_to_delete) << endl;
		vector<VertexId> vertices_to_delete = VerticesToDelete(corrected_path);
		//		cerr << "To delete " << PrintVertices(vertices_to_delete) << endl;

		FireDeletePath(edges_to_delete, vertices_to_delete);
		FireAddEdge(new_edge);
		HiddenDeletePath(edges_to_delete, vertices_to_delete);
		TRACE("Path merged. Corrected path " << SimplePrint(corrected_path)
				<< "merged into " << PrintEdge(new_edge));
		return new_edge;
	}

	pair<EdgeId, EdgeId> SplitEdge(EdgeId edge, size_t position) {
		TRACE("Splitting edge " << PrintEdge(edge));
		pair<VertexData, pair<EdgeData, EdgeData>> newData = master_.SplitData(
				data(edge), position);
		VertexId splitVertex = HiddenAddVertex(newData.first);
		EdgeId new_edge1 = HiddenAddEdge(EdgeStart(edge), splitVertex,
				newData.second.first);
		EdgeId new_edge2 = HiddenAddEdge(splitVertex, EdgeEnd(edge),
				newData.second.second);
		FireSplit(edge, new_edge1, new_edge2);
		FireDeleteEdge(edge);
		FireAddVertex(splitVertex);
		FireAddEdge(new_edge1);

		FireAddEdge(new_edge2);
		HiddenDeleteEdge(edge);
		TRACE("Edge " << edge << " split into " << PrintEdge(new_edge1) <<" and " << PrintEdge(new_edge2));
		return make_pair(new_edge1, new_edge2);
	}

				void GlueEdges(EdgeId edge1, EdgeId edge2) {
					TRACE("Gluing edges " << PrintEdge(edge1) << " and " <<PrintEdge(edge2));
	EdgeId new_edge = HiddenAddEdge(EdgeStart(edge2), EdgeEnd(edge2),
			master_.GlueData(data(edge1), data(edge2)));
	FireGlue(new_edge, edge1, edge2);
	FireDeleteEdge(edge1);
	FireDeleteEdge(edge2);
	FireAddEdge(new_edge);
	VertexId start = EdgeStart(edge1);
	VertexId end = EdgeEnd(edge1);
	HiddenDeleteEdge(edge1);
	HiddenDeleteEdge(edge2);
	if (IsDeadStart(start) && IsDeadEnd(start)) {
		DeleteVertex(start);
	}
	if (IsDeadStart(end) && IsDeadEnd(end)) {
		DeleteVertex(end);
	}
	TRACE("Edges " << edge1 << " and " << edge2 << " glued into " << PrintEdge(new_edge));
}

private:
DECL_LOGGER("AbstractEditableGraph")
};

}

