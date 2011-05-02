#ifndef EDGE_GRAPH_HPP_
#define EDGE_GRAPH_HPP_

#include <vector>
#include <set>
//#include <ext/hash_map>
#include <tr1/unordered_map>
#include <cstring>
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "sequence.hpp"
#include "logging.hpp"
#include "nucl.hpp"
#include "debruijn.hpp"
#include "strobe_read.hpp"
#include "utils.hpp"

using namespace std;

namespace edge_graph {
LOGGER("d.edge_graph");

class Vertex;

class Edge {
//public:
private:
	const Sequence& nucls() const {
		return nucls_;
	}
	friend class EdgeGraph;
	Sequence nucls_;
	Vertex* end_;
	size_t coverage_;
	Edge *complement_;

	Edge(const Sequence& nucls, Vertex* end, size_t coverage) :
		nucls_(nucls), end_(end), coverage_(coverage) {
	}

	Vertex* end() const {
		return end_;
	}

	size_t size() const {
		return nucls_.size();
	}

	Edge *Complement() {
		return complement_;
	}

	void SetComplement(Edge* complement) {
		complement_ = complement;
	}

	~Edge() {
	}

};

class Vertex {
public:
	typedef vector<Edge*>::const_iterator EdgeIterator;
private:
	friend class EdgeGraph;

	vector<Edge*> outgoing_edges_;

	Vertex* complement_;

	void SetComplement(Vertex* complement) {
		complement_ = complement;
	}

	EdgeIterator begin() const {
		return outgoing_edges_.begin();
	}

	EdgeIterator end() const {
		return outgoing_edges_.end();
	}

	Vertex() {
	}

	size_t OutgoingEdgeCount() const {
		return outgoing_edges_.size();
	}

	bool IsDeadend() {
		return outgoing_edges_.size() == 0;
	}

	void AddOutgoingEdge(Edge* e) {
		outgoing_edges_.push_back(e);
	}

	bool RemoveOutgoingEdge(const Edge* e) {
		vector<Edge *>::iterator it = outgoing_edges_.begin();
		while (it != outgoing_edges_.end() && *it != e) {
			++it;
		}
		if (it == outgoing_edges_.end()) {
			return false;
		}
		outgoing_edges_.erase(it);
		return true;
	}

	Vertex* Complement() const {
		return complement_;
	}

	~Vertex() {
		assert(outgoing_edges_.size() == 0);
	}

};

using de_bruijn::SmartVertexIterator;
using de_bruijn::SmartEdgeIterator;
using de_bruijn::HandlerApplier;
using de_bruijn::PairedHandlerApplier;

class EdgeGraph {
public:
	typedef Edge* EdgeId;
	typedef Edge EdgeData;
	typedef Vertex* VertexId;

	typedef set<Vertex*>::const_iterator VertexIterator;
	typedef Vertex::EdgeIterator EdgeIterator;
	typedef de_bruijn::GraphActionHandler<EdgeGraph> ActionHandler;
//	typedef de_bruijn::PairedActionHandler<EdgeGraph> PairedActionHandler;
	//	typedef de_bruijn::SmartVertexIterator<EdgeGraph> SmartVertexIterator;
	//	typedef de_bruijn::SmartEdgeIterator<EdgeGraph> SmartEdgeIterator;

private:
	const size_t k_;

	const HandlerApplier<EdgeGraph> *applier_;

	vector<ActionHandler*> action_handler_list_;

	set<Vertex*> vertices_;

	VertexId HiddenAddVertex();

	EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const Sequence &nucls,
			size_t coverage = 0);

	EdgeId AddSingleEdge(VertexId v1, VertexId v2, const Sequence& s,
			size_t coverage);

	void DeleteAllOutgoing(Vertex* v);

	bool GoUniqueWay(VertexId &v);

	void FireAddVertex(VertexId v);
	void FireAddEdge(EdgeId edge);
	void FireDeleteVertex(VertexId v);
	void FireDeleteEdge(EdgeId edge);
	void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge);
	void FireGlue(EdgeId edge1, EdgeId edge2);
	void FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2);

public:
	VertexIterator begin() const {
		return vertices_.begin();
	}

	VertexIterator end() const {
		return vertices_.end();
	}

	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<EdgeGraph, Comparator> SmartVertexBegin(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<EdgeGraph, Comparator> (*this, true,
				comparator);
	}

	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<EdgeGraph, Comparator> SmartVertexEnd(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<EdgeGraph, Comparator> (*this, false,
				comparator);
	}

	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<EdgeGraph, Comparator> SmartEdgeBegin(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<EdgeGraph, Comparator> (*this, true,
				comparator);
	}

	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<EdgeGraph, Comparator> SmartEdgeEnd(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<EdgeGraph, Comparator> (*this, false,
				comparator);
	}

	size_t size() {
		return vertices_.size();
	}

	/**
	 * Constructs empty graph to work with k-mers.
	 *
	 * @param k Main parameter that defines the size of k-mers
	 * //@param action_handler Graph actions handler
	 */
	EdgeGraph(size_t k) : k_(k), applier_(new PairedHandlerApplier<EdgeGraph>(*this)) {
		assert(k % 2 == 1);
	}

	/**
	 * Deletes action_handler.
	 */
	~EdgeGraph() {
		while (!vertices_.empty()) {
			ForceDeleteVertex(*vertices_.begin());
		}
		delete applier_;
	}

	size_t k() {
		return k_;
	}

	void AddActionHandler(ActionHandler* action_handler) {
		DEBUG("Action handler added");
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(ActionHandler* action_handler) {
		DEBUG("Trying to remove action handler");
		for (vector<ActionHandler*>::iterator it =
				action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				delete *it;
				action_handler_list_.erase(it);
				DEBUG("Action handler removed");
				return true;
			}
		}
		//		assert(false);
		return false;
	}

	//	//todo remove
	//	const vector<ActionHandler*> GetHandlers() {
	//		return action_handler_list_;
	//	}

public:
	void OutgoingEdges(VertexId v, EdgeIterator& begin, EdgeIterator& end) const;

	const vector<EdgeId> OutgoingEdges(VertexId v) const;

	const vector<EdgeId> IncomingEdges(VertexId v) const;

	EdgeId OutgoingEdge(VertexId v, char nucl) const;

	size_t OutgoingEdgeCount(VertexId v) const {
		return v->OutgoingEdgeCount();
	}

	size_t IncomingEdgeCount(VertexId v) const {
		return v->Complement()->OutgoingEdgeCount();
	}

	bool CheckUniqueOutgiongEdge(VertexId v) const {
		return v->OutgoingEdgeCount() == 1;
	}

	EdgeId GetUniqueOutgoingEdge(VertexId v) const {
		assert(CheckUniqueOutgiongEdge(v));
		return *(v->begin());
	}

	bool CheckUniqueIncomingEdge(const VertexId v) const {
		return CheckUniqueOutgiongEdge(v->Complement());
	}

	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		return Complement(GetUniqueOutgoingEdge(v->Complement()));
	}

	//	Edge* ComplementEdge(const Edge* edge) const;

	const Sequence& EdgeNucls(EdgeId edge) const {
		return edge->nucls();
	}

	void SetCoverage(EdgeId edge, size_t cov) {
		edge->coverage_ = cov;
	}

	double coverage(EdgeId edge) const {
		return (double) edge->coverage_ / length(edge);
	}

//	size_t KPlusOneMerCoverage(EdgeId edge) const {
//		return edge->coverage_;
//	}

	void IncCoverage(EdgeId edge, int toAdd) {
		edge->coverage_ += toAdd;
		EdgeId rc = Complement(edge);
		if (edge != rc) {
			rc->coverage_ += toAdd;
		}
	}

	void IncCoverage(EdgeId edge) {
		edge->coverage_++;
		EdgeId rc = Complement(edge);
		if (edge != rc) {
			rc->coverage_++;
		}
	}

	/**
	 * adds vertex and its complement
	 */
	VertexId AddVertex();

	Sequence VertexNucls(VertexId v) const;

	/**
	 * deletes vertex and its complement
	 */
	void DeleteVertex(VertexId v);

	void ForceDeleteVertex(VertexId v);

	Edge* AddEdge(VertexId v1, VertexId v2, const Sequence &nucls,
			size_t coverage = 0);

	void DeleteEdge(EdgeId edge);

	size_t length(EdgeId edge) const {
		return edge->nucls_.size() - k_;
	}

	bool AreLinkable(VertexId v1, VertexId v2, const Sequence &nucls) const;

	bool IsDeadEnd(VertexId v) const {
		return v->IsDeadend();
	}

	bool IsDeadStart(VertexId v) const {
		return IsDeadEnd(v->Complement());
	}

	VertexId EdgeStart(EdgeId edge) const;

	VertexId EdgeEnd(EdgeId edge) const;

	VertexId Complement(VertexId v) const {
		return v->Complement();
	}

	EdgeId Complement(EdgeId e) const;

	const EdgeData& GetData(EdgeId e) const {
		return *e;
	}

	bool CanCompressVertex(VertexId v) const;

	void CompressVertex(VertexId v);

	EdgeId CompressPath(const vector<VertexId>& path);

	void CompressAllVertices();

	pair<EdgeId, EdgeId> SplitEdge(EdgeId edge, size_t position);

	void GlueEdges(EdgeId edge1, EdgeId edge2);

};

typedef EdgeGraph::EdgeId EdgeId;
typedef EdgeGraph::EdgeData EdgeData;
typedef EdgeGraph::VertexId VertexId;
typedef EdgeGraph::VertexIterator VertexIterator;
typedef EdgeGraph::EdgeIterator EdgeIterator;
typedef EdgeGraph::ActionHandler ActionHandler;
//typedef EdgeGraph::SmartVertexIterator SmartVertexIterator;
//typedef EdgeGraph::SmartEdgeIterator SmartEdgeIterator;

typedef de_bruijn::TraversalHandler<EdgeGraph> TraversalHandler;

}
#endif /* EDGE_GRAPH_HPP_ */

