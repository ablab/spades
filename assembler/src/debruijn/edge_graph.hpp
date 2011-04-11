/*
 ReadOnlyGraph -> Graph -> EdgeGraph/VertexGraph
 ComplementGraph
 PairedGraph/Debruijn



 interface ReadonlyGraph {
 NodeIt begin_node();
 NodeIt end_node();

 void outgoing_edges(const NodeId &n, EdgeIt &begin, EdgeIt &end);

 bool IsLast(const NodeId &n);

 bool IsFirst(const NodeId &n);

 NodeId start(const EdgeId &e);
 NodeId end(const EdgeId &e);

 };

 interface Graph {

 };

 interface EdgeGraph {
 //	const NodeId addNode(const NodeData& n_d);

 const NodeId addNode();

 const EdgeId addEdge(const EdgeData& e_d, const NodeId& n1, const NodeId& n2);
 const EdgeData& getEdgeData(const EdgeId &e);
 }

 interface ComplementGraph {
 const NodeId complement_node(const& NodeId);
 const EdgeId complement_edge(const& EdgeId);
 }

 interface NodeGraph {
 const NodeData& getNodeData(const NodeId &n);
 }

 interface DeBruijnGraph {
 const EdgeData* getNodeData(const NodeId &n, char c);
 }

 */

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
public:
	const Sequence& nucls() const {
		return nucls_;
	}
private:
	friend class EdgeGraph;
	Sequence nucls_;
	Vertex* end_;
	size_t coverage_;

	Edge(const Sequence& nucls, Vertex* end) :
		nucls_(nucls), end_(end), coverage_(0) {
	}

	Vertex* end() const {
		return end_;
	}

	size_t size() const {
		return nucls_.size();
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

	void set_complement(Vertex* complement) {
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

	Vertex* complement() const {
		return complement_;
	}

	~Vertex() {
		assert(outgoing_edges_.size() == 0);
	}

};

class EdgeGraph {
public:
	typedef Edge* EdgeId;
	typedef Edge EdgeData;
	typedef Vertex* VertexId;

	typedef set<Vertex*>::const_iterator VertexIterator;
	typedef Vertex::EdgeIterator EdgeIterator;
	typedef de_bruijn::GraphActionHandler<EdgeGraph> ActionHandler;
	typedef de_bruijn::SmartVertexIterator<EdgeGraph> SmartVertexIterator;
	typedef de_bruijn::SmartEdgeIterator<EdgeGraph> SmartEdgeIterator;

	VertexIterator begin() const {
		return vertices_.begin();
	}

	VertexIterator end() const {
		return vertices_.end();
	}

	SmartVertexIterator SmartVertexBegin() {
		return SmartVertexIterator(*this);
	}

	SmartVertexIterator SmartVertexEnd() const {
		return SmartVertexIterator();
	}

	SmartEdgeIterator SmartEdgeBegin() {
		return SmartEdgeIterator(*this);
	}

	SmartEdgeIterator SmartEdgeEnd() const {
		return SmartEdgeIterator();
	}

	size_t size() {
		return vertices_.size();
	}

	/**
	 * Constructs empty graph to work with k-mers.
	 *
	 * @param k Main parameter that defines the size of k-mers
	 * @param action_handler Graph actions handler
	 */
	EdgeGraph(size_t k) {
		assert(k % 2 == 1);
		k_ = k;
	}

	/**
	 * Deletes action_handler.
	 */
	~EdgeGraph() {
		while (!vertices_.empty()) {
			ForceDeleteVertex(*vertices_.begin());
		}
	}

	size_t k() {
		return k_;
	}

	void AddActionHandler(ActionHandler* action_handler) {
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(ActionHandler* action_handler) {
		for (vector<ActionHandler*>::iterator it =
				action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
			if(*it == action_handler) {
				action_handler_list_.erase(it);
				return true;
			}
		}
		assert(false);
		return false;
	}

	//todo remove
	const vector<ActionHandler*> GetHandlers() {
		return action_handler_list_;
	}

	void OutgoingEdges(VertexId v, EdgeIterator& begin,
			EdgeIterator& end) const;

	const vector<EdgeId> OutgoingEdges(VertexId v) const;

	const vector<EdgeId> IncomingEdges(VertexId v) const;

	EdgeId OutgoingEdge(VertexId v, char nucl) const;

	size_t OutgoingEdgeCount(VertexId v) const {
		return v->OutgoingEdgeCount();
	}

	size_t IncomingEdgeCount(VertexId v) const {
		return v->complement()->OutgoingEdgeCount();
	}

	bool CheckUniqueOutgiongEdge(VertexId v) const {
		return v->OutgoingEdgeCount() == 1;
	}

	EdgeId GetUniqueOutgoingEdge(VertexId v) const {
		assert(CheckUniqueOutgiongEdge(v));
		return *(v->begin());
	}

	bool CheckUniqueIncomingEdge(const VertexId v) const {
		return CheckUniqueOutgiongEdge(v->complement());
	}

	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		return Complement(GetUniqueOutgoingEdge(v->complement()));
	}

//	Edge* ComplementEdge(const Edge* edge) const;

	const Sequence& EdgeNucls(EdgeId edge) const {
		return edge->nucls();
	}

	void set_coverage(EdgeId edge, size_t cov) {
		edge->coverage_ = cov;
	}

	size_t coverage(EdgeId edge) const {
		return edge->coverage_;
	}

	void inc_coverage(EdgeId edge, int toAdd) {
		edge->coverage_ += toAdd;
	}

	void inc_coverage(EdgeId edge) {
		edge->coverage_++;
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

	Edge* AddEdge(VertexId v1, VertexId v2, const Sequence &nucls);

	void DeleteEdge(EdgeId edge);

	size_t length(EdgeId edge) {
		return edge->nucls_.size() - k_;
	}

	bool AreLinkable(VertexId v1, VertexId v2, const Sequence &nucls) const;

	bool IsDeadEnd(VertexId v) const {
		return v->IsDeadend();
	}

	bool IsDeadStart(VertexId v) const {
		return IsDeadEnd(v->complement());
	}

	VertexId EdgeStart(EdgeId edge) const;

	VertexId EdgeEnd(EdgeId edge) const;

	VertexId Complement(VertexId v) const {
		return v->complement();
	}

	EdgeId Complement(EdgeId e) const;

	const EdgeData& GetData(EdgeId e) const {
		return *e;
	}

	bool CanCompressVertex(VertexId v) const;

	EdgeId CompressVertex(VertexId v);

	EdgeId CompressPath(const vector<VertexId>& path);

	void CompressAllVertices();

private:
	size_t k_;

	EdgeId AddSingleEdge(VertexId v1, VertexId v2, const Sequence& s);

	vector<ActionHandler*> action_handler_list_;

	set<Vertex*> vertices_;

	void DeleteAllOutgoing(Vertex* v);

};

typedef EdgeGraph::EdgeId EdgeId;
typedef EdgeGraph::EdgeData EdgeData;
typedef EdgeGraph::VertexId VertexId;
typedef EdgeGraph::VertexIterator VertexIterator;
typedef EdgeGraph::EdgeIterator EdgeIterator;
typedef EdgeGraph::ActionHandler ActionHandler;
typedef EdgeGraph::SmartVertexIterator SmartVertexIterator;
typedef EdgeGraph::SmartEdgeIterator SmartEdgeIterator;

typedef de_bruijn::TraversalHandler<EdgeGraph> TraversalHandler;
//////////////////////////////////////////////////////////////////

class VisHandler: public TraversalHandler {
	const EdgeGraph& g_;
	gvis::GraphPrinter<VertexId>& pr_;
public:

	VisHandler(const EdgeGraph& g, gvis::GraphPrinter<VertexId>& pr) :
		g_(g), pr_(pr) {
	}

	virtual void HandleVertex(VertexId v) {
		pr_.addVertex(v, "");
	}

	virtual void HandleEdge(EdgeId e) {
		stringstream ss;
		ss << e->nucls().size();

		pr_.addEdge(g_.EdgeStart(e), g_.EdgeStart(e), ss.str());
	}

};

class ComplementVisHandler: public TraversalHandler {
	const EdgeGraph& g_;
	gvis::PairedGraphPrinter<VertexId>& pr_;
public:

	ComplementVisHandler(const EdgeGraph& g, gvis::PairedGraphPrinter<VertexId>& pr) :
		g_(g), pr_(pr) {
	}

	virtual void HandleVertex(VertexId v) {
		pr_.addVertex(v, "", g_.Complement(v), "");
	}

	virtual void HandleEdge(EdgeId e) {
		stringstream ss;
		ss << e->nucls().size();
		VertexId v1 = g_.EdgeStart(e);
		VertexId v2 = g_.EdgeEnd(e);
		pr_.addEdge(make_pair(v1, g_.Complement(v1)), make_pair(v2, g_.Complement(v2)), ss.str());
	}

};


class GraphVisualizer {
public:
	virtual void Visualize(const EdgeGraph& g) = 0;
};

class SimpleGraphVisualizer: public GraphVisualizer {
	gvis::GraphPrinter<VertexId>& gp_;
public:
	SimpleGraphVisualizer(gvis::GraphPrinter<VertexId>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const EdgeGraph& g);
};

class ComplementGraphVisualizer: public GraphVisualizer {
	gvis::PairedGraphPrinter<VertexId>& gp_;
public:
	ComplementGraphVisualizer(gvis::PairedGraphPrinter<VertexId>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const EdgeGraph& g);
};

}
#endif /* EDGE_GRAPH_HPP_ */

