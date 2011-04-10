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
using de_bruijn::SmartEdgeIterator;
using de_bruijn::SmartVertexIterator;

namespace edge_graph {

using de_bruijn::GraphActionHandler;

LOGGER("d.edge_graph");

class Vertex;

class Edge {
	Sequence nucls_;
	Vertex* end_;
	size_t coverage_;
	size_t incoming_coverage_;
	size_t outgoing_coverage_;

	friend class EdgeGraph;
	Edge(const Sequence& nucls, Vertex* end) :
		nucls_(nucls), end_(end), coverage_(0) {
	}

	Vertex* end() const {
		return end_;
	}

	size_t size() const {
		return nucls_.size();
	}

public:
	const Sequence& nucls() const {
		return nucls_;
	}

};

class Vertex {
public:
	typedef vector<Edge *>::const_iterator EdgeIterator;
private:
	vector<Edge *> outgoing_edges_;

	Vertex* complement_;

	friend class EdgeGraph;

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

	//
	//	void RemoveOutgoingEdge(vector<Edge *>::iterator iter) {
	//		assert(iter != outgoing_edges_.end());
	//		outgoing_edges_.erase(iter);
	//	}

	Vertex* complement() const {
		return complement_;
	}

	~Vertex() {
		assert(outgoing_edges_.size() == 0);
	}

};
/*
class GraphActionHandler {
public:

	virtual void HandleAdd(Vertex* v) {
	}

	virtual void HandleAdd(Edge* e) {
	}

	virtual void HandleDelete(Vertex* v) {
	}

	virtual void HandleDelete(Edge* e) {
	}

};*/

class EdgeGraph {

	size_t k_;

	Edge* AddSingleEdge(Vertex* v1, Vertex* v2, const Sequence& s);

	//	void DeleteSingleEdge(const Edge* edge);

	vector<GraphActionHandler<EdgeGraph> *> action_handler_list_;

	set<Vertex*> vertices_;

	void DeleteAllOutgoing(Vertex *v);

//	const set<Vertex*>& vertices() const {
//		return vertices_;
//	}

public:

	typedef Edge* EdgeId;
	typedef Vertex* VertexId;

	typedef set<Vertex *>::const_iterator VertexIterator;

	VertexIterator begin() {
		return vertices_.begin();
	}

	VertexIterator end() {
		return vertices_.end();
	}

	SmartVertexIterator<EdgeGraph> SmartVertexBegin() {
		return de_bruijn::SmartVertexIterator<EdgeGraph>(*this);
	}

	SmartVertexIterator<EdgeGraph> SmartVertexEnd() const {
		return de_bruijn::SmartVertexIterator<EdgeGraph>();
	}

	SmartEdgeIterator<EdgeGraph> SmartEdgeBegin() {
		return de_bruijn::SmartEdgeIterator<EdgeGraph>(*this);
	}

	SmartEdgeIterator<EdgeGraph> SmartEdgeEnd() const {
		return de_bruijn::SmartEdgeIterator<EdgeGraph>();
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

	void AddActionHandler(GraphActionHandler<EdgeGraph>* action_handler) {
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(GraphActionHandler<EdgeGraph>* action_handler) {
		for (vector<GraphActionHandler<EdgeGraph> *>::iterator it =
				action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
			if(*it == action_handler) {
				action_handler_list_.erase(it);
				return true;
			}
		}
		return false;
	}

	const vector<GraphActionHandler<EdgeGraph> *> GetHandlers() {
		return action_handler_list_;
	}

	void OutgoingEdges(const Vertex* v, Vertex::EdgeIterator &begin,
			Vertex::EdgeIterator &end) const;

	const vector<Edge *> OutgoingEdges(const Vertex* v) const;

	const vector<Edge *> IncomingEdges(const Vertex* v) const;

	Edge* OutgoingEdge(const Vertex* v, char nucl) const;

	size_t OutgoingEdgeCount(Vertex *v) const {
		return v->OutgoingEdgeCount();
	}

	size_t IncomingEdgeCount(Vertex *v) const {
		return v->complement()->OutgoingEdgeCount();
	}

	bool CheckUniqueOutgiongEdge(const Vertex *v) const {
		return v->OutgoingEdgeCount() == 1;
	}

	Edge *GetUniqueOutgoingEdge(const Vertex *v) const {
		assert(CheckUniqueOutgiongEdge(v));
		return *(v->begin());
	}

	bool CheckUniqueIncomingEdge(const Vertex *v) const {
		return CheckUniqueOutgiongEdge(v->complement());
	}

	Edge *GetUniqueIncomingEdge(const Vertex *v) const {
		return ComplementEdge(GetUniqueOutgoingEdge(v->complement()));
	}

	Edge *ComplementEdge(const Edge* edge) const;

	const Sequence &EdgeNucls(const Edge *edge) const {
		return edge->nucls();
	}

	void set_coverage(Edge *edge, size_t cov) {
		edge->coverage_ = cov;
	}

	size_t coverage(Edge *edge) {
		return edge->coverage_;
	}

	void inc_coverage(Edge *edge, int toAdd) {
		edge->coverage_ += toAdd;
	}

	void inc_coverage(Edge *edge) {
		edge->coverage_++;
	}

	/**
	 * adds vertex and its complement
	 */
	Vertex* AddVertex();

	Sequence vertexNucls(const Vertex *v) const;

	/**
	 * deletes vertex and its complement
	 */
	void DeleteVertex(Vertex* v);

	void ForceDeleteVertex(Vertex* v);

	Edge* AddEdge(Vertex* v1, Vertex* v2, const Sequence &nucls);

	void DeleteEdge(Edge* edge);

	size_t length(Edge *edge) {
		return edge->nucls_.size() - k_;
	}

	bool AreLinkable(Vertex* v1, Vertex* v2, const Sequence &nucls) const;

	bool IsDeadEnd(Vertex* v) const {
		return v->IsDeadend();
	}

	bool IsDeadStart(Vertex* v) const {
		return IsDeadEnd(v->complement());
	}

	Vertex *EdgeStart(const Edge *edge) const;

	Vertex *EdgeEnd(const Edge *edge) const;

	Vertex *ComplementVertex(const Vertex* v) const {
		return v->complement();
	}

	Vertex *Complement(const Vertex* v) const {
		return ComplementVertex(v);
	}

	Edge *Complement(const Edge* e) const {
		return ComplementEdge(e);
	}

	const Edge& GetData(Edge* e) const {
		return *e;
	}

	bool CanCompressVertex(Vertex *v) const;

	Edge *CompressVertex(Vertex *v);

	Edge *CompressPath(const vector<Vertex *> path);

	void CompressAllVertices();
};

//////////////////////////////////////////////////////////////////

class GraphVisualizer {
public:
	virtual void Visualize(const EdgeGraph& g) = 0;
};

class SimpleGraphVisualizer: public GraphVisualizer {
	gvis::GraphPrinter<const Vertex*>& gp_;
public:
	SimpleGraphVisualizer(gvis::GraphPrinter<const Vertex*>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const EdgeGraph& g);
};

class ComplementGraphVisualizer: public GraphVisualizer {
	gvis::PairedGraphPrinter<const Vertex*>& gp_;
public:
	ComplementGraphVisualizer(gvis::PairedGraphPrinter<const Vertex*>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const EdgeGraph& g);
};

}
#endif /* EDGE_GRAPH_HPP_ */

