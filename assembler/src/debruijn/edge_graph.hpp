#ifndef EDGE_GRAPH_HPP_
#define EDGE_GRAPH_HPP_

#include <vector>
#include <set>
#include <cstring>
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "sequence.hpp"
#include "logging.hpp"
#include "nucl.hpp"
#include "strobe_read.hpp"
#include "utils.hpp"
#include "omni_utils.hpp"

using namespace std;

namespace edge_graph {
LOGGER("d.edge_graph");

using omnigraph::GraphActionHandler;
using omnigraph::HandlerApplier;
using omnigraph::PairedHandlerApplier;
using omnigraph::SmartVertexIterator;
using omnigraph::SmartEdgeIterator;


class Vertex;

class Edge {
private:
	friend class EdgeGraph;
	const Sequence& nucls() const {
		return nucls_;
	}
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

	const vector<Edge*> OutgoingEdges() const {
		return outgoing_edges_;
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

class EdgeGraph {
public:
	typedef Edge* EdgeId;
	typedef Edge EdgeData;
	typedef Vertex* VertexId;

	typedef set<Vertex*>::const_iterator VertexIterator;
	typedef Vertex::EdgeIterator EdgeIterator;
	typedef de_bruijn::GraphActionHandler<EdgeGraph> ActionHandler;

private:
	const size_t k_;

	const PairedHandlerApplier<EdgeGraph> applier_;

	vector<ActionHandler*> action_handler_list_;

	set<Vertex*> vertices_;


	VertexId HiddenAddVertex();

	EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const Sequence &nucls,
			size_t coverage = 0);

	EdgeId AddSingleEdge(VertexId v1, VertexId v2, const Sequence& s,
			size_t coverage);

	void DeleteAllOutgoing(Vertex* v);

	void FireAddVertex(VertexId v);
	void FireAddEdge(EdgeId edge);
	void FireDeleteVertex(VertexId v);
	void FireDeleteEdge(EdgeId edge);
	void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge);
	void FireGlue(EdgeId edge1, EdgeId edge2);
	void FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2);

public:

	/**
	 * @return const iterator pointing to the beginning of collection of vertices
	 */
	VertexIterator begin() const {
		return vertices_.begin();
	}

	/**
	 * @return const iterator pointing to the end of collection of vertices
	 */
	VertexIterator end() const {
		return vertices_.end();
	}

	/**
	 * Method returns smart iterator over vertices of graph pointing to the beginning of vertex collection
	 * @param comparator comparator which defines order in which vertices would be iterated.
	 */
	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<EdgeGraph, Comparator> SmartVertexBegin(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<EdgeGraph, Comparator> (*this, true,
				comparator);
	}

	/**
	 * Method returns smart iterator over vertices of graph pointing to the end of vertex collection
	 * @param comparator comparator which defines order in which vertices would be iterated.
	 */
	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<EdgeGraph, Comparator> SmartVertexEnd(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<EdgeGraph, Comparator> (*this, false,
				comparator);
	}

	/**
	 * Method returns smart iterator over vertices of graph pointing to the beginning of edge collection
	 * @param comparator comparator which defines order in which edges would be iterated.
	 */
	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<EdgeGraph, Comparator> SmartEdgeBegin(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<EdgeGraph, Comparator> (*this, true,
				comparator);
	}

	/**
	 * Method returns smart iterator over vertices of graph pointing to the end of edge collection
	 * @param comparator comparator which defines order in which edges would be iterated.
	 */
	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<EdgeGraph, Comparator> SmartEdgeEnd(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<EdgeGraph, Comparator> (*this, false,
				comparator);
	}

	/**
	 * @return number of vertices
	 */
	size_t size() {
		return vertices_.size();
	}

	/**
	 * Constructs empty graph to work with k-mers.
	 *
	 * @param k Main parameter that defines the size of k-mers
	 * //@param action_handler Graph actions handler
	 */
	EdgeGraph(size_t k) : k_(k) , applier_(*this) {
		assert(k % 2 == 1);
	}

	/**
	 * Deletes action_handler.
	 */
	~EdgeGraph() {
		while (!vertices_.empty()) {
			ForceDeleteVertex(*vertices_.begin());
		}
	}

	/**
	 *
	 * @return value of k, which is number of nucleotides stored in vertices
	 */
	size_t k() {
		return k_;
	}

	/**
	 * Method adds new action handler to graph
	 * @param action handler to add
	 */
	void AddActionHandler(ActionHandler* action_handler) {
		TRACE("Action handler added");
		action_handler_list_.push_back(action_handler);
	}

	/**
	 * Method removes action handler from graph
	 * @param action handler to delete
	 * @return true if given action handler was among graph handlers and false otherwise
	 */
	bool RemoveActionHandler(ActionHandler* action_handler) {
		TRACE("Trying to remove action handler");
		for (vector<ActionHandler*>::iterator it =
				action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				action_handler_list_.erase(it);
				TRACE("Action handler removed");
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

//	void OutgoingEdges(VertexId v, EdgeIterator& begin, EdgeIterator& end) const;

	/**
	 * Method returnes vector of all outgoing edges of given vertex
	 * @param v vertex to get outgoing edges from
	 */
	const vector<EdgeId> OutgoingEdges(VertexId v) const;

	/**
	 * Method returnes vector of all incoming edges of given vertex
	 * @param v vertex to get incoming edges from
	 */
	const vector<EdgeId> IncomingEdges(VertexId v) const;

	/**
	 * @depricated
	 */
	const vector<EdgeId> IncidentEdges(VertexId v) const;

	/**
	 * @depricated
	 */
	const vector<EdgeId> NeighbouringEdges(EdgeId e) const;

	/**
	 * Method returns outgoing edge with given nucleotide at k-th position in sequence of
	 * edge returned, which is the first position outgoing edges have different nucleotides in.
	 * @param v vertex to find outgoing edge for
	 * @param nucl nucleotide to be found at k-th position of edge
	 */
	EdgeId OutgoingEdge(VertexId v, char nucl) const;

	/**
	 * Method returns the number of outgoing edges.
	 * @param v vertex to count outgoing edges for
	 */
	size_t OutgoingEdgeCount(VertexId v) const {
		return v->OutgoingEdgeCount();
	}

	/**
	 * Method returns the number of incoming edges.
	 * @param v vertex to count incoming edges for
	 */
	size_t IncomingEdgeCount(VertexId v) const {
		return v->Complement()->OutgoingEdgeCount();
	}

	/**
	 * Method returns true if vertex has only one outgoing edge and false otherwise.
	 * @param v vertex to check
	 */
	bool CheckUniqueOutgiongEdge(VertexId v) const {
		return v->OutgoingEdgeCount() == 1;
	}

	/**
	 * Method returns unique outgoing edge. Asserts if outgoing edge is not unique
	 * @param v vertex to find unique outgoing edge for
	 */
	EdgeId GetUniqueOutgoingEdge(VertexId v) const {
		assert(CheckUniqueOutgiongEdge(v));
		return (v->OutgoingEdges())[0];
	}

	/**
	 * Method returns true if vertex has only one incoming edge and false otherwise.
	 * @param v vertex to check
	 */
	bool CheckUniqueIncomingEdge(VertexId v) const {
		return CheckUniqueOutgiongEdge(v->Complement());
	}

	/**
	 * Method returns unique incoming edge. Asserts if incoming edge is not unique
	 * @param v vertex to find unique incoming edge for
	 */
	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		return Complement(GetUniqueOutgoingEdge(v->Complement()));
	}

	//	Edge* ComplementEdge(const Edge* edge) const;

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return edge->nucls();
	}

	/**
	 * Method sets coverage value for the edge
	 */
	void SetCoverage(EdgeId edge, size_t cov) {
		edge->coverage_ = cov;
	}

	/**
	 * Method returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return (double) edge->coverage_ / length(edge);
	}

	/**
	 * Method increases coverage value
	 */
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

	void Merge(EdgeId edge1, EdgeId edg2);

	EdgeId MergePath(const vector<EdgeId>& path);

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

