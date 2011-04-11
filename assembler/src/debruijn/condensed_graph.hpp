/**
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

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

using namespace std;

namespace condensed_graph {
LOGGER("d.condensed_graph");

/**
 * @brief Vertex of a condensed graph.
 *
 * Class representing vertex of condensed graph.
 * Contains the Sequence of nucleotides, average coverage of this sequence,
 * pointers to children and coverage of corresponding edges.
 *
 */
class Vertex {
private:
	Sequence nucls_;
	Vertex* right_neighbours_[4];
	Vertex* complement_;

//	int coverage_;
	size_t edge_coverage_[4];

	friend class CondensedGraph;
public:
	//bool deleted;
	Vertex(const Sequence &nucls) :
		nucls_(nucls) {
		fill_n(right_neighbours_, 4, (Vertex*) NULL);
		fill_n(edge_coverage_, 4, 0);

	}

	int RightNeighbourCount() {
		int c = 0;
		for (int i = 0; i < 4; ++i)
			if (right_neighbours_[i] != NULL)
				c++;
		return c;
	}

	bool IsDeadend() {
		return RightNeighbourCount() == 0;
	}

	//	Vertex* const * right_neighbours() const {
	//		return right_neighbours_;
	//	}

	Vertex* right_neighbour(char nucl) const {
		return right_neighbours_[(int) nucl];
	}

	size_t coverage(char nucl) {
		return edge_coverage_[(int) nucl];
	}

	size_t size() const {
		return nucls_.size();
	}

	const Sequence& nucls() const {
		return nucls_;
	}

	void set_right_neigbour(Vertex* v, char nucl) {
		right_neighbours_[(int) nucl] = v;
	}

	void set_coverage(size_t coverage, char nucl) {
		edge_coverage_[(int) nucl] = coverage;
	}

	void inc_coverage(char nucl) {
		++edge_coverage_[(int) nucl];
	}

	Vertex* complement() const {
		return complement_;
	}
	void set_complement(Vertex* complement) {
		complement_ = complement;
	}
//	int coverage() {
//		return coverage_;
//	}
//	void set_coverage(int coverage) {
//		coverage_ = coverage;
//	}
};

/**
 * @brief Condensed graph action listener.
 *
 * Listener of condensed graph public actions.
 *
 */
class GraphActionHandler {
public:

	/**
	 * Handle addition of vertex and complement.
	 *
	 * @see CondensedGraph::AddVertex()
	 * @param v Pointer to vertex that has been added
	 */
	virtual void HandleAdd(Vertex* v) {
	}

	/**
	 * Handle deletion of vertex and complement.
	 *
	 * @see CondensedGraph::DeleteVertex()
	 * @param v Pointer to vertex that has been deleted
	 */
	virtual void HandleDelete(Vertex* v) {
	}

	/**
	 * Handle merge of two vertices and complement.
	 *
	 * @see CondensedGraph::Merge()
	 * @param v1 Pointer to left vertex that was merged
	 * @param v2 Pointer to right vertex that was merged
	 * @param v Pointer to vertex that is merge result
	 */
	virtual void HandleMerge(Vertex* v1, Vertex* v2, Vertex* v) {
	}

	/**
	 * Handle split of vertex and its complement into two.
	 *
	 * @see CondensedGraph::SplitVertex()
	 * @param v Pointer to vertex that was split
	 * @param pos Position by which v was split (this position corresponds to the k-1 position of v2)
	 * @param v1 Pointer to left vertex that is split result
	 * @param v2 Pointer to right vertex that is split result
	 */
	virtual void HandleSplit(Vertex* v, size_t pos, Vertex* v1, Vertex* v2) {
	}
};

/**
 * @brief Condensed DeBruijn Graph with empty edges and sequences in vertices
 *
 * Condensed DeBruijn Graph with empty edges and sequences in vertices.
 *
 * When some graph action happens (Add, Delete, Split, Merge),
 * corresponding method of handler with proper arguments is called.
 *
 * Passing pointer on action handler to graph, forget about allocated memory,\
 * cause graph will delete it when no longer needed.
 */
class CondensedGraph {

	/**
	 * Fixes incoming edges during split.
	 *
	 * After split of vertex v edges that was pointing to it should point to v1.
	 *
	 * @param v Vertex that was split.
	 * @param v1 Left vertex that was obtained during split. Incoming target for new edges.
	 * @param v2 Right vertex that was obtained during split. Need it to handle one of special cases.
	 */
	void FixIncomingOnSplit(Vertex* v, Vertex* v1, Vertex* v2);

	/**
	 * Fixes incoming edges during merge.
	 *
	 * After merge of vertices v1 and v2 edges that was pointing to v1 should point to merge result v.
	 *
	 * @param v1 Left vertex that was merged.
	 * @param v2 Right vertex that was merged. Need it to handle one of special cases.
	 * @param v Merge result. Incoming target for new edges.
	 */
	void FixIncomingOnMerge(Vertex* v1, Vertex* v2, Vertex* v);

	bool CheckIfNoIncoming(Vertex* v) const;

	/**
	 * Checks that neighbours (its right neighbours and right neighbours of complement)
	 * of the v don't have v as neighbour.
	 *
	 * @param v Vertex to check condition for.
	 * @return Check result.
	 */
	bool CanBeDeleted(Vertex* v) const;

	/**
	 * Adds v2 as right neighbour to v1.
	 */
	void AddRightNeighbour(Vertex* v1, Vertex* v2);

	size_t k_;
	GraphActionHandler* action_handler_;
	set<Vertex*> vertices_;

public:

	/**
	 * Constructs empty graph to work with k-mers.
	 *
	 * @param k Main parameter that defines the size of k-mers
	 * @param action_handler Graph actions handler
	 */
	CondensedGraph(size_t k, GraphActionHandler* action_handler) :
		k_(k), action_handler_(action_handler) {
	}

	/**
	 * Deletes action_handler.
	 */
	~CondensedGraph() {
		delete action_handler_;
	}

	const set<Vertex*>& vertices() const {
		return vertices_;
	}

	size_t k() {
		return k_;
	}

	void set_action_handler(GraphActionHandler* action_handler) {
		delete action_handler_;

		action_handler_ = action_handler;
	}

	vector<Vertex*> LeftNeighbours(const Vertex* v) const;

	vector<Vertex*> RightNeighbours(const Vertex* v) const;

	//todo make private
	/**
	 * adds vertex and its complement
	 */
	Vertex* AddVertex(const Sequence &nucls);
	//	Vertex* AddVertex(Sequence nucls);

	/**
	 * deletes vertex and its complement
	 */
	void DeleteVertex(Vertex* v);

	//pos exclusive! (goes into second vertex)
	//deletes vertex if actual split happens
	//returns first of the new vertices
	Vertex* SplitVertex(Vertex* v, size_t pos);

	Vertex* Merge(Vertex* v1, Vertex* v2);

	bool IsMergePossible(Vertex* v1, Vertex* v2) const {
		return IsLast(v1) && IsFirst(v2) && v1->complement() != v2 && v1 != v2
				&& AreLinkable(v1, v2);
	}

	void LinkVertices(Vertex* v1, Vertex* v2);

	void UnLinkVertices(Vertex* v1, Vertex* v2);

	void UnLinkAll(Vertex* v);

	bool AreLinkable(Vertex* v1, Vertex* v2) const {
		return v2->nucls().Subseq(0, k_ - 1) == v1->nucls().Subseq(
				v1->size() - (k_ - 1));
	}

	bool IsLast(Vertex* v) const {
		return v->IsDeadend();
	}

	bool IsFirst(Vertex* v) const {
		return IsLast(v->complement());
	}

	bool IsLastKmer(Vertex* v, size_t pos) const {
		return pos + k_ == v->size();
	}

	bool IsFirstKmer(Vertex* v, size_t pos) const {
		return pos == 0;
	}
};

//////////////////////////////////////////////////////////////////

/**
 * @brief Base class for condensed graph traversals.
 */
class Traversal {
public:

	/**
	 * Stub base class for handling graph primitives during traversal.
	 */
	class Handler {
	public:
		virtual void HandleStartVertex(const Vertex* v) {
		}
		virtual void HandleEndVertex(const Vertex* v) {
		}
		virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		}
	};

	Traversal(const CondensedGraph& g) :
		g_(g) {
	}

	/**
	 *
	 */
	virtual void Traverse(Handler& h) =0;

protected:
	const CondensedGraph& g_;
};

class DFS: public Traversal {
	set<Vertex*> visited_;
	void ProcessVertex(Vertex* v, vector<Vertex*>& stack, Handler& h);
public:
	DFS(const CondensedGraph& g) :
		Traversal(g) {

	}
	virtual void Traverse(Handler& h);
};

class GraphVisualizer {
public:
	virtual void Visualize(const CondensedGraph& g) = 0;
};

class SimpleGraphVisualizer: public GraphVisualizer {
	gvis::GraphPrinter<const Vertex*>& gp_;
public:
	SimpleGraphVisualizer(gvis::GraphPrinter<const Vertex*>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const CondensedGraph& g);
};

class ComplementGraphVisualizer: public GraphVisualizer {
	gvis::PairedGraphPrinter<const Vertex*>& gp_;
public:
	ComplementGraphVisualizer(gvis::PairedGraphPrinter<const Vertex*>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const CondensedGraph& g);
};

class SimpleStatCounter: public Traversal::Handler {
	size_t v_count_;
	size_t e_count_;
public:
	SimpleStatCounter() :
		v_count_(0), e_count_(0) {
	}
	virtual void HandleStartVertex(const Vertex* v) {
		v_count_++;
	}
	virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		e_count_++;
	}

	size_t v_count() const {
		return v_count_;
	}

	size_t e_count() const {
		return e_count_;
	}
};

//class CountHandler: public Traversal::Handler {
//	tr1::unordered_map<const Vertex*, size_t>& map_;
//	size_t count_;
//public:
//
//	CountHandler(tr1::unordered_map<const Vertex*, size_t>& map) :
//		map_(map), count_(0) {
//	}
//
//	virtual void HandleStartVertex(const Vertex* v) {
//		map_.insert(make_pair(v, count_++));
//	}
//};
//
class VisHandler: public Traversal::Handler {
	gvis::GraphPrinter<const Vertex*>& pr_;
public:

	VisHandler(gvis::GraphPrinter<const Vertex*>& pr) :
		pr_(pr) {
	}

	virtual void HandleStartVertex(const Vertex* v) {
		stringstream ss;
		ss << v->nucls().size();
		pr_.addVertex(v, ss.str());
	}

	virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		pr_.addEdge(v1, v2, "");
	}

};

class ComplementVisHandler: public Traversal::Handler {
	gvis::PairedGraphPrinter<const Vertex*>& pr_;
public:

	ComplementVisHandler(gvis::PairedGraphPrinter<const Vertex*>& pr) :
		pr_(pr) {
	}

	virtual void HandleStartVertex(const Vertex* v) {
		stringstream ss;
		ss << v->nucls().size();

		//todo delete after debug
		stringstream ss2;
		ss2 << v->complement()->nucls().size();
		pr_.addVertex(v, ss.str(), v->complement(), ss2.str());
	}

	virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		pr_.addEdge(make_pair(v1, v1->complement()),
				make_pair(v2, v2->complement()), "");
	}

};
}

#endif /* CONDENSED_GRAPH_H_ */
