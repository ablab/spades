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
//typedef Seq<K> Kmer;
//typedef Seq<K - 1> KMinusOneMer;
//typedef Seq<N> Read;
LOGGER("d.condensed_graph");

class Vertex;

class Vertex {
private:
	Sequence nucls_;
	Vertex* desc_[4];
	Vertex* complement_;

	int coverage_;
	int arc_coverage_[4];

public:
	//bool deleted;
	Vertex(const Sequence &nucls) :
		nucls_(nucls) {
		fill_n(desc_, 4, (Vertex*) NULL);
		fill_n(arc_coverage_, 4, 0);

	}

	Vertex(const Sequence &nucls, Vertex** desc) :
		nucls_(nucls) {
		memcpy(desc, desc_, 4 * sizeof(Vertex*));
		fill_n(arc_coverage_, 4, 0);
	}
	~Vertex() {
	}
	int DescCount() {
		int c = 0;
		for (int i = 0; i < 4; ++i)
			if (desc_[i] != NULL)
				c++;
		return c;
	}
	bool IsDeadend() {
		return DescCount() == 0;
	}
	Vertex* const * desc() const {
		return desc_;
	}
	Vertex* desc(char nucl) {
		return desc_[(int) nucl];
	}
	size_t size() {
		return nucls_.size();
	}
	const Sequence& nucls() const {
		return nucls_;
	}
	void AddDesc(Vertex* v, char nucl) {
		desc_[(int) nucl] = v;
	}
	Vertex* complement() const {
		return complement_;
	}
	void set_complement(Vertex* complement) {
		complement_ = complement;
	}
	int coverage() {
		return coverage_;
	}
	void set_coverage(int coverage) {
		coverage_ = coverage;
	}
};

//template<size_t kmer_size_> class GraphConstructor;

class ActionHandler {
public:
	virtual void HandleAdd(Vertex* v) {
	}
	virtual void HandleDelete(Vertex* v) {
	}
	virtual void HandleMerge(Vertex* v1, Vertex* v2, Vertex* v) {
	}
	virtual void HandleSplit(Vertex* v, size_t pos, Vertex* v1, Vertex* v2) {
	}
};

//////////////////////////////////////////////////

class Graph {
	set<Vertex*> vertices_;

	/**
	 * deals with incoming links and their complement only!!!
	 */
	void FixIncomingOnSplit(Vertex* v, Vertex* v1, Vertex* v2);

	void FixIncomingOnMerge(Vertex* v1, Vertex* v2, Vertex* v);

	bool CanBeDeleted(Vertex* v) const;

	void AddDesc(Vertex* anc, Vertex* desc);

	size_t k_;
	ActionHandler* action_handler_;
	//	template<size_t kmer_size_> friend class GraphConstructor;
public:

	Graph(size_t k, ActionHandler* action_handler) :
		k_(k), action_handler_(action_handler) {
	}

	~Graph() {
		delete action_handler_;
	}

	const set<Vertex*>& vertices() const {
		return vertices_;
	}

	void set_action_handler(ActionHandler* action_handler) {
		delete action_handler_;

		action_handler_ = action_handler;
	}

	vector<Vertex*> Anc(const Vertex* v) const;

	vector<Vertex*> Desc(const Vertex* v) const;

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

	void LinkVertices(Vertex* anc, Vertex* desc);

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

class Traversal {
public:
	class Handler {
	public:
		virtual void HandleStartVertex(const Vertex* v) {
		}
		virtual void HandleEndVertex(const Vertex* v) {
		}
		virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		}
	};

	Traversal(const Graph* g) :
		g_(g) {
	}
	virtual void Traverse(Handler& h) {
	}
protected:
	const Graph* g_;
};

class DFS: public Traversal {
	set<Vertex*> visited_;
	void go(Vertex* v, vector<Vertex*>& stack, Handler& h);
public:
	DFS(const Graph* g) :
		Traversal(g) {

	}
	virtual void Traverse(Handler& h);
};

class GraphVisualizer {
public:
	virtual void Visualize(const Graph& g) = 0;
};

class SimpleGraphVisualizer: public GraphVisualizer {
	gvis::GraphPrinter<const Vertex*>& gp_;
public:
	SimpleGraphVisualizer(gvis::GraphPrinter<const Vertex*>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const Graph& g);
};

class ComplementGraphVisualizer: public GraphVisualizer {
	gvis::PairedGraphPrinter<const Vertex*>& gp_;
public:
	ComplementGraphVisualizer(gvis::PairedGraphPrinter<const Vertex*>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const Graph& g);
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

class CountHandler: public Traversal::Handler {
	tr1::unordered_map<const Vertex*, size_t>& map_;
	size_t count_;
public:

	CountHandler(tr1::unordered_map<const Vertex*, size_t>& map) :
		map_(map), count_(0) {
	}

	virtual void HandleStartVertex(const Vertex* v) {
		map_.insert(make_pair(v, count_++));
	}
};

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
		pr_.addEdge(make_pair(v1, v1->complement()), make_pair(v2, v2->complement()), "");
	}

};
}

#endif /* CONDENSED_GRAPH_H_ */
