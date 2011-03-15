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
#include <iostream>

using namespace std;

namespace condensed_graph {
#define HASH_SEED 128
typedef Seq<K> Kmer;
typedef Seq<K - 1> KMinusOneMer;
typedef Seq<N> Read;

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
	void AddDesc(Vertex* v) {
		int k_th_nucl = v->nucls()[K - 1];
		desc_[k_th_nucl] = v;
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

class SimpleHashTable { // To Sergey: it's C++, not Java Collections ;)
private:
	//typedef __gnu_cxx::hash_map<const Kmer, pair<Vertex*, size_t> , Kmer::hash<HASH_SEED>, Kmer::equal_to> hmap;
	typedef tr1::unordered_map<const Kmer, pair<Vertex*, size_t> , Kmer::hash<HASH_SEED>, Kmer::equal_to> hmap;
	hmap h_;
public:
	void put(Kmer k, Vertex* v, size_t s) {
		//DEBUG("Putting position for k-mer '" << k.str() <<  "' : position " << v.second)
		hmap::iterator hi = h_.find(k);
		if (hi == h_.end()) { // put new element
			h_[k] = make_pair(v, s);
		} else { // change existing element
			hi->second = make_pair(v, s);
		}
		/*if (contains(k)) {
		 h_.erase(k);
		 }
		 h_.insert(make_pair(k, v));*/
	}

	bool contains(Kmer k) {
		return h_.find(k) != h_.end();
	}

	const pair<Vertex*, size_t> get(const Kmer &k) {
		hmap::iterator hi = h_.find(k);
		assert(hi != h_.end()); // contains
		//DEBUG("Getting position of k-mer '" + k.str() + "' Position is " << hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}
};

class Graph {
	set<Vertex*> vertices_;
	SimpleHashTable h_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(Vertex* v) {
		assert(v->nucls().size() >= K);
		Kmer k(v->nucls());
		h_.put(k, v, 0);
		for (size_t i = K, n = v->nucls().size(); i < n; ++i) {
			k = k << v->nucls()[i];
			h_.put(k, v, i - K + 1);
		}
	}

	/**
	 * deals with incoming links and their complement only!!!
	 */
	void FixIncomingOnSplit(Vertex* v, Vertex* v1, Vertex* v2);

	void FixIncomingOnMerge(Vertex* v1, Vertex* v2, Vertex* v);

	bool CanBeDeleted(Vertex* v) const;

public:
	const set<Vertex*>& vertices() const;

	vector<Vertex*> Anc(const Vertex* v) const;

	vector<Vertex*> Desc(const Vertex* v) const;

	bool IsLast(Vertex* v) const {
		return v->IsDeadend();
	}

	bool IsFirst(Vertex* v) const {
		return IsLast(v->complement());
	}

	/**
	 * adds vertex and its complement
	 */
	Vertex* AddVertex(const Sequence &nucls);

	bool IsMergePossible(Vertex* v1, Vertex* v2) const;

	Vertex* Merge(Vertex* v1, Vertex* v2);

	/**
	 * deletes vertex and its complement
	 */
	void DeleteVertex(Vertex* v);

	bool AreLinkable(Vertex* v1, Vertex* v2) const;

	void LinkVertices(Vertex* anc, Vertex* desc);

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewHashForVertexKmers(Vertex* v);

	//pos exclusive! (goes into second vertex)
	//deletes vertex if actual split happens
	//returns first of the new vertices
	Vertex* SplitVertex(Vertex* v, size_t pos);

	pair<Vertex*, int> GetPosition(Kmer k) {
		assert(h_.contains(k));
		return h_.get(k);
	}

	bool Contains(Kmer k) {
		return h_.contains(k);
	}

	pair<Vertex*, int> GetPosMaybeMissing(Kmer k) {
		if (!h_.contains(k)) {
			AddVertex(Sequence(k));
		}
		return h_.get(k);
	}

	void ThreadRead(const Read &r);

	bool IsLastKmer(Vertex* v, size_t pos) const {
		return pos + K == v->size();
	}

	bool IsFirstKmer(Vertex* v, size_t pos) const {
		return pos == 0;
	}
};

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

	Traversal(const Graph& g);
	virtual void Traverse(Handler& h) {
	}
protected:
	const Graph& g_;
};

class DFS: public Traversal {
	set<Vertex*> visited_;
	void go(Vertex* v, vector<Vertex*>& stack, Handler& h);
public:
	DFS(const Graph& g);
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

void CondenseGraph(DeBruijn<K>& origin, Graph& g);

}
#endif /* CONDENSED_GRAPH_H_ */
