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
#include <iostream>

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

template<size_t kmer_size_>
class SimpleHashTable {
private:
	typedef Seq<kmer_size_> Kmer;
	typedef tr1::unordered_map<Kmer, pair<Vertex*, size_t> ,
			typename Kmer::hash, typename Kmer::equal_to> hmap;
	typedef typename hmap::iterator map_iter;
	//	typedef __gnu_cxx::hash_map<const Kmer, pair<Vertex*, size_t> , myhash, Kmer::equal_to> hmap;
	hmap h_;
public:
	void put(Kmer k, Vertex* v, size_t s) {
		//DEBUG("Putting position for k-mer '" << k.str() <<  "' : position " << v.second)
		map_iter hi = h_.find(k);
		if (hi == h_.end()) { // put new element
			h_[k] = make_pair(v, s);
		} else { // change existing element
			hi->second = make_pair(v, s);
		}
	}

	bool contains(Kmer k) {
		return h_.find(k) != h_.end();
	}

	const pair<Vertex*, size_t> get(const Kmer &k) {
		map_iter hi = h_.find(k);
		assert(hi != h_.end()); // contains
		//DEBUG("Getting position of k-mer '" + k.str() + "' Position is " << hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}

	bool deleteIfVertex(const Kmer &k, Vertex* v) {
		map_iter hi = h_.find(k);
		if (hi != h_.end() && (*hi).second.first == v) {
			h_.erase(k);
			return true;
		}
		return false;
	}
};

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

template<size_t kmer_size_>
class HashRenewer: public ActionHandler {
	typedef Seq<kmer_size_> Kmer;

	SimpleHashTable<kmer_size_>& h_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(Vertex* v) {
		assert(v->nucls().size() >= kmer_size_);
		Kmer k(v->nucls());
		h_.put(k, v, 0);
		for (size_t i = kmer_size_, n = v->nucls().size(); i < n; ++i) {
			k = k << v->nucls()[i];
			h_.put(k, v, i - kmer_size_ + 1);
		}
	}

	void DeleteKmersHash(Vertex* v) {
		assert(v->nucls().size() >= kmer_size_);
		Kmer k(v->nucls());
		h_.deleteIfVertex(k, v);
		for (size_t i = kmer_size_, n = v->nucls().size(); i < n; ++i) {
			k = k << v->nucls()[i];
			h_.deleteIfVertex(k, v);
		}
	}

public:
	HashRenewer(SimpleHashTable<kmer_size_>& h) :
		h_(h) {
	}
	;

	virtual void HandleAdd(Vertex* v) {
		RenewKmersHash(v);
		RenewKmersHash(v->complement());
	}

	virtual void HandleDelete(Vertex* v) {
		DeleteKmersHash(v);
		DeleteKmersHash(v->complement());
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
	ActionHandler& action_handler_;
	//	template<size_t kmer_size_> friend class GraphConstructor;
public:

	Graph(size_t k, ActionHandler& action_handler) :
		k_(k), action_handler_(action_handler) {
	}

	const set<Vertex*>& vertices() const {
		return vertices_;
	}

	void set_action_handler(ActionHandler& action_handler) {
		action_handler_ = action_handler;
	}

	vector<Vertex*> Anc(const Vertex* v) const;

	vector<Vertex*> Desc(const Vertex* v) const;

	//todo make private
	/**
	 * adds vertex and its complement
	 */
	Vertex* AddVertex(const Sequence &nucls);

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
		return v2->nucls().Subseq(0, k_) == v1->nucls().Subseq(
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

template<size_t kmer_size_>
class GraphConstructor {
protected:
	typedef Seq<kmer_size_> Kmer;
	Graph *g_;
	SimpleHashTable<kmer_size_> *h_;
	HashRenewer<kmer_size_> action_handler_;

	pair<Vertex*, int> GetPosition(Kmer k) {
		assert(h_->contains(k));
		return h_->get(k);
	}

	bool Contains(Kmer k) {
		return h_->contains(k);
	}

	pair<Vertex*, int> GetPosMaybeMissing(Kmer k) {
		if (!h_->contains(k)) {
			g_->AddVertex(Sequence(k));
		}
		return h_->get(k);
	}

	GraphConstructor() :
		h_(new SimpleHashTable<kmer_size_> ()), action_handler_(*h_),
				g_(new Graph(kmer_size_, action_handler_)) {
	}

public:
	virtual void ConstructGraph(Graph* &g, SimpleHashTable<kmer_size_>* &h) {
		g = g_;
		h = h_;
	}
};

//for tests!!!
template<size_t kmer_size_, size_t read_size_>
class DirectConstructor: public GraphConstructor<kmer_size_> {
	typedef Seq<kmer_size_> Kmer;
	typedef GraphConstructor<kmer_size_> super;

	vector<strobe_read<read_size_, 1>>& reads_;

	//todo extract from class definition!!!
	void ThreadRead(const Seq<read_size_> &r);

public:
	DirectConstructor(vector<strobe_read<read_size_, 1>>& reads) : super(),
		reads_(reads) {
	}

	virtual void ConstructGraph(Graph* &g, SimpleHashTable<kmer_size_>* &h) {
		for (size_t i = 0; i < reads_.size(); ++i) {
			ThreadRead(reads_[i][0]);
		}
		super::ConstructGraph(g, h);
	}
};

template<size_t kmer_size_, size_t read_size_>
void DirectConstructor<kmer_size_, read_size_>::ThreadRead(const Seq<read_size_> &r) {
	Kmer k(r);
	DEBUG("Threading k-mer: " + k.str())
	for (size_t i = kmer_size_; i < N; ++i) {
		pair<Vertex*, int> prev_pos = GetPosMaybeMissing(k);
		Kmer old_k = k;
		k = k << r[i];
		DEBUG("Threading k-mer: " + k.str())
		pair<Vertex*, int> curr_pos = GetPosMaybeMissing(k);

		Vertex* prev_v = prev_pos.first;
		Vertex* curr_v = curr_pos.first;
		size_t prev_offset = prev_pos.second;
		size_t curr_offset = curr_pos.second;

		if (super::g_->IsLastKmer(prev_v, prev_offset)
				&& super::g_->IsFirstKmer(curr_v, curr_offset)
				&& super::g_->IsMergePossible(prev_v, curr_v)) {
			super::g_->Merge(prev_v, curr_v);
		} else if (prev_v == curr_v && prev_offset + 1 == curr_offset) {
			//todo check links here to optimize???
			//do nothing
		} else {
			super::g_->SplitVertex(prev_v, prev_offset + kmer_size_);
			//need if k-mers were on same or complementary vertices
			curr_pos = GetPosition(k);
			Vertex* curr_v = curr_pos.first;
			size_t curr_offset = curr_pos.second;
			Vertex* v2 = super::g_->SplitVertex(curr_v->complement(),
					curr_v->size() - curr_offset)->complement();
			Vertex* v1 = GetPosition(old_k).first;
			super::g_->LinkVertices(v1, v2);
		}
	}
}

template<size_t kmer_size_>
class CondenseConstructor: public GraphConstructor<kmer_size_> {
	typedef Seq<kmer_size_> Kmer;
	typedef GraphConstructor<kmer_size_> super;
	typedef DeBruijn<kmer_size_> debruijn;

	DeBruijn<kmer_size_>& origin_;
public:
	CondenseConstructor(DeBruijn<kmer_size_>& origin) : origin_(origin) {

	}
	virtual void ConstructGraph(Graph* &g, SimpleHashTable<kmer_size_>* &h) {
		for (typename debruijn::kmer_iterator it = origin_.kmer_begin(), end =
				origin_.kmer_end(); it != end; it++) {
			Kmer kmer = *it;
			if (!super::Contains(kmer)) {
				//DEBUG("Starting proces for " << kmer);
				Kmer initial_kmer = kmer;
				//			DeBruijn<K>::Data d = origin.get(kmer);
				//DEBUG("Prev Count " << d.PrevCount());
				//go left while can
				while (origin_.PrevCount(kmer) == 1) {
					Kmer prev_kmer = *(origin_.begin_prev(kmer));
					//DEBUG("Prev kmer " << prev_kmer);
					//				DeBruijn<K>::Data prev_d = origin.get(prev_kmer);
					//DEBUG("Next Count of prev " << prev_d.NextCount());
					if (origin_.NextCount(prev_kmer) == 1 && kmer != !prev_kmer
							&& prev_kmer != initial_kmer) {
						kmer = prev_kmer;
					} else {
						break;
					}
				}
				//DEBUG("Stopped going left at " << kmer);
				//go right, appending sequence
				SequenceBuilder s;
				initial_kmer = kmer;
				s.append(kmer);

				//DEBUG("Next Count " << d.NextCount());

				while (origin_.NextCount(kmer) == 1) {
					Kmer next_kmer = *(origin_.begin_next(kmer));
					//DEBUG("Next kmer " << next_kmer);
					//				DeBruijn<K>::Data next_d = origin.get(next_kmer);
					//DEBUG("Prev Count of next " << next_d.PrevCount());
					if (origin_.PrevCount(next_kmer) == 1 && kmer != !next_kmer
							&& next_kmer != initial_kmer) {
						kmer = next_kmer;
						s.append(kmer[kmer_size_ - 1]);
					} else {
						break;
					}
				}
				//DEBUG("Stopped going right at " << kmer);
				super::g_->AddVertex(s.BuildSequence());
			}
		}
		for (set<Vertex*>::iterator it = super::g_->vertices().begin(), end =
				super::g_->vertices().end(); it != end; it++) {
			Vertex* v = *it;
			Kmer kmer = v->nucls().end<kmer_size_> ();

			typename debruijn::neighbour_iterator n_it = origin_.begin_next(kmer);
			for (size_t i = 0, n = origin_.NextCount(kmer); i < n; ++i, ++n_it) {
				super::g_->LinkVertices(v, super::g_->GetPosition(*n_it).first);
			}
			//todo now linking twice!!!
		}
	}
};


class Handler {
public:
	virtual void HandleStartVertex(const Vertex* v) {
	}
	virtual void HandleEndVertex(const Vertex* v) {
	}
	virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
	}
};

class Traversal {
public:
	Traversal(const Graph& g) :
		g_(g) {
	}
	virtual void Traverse(Handler& h) {
	}
protected:
	const Graph& g_;
};

class DFS: public Traversal {
	set<Vertex*> visited_;
	void go(Vertex* v, vector<Vertex*>& stack, Handler& h);
public:
	DFS(const Graph& g) :
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

class SimpleStatCounter: public Handler {
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

class CountHandler: public Handler {
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

class VisHandler: public Handler {
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
}

#endif /* CONDENSED_GRAPH_H_ */
