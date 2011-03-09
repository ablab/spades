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
#include "sequence.hpp"
#include "logging.hpp"
#include "nucl.hpp"
#include <iostream>

using namespace std;

namespace condensed_graph {
#define K 11//5//25
#define N 100//11//100
#define HASH_SEED 1845724623

typedef Seq<K> Kmer;
typedef Seq<K - 1> KMinusOneMer;
typedef Seq<N> Read;

LOGGER("debruijn.condensed_graph")

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
	Vertex(const Sequence &nucls);
	Vertex(const Sequence &nucls, Vertex** desc);
	~Vertex();
	int DescCount();
	bool IsDeadend();
	Vertex** desc();
	Vertex* desc(char nucl);
	size_t size();
	Sequence nucls();
	void AddDesc(Vertex* v);
	Vertex* complement();
	void set_complement(Vertex* complement);
	int coverage();
	void set_coverage(int coverage);
};


/*
 * To Sergey: please use Kmer::hash instead
 */
/*class SimpleHash {
public:
	unsigned int operator()(const Kmer& seq) const {
		unsigned int h = HASH_SEED;
		for (size_t i = 0; i < seq.size(); i++) {
			h = ((h << 5) - h) + seq[i];
		}
		return h;
	}
};*/

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
			h_[k] = make_pair(v,s);
		}
		else { // change existing element
			hi->second = make_pair(v,s);
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
		DEBUG("Getting position of k-mer '" + k.str() +  "' Position is " <<  hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}
};

class Graph {
	set<Vertex*> component_roots_;
	SimpleHashTable h_;

	void RenewKmersHash(Vertex* v);

	/**
	 * deals with incoming links and their complement only!!!
	 */
	void FixIncomingOnSplit(Vertex* v, Vertex* v1, Vertex* v2);

	void FixIncomingOnMerge(Vertex* v1, Vertex* v2, Vertex* v);

	bool CanBeDeleted(Vertex* v) const;

public:
	const set<Vertex*>& component_roots() const;

	vector<Vertex*> Anc(Vertex* v) const;

	vector<Vertex*> Desc(Vertex* v) const;

	bool IsLast(Vertex* v) const;

	bool IsFirst(Vertex* v) const;

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

	pair<Vertex*, int> GetPosMaybeMissing(Kmer k);

	void ThreadRead(const Read &r);

	bool IsLastKmer(Vertex* v, size_t pos) const;
	bool IsFirstKmer(Vertex* v, size_t pos) const;
};

/*class VertexPool {
 _v_idx _size;
 bool* _free;
 Vertex* _vertices;
 _v_idx _max_idx;
 public:
 VertexPool(_v_idx size);
 ~VertexPool();
 Vertex& operator[](_v_idx index) const;
 bool IsFree(_v_idx index) const;
 void Free(_v_idx index);
 _v_idx AllocatePair();
 };*/

}
#endif /* CONDENSED_GRAPH_H_ */
