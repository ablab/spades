/**
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include <vector>
#include <set>
#include <ext/hash_map>
#include <cstring>
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"

using namespace std;
using namespace __gnu_cxx;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

namespace condensed_graph {

#define K 25
#define N 100
#define HASH_SEED 1845724623

typedef Seq<K> Kmer;
typedef Seq<K - 1> KMinusOneMer;
typedef Seq<N> Read;

class Vertex;

class Vertex {

	Sequence nucls_;
	Vertex* desc_[4];
	Vertex* complement_;

	int coverage_;
	int arc_coverage_[4];
public:
	Vertex(Sequence nucls);
	Vertex(Sequence nucls, Vertex** desc);
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

class SimpleHash {
public:
	unsigned int operator()(const Kmer& seq) const {
		unsigned int h = HASH_SEED;
		for (int i = 0; i < seq.size(); i++) {
			h = ((h << 5) - h) + seq[i];
		}
		return h;
	}
};

class MySimpleHashTable {
	hash_map<const Kmer, pair<Vertex*, size_t> , SimpleHash, Kmer::equal_to> h_;
public:
	void put(Kmer k, pair<Vertex*, size_t> v) {
		h_.insert(make_pair(k, v));
	}
	bool contains(Kmer k) {
		return h_.count(k) == 1;
	}
	const pair<Vertex*, size_t> get(Kmer k) {
		assert(contains(k));
		return h_[k];
	}
};

class Graph {
	set<Vertex*> component_roots_;
	MySimpleHashTable h_;

	void RenewHashForSingleVertexKmers(Vertex* v);

public:
	const set<Vertex*>& component_roots();

	vector<Vertex*> Anc(Vertex* v);

	vector<Vertex*> Desc(Vertex* v);

	bool IsLast(Vertex* v);

	bool IsFirst(Vertex* v);

	bool AddIfRoot(Vertex* v);

	/**
	 * adds vertex and its complement
	 */
	Vertex* AddVertex(Sequence nucls);

	bool IsMergePossible(Vertex* v1, Vertex* v2);

	bool CanBeDeleted(Vertex* v);

	Vertex* Merge(Vertex* v1, Vertex* v2);

	/**
	 * deletes vertex and its complement
	 */
	void DeleteVertex(Vertex* v);

	bool AreLinkable(Vertex* v1, Vertex* v2);

	void LinkVertices(Vertex* anc, Vertex* desc);

	/**
	 * deals with incoming links and their complement only!!!
	 */
	void FixIncoming(Vertex* v, Vertex* new_v);

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewHashForVertexKmers(Vertex* v);

	//pos exclusive! (goes into second vertex)
	//deletes vertex if actual split happens
	Vertex* SplitVertex(Vertex* v, size_t pos, bool return_second);

	pair<Vertex*, int> GetPosMaybeMissing(Kmer k);

	void ThreadRead(Read r);

	bool IsLastKmer(Vertex* v, size_t pos);
	bool IsFirstKmer(Vertex* v, size_t pos);
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
