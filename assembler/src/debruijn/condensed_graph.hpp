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
#include "../seq.hpp"

using namespace std;
using namespace __gnu_cxx;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

//typedef int _v_idx;

namespace assembler {

#define K 25
#define HASH_SEED 1845724623

typedef Seq<K> Kmer;

static const int k = 25;

char toIndex(char c) {
	switch (c) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		break;//throw exception
	}
}

class Vertex;

class Vertex {

	Sequence _nucls;
	Vertex* _desc[4];
	Vertex* _complement;

	int _coverage;
	int _arc_coverage[4];
public:

	Vertex(Sequence nucls) :
		_nucls(nucls) {
		fill_n(_desc, 4, (Vertex*) NULL);
	}

	//todo talk with Kolya about problem with Sequence copying!!!
	Vertex(Sequence nucls, Vertex** desc) :
		_nucls(nucls) {
		memcpy(desc, _desc, 4 * sizeof(Vertex*));
		//		_arcs = arcs;
	}
	;

	~Vertex() {
		delete [] _desc;
		delete [] _arc_coverage;
 	}

	//static Vertex AbsentVertex = Vertex(0, 0, NULL, true, 0, NULL);

	size_t size() {
		return _nucls.size();
	}
	;

	char operator[](const int &index) const {
		return _nucls[index];
	}

	int DescCount() {
		int c = 0;
		for (int i = 0; i < 4; ++i)
			if (_desc[i] != NULL)
				c++;
		return c;
	}
	;

	Vertex** desc() {
		return _desc;
	}
	;

	Vertex* desc(char nucl) {
		return _desc[(int)nucl];
	}

	Sequence nucls() {
		return _nucls;
	}

	void AddDesc(Vertex* v) {
		//check if k-1 mers differ
		int k_th_nucl = (*v)[K];
		_desc[k_th_nucl] = v;
	}

	Vertex* complement() {
		return _complement;
	}
	;

	void set_complement(Vertex* complement) {
		_complement = complement;
	}
	;

	int coverage() {
		return _coverage;
	}
	;

	void set_coverage(int coverage) {
		_coverage = coverage;
	}
	;
};

class Graph {
	set<Vertex*> _component_roots;

	bool Empty(Vertex** vs) {
		for (int i = 0; i < 4; ++i) {
			if (vs[i] != NULL)
				return false;
		}
		return true;
	}

public:
	const set<Vertex*>& component_roots() {
		return _component_roots;
	}

	vector<Vertex*> Anc(Vertex* v) {
		vector<Vertex*> ans;
		Vertex** compl_desc = v->complement()->desc();
		for (int i = 3; i >= 0; --i) {
			if (compl_desc[i] != NULL) {
				ans.push_back(compl_desc[i]->complement());
			}
		}
		return ans;
	}

	vector<Vertex*> Desc(Vertex* v) {
		vector<Vertex*> ans;
		Vertex** desc = v->desc();
		for (int i = 0; i < 4; ++i) {
			if (desc[i] != NULL) {
				ans.push_back(desc[i]);
			}
		}
		return ans;
	}
	/**
	 * adds two complement vertices
	 */
	Vertex* AddVertices(Sequence nucls, Vertex** outgoing_vert, Vertex** outgoing_vert_for_compl) {
		Vertex* v1 = new Vertex(nucls, outgoing_vert);
		Vertex* v2 = new Vertex(!nucls, outgoing_vert_for_compl);
		v1->set_complement(v2);
		v2->set_complement(v1);
		if (Empty(outgoing_vert)) {
			_component_roots.insert(v1);
		}
		return v1;
	}

	Vertex* AddVertices(Sequence nucls) {
		Vertex* v1 = new Vertex(nucls);
		Vertex* v2 = new Vertex(!nucls);
		v1->set_complement(v2);
		v2->set_complement(v1);
		return v1;
	}

	/**
	 * deletes vertex and its complement
	 */
	void DeleteVertices(Vertex* v) {
		Vertex* complement = v->complement();
		delete v;
		delete complement;
	}

	void LinkVertices(Vertex* anc, Vertex* desc) {
		anc->AddDesc(desc);
		desc->complement()->AddDesc(anc->complement());
	}

	/**
	 * deals with incoming links and their complement only!!!
	 */
	void FixIncoming(Vertex* v, Vertex* new_v) {
		//todo throw exception if start k-mers differ
		vector<Vertex*> anc = Anc(v);
		for (size_t i = 0; i < anc.size(); ++i) {
			LinkVertices(anc[i], new_v);
		}
	}

	//pos exclusive! (goes into second vertex)
	Vertex* SplitVertex(Vertex* v, size_t pos) {
		if (pos == v->size() - 1) {
			return v;
		};

		Sequence nucls = v->nucls();

		Vertex* v1 = AddVertices(nucls.Subseq(0, pos));
		Vertex* v2 = AddVertices(nucls.Subseq(pos - (k - 1), nucls.size()));

		LinkVertices(v1, v2);

		FixIncoming(v, v1);

		FixIncoming(v->complement(), v2->complement());

		return v2;
	}

	bool SplitVertexRenewHash(Vertex* v, size_t pos) {
		return false;
	}
};

/*class KmerPos {
	//vertex, offset and direction of the canonical
	Vertex* _v;
	int _offset;
public:
	KmerPos(Vertex* v, int offset) :
		_v(v), _offset(offset) {
	}
};*/

class SimpleHash {
public:
	unsigned int operator() (const Kmer& seq) const {
	    unsigned int h = HASH_SEED;
		for (int i = 0; i < seq.size(); i++) {
			h = ((h << 5) - h) + seq[i];
		}
		return h;
	}
};

class MySimpleHashTable {
	//pair<Vertex*, int> - vertex and offset
//	hash_map<const Kmer, pair <Vertex*, int> , SimpleHash> h;
	hash_map<const Kmer, pair<Vertex*, size_t>, SimpleHash, Kmer::equal_to > h;
	//vector<V> array[size];
public:
	//todo think of using references
	void put(Kmer k, pair <Vertex*, size_t> v) {
		h.insert(make_pair(k, v));
	}

	const pair<Vertex*, size_t> get(Kmer k) {
		return h[k];
	}

//	void remove(Kmer k) {
//		h.
//	}
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
