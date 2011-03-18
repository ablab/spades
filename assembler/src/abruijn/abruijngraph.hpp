#ifndef _abruijngraph_h_
#define _abruijngraph_h_

#include <cassert>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <ext/hash_map>
#include "seq.hpp"
#include "sequence.hpp"
#include "parameters.hpp"
#include "hash.hpp"
#include "logging.hpp"

using namespace std;
using namespace __gnu_cxx;

namespace abruijn {

class Edge {
public:
	map<size_t, int> lengths_;

	void addLength(int len) {
		lengths_[len]++;
	}

	string toString() {
		stringstream ss;
		for (map<size_t, int>::iterator it = lengths_.begin(); it != lengths_.end(); ++it) {
			ss << " " << it->first << " {" << it->second << "};";
		}
		return ss.str();
	}
};

class Vertex {
	const Sequence* data_;
	int size_;
public:
	Vertex* complement_;
	typedef map<Vertex*, Edge> Edges;
	Edges edges_;
	Vertex(const Sequence* kmer) : data_(kmer), size_(kmer->size()) {};

	int size() {
		return size_;
	}
	//TODO trash
	void setSize(int size) {
		size_ = size;
	}

	void addEdge(Vertex* to, int len) {
		edges_[to].addLength(len);
	}

	int degree() {
		return edges_.size();
	}

	bool is(const Sequence* kmer) {
		return (*data_) == (*kmer);
	}

	string str() {
		return data_->str();
	}

	Sequence* concat(Vertex* u) {
		SequenceBuilder sb;
//		size_t sum = data_->size() + u->data_->size();
//		sb.append(*(v->kmer_));
//		if (len > sum) {
//			TRACE("Should've been filled correctly."); // TODO
//			for (size_t i = sum; i < len; i++) {
//				sb.append('A');
//			}
//			sb.append(*(u->kmer_));
//		} else {
//			sb.append(u->kmer_->Subseq(sum - len));
//		}
//		assert(sb.size() == len);
//		return new Sequence(sb.BuildSequence());

		sb.append(data_->Subseq(0, K / 2));
		sb.append(u->data_->Subseq(u->data_->size() - K / 2));
		return new Sequence(sb.BuildSequence());

//		for (int i = 0; i < K / 2; i++) {
//			sb.append(data_->operator [](i));
//			DEBUG(i);
//		}
//		return sb.BuildSequence();
//		size_t start = u->data_->size() - K / 2;
//		for (int i = 0; i < K / 2; i++) {
//			sb.append(u->data_[start + i]);
//		}
//		return sb.BuildSequence();
	}



	string toString() {
		#ifdef OUTPUT_PAIRED
			return data_->Subseq(0, LABEL).str() + "_" + itoa(size());
		#endif
		#ifndef OUTPUT_PAIRED
			return data_->Subseq(0, LABEL).str() + "_" + itoa(size()) + "_"+ data_->Subseq(data_->size() - LABEL).str();
		#endif
	}
};

class Graph {
public:
	typedef hash_map < Sequence, Vertex*, HashSym<Sequence>, EqSym<Sequence> > SeqVertice;
	SeqVertice seqVertice;

	typedef set<Vertex*> Vertices;
	Vertices vertices;

	Graph() {}
	Vertex* createVertex(const Sequence* kmer);
	Vertex* createVertex(const Sequence* kmer, size_t size);
	void addEdge(Vertex* from, Vertex* to, int len);
	void removeVertex(Vertex* v);
	void removeVertex_single(Vertex* v);
	bool hasVertex(const Sequence* kmer);
	Vertex* getVertex(const Sequence kmer);

	Vertex* condense(Vertex* v);

	void output(std::ofstream &out);
	void output(string filename);
};

}

#endif // _abruijngraph_h_
