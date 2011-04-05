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
#include "filter_iterator.hpp"

using namespace std;
using namespace __gnu_cxx;

namespace abruijn {

//class Profile {
//protected:
//	size_t len_;
//public:
//	Profile(size_t len) : len_(len) {}
//	virtual ~Profile() {}
//	virtual void addSequence(const Sequence& seq);
//};

class FrequencyProfile {
	vector<size_t> freq_;
public:
	const static size_t SURE = -1;
	FrequencyProfile(size_t len) {
		assert(len > 0);
		for (size_t i = 0; i < (len * 4); ++i) {
			freq_.push_back(0);
		}
	}

	size_t size() const {
		return freq_.size() / 4;
	}

	void addSequence(const Sequence& seq) {
		assert(seq.size() == size());
		for (size_t i = 0; i < seq.size(); ++i) {
			freq_[i * 4 + seq[i]]++;
		}
	}

	void appendNucl(char nucl) {
		for (int i = 0; i < 4; i++) {
			if (i == nucl) {
				freq_.push_back(SURE);
			} else {
				freq_.push_back(0);
			}
		}
	}
	void append(FrequencyProfile p) {
		freq_.insert(freq_.end(), p.freq_.begin(), p.freq_.end());
	}

	ostream& output(ostream& os) const {
		for (size_t i = 0; i < freq_.size() / 4; ++i) {
			int num = 0;
			int jj = -1;
			for (int j = 0; j < 4; ++j) {
				if (freq_[(i << 2) + j]) {
					num++;
					jj = j;
				}
			}
			if (num == 1) {
				os << nucl(jj);
			} else {
				os << "(";
				for (int j = 0; j < 4; ++j) {
					if (freq_[(i << 2) + j]) {
						os << (nucl(j) - 'A' + 'a');
					}
				}
				os << ")";
			}
		}
		return os;
	}
};

typedef FrequencyProfile Profile;

ostream& operator<< (ostream& os, const Profile& p);

class Edge {
	map<size_t, Profile> long_;
	// minus keys
	map<size_t, size_t> short_;
public:
	ostream& output(ostream& os) const {
		for (map<size_t, size_t>::const_iterator it = short_.begin(); it != short_.end(); it++) {
			os << (int) -(it->first) << ";";
		}
		for (map<size_t, Profile>::const_iterator it = long_.begin(); it != long_.end(); it++) {
			os << (it->second) << ";";
		}
		return os;
	}

	void addSequence(size_t fromSize, size_t toSize, const Sequence& seq) {
		if (seq.size() > fromSize + toSize) {
			size_t len = seq.size() - fromSize - toSize;
			map<size_t, Profile>::iterator it = long_.find(len);
			if (it == long_.end()) {
				Profile p(len);
				p.addSequence(seq.Subseq(fromSize, seq.size() - toSize));
				long_.insert(make_pair(len, p));
			} else {
				it->second.addSequence(seq.Subseq(fromSize, seq.size() - toSize));
			}
		} else {
			short_[fromSize + toSize - seq.size()]++;
		}
	}
};

ostream& operator<< (ostream& os, const Edge& e);

class Vertex;

typedef map<Vertex*, Edge> Edges;

class Vertex {
	Vertex* complement_;
	const Sequence data_;
	friend ostream& operator<<(ostream&, const Vertex&);
	Edges edges_;
public:
	Vertex(const Sequence& kmer) : data_(kmer) {};
	Vertex(const Sequence& kmer, bool withComplement) : data_(kmer) {
		complement_ = new Vertex(!kmer);
		complement_->complement_ = this;
	};

	Vertex* complement() const {
		return complement_;
	}

	int size() const {
		return data_.size();
	}

	Edges& edges() {
		return edges_;
	}

	void addEdge(Vertex* to, const Sequence& seq) {
		edges_[to].addSequence(size(), to->size(), seq);
	}

	int degree() const {
		return edges_.size();
	}

	const Sequence data() const {
		return data_;
	}

	string str() const {
		return data_.str();
	}

	ostream& output(ostream& os) const {
		os << str();
		return os;
	}
};

ostream& operator<< (ostream& os, const Vertex& e);

typedef set<Vertex*> Vertices;

//class GraphIterator {
//	const Vertices& vertices_;
//	Vertices::iterator it_;
//public:
//	GraphIterator(Vertices vertices) {
//		vertices_ = vertices;
//		it_ = vertices_.begin();
//	}
//};

class Graph {
public:
	typedef hash_map < Sequence, Vertex*, HashSym<Sequence>, EqSym<Sequence> > SeqVertice;
	SeqVertice seqVertice;

	Vertices vertices;
	Vertices removedVertices;

	Graph() {}
	Vertex* createVertex(const Sequence& kmer);
	void addEdge(Vertex* from, Vertex* to, const Sequence& seq);
	void removeVertex(Vertex* v);
	void removeVertex_single(Vertex* v);
	bool hasVertex(const Sequence& kmer);
	Vertex* getVertex(const Sequence& kmer);

	bool tryCondenseA(Vertex* v);
	void cleanup();

	void output(std::ofstream &out);
	void output(string filename);

	class VertexIsAlive {
		const Graph& graph_;
	public:
		VertexIsAlive(Graph& graph) : graph_(graph) {}
		bool operator() (Vertex* v) const {
			return (graph_.removedVertices.count(v) == 0) && (graph_.vertices.count(v) > 0);
		}
	};
	typedef filter_iterator<Vertices::iterator, VertexIsAlive> iterator;
	iterator begin() {
		return iterator(vertices.begin(), vertices.end(), VertexIsAlive(*this));
	}
	iterator end() {
		return iterator(vertices.end(), vertices.end(), VertexIsAlive(*this));
	}
};

}

#endif // _abruijngraph_h_
