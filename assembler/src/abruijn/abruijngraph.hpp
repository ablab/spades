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

namespace abruijn {

using namespace std;
using namespace __gnu_cxx;

LOGGER("a.graph");

/**
* @name FrequencyProfile 
*/
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

/**
 * @see Sequence
 */
class Edge {
	typedef map<size_t, size_t> Lengths;
	Lengths lengths_;
public:
	ostream& output(ostream& os) const {
		for (map<size_t, size_t>::const_iterator it = lengths_.begin(); it != lengths_.end(); it++) {
			os << (int) (it->first /*- K*/) << "(*" << it->second << ");";
		}
		return os;
	}

	void addSequence(const Sequence& seq, size_t quantity = 1) {
		size_t length = seq.size() - K;
		lengths_[length] += length * quantity;
	}

//	void touch() {
//		lengths_[7] += 7;
//	}

	Edge& operator+=(Edge& edge) {
//		for (Lengths::const_iterator it = edge.lengths_.begin(); it != edge.lengths_.end(); ++it) {
//			lengths_[it->first] += it->second;
//		}
		return *this;
	}

	Edge add(const Edge& edge) {
		Edge r;
		for (Lengths::const_iterator it = lengths_.begin(); it != lengths_.end(); ++it) {
			r.lengths_[it->first] += it->second;
		}
		for (Lengths::const_iterator it = edge.lengths_.begin(); it != edge.lengths_.end(); ++it) {
			r.lengths_[it->first] += it->second;
		}
		return r;
	}

	Edge concat(const Edge& edge) {
		Edge r;
		for (Lengths::const_iterator it = lengths_.begin(); it != lengths_.end(); ++it) {
			for (Lengths::const_iterator it2 = edge.lengths_.begin(); it2 != edge.lengths_.end(); ++it2) {
				r.lengths_[it->first + it2->first] += it->second + it2->second;
			}
		}
		return r;
	}
};

/**
 * @see Profile
 * @see Sequence
 */
class ProfileEdge {
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

/**
* @see Edges
*/
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
		edges_[to].addSequence(seq);
	}

	int degree() const {
		return edges_.size();
	}

	/**
	 * Returns ANY right neighbour
	 */
	Vertex* forward() const {
		return edges_.begin()->first;
	}

	/**
	 * Returns ANY left neighbour
	 */
	Vertex* backward() const {
		return complement_->forward()->complement_;
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

/**
* @see Sequence
* @see SeqVertice
* @see Vertices
*/
class Graph {
public:
	typedef hash_map < Sequence, Vertex*, hashing::HashSym<Sequence>, hashing::EqSym<Sequence> > SeqVertice;
	SeqVertice seqVertice;

	Vertices vertices;
	Vertices removedVertices;

	Graph() : vertex_is_alive(*this) {}
	Vertex* createVertex(const Sequence& kmer);
	void addEdge(Vertex* from, Vertex* to, const Sequence& seq);
	void removeVertex(Vertex* v);
	void removeVertex_single(Vertex* v);
	bool hasVertex(const Sequence& kmer);
	Vertex* getVertex(const Sequence& kmer);

	void Condense();
	bool Condense(Vertex* v);
	void Condense_single(Vertex* v);
	void cleanup();
	void stats();

	void output(std::ofstream &out, bool paired);
	void output(string filename, bool paired);

	class VertexIsAlive {
		const Graph& graph_;
	public:
		VertexIsAlive(Graph& graph) : graph_(graph) {}
		bool operator() (Vertex* v) const {
			return (graph_.removedVertices.count(v) == 0) && (graph_.vertices.count(v) > 0);
		}
	};
	VertexIsAlive vertex_is_alive;

	typedef filter_iterator<Vertices::iterator, VertexIsAlive> iterator;
	iterator begin() {
		return iterator(vertices.begin(), vertices.end(), vertex_is_alive);
	}
	iterator end() {
		return iterator(vertices.end(), vertices.end(), vertex_is_alive);
	}
};

}

#endif // _abruijngraph_h_
