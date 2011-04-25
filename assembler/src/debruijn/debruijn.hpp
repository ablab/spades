/*
 * debruijn.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 */

#ifndef DEBRUIJN_HPP_
#define DEBRUIJN_HPP_

//#include "seq.hpp"
//#include "strobe_read.hpp"
////#include <google/sparse_hash_map> // ./configure, make and sudo make install from libs/sparsehash-1.10
//#include <map>
//#include <vector>
//#include <tr1/unordered_map>
//#include <numeric>
//#include <bitset>
//#include "read.hpp"
//#include "sequence.hpp"
//
//namespace de_bruijn {
//
//template<size_t size_>
//class DeBruijn {
//private:
//	typedef Seq<size_> Kmer;
//	typedef Seq<size_ + 1> KPlusOneMer;
//
//	class Data {
//	private:
//		unsigned char edges_; // 0123 bits are outgoing acgt edges, 4567 bits are incoming acgt edges;
//	public:
//		Data() : edges_(0) {
//		}
//		void addEdgeTo(const Kmer &kmer) {
//			edges_ |=  1 << kmer.last();
//		}
//		void addEdgeFrom(const Kmer &kmer) {
//			edges_ |= 16 << kmer.first();
//		}
//		bool hasEdgeTo(char nucl) const { // nucl = 0123
//			return edges_ & (1 << nucl);
//		}
//		bool hasEdgeFrom(char nucl) const { // nucl = 0123
//			return edges_ & (16 << nucl);
//		}
//		char outEdgesCount() const {
//			return __builtin_popcountl(edges_ & 15); // count number of bits = 1 on positions 0123
//		}
//		char inEdgesCount() const {
//			return __builtin_popcountl(edges_ & 240); // count number of bits = 1 on positions 4567
//		}
//		unsigned char getInEdges() const {
//			return edges_ >> 4;
//		}
//		unsigned char getOutEdges() const {
//			return edges_ & 15;
//		}
//	};
//
//	size_t CountPositive(const bool* a) const {
//		size_t c = 0;
//		for (size_t i = 0; i < 4; ++i) {
//			if (a[i]) {
//				c++;
//			}
//		}
//		return c;
//	}
//
//	//	typedef Data value;
//	//typedef google::sparse_hash_map<key, value,	typename key::hash, typename key::equal_to> hash_map;
//	//	typedef std::map<key, value, typename key::less> map_type;
//	typedef std::tr1::unordered_map<Kmer, Data, typename Kmer::hash> map_type;
//	typedef typename map_type::iterator map_iterator;
//	typedef typename map_type::const_iterator map_const_iterator;
//	map_type nodes_;
//
//	void CountSequence(const Sequence& s) {
//		Seq<size_> head = s.start<size_> ();
//		for (size_t j = size_; j < s.size(); ++j) {
//			Seq<size_> tail = head << s[j];
//			addEdge(head, tail);
//			head = tail;
//		}
//	}
//
//	const Data& get(const Kmer &kmer) const {
//		map_const_iterator mci = nodes_.find(kmer);
//		assert(mci != nodes_.end());
//		return mci->second;
//	}
//
//	Data& addNode(const Kmer &seq) {
//		std::pair<const Kmer, Data> p = make_pair(seq, Data());
//		std::pair<typename map_type::iterator, bool> node = nodes_.insert(p);
//		return node.first->second; // return node's data
//	}
//
//	void CountRead(const Read &read) {
//		if (read.isValid()) {
//			Sequence s = read.getSequence();
//			CountSequence(s);
////			CountSequence(!s);
//		}
//	}
//
//public:
//
//	DeBruijn() {
//		// empty debruijn graph
//	}
//
//	void addEdge(const Kmer &from, const Kmer &to) {
//		Data &d_from = addNode(from);
//		Data &d_to = addNode(to);
//		d_from.addEdgeTo(to);
//		d_to.addEdgeFrom(from);
//		//d_from.out_edges_[(size_t) to[size_ - 1]] = true; // was ++
//		//d_to.in_edges_[(size_t) from[0]] = true; // was ++
//	}
//
//	//void addEdge(const KPlusOneMer &edge) {
//	//	addEdge(edge.start(), edge.end());
//	//}
//
//	//	size_t size() const {
//	//		return nodes_.size();
//	//	}
//
//	class edge_iterator {
//	private:
//		Kmer kmer_;
//		char pos_;
//		unsigned char neighbours_;
//		bool right_neighbours_;
//
//		void ShiftPos() {
//			while (pos_ < 4 && !(neighbours_ & (1 << pos_)) ) {
//				pos_++;
//			}
//		}
//	public:
//
//		bool isEnd() const {
//			return (pos_ >= 4);
//		}
//
//		edge_iterator(Kmer kmer, char pos, unsigned char neighbours,
//				bool right_neighbours) :
//			kmer_(kmer), pos_(pos), neighbours_(neighbours),
//					right_neighbours_(right_neighbours) {
//			ShiftPos();
//		}
//
//		edge_iterator& operator++() {
//			if (pos_ < 4) {
//				pos_++;
//				ShiftPos();
//			}
//			return *this;
//		}
//
//		KPlusOneMer operator *() const {
//			return right_neighbours_ ? kmer_.pushBack(pos_) : kmer_.pushFront(pos_);
//		}
//	};
//
//	size_t OutgoingEdgeCount(const Kmer &kmer) const {
//		return get(kmer).outEdgesCount();
//	}
//
//	size_t IncomingEdgeCount(const Kmer &kmer) const {
//		return get(kmer).inEdgesCount();
//	}
//
//	edge_iterator OutgoingEdges(const Kmer &kmer) const {
//		unsigned char out_edges = get(kmer).getOutEdges();
//		return edge_iterator(kmer, 0, out_edges, true);
//	}
//
//	edge_iterator IncomingEdges(const Kmer &kmer) const {
//		unsigned char in_edges = get(kmer).getInEdges();
//		return edge_iterator(kmer, 0, in_edges, false);
//	}
//
//	class kmer_iterator: public map_type::const_iterator {
//	public:
//		typedef typename map_type::const_iterator map_iterator;
//
//		kmer_iterator(const map_iterator& other) :
//			map_type::const_iterator(other) {
//		}
//		;
//
//		const Kmer& operator *() const {
//			return map_iterator::operator*().first;
//		}
//	};
//
//	kmer_iterator begin() const {
//		return kmer_iterator(nodes_.begin());
//	}
//
//	kmer_iterator end() const {
//		return kmer_iterator(nodes_.end());
//	}
//
//	typedef edge_iterator EdgeIterator;
//
//	typedef kmer_iterator VertexIterator;
//
//	typedef Seq<size_> VertexId;
//	typedef Seq<size_ + 1> EdgeId;
//
//	//template<size_t size2_, size_t count_>
//	void ConstructGraph(const vector<Read> &v) {
//		for (size_t i = 0; i < v.size(); ++i) {
//			CountRead(v[i]);
//		}
//	}
//
//	template<class ReadStream>
//	void ConstructGraphFromStream(ReadStream& stream) {
//		Read r;
//		while (!stream.eof()) {
//			stream >> r;
//			CountRead(r);
//		}
//	}
//
//	//	void show(string genome) {
//	//		int arr[100];
//	//		for (int i = 0; i < 100; i++)
//	//			arr[i] = 0;
//	//		Kmer oppa(genome);
//	//		set<KPlusOneMer, typename KPlusOneMer::less> s;
//	//		for (int i = size_; i < genome.length(); i++) {
//	//			typename map_type::iterator it = nodes_.find(oppa);
//	//			int c = dignucl(genome[i]);
//	//			if (it == nodes_.end()) {
//	//				cout << "gopa" << endl;
//	//				if(s.find(oppa.pushBack(c)) == s.end()) {
//	//					arr[0]++;
//	//					s.insert(oppa.pushBack(c));
//	//				}
//	//			} else {
//	//				Data& d = it->second;
//	//				cout << d.out_edges_[(int) c] << endl;
//	////				if (d.out_edges_[(int) c] > 0) {
//	//					if (d.out_edges_[(int) c] < 100)
//	//						if(s.find(oppa.pushBack(c)) == s.end()) {
//	//							arr[d.out_edges_[(int) c]]++;
//	//							s.insert(oppa.pushBack(c));
//	//						}
//	////				}
//	//			}
//	//			oppa = oppa << genome[i];
//	//		}
//	//		for (int i = 0; i < 100; i++) {
//	//			cout << i << " " << arr[i] << endl;
//	//		}
//	//	}
//
//};
//
//}

#endif /* DEBRUIJN_HPP_ */
