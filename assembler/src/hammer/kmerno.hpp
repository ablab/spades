#ifndef KMERNO_HPP_
#define KMERNO_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include <boost/unordered_map.hpp>
#include "google/sparse_hash_map"
#include "google/dense_hash_map"

#define GOOGLE_SPARSE_MAP
//#define BOOST_UNORDERED_MAP
//#define GOOGLE_DENSE_MAP

#include "kmer_stat.hpp"

const int KMERNO_HASH_MODULUS = 32451233;
const int KMERNO_HASH_Q = 1299673;
const int KMERNO_HASH_Q_INV = 31471908;
const int KMERNO_HASH_Q_POW_K_MINUS_ONE = 12533099;

struct KMerNo {
	hint_t index;
	double errprob;

	KMerNo( hint_t no, double qual ) : index(no), errprob(qual) { }
	KMerNo( hint_t no ) : index(no), errprob(1) { }
	KMerNo( ) : index(-1), errprob(1) { }

	bool equal(const KMerNo & kmerno) const;
	std::string str() const;
	bool less(const KMerNo &r) const;
	bool greater(const KMerNo &r) const;

	static hint_t new_hash( hint_t index );
	static hint_t next_hash( hint_t old_hash, hint_t new_index );
	static void precomputeHashes();

	struct hash {
		uint32_t operator() (const KMerNo &kn) const;
	};

	struct are_equal {
		bool operator() (const KMerNo &l, const KMerNo &r) const;
	};
};

#ifdef GOOGLE_SPARSE_MAP
	typedef google::sparse_hash_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;
#endif
#ifdef BOOST_UNORDERED_MAP
	typedef boost::unordered_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;
#endif
#ifdef GOOGLE_DENSE_MAP
	typedef google::dense_hash_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;
#endif


#endif

