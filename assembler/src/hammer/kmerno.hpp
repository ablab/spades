#ifndef KMERNO_HPP_
#define KMERNO_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include <boost/unordered_map.hpp>
#include "google/sparse_hash_map"
#include "google/dense_hash_map"

//#define GOOGLE_SPARSE_MAP
#define BOOST_UNORDERED_MAP
//#define GOOGLE_DENSE_MAP

#include "kmer_stat.hpp"

const uint64_t KMERNO_HASH_MODULUS = 2305843009213693951;
const uint64_t KMERNO_HASH_Q = 3712758430079221;
const uint64_t KMERNO_HASH_Q_INV = 2250585152990002931;
const uint64_t KMERNO_HASH_Q_POW_K_MINUS_ONE = 412252044596125152;

struct KMerNo {
	hint_t index;
	double errprob;

	KMerNo( hint_t no, double qual ) : index(no), errprob(qual) { }
	KMerNo( hint_t no ) : index(no), errprob(1) { }
	KMerNo( ) : index(-1), errprob(1) { }

	bool equal(const KMerNo & kmerno) const;
	bool equal(const KMerCount & kmc) const;
	std::string str() const;
	bool less(const KMerNo &r) const;
	bool greater(const KMerNo &r) const;

	static uint64_t new_hash( hint_t index );
	static uint64_t next_hash( uint64_t old_hash, hint_t new_index );
	//static void precomputeHashes();

	struct hash {
		uint64_t operator() (const KMerNo &kn) const;
	};

	struct string_hash {
		uint64_t operator() (const std::string &kn) const;
	};

	struct are_equal {
		bool operator() (const KMerNo &l, const KMerNo &r) const;
	};

	struct is_less {
		bool operator() (const KMerNo &l, const KMerNo &r) const;
	};

	struct is_less_kmercount {
		bool operator() (const KMerCount &l, const KMerCount &r) const;
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

