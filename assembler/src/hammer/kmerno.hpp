#ifndef KMERNO_HPP_
#define KMERNO_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include <boost/unordered_map.hpp>

#include "kmer_stat.hpp"
#include "globals.hpp"

struct KMerNo {
	hint_t index;
	double errprob;

	KMerNo( hint_t no, double qual ) : index(no), errprob(qual) { }
	KMerNo( hint_t no ) : index(no), errprob(1) { }

	bool equal(const KMerNo & kmerno) const;
	string str() const;
	bool less(const KMerNo &r) const;
	bool greater(const KMerNo &r) const;

	struct hash {
		size_t operator() (const KMerNo &kn) const;
	};

	struct are_equal {
		bool operator() (const KMerNo &l, const KMerNo &r) const;
	};

};

typedef boost::unordered_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;

#endif

