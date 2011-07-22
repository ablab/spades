/*
 * config.hpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#ifndef HAMMER_CONFIG_HPP
#define HAMMER_CONFIG_HPP


#include <unordered_map>
#include <algorithm>
#include <map>
#include "sequence/seq.hpp"
#include "kmer_stat.hpp"
#include "read/read.hpp"

#include "uf.hpp"

using namespace std;
using namespace __gnu_cxx;

#define K 55

#define GOOD_SINGLETON_THRESHOLD 1 

typedef Seq<K> KMer;
typedef iufstream<K> UFStream;
typedef UFCluster<K> MyUFC;
typedef unordered_map<KMer, KMer, KMer::hash> KMerHashMap;

//typedef pair<KMer, KMerStat> KMerCount;
//typedef vector<KMerCount> KMerStatVector;

struct StringKMer{
	uint32_t start;
	uint64_t kmerno;
};
typedef vector<StringKMer> StringKMerVector;

/*inline bool SKgreater ( const StringKMer & elem1, const StringKMer & elem2 ) {
   return lexicographical_compare(elem1.sub.begin(), elem1.sub.end(), elem2.sub.begin(), elem2.sub.end());
}*/

/*inline bool SKequal ( const StringKMer & elem1, const StringKMer & elem2, const vector<KMerCount> & km, const int tau) {
	size_t j = elem2.start;
	for (size_t i = elem1.start; i < K; i += tau+1) {
		if (km[elem1.kmerno].first[i] != km[elem2.kmerno].first[j]) {
			return false;
		}
		j += tau+1;
	}
	return true;
}

inline bool SKgreater2 ( const StringKMer & elem1, const StringKMer & elem2, const vector<KMerCount> & km, const int tau ) {
	size_t j = elem2.start;
	for (size_t i = elem1.start; i < K; i += tau+1) {
		if (km[elem1.kmerno].first[i] != km[elem2.kmerno].first[j]) {
			return (km[elem1.kmerno].first[i] < km[elem2.kmerno].first[j]);
		}
		j += tau+1;
	}
	return false;
}*/


/*inline bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	for (size_t i = 0; i < K; ++i) {
		if (l.first[i] != r.first[i]) {
			return (l.first[i] < r.first[i]);
		}
	}
	return false;
}*/

inline bool KMerLess(const KMer &l, const KMer &r) {
	for (size_t i = 0; i < l.size(); ++i) {
		if (l[i] != r[i]) {
			return (l[i] < r[i]);
		}
	}
	return false;
}

// typedef map<KMer, KMerStat, KMer::less2> KMerStatMap;

#endif



