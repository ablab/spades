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
#include "uf.hpp"

using namespace std;
using namespace __gnu_cxx;

#define K 55

typedef Seq<K> KMer;
typedef iufstream<K> UFStream;
typedef UFCluster<K> MyUFC;
typedef unordered_map<KMer, KMer, KMer::hash> KMerHashMap;

struct KMerStat {
	size_t count;
	float freq;
};
typedef pair<KMer, KMerStat> KMerCount;
typedef vector<KMerCount> KMerStatVector;
typedef map<string, size_t> StringCountMap;
typedef pair<string, size_t> StringCount;
typedef vector<StringCount> StringCountVector;

inline bool SCgreater ( const StringCount & elem1, const StringCount & elem2 ) {
   return lexicographical_compare(elem1.first.begin(), elem1.first.end(), elem2.first.begin(), elem2.first.end());
}

struct StringKMer{
	string sub;
	size_t count;
	KMer kmer;
};
typedef vector<StringKMer> StringKMerVector;

inline bool SKgreater ( const StringKMer & elem1, const StringKMer & elem2 ) {
   return lexicographical_compare(elem1.sub.begin(), elem1.sub.end(), elem2.sub.begin(), elem2.sub.end());
}

inline bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	for (size_t i = 0; i < K; ++i) {
		if (l.first[i] != r.first[i]) {
			return (l.first[i] < r.first[i]);
		}
	}
	return false;
}

inline bool KMerLess(const KMer &l, const KMer &r) {
	for (size_t i = 0; i < l.size(); ++i) {
		if (l[i] != r[i]) {
			return (l[i] < r[i]);
		}
	}
	return false;
}

typedef map<KMer, KMerStat, KMer::less2> KMerStatMap;

#endif



