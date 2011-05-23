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
typedef map<KMer, KMerStat, KMer::less2> KMerStatMap;
typedef pair<KMer, KMerStat> KMerCount;
typedef vector<KMerCount> KMerStatVector;
typedef map<string, size_t> StringCountMap;
typedef pair<string, size_t> StringCount;
typedef vector<StringCount> StringCountVector;

bool SCgreater ( const StringCount & elem1, const StringCount & elem2 )
{
   return lexicographical_compare(elem1.first.begin(), elem1.first.end(), elem2.first.begin(), elem2.first.end());
}

bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	for (size_t i = 0; i < K; ++i) {
		if (l.first[i] != r.first[i]) {
			return (l.first[i] < r.first[i]);
		}
	}
	return false;
}


#endif



