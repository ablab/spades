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

#endif



