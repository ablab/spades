#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <ext/hash_map>
#include <algorithm>
#include <string>
#include <set>

#include <tr1/unordered_map>
#include "logging.hpp"

#define forn(i, n) for(size_t i = 0; i < (size_t) n; i++)
#define ll long long
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define edgesMap  std::tr1::unordered_map<ll, vector<EdgePrototype *> >
#define verticesMap  std::tr1::unordered_map<ll, vector<VertexPrototype *> >
#define longEdgesMap  std::tr1::unordered_map<int, Edge*>

#define otherDirection(direction) (direction == LEFT ? RIGHT : LEFT)
#define RIGHT 1
#define LEFT -1

#define IN_EDGE 0
#define OUT_EDGE 1

//LOGGER("paireddebruijn.common");

#define MAX_VERT_NUMBER 50000
#define MAX_DEGREE 50
#define suffix "_const_d"


using namespace std;
//const string parsed_reads = "data/reads_const_d.txt";
//const string parsed_reads = "data/filtered_reads";
/*
const string parsed_k_l_mers = "data/klmers_const_d.txt";
const string parsed_k_sequence = "data/vertices_const_d.txt";
const string error_log = "data/error.log";
const string parsed_l_mers = "data/lmers_const_d.txt";
const string graph_file = "data/graph_const_d.dot";
const string graph2 = "data/graph2_const_d.dot";
const string threaded_graph = "data/threaded_graph_const_d.dot";
const string auxilary_lmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
*/
//const string parsed_reads = string("data/reads") + suffix + ".txt";


/*const string parsed_reads = "data/filtered_reads";
const string parsed_k_l_mers = string("data/klmers") + suffix  + ".txt";
const string parsed_k_sequence = string("data/vertices") + suffix  + ".txt";
const string error_log = "data/error.log";
const string parsed_l_mers = string("data/lmers")  + suffix  + ".txt";
const string graph_file = string("data/grapht")  + suffix  + ".dot";
const string graph2 = string("data/graph2") + suffix  + ".dot";
const string threaded_graph = string("data/threaded_graph") +  + ".dot";
*/
const string auxilary_lmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

const string ini_file = "data/paireddebruijn/paired.ini";
extern string distance_type;
extern string parsed_reads;
extern string parsed_k_l_mers;
extern string parsed_l_mers;
extern string parsed_l_k_mers;
extern string parsed_k_mers;

extern string parsed_k_sequence;
extern string error_log;
extern string graph_file;
extern string graph2;
extern string threaded_graph;
extern string folder;

extern int k;
extern int l;
extern int readLength;
const int maxSeqLength = 200;
extern int insertLength;
extern int minIntersect;
extern int inClusterMaxShift;
extern int useKmersVertices;
extern int useRevertedPairs;


extern int fictiveSecondReads;
extern int needPairs;
extern int needLmers;
extern int needRevertedPairs;
extern int needSequences;
extern int needGraph;
extern int useExpandDefinite;
extern int useExtractDefinite;
extern int useTraceReads;
extern int useProcessLower;


void initConstants(string ini_file);
ll pushNucleotide(ll kMer, int length, int direction, int nucl);
ll popNucleotide(ll kMer, int length, int direction);

using namespace __gnu_cxx;

namespace __gnu_cxx

{

	template<> struct hash< std::string > {

		size_t operator()( const std::string& x ) const {

         return hash< const char* >()( x.c_str() );

		}

	};

	template<> struct hash< ll > {

		size_t operator()( const ll& x ) const {

         return (x >> 32L) ^ hash< int >()( x & 0xFFFFFFFF );

		}

	};

}


#endif /*COMMON_HPP_*/
