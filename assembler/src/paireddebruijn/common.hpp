#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <math.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <ext/hash_map>
#include <algorithm>
#include <string>
#include <set>
#include "logging.hpp"

#define forn(i, n) for(int i = 0; i < (int) n; i++)
#define ll long long
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define edgesMap  map<ll, vector<VertexPrototype *> >
#define verticesMap  map<ll, vector<VertexPrototype *> >
#define longEdgesMap  map<int, Edge*>

LOGGER("paireddebruijn.common");

#define MAX_VERT_NUMBER 100000
#define MAX_DEGREE 30

using namespace std;
const string parsed_reads = "data/reads_const_d.txt";
//const string parsed_reads = "data/filtered_reads";
const string parsed_k_l_mers = "data/klmers_const_d.txt";
const string parsed_k_sequence = "data/vertices_const_d.txt";
const string error_log = "data/error.log";
const string parsed_l_mers = "data/lmers_const_d.txt";
const string graph_file = "data/graph_const_d.dot";
const string threaded_graph = "data/threaded_graph_const_d.dot";
const string auxilary_lmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

const int k = 25;
const int l = 31;
const int readLength = 100;
const int maxSeqLength = 200;
#endif /*COMMON_HPP_*/
string decompress(ll a, int l);
