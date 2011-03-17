#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>
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
#define RIGHT 1
#define LEFT -1

LOGGER("paireddebruijn.common");

#define MAX_VERT_NUMBER 100000
#define MAX_DEGREE 30
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

const string ini_file = "data/paired.ini";
extern string parsed_reads;
extern string parsed_k_l_mers;
extern string parsed_k_sequence;
extern string error_log;
extern string parsed_l_mers;
extern string graph_file;
extern string graph2;
extern string threaded_graph;

extern int k;
extern int l;
extern int readLength;
extern int maxSeqLength;
void initConstants(string ini_file);
#endif /*COMMON_HPP_*/
