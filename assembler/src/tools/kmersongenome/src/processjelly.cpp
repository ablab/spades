//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//============================================================================
// Name        : processjelly.cpp
// Author      : Sergey Nikolenko
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include "libtrie/Trie.hxx"

using namespace std;

typedef ToolBox::Trie<vector<long> > TTrie;

const char nucl_map[4] = {'A', 'C', 'G', 'T'};

const char nucl_complement_map['T' + 1] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A'};
inline char nucl_complement(signed char c){    return nucl_complement_map[c]; }

inline const std::string ReverseComplement(const std::string &s) {
    std::string res(s.size(), 0);
    transform(s.begin(), s.end(), res.rbegin(), nucl_complement);
    return res;
}

string getKmer(const string & g, size_t i, size_t k) {
    string s = g.substr(i, k);
    string r = ReverseComplement(s);
    if ( s.compare(r) > 0 ) return r;
    else return s;
}

const size_t k_num = 11;
//const size_t kfilter[k_num] = {17, 19, 21, 25, 29, 33, 37, 45, 55, 65, 75};
const size_t kfilter[k_num] = {13, 15, 17, 19, 21, 25, 29, 33};

void read_jelly_dump( const char * fname, size_t k, const TTrie & t, vector< size_t > & v ) {
    FILE * pFile;
    char buf[k+1];
    int count = 0;
    pFile = fopen (fname,"r");
    while (!feof(pFile)) {
        if (!fscanf(pFile, "%s", buf) || !fscanf(pFile, "%i", &count)) {
            cout << "Error reading dump!\n";
        }
        if (feof(pFile)) break;
        const vector<long> l = t.getEntry(buf, k);
        for ( size_t i=0; i<l.size(); ++i ) {
            v[l[i]] = count;
        }
    }
    fclose(pFile);
}

void read_genome( const char * fname, string & s ) {
    s = "";
    FILE * pFile;
    char buf[1024];
    pFile = fopen (fname,"r");
    while (!feof(pFile)) {
        if (!fscanf(pFile, "%s", buf)) {
            cout << "Error reading dump!\n";
        }
        if (feof(pFile)) break;
        if (!( (buf[0] == 'A') || (buf[0] == 'C') || (buf[0] == 'G') || (buf[0] == 'T') )) continue;
        s.append(buf);
    }
    fclose(pFile);
}

void read_genome_trie( const string & g, size_t k, TTrie & t ) {
    for ( size_t i=0; i < g.size() - k; ++i ) {
        string s = getKmer(g, i, k);
        vector<long> l = t.getEntry(s.c_str(), k);
        l.push_back(i);
        t.setEntry(s.c_str(), k, l);
    }
}

void write_coverage( const char * fname, vector<size_t> v ) {
    FILE * pFile = fopen (fname, "w");
    for ( size_t i=0; i < v.size(); ++i ) {
        fprintf(pFile, "(%li, %li)\n", i, v[i]);
    }
    fclose(pFile);
}

int main(int argc, char* argv[]) {
    size_t k_min = 17;
    size_t k_max = 19;
    int nthreads = 1;
    string pref = "";
    string outpref = "res";
    string genomefname = "../MG1655-K12.fasta";
    for(int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-kmin")) {
            k_min = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i], "-kmax")) {
            k_max = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i], "-nthreads")) {
            nthreads = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i], "-dump")) {
            pref = argv[i+1];
        }
        if (!strcmp(argv[i], "-out")) {
            outpref = argv[i+1];
        }
        if (!strcmp(argv[i], "-genome")) {
            genomefname = argv[i+1];
        }
    }

    cout << "Reading genome..." << endl;
    string g;
    read_genome( genomefname.c_str(), g );
    cout << "  ...genome size: " << g.size() << endl;

    //size_t i = 1603380;
    //cout << g.substr(1603380, 29) << endl << ReverseComplement(g.substr(1603380, 29)) << endl;
    //return 0;

    #pragma omp parallel for shared(g) num_threads(nthreads)
    for ( size_t j = 0; j <= k_num; ++j ) {
        size_t k = kfilter[j]; if ( (k < k_min) || (k > k_max) ) continue;
        cout << endl;
        TTrie trie(vector<long>(0));
        cout << "Creating genome trie for k=" << k << "..."; flush(cout);
        read_genome_trie( g, k, trie );
        cout << "  OK!" << endl;

        vector< size_t > cov(g.size(), 0);
        ostringstream dfname; dfname << pref << "/" << k << "/dump.tsv";
        cout << "Reading jelly dump for k=" << k << "..."; flush(cout);
        cout << dfname.str() << endl;
        read_jelly_dump( dfname.str().c_str(), k, trie, cov );
        cout << " OK!" << endl;

        size_t sum = accumulate(cov.begin(), cov.end(), 0);
        cout << "Sum of genome kmer counts for k=" << k << " is " << sum << endl;

        ostringstream resfname; resfname << outpref << "." << k << ".tex";
        cout << "Writing coverage for k=" << k << " to " << resfname.str() << "..."; flush(cout);
        write_coverage( resfname.str().c_str(), cov );
        cout << " OK!" << endl;
    }

    return 0;
}
