//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef BITHASH_H
#define BITHASH_H

const int kK = 55; //k-mer length
#include <string>
#include <vector>
#include <cmath>
#include <boost/unordered_set.hpp>
#include "sequence/seq.hpp"
using namespace::std;

class bithash {
 public:
  bithash(int _k);
  ~bithash();
  void add(Seq<kK> kmer);
  bool check(unsigned kmer[]);
  bool check(unsigned kmer[], Seq<kK> &kmermap);
  bool check(Seq<kK> &kmermap, unsigned last, unsigned next);
  bool check(Seq<kK> kmermap);
  void meryl_file_load(const char* merf, const double boundary);
  void tab_file_load(istream & mer_in, const double boundary, unsigned long long atgc[]);
  void tab_file_load(istream & mer_in, const vector<double> boundary, unsigned long long atgc[]);
  void hammer_file_load(istream & hammer_in, unsigned long long atgc[]);
  Seq<kK> binary_kmer(const string &s);
  Seq<kK> binary_rckmer(const string &s);
  void binary_file_output(char* outf);

  void binary_file_input(char* inf, unsigned long long atgc[]);
  unsigned int num_kmers();

  static int k;
 private:  
  unsigned binary_nt(char ch);
  int count_at(string seq);
  int count_at(Seq<kK> seq);

  boost::unordered_set< Seq<kK>, Seq<kK>::hash > bits;
};

#endif
