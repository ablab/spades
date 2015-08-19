//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef BITHASH_H
#define BITHASH_H

#include <string>
#include <vector>
#include <cmath>
#include <boost/unordered_set.hpp>
using namespace::std;

class bithash {
 public:
  bithash(int _k);
  ~bithash();
  void add(long long unsigned kmer);
  bool check(unsigned kmer[]);
  bool check(unsigned kmer[], long long unsigned & kmermap);
  bool check(long long unsigned & kmermap, unsigned last, unsigned next);
  bool check(long long unsigned kmermap);
  void meryl_file_load(const char* merf, const double boundary);
  void tab_file_load(istream & mer_in, const double boundary, unsigned long long atgc[]);
  void tab_file_load(istream & mer_in, const vector<double> boundary, unsigned long long atgc[]);
  void hammer_file_load(istream & hammer_in, unsigned long long atgc[]);
  long long unsigned binary_kmer(const string &s);
  long long unsigned binary_rckmer(const string &s);
  void binary_file_output(char* outf);

  void binary_file_input(char* inf, unsigned long long atgc[]);
  unsigned int num_kmers();

  static int k;
 private:  
  unsigned binary_nt(char ch);
  int count_at(string seq);
  int count_at(unsigned long long seq);

  boost::unordered_set<unsigned long long> bits;
  unsigned long long mask;
};

#endif
