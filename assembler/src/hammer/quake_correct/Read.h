//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef READ_H
#define READ_H

#include "bithash.h"
#include <string>
#include <vector>
#include <fstream>
#include <bitset>

using namespace::std;

//const int bitsize = 22;
const int bitsize = 200;

////////////////////////////////////////////////////////////
// correction
//
// Simple structure for corrections to reads
////////////////////////////////////////////////////////////
class correction {
public:
  correction(short i, short t) {
    index = i;
    to = t;
  };
  // default copy constructor should suffice

  short index;
  //int from;
  short to;
};

////////////////////////////////////////////////////////////
// corrected_read
//
// Simple structure for corrected reads
////////////////////////////////////////////////////////////
class corrected_read {
public:

 corrected_read(vector<correction> & c, bitset<bitsize> & u, float l, short re)
    :untrusted(u) {
    likelihood = l;
    region_edits = re;
    for(int i = 0; i < c.size(); i++)
      corrections.push_back(correction(c[i]));
  };
 corrected_read(bitset<bitsize> & u, float l, short re)
    :untrusted(u) {
    likelihood = l;
    region_edits = re;
  };
  ~corrected_read() {
    /*
    while(corrections.size() > 0) {
      delete corrections.back();
      corrections.pop_back();
    }
    */
  }
  // default destructor should call the correction vector
  // destructor which should call the correction destructor

  vector<correction> corrections; 
  bitset<bitsize> untrusted; // inaccurate until pop'd off queue and processed
  float likelihood;
  short region_edits;
};

////////////////////////////////////////////////////////////
// Read
////////////////////////////////////////////////////////////
class Read {
 public:
  Read(const string & h, const unsigned int* s, const string & q, vector<int> & u, const int read_length);
  ~Read();

  string trim(int t);
  string correct(bithash *trusted, double ntnt_prob[][4][4], double prior_prob[4], bool learning = false);
  int correct_cc(vector<short>, vector<int> untrusted_subset, bithash* trusted, double ntnt_prob[][4][4], double prior_prob[4], bool learning);
  vector<short> error_region(vector<int> untrusted_subset);
  vector<short> error_region_chop(vector<int> untrusted_subset);
  bool check_trust(corrected_read *cr, bithash *trusted, unsigned int & check_count);
  string print_seq();
  string print_corrected(vector<correction> & cor);
  string print_corrected(vector<correction> & cor, int print_nt);


  string header;
  int read_length;
  int trim_length;
  unsigned int* seq;
  unsigned int* quals;
  float* prob;
  vector<int> untrusted;
  corrected_read *trusted_read;

  const static float trust_spread_t;
  const static float correct_min_t;
  const static float learning_min_t;
  const static unsigned int max_qual = 60;
  static int quality_scale;

 private:
  bool untrusted_intersect(vector<int> untrusted_subset, vector<short> & region);
  void untrusted_union(vector<int> untrusted_subset, vector<short> & region);
  void quality_quicksort(vector<short> & indexes, int left, int right);

  float global_like;  // to track likelihood across components

  const static bool aggressive = true;
  const static short expand_region = 1;
};

#endif
