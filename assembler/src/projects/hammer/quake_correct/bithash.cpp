//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bithash.h"
#include "sequence/nucl.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

using namespace::std;

bithash::bithash(int _k)
   :bits()
{
  k = _k;
  assert(_k == kK && "kK and k passed from he programm does not match");
}

bithash::~bithash() {
}

////////////////////////////////////////////////////////////
// add
//
// Add a single sequence to the bitmap
////////////////////////////////////////////////////////////
void bithash::add(Seq<kK> kmer) {
  bits.insert(kmer);
}


////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree
//
// Can handle N's!  Returns False!
////////////////////////////////////////////////////////////
bool bithash::check(unsigned kmer[]) {
  for(int i = 0; i < k; i++) { // ToDo: if we add constructor which
                               // can soft fail if we pass N's in seq
                               // we can optimize this code.
    if (!is_dignucl(kmer[i])) {
      return false;
    }
  }
  return bits.count(Seq<kK>(kmer)) != 0;
}

////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree.
// Pass the kmer map value back by reference to be re-used
//
// Can't handle N's!
////////////////////////////////////////////////////////////
bool bithash::check(unsigned kmer[], Seq<kK> &kmermap) {
  kmermap = Seq<kK>(kmer);
  return bits.count(kmermap) != 0;
}

////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree.
////////////////////////////////////////////////////////////
bool bithash::check(Seq<kK> kmermap) {
  return bits.count(kmermap) != 0;
}

////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree.
// Pass the kmer map value back by reference to be re-used
//
// Can't handle N's!
////////////////////////////////////////////////////////////
bool bithash::check(Seq<kK> &kmermap, unsigned last, unsigned next) {
  kmermap = kmermap << next; 
  // ToDo we can optimize this if Seq will
  // have << operator
  return bits.count(kmermap) != 0;
}

////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the FASTA-file 
////////////////////////////////////////////////////////////
void bithash::hammer_file_load(istream & hammer_in, unsigned long long atgc[]) {
  string line;
  while(getline(hammer_in, line)) {
    if (line[0] != '>') {
      // add to tree
      string kmer = line.substr(0,k);
      add(binary_kmer(kmer));

      // add reverse to tree
      add(binary_rckmer(kmer));

      // count gc
      if(atgc != NULL) {
    unsigned int at = count_at(kmer);
    atgc[0] += at;
    atgc[1] += (k-at);
      }  
    }
  }
}


////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the file given that
// occur >= "boundary" times
////////////////////////////////////////////////////////////
void bithash::meryl_file_load(const char* merf, const double boundary) {
  ifstream mer_in(merf);
  string line;
  double count;
  bool add_kmer = false;

  while(getline(mer_in, line)) {
    if(line[0] == '>') {
      // get count
      count = atof(line.substr(1).c_str());
      //cout << count << endl;
      
      // compare to boundary
      if(count >= boundary) {
    add_kmer = true;
      } else {
    add_kmer = false;
      }

    } else if(add_kmer) {
      // add to tree
      add(binary_kmer(line));

      // add reverse to tree
      add(binary_rckmer(line));
    }
  }
}

////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the file given that
// occur >= "boundary" times
////////////////////////////////////////////////////////////
void bithash::tab_file_load(istream & mer_in, const double boundary, unsigned long long atgc[]) {
  string line;
  double count;

  while(getline(mer_in, line)) {
    if(line[k] != ' ' && line[k] != '\t') {
      cerr << "Kmers are not of expected length " << k << endl;
      exit(EXIT_FAILURE);
    }

    // get count
    count = atof(line.substr(k+1).c_str());
    //cout << count << endl;
      
    // compare to boundary
    if(count >= boundary) {
      // add to tree
      add(binary_kmer(line.substr(0,k)));

      // add reverse to tree
      add(binary_rckmer(line.substr(0,k)));

      // count gc
      if(atgc != NULL) {
    unsigned int at = count_at(line.substr(0,k));
    atgc[0] += at;
    atgc[1] += (k-at);
      }
    }
  }
}

////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the file given that
// occur >= "boundary" times
////////////////////////////////////////////////////////////
void bithash::tab_file_load(istream & mer_in, const vector<double> boundary, unsigned long long atgc[]) {
  string line;
  double count;
  int at;

  while(getline(mer_in, line)) {
    if(line[k] != '\t') {
      cerr << "Kmers are not of expected length " << k << endl;
      exit(EXIT_FAILURE);
    }

    at = count_at(line.substr(0,k));

    // get count
    count = atof(line.substr(k+1).c_str());
    //cout << count << endl;
      
    // compare to boundary
    if(count >= boundary[at]) {
      // add to tree
      add(binary_kmer(line.substr(0,k)));

      // add reverse to tree
      add(binary_rckmer(line.substr(0,k)));

      // count gc
      if(atgc != NULL) {
    unsigned int at = count_at(line.substr(0,k));
    atgc[0] += at;
    atgc[1] += (k-at);
      }
    }
  }
}

////////////////////////////////////////////////////////////
// binary_file_output
//
// Write bithash to file in binary format
////////////////////////////////////////////////////////////
void bithash::binary_file_output(char* outf) {
  /*  unsigned long long mysize = (unsigned long long)bits.size() / 8ULL;
  char* buffer = new char[mysize];
  unsigned int flag = 1;
  for(unsigned long long i = 0; i < mysize; i++) {
    unsigned int temp = 0;
    for(unsigned int j = 0; j < 8; j++) { // read 8 bits from the bitset
      temp <<= 1;
      //unsigned int tmp = i*8 + j;
      //cout << tmp << ",";
      if(bits.count(i*8 + j) != 0)
    temp |= flag;
    }
    buffer[i] = (char)temp;
  }
  ofstream ofs(outf, ios::out | ios::binary);
  ofs.write(buffer, mysize);
  ofs.close();*/
}

////////////////////////////////////////////////////////////
// binary_file_input
//
// Read bithash from file in binary format
////////////////////////////////////////////////////////////
/*
void bithash::binary_file_input(char* inf) {
  ifstream ifs(inf, ios::binary);

  // get size of file
  ifs.seekg(0,ifstream::end);
  unsigned long long mysize = ifs.tellg();
  ifs.seekg(0);

  // allocate memory for file content
  char* buffer = new char[mysize];

  // read content of ifs
  ifs.read (buffer, mysize);

  // parse bits
  unsigned int flag = 128;
  unsigned int temp;
  for(unsigned long i = 0; i < mysize; i++) {
    temp = (unsigned int)buffer[i];
    for(unsigned int j = 0; j < 8; j++) {
      if((temp & flag) == flag)
    bits.set(i*8 + j);
      temp <<= 1;
    }
  }

  delete[] buffer;
}
*/

////////////////////////////////////////////////////////////
// binary_file_input
//
// Read bithash from file in binary format
////////////////////////////////////////////////////////////
void bithash::binary_file_input(char* inf, unsigned long long atgc[]) {
  /*unsigned int flag = 128;
  unsigned int temp;

  ifstream ifs(inf, ios::binary);

  // get size of file
  ifs.seekg(0,ifstream::end);
  unsigned long long mysize = ifs.tellg();
  ifs.seekg(0);

  // allocate memory for file content
  unsigned long long buffersize = 134217728;  // i.e. 4^15 / 8, 16 MB
  if(mysize < buffersize)
       buffersize = mysize;
  char* buffer = new char[buffersize];

  for(unsigned long long b = 0; b < mysize/buffersize; b++) {

    // read content of ifs
    ifs.read (buffer, buffersize);

    // parse bits
    for(unsigned long long i = 0; i < buffersize; i++) {
      temp = (unsigned int)buffer[i];
      for(int j = 0; j < 8; j++) {
    if((temp & flag) == flag) {
      bits.set((buffersize*b + i)*8 + j);
      
      // count gc
      unsigned int at = count_at((buffersize*b + i)*8 + j);
      atgc[0] += at;
      atgc[1] += (k-at);
    }
    temp <<= 1;
      }
    }
  }

  delete[] buffer;*/
}

////////////////////////////////////////////////////////////
// count_at
//
// Count the A's and T's in the sequence given
////////////////////////////////////////////////////////////
int bithash::count_at(string seq) {
  int at = 0;
  for(int i = 0; i < seq.size(); i++)
    if(seq[i] == 'A' || seq[i] == 'T')
      at +=  1;
  return at;
}

int bithash::count_at(Seq<kK> seq) {
  int at = 0;
  unsigned long long mask = 3;
  unsigned long long nt;
  for(int i = 0; i < k; i++) {
    if(seq[i] == 0 || seq[i] == 3)
      at++;
  }
  return at;
}

//  Convert string  s  to its binary equivalent in  mer .
Seq<kK> bithash::binary_kmer(const string & s) {
  return Seq<kK>(s);
}

//  Convert string s to its binary equivalent in mer .
Seq<kK> bithash::binary_rckmer(const string & s) {
  return !Seq<kK>(s); //ToDo: optimize
}

//  Return the binary equivalent of  ch .
unsigned bithash::binary_nt(char ch) {
  switch  (tolower (ch)) {
  case  'a' : return  0;
  case  'c' : return  1;
  case  'g' : return  2;
  case  't' : return  3;
  }
}


unsigned int bithash::num_kmers() {
  return (unsigned int)bits.size();
}
