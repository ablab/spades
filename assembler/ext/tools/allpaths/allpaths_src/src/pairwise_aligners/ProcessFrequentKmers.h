// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef PROCESS_FREQUENT_KMERS
#define PROCESS_FREQUENT_KMERS

#include "Basevector.h"
#include "kmers/KmerRecord.h"
#include "String.h"
#include "Vec.h"

template<int I, int k>
void ProcessFrequentKmers( vec<bvec const*> const& EE,
                           vec< kmer_record<k,I> > R,
                           unsigned int S,
                           int max_clique,
                           String run_dir );

template<int I, int k>
inline void ProcessFrequentKmers( const vecbasevector& EE,
                                  vec< kmer_record<k,I> > R,
                                  unsigned int S,
                                  int max_clique,
                                  String run_dir )
{ vec<bvec const*> vvv; vvv.reserve(EE.size());
  vecbvec::const_iterator end(EE.end());
  for ( vecbvec::const_iterator itr(EE.begin()); itr != end; ++itr )
      vvv.push_back(&*itr);
  ProcessFrequentKmers<I,k>(vvv,R,S,max_clique,run_dir); }

#endif
