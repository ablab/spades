////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef LOOKUP__FIRST_LOOKUP_PERFECT_H
#define LOOKUP__FIRST_LOOKUP_PERFECT_H

#include "Basevector.h"
#include "Vec.h"
#include "lookup/FirstLookAlign.h"
#include "lookup/LookupTab.h"
#include <algorithm>

/**
 * class FirstLookupPerfect.h
 *
 * Adaptated from FirstLookupFinderECJ. It finds all and only the
 * perfect alignments of query sequences onto the LookupTab of a
 * target set of unibases. These are assumed to come from unipaths
 * (generated with kmer size UNIPATH_K), and hence to overlap by
 * UNIPATH_K - 1 bases.
 */
class FirstLookupPerfect 
{
private:
  
  const LookupTab & _lookup_tab;
  const BaseVecVec & _contigs;
  
  static size_t const MIN_SIZE = 96;         // to extend ec reads
  static size_t const MAX_HITS = 2;          // max aligns per query
  static size_t const MAX_FREQ = 10000;      // max attempts per query
  static size_t const CHUNK_SIZE = 1000;     // reads in each worklist

public:
  
  FirstLookupPerfect(const LookupTab & lookup_tab,
                     const BaseVecVec & contigs)
    : _lookup_tab(lookup_tab),
      _contigs(contigs)
  {}
  
  // Find all alignments for a specified query.
  void getAlignments(const BaseVec & query,
                     const size_t query_ID,
                     vec<first_look_align> * result) const;
  
  // Find all alignments for all queries (runs getAlignments on all queries).
  void getAllAlignments(const BaseVecVec & queries,
                        vec<first_look_align> * alignments,
                        const size_t n_threads) const;
  
};

#endif
