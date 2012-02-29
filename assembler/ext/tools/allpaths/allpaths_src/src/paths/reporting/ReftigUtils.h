///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REFTIG_UTILS_H
#define REFTIG_UTILS_H

#include "Bitvector.h"
#include "String.h"
#include "graph/Digraph.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "lookup/LookAlign.h"

/**
 * GetAlignsFast
 *
 * WARNING! This code is parallelized with omp.
 *
 * USE_CACHE: if true, load cached aligns (if found)
 */
void GetAlignsFast( const int K, const String& fastb_file,
		    const String& lookup_file, const String& aligns_file,
		    vec<look_align>& aligns, const Bool USE_CACHE,
		    const String& tmp_dir);

/**
 * HyperToReftigsCore
 *
 * Find reftigs by generating an acyclic graph in which vertices are
 * alignments, and edges correspond do adjacent alignments. Store the
 * alignment graph in alignsG, if the pointer is not NULL.
 */
template<class HYPER_T>
void HyperToReftigsCore( const int K, const HYPER_T& h,
			 const vec<look_align>& aligns,
			 vec< pair<int,ho_interval> >& reftigs,
			 digraph *alignsG = 0 );

/**
 * GenerateDot
 *
 * Generate dot file. Output file saved as <dot_base>.dot
 *
 * VERBOSITY can be 0, 1, or 2, and it affects how edges are
 * decorated:
 *     0: show base alpha (as in PrettyDOT's edge_labels_base_alpha)
 *     1: specify wich reftig this edge belongs to
 *     2: use integer ids, and show interval of alignment for this edge
 */
template<class HYPER_T>
void GenerateDot( const int ORIGIN, const String &dot_base,
		  const HYPER_T &hyper, const digraph &agraph,
		  const vec<look_align> &aligns,
		  const vec< pair<int,ho_interval> >& reftigs,
		  int VERBOSITY = 0 );
  
/**
 * PrintReftigs
 *
 * Print given reftigs.
 *
 * K: the kmer used to align (as in HyperToReftigsCore above)
 * ORIGIN: offset the start point of the reftigs
 * out: where output is sent
 * amb: if given, use it to disply gaps
 * min_len: report only intervals >= min_len kmers
 */
void PrintReftigs ( ostream &out,
		    const int K,
		    const int ORIGIN,
		    const vec< pair<int,ho_interval> > &reftigs,
		    const vecbitvector *amb = 0,
		    const int min_len = 0 );

#endif
