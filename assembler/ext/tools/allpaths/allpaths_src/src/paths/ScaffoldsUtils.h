///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SCAFFOLDS_UTILS_H
#define SCAFFOLDS_UTILS_H

#include "Basevector.h"
#include "Fastavector.h"
#include "PairsManager.h"
#include "paths/Alignlet.h"
#include "paths/Alignlets2ReadLocs.h"
#include "util/SearchFastb2Core.h"

/**
 * AlignReadsToContigs
 *
 * Align reads to contigs. Output is saved in aligns, index.
 */
void AlignReadsToContigs( const int K,
			  const String &tmp_dir,
			  const String &reads_fastb_file,
			  const vec<fastavector> &contigs,
			  vec<alignlet> &aligns,
			  vec<int> &index,
			  ostream &log,
			  bool VERBOSE = false );

/**
 * ReportScaffoldsN50
 *
 * Print N-stats for contigs and scaffolds.
 */
void ReportScaffoldsN50( const vec<superb> &supers,
			 ostream &out );

/**
 * ReportScaffoldsBrief
 *
 * More compact version (one liner)
 */
void ReportScaffoldsBrief( const vec<superb> &supers,
			   const int min_links,
			   const int step,
			   ostream &out );

/**
 * SaveInterimScaffolds
 *
 * Used to help debugging and testing the iterative MakeScaffoldsLG
 * (it dumps a lot of data). It also runs EvalScaffolds, if the proper
 * genome.lookup file is found in data_dir.
 */
void SaveInterimScaffolds( const String &data_dir,
			   const String &out_dir,
			   const PairsManager &pairs,
			   const vec<fastavector> &contigs,
			   const vec<superb> &supers,
			   const vec<alignlet> *aligns = 0,
			   const vec<int> *index = 0);

/**
 * UpdateIndexFile
 *
 * Update an index file after the changes to a filtered_index file.
 * The two files are assumed to be connected, where the filtered index
 * was generated after resetting some entries in index to a <0 value
 * (by RemoveHighCNAligns, for example). Run UpdateIndexFile if some
 * reads have been removed in the original set, to reflect the change
 * in the filtered set. It returns the number of updates.
 */
size_t UpdateIndexFile( const vec<int> &orig_index, vec<int> &filt_index );

/**
 * CheckScaffoldsIntegrity
 *
 * Make sure all contigs are accounted for (each contig appears in
 * exactly one super).
 */
bool CheckScaffoldsIntegrity( const size_t n_contigs,
			      const vec<superb> &supers );

#endif
