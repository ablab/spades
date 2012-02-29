/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef CONTIGS_MANAGER_H
#define CONTIGS_MANAGER_H

#include "MainTools.h"
#include "Alignment.h"
#include "Fastavector.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"
#include "paths/Alignlet.h"

/**
 * class ContigsManager
 *
 * Class that allows merging and breaking of contigs, with bookkeeping
 * for existing alignments of reads. Warning: read alignments are
 * removed, by resetting their index to -1.
 */
class ContigsManager  {
  
public:
  
  ContigsManager( vec<fastavector> &contigs,
		  vec<alignlet> &aligns0, vec<int> &aligns0_index, 
		  vec<alignlet> &ualigns0, vec< vec<int> > &ualigns0_index);

  const vec<fastavector> &Contigs( ) const { return contigs_; }
  
  // Split contig cg_id at the given pos. It returns the ids of all chunks.
  vec<size_t> SplitContig( size_t cg_id, size_t pos1, size_t pos2 );
  vec<size_t> SplitContig( size_t cg_id, size_t pos );
  
  // Split contig in such a way that resulting contigs overlap in the region [pos1,pos2]
  vec<size_t> SlideSplitContig( size_t cg_id, size_t pos1, size_t pos2 );

  // Cut away bases [0, pos) from cg_id.
  void CutHead( size_t cg_id, size_t pos );

  // Cut away bases [pos, end_of_contig) from cg_id.
  void CutTail( size_t cg_id, size_t pos );

  // Ids of contigs (sig_cg1 and sig_cg2) are signed, to capture orientation.
  void MergeContigs( int sig_cg1, int sig_cg2, const alignment &al );
  
  // Remove all aligns on this contig (by resetting index to -1).
  void IsolateContig( size_t cg_id );

  
private:
  
  // Generate the map cg2reads_.
  void Setup( );
  
  // Reverse complement the given contig, and update aligns accordingly.
  void ReverseComplement( size_t cg_id );
  
  // Romove (ie reset index to -1) all aligns overlapping region in contig.
  void RemoveAligns( size_t cg_id, const ho_interval *win = 0 );
  
  
private:

  vec<fastavector> &contigs_;  // contigs
  vec<alignlet>& aligns0_; 
  vec<int>& aligns0_index_;
  vec<alignlet>& ualigns0_; 
  vec< vec<int> >& ualigns0_index_;
 
  vec< vec<int> > cg2seqs_;   // a map contig id to ids of reads on contig
  vec< vec<int> > cg2useqs_; 
  
};

#endif
