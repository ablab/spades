///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SHORT_LOCS_HANDLER_H
#define SHORT_LOCS_HANDLER_H

#include "ReadLocationLG.h"
#include "String.h"
#include "paths/KmerPath.h"
#include "Intvector.h"

/**
 * class LocsHandlerLG
 *
 * Generate indices for ReadLocationLG's.
 */
class LocsHandlerLG {

public:

  // Empty constructor.
  LocsHandlerLG( );

  // Constructor that loads the data.
  LocsHandlerLG( const int K, const String &reads_head );
  
  // Constructor that is passed the data as pointers.
  LocsHandlerLG( const vec<int> *rlens,
	     const vec<int> *clens,
	     const vec<int> *to_rc,
	     const vec<ReadLocationLG> *locs );
  
  // Load from KmerPaths and unilocs.
  void LoadFromKmers( const int K, const String &reads_head );
  
  // Bracket on read for given loc.
  pair<int,int> BracketOnRead( longlong loc_id ) const;
  
  // Bracket on contig for given loc.
  pair<int,int> BracketOnContig( longlong loc_id ) const;
  
  // Return true if chain is valid. See .cc for details.
  bool GetFwChain( longlong read_id, vec<longlong> &locpos ) const;

  // Print given chain.
  void PrintChain( vec<longlong> &locpos, ostream &out ) const;

  // Returns all possible placements (as loc ids) for read_id.
  vec<longlong> GetAllPlacements( longlong read_id ) const;
  
  // Returns all fw placements (as loc ids) for read_id.
  vec<longlong> GetFwPlacements( longlong read_id ) const;

  // Returns all rc placements (as loc ids) for read_id.
  vec<longlong> GetRcPlacements( longlong read_id ) const;
  
  // Number of contigs.
  longlong NContigs( ) const { return clens_->size( ); }

  // Number of reads.
  longlong NReads( ) const { return rlens_->size( ); }

  // Length of contig.
  int ContigLength( int contig_id ) { return (*clens_)[contig_id]; }

  // Length of read.
  int ReadLength( longlong read_id ) { return (*rlens_)[read_id]; }

  // Find rc of this contig.
  int ToRc( int contig_id ) const { return (*to_rc_)[contig_id]; }
  
  // For backward compatibility: returns vector with all locs.
  const vec<ReadLocationLG> &Locs( ) const { return (*locs_); }
  
  // Returns location at ii.
  const ReadLocationLG &operator[]( longlong ii ) const { return (*locs_)[ii]; }
  
  
private:
  
  // Set indices (to_loc_).
  void SetIndices( );
  
  
private:
  
  // Pointers to core data.
  const vec<int> *rlens_;                 // lengths of reads (in kmers)
  const vec<int> *clens_;                 // lengths of contigs (in kmers)
  const vec<int> *to_rc_;                 // which contig is rc of contig_id
  const vec<ReadLocationLG> *locs_;  // locs (sorted)

  // These may be empty (filled only if data are loaded rather than passed in).
  vec<int> core_rlens_;
  vec<int> core_clens_;
  vec<int> core_to_rc_;
  vec<ReadLocationLG> core_locs_;
  
  // Various maps.
  Int64VecVec to_loc_;       // all placements for given read_id

};

#endif
