///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_LIN_FASTAVEC_H
#define C_LIN_FASTAVEC_H

#include "Fastavector.h"
#include "Vec.h"
#include "paths/HyperFastavector.h"

/**
 * CLinFastavec
 *
 * Class with tools to incrementally linearize (flatten) a given
 * HyperFastavector.
 */
class CLinFastavec {
  
public:
  
  CLinFastavec( HyperFastavector &hyper, ostream *plog = 0 );

  void SetNewAlgorithm( ) { new_algorithm_ = True; }

  Bool NewAlgorithm( ) const { return new_algorithm_; }
  
  // Save all intermediate steps as base_name.<iter_>
  void SetBaseInterSave( const String &base_name );

  // Find and collapse cells (see graph/FindCells.h for details).
  int FlattenCells( const int MAX_CELL_SIZE );
  
  // Flatten frayed ends (multiple sources or sinks).
  int FlattenFrayedEnds( );
  
  // Flatten bubbles.
  int FlattenBubbles( );

  // Report info on existing bubbles.
  void ReportBubbles( );
  
  
private:

  // Save intermediate assembly.
  void SaveIntermediate( ) const;

  // Return fastavector of flattened region (empty if failed).
  fastavector ConsensusPath( vec< vec<int> > &paths, bool same_size ) const;

  // Similar to Combine in Fastavector, but it allows for some slack (see .cc).
  fastavector CombineSlack( const vec<fastavector> &fvec ) const;
  
  // Find all bubbles (pairs of vertices with exactly two edges between them).
  void FindBubbles( vec< pair<int,int> > &bubbles) const;
  
  // Clean up graph, and increment iter_ by 1.
  void CleanUp( bool alt = false );
  
  
private:

  HyperFastavector &hyper_;   // HyperFastavector
  String base_inter_save_;    // base name for intermediate saves (may be empty)
  ostream *log_;              // log stream (may be null)
  int iter_;                  // iteration counter
  Bool new_algorithm_;        // use new algorithm?
  
};

#endif
