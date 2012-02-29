///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Superb.h"
#include "STLExtensions.h"
#include "lookup/LookAlign.h"
#include "lookup/QueryLookupTableCore.h"
#include "paths/SquashHyperFastavector.h"
#include "system/System.h"
#include "util/RunCommand.h"

/**
 * SquashHyperFastavector
 */
void SquashHyperFastavector( const String &work_dir,
			     const int num_threads,
			     const float max_error_rate,
			     HyperFastavector &hfv,
			     ostream *plog )
{
  // HEURISTIC - this could be made into arguments. Accept an
  // alignment between two scaffolds even if one is not embedded in
  // the other, if:
  //    1. the two scaffolds overlap by at least slack_min_overlap, and
  //    2. no hanging end exceeds slack_max_hang bases
  const int slack_max_hang = 50;
  const int slack_min_overlap = 1000;

  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &log = plog ? *plog : devnull;

  // File names.
  ForceAssert( IsDirectory( work_dir ) );
  String head = work_dir + "/squash";
  String fasta_file = head + ".fasta";
  String fastb_file = head + ".fastb";
  String qlt_file = head + ".qlt";
  
  // Generating superbs.
  log << Date( ) << ": generating supers structure" << endl;
  vec<superb> supers;
  const int MIN_EDGE_TO_SAVE = 1000; // should not be hardcoded here!
  hfv.ConvertToSuperbs( supers, num_threads, MIN_EDGE_TO_SAVE );
  
  // Make a lookup table out of the initial linearized assembly.
  log << Date( ) << ": saving linearized structure" << endl;
  hfv.WriteScaffoldsFasta( fasta_file, supers );
  
  log << Date( ) << ": generating lookup table" << endl;
  String theCommand
    = "MakeLookupTable SOURCE=" + fasta_file
    + " OUT_HEAD=" + head
    + " QUIET=True";
  RunCommandWithLog( theCommand, "/dev/null" );
  
  // Align assembly.
  log << Date( ) << ": running QLT" << endl;
  theCommand
    = "QueryLookupTable L=" + head + ".lookup"
    + " SEQS=" + fastb_file
    + " PARSEABLE=True"
    + " SMITH_WAT=True"
    + " MM=96"
    + " K=12";
  RunCommandWithLog( theCommand, qlt_file );
  
  // Tag supers for removal.
  log << Date( ) << ": tag supers for removal" << endl;
  vec<bool> duplicate( supers.size( ), false );
  {
    vec<look_align_plus> all_hits;
    LoadLookAlignPlus( qlt_file, all_hits );
    for (uint ii=0; ii<all_hits.size( ); ii++) {
      const look_align_plus &hit = all_hits[ii];

      // Discard query == target, and (for unicity) query > target.
      if ( hit.query_id >= hit.target_id ) continue;

      // Not a proper align.
      int qlen = hit.query_length;
      int tlen = hit.target_length;
      bool is_short = ( hit.Pos2( ) - hit.pos2( ) < slack_min_overlap );
      int beg = is_short ? 0 : slack_max_hang;
      int end = is_short ? qlen : qlen - slack_max_hang;
      bool QInT = ( hit.pos1( ) <= beg && hit.Pos1( ) >= end );
      bool TInQ = ( hit.pos2( ) <= beg && hit.Pos2( ) >= end );
      if ( ! ( QInT || TInQ ) ) continue;
      
      // Too many error rates.
      if ( hit.ErrorRate( ) > max_error_rate ) continue;
      
      // Ok.
      int dup_id = ( QInT ) ? hit.query_id : hit.target_id;
      int stays_id = ( QInT ) ? hit.target_id : hit.query_id;
      duplicate[dup_id] = true;
      
      // Log info.
      log << " removing s_" << dup_id << " (it matches s_" << stays_id << ")\n";
    }
  }
  
  // Remove tagged supers.
  log << Date( ) << ": removing edges of tagged supers\n" << endl;
  vec<int> to_delete;
  for (uint super_id=0; super_id<supers.size( ); super_id++) {
    if ( ! duplicate[super_id] ) continue;
    const superb &sup = supers[super_id];
    for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
      log << " removing e_" << sup.Tig( cgpos ) << " from fastavector\n";
      to_delete.push_back( sup.Tig( cgpos ) );
    }
  }
  hfv.DeleteEdges( to_delete );
  hfv.RemoveDeadEdgeObjects( );
  
}

