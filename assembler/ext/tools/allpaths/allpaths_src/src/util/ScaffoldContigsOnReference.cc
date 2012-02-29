/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "FetchReads.h"
#include "lookup/LookAlign.h"
#include "paths/SaveScaffoldGraph.h"
#include "SupersHandler.h"
#include "util/RunCommand.h"
// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency EvalScaffolds

/**
 * ScaffoldContigsOnReference
 *
 * Align contigs to a reference, and generate scaffolds.
 *
 * HEAD_REF: head of reference genome
 * FORCE: do not use cached aligns
 * */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( ASSEMBLY_IN );
  CommandArgument_String( ASSEMBLY_OUT );
  CommandArgument_String( HEAD_REF );
  CommandArgument_Int_OrDefault( MAX_GAP, 2000 );
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 1000 );
  CommandArgument_Int_OrDefault( MIN_GAP_DEV, 20 );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir and file names.
  String tmp_dir = ASSEMBLY_OUT + ".tmp";

  String contigs_file = ASSEMBLY_IN + ".contigs.fasta";
  String supers_file = ASSEMBLY_IN + ".superb";
  String lookup_file = HEAD_REF + ".lookup";
  String aligns_file = tmp_dir + "/aligns.qlt";

  Mkpath( tmp_dir );

  // Load aligns.
  vec<look_align_plus> aligns;
  if ( FORCE || ! IsRegularFile( aligns_file ) ) {
    String comm
      = "EvalScaffolds LOOKUP=" + lookup_file
      + " SCAFFOLDS=" + ASSEMBLY_IN
      + " OUT_DIR=" + tmp_dir;
    RunCommand( comm );
  }
  LoadLookAlignPlus( aligns_file, aligns );
  
  // Load contigs.
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  cout << Date( ) << ": loading supers" << endl;
  shandler initial_supers( -1, supers_file );

  // Flag which contigs are actually in the inital supers
  vec<bool> initial( contigs.size( ), false );
  for (int s = 0; s < initial_supers.Size(); ++s) {
    for (int t = 0; t < initial_supers[s].Ntigs(); ++t) {
      initial[initial_supers[s].Tig(t)] = True;
    }      
  }

  // Place contigs on reference (rc if needed).
  vec<bool> placed( contigs.size( ), false );
  vec<seq_interval> placs;
  for (size_t ii=0; ii<aligns.size( ); ii++) {
    const look_align_plus &al = aligns[ii];
    int cid = al.query_id;
    int tid = al.target_id;
    int clen = al.query_length;
    int tlen = al.target_length;

    // Skip contigs we don't care about
    if ( !initial[cid] ) continue;

    // Skip already placed contigs.
    if ( placed[cid] ) continue;

    // Only embedded contigs are allowed.
    if ( al.a.pos1( ) != 0 || al.a.Pos1( ) != clen ) continue;

    // This should not happen.
    if ( al.a.pos2( ) < 0 || al.a.Pos2( ) > tlen ) continue;
    
    // Place contig (eventually rc it).
    placed[cid] = true;
    placs.push_back( seq_interval( cid, tid, al.a.pos2( ), al.a.Pos2( ) ) );
    if ( al.Rc1( ) ) contigs[cid].ReverseComplement( );
  }

  sort( placs.begin( ), placs.end( ) );
  
  // Create scaffolds (placed contigs).
  cout << Date( ) << ": Create scaffolds from alignments" << endl;
  vec<superb> supers;
  for (size_t ii=0; ii<placs.size( ); ii++) {
    int cid = placs[ii].IntervalId( );
    int clen = contigs[cid].size( );

    // New super.
    if ( ii == 0 || placs[ii].SeqId( ) != placs[ii-1].SeqId( )) {
      superb sup;
      sup.PlaceFirstTig( cid, clen );
      supers.push_back( sup );
    } else {
      int gap = placs[ii].Begin( ) - placs[ii-1].End( );
      int dev = Max( MIN_GAP_DEV, gap / 4 );
      if (gap >= -MAX_OVERLAP && gap <= MAX_GAP) {
	// Append to existing super.
	supers[supers.size( )-1].AppendTig( cid, clen, gap, dev );
      } else {
	// New super.
	superb sup;
	sup.PlaceFirstTig( cid, clen );
	supers.push_back( sup );
      }
    }
  }

  cout << Date( ) << ": Add unplaced contigs" << endl;
  // Add unplaced contigs.
  for (int cid=0; cid<(int)contigs.size( ); cid++) {
    if ( !initial[cid] ) continue;
    if ( placed[cid]  ) continue;

    int clen = contigs[cid].size( );
    superb sup;
    sup.PlaceFirstTig( cid, clen );
    supers.push_back( sup );    
  }


  // Save.
  SaveScaffoldAssembly( ASSEMBLY_OUT, supers, contigs, &cout, true );
  
  // Done.
  cout << Date( ) << ": done" << endl;

}
