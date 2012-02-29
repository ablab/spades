///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/ContigsManager.h"

/**
 * ContigsManager
 * Constructor
 */
ContigsManager::ContigsManager( vec<fastavector> &contigs,
		  vec<alignlet> &aligns0, vec<int> &aligns0_index, 
		  vec<alignlet> &ualigns0, vec< vec<int> > &ualigns0_index ) :
  contigs_ ( contigs ),
  aligns0_ ( aligns0 ),
  aligns0_index_ ( aligns0_index ),
  ualigns0_ ( ualigns0 ),
  ualigns0_index_ ( ualigns0_index )
{
  this->Setup( );
}



/**
 * ContigsManager
 * SplitContig
 */
vec<size_t> ContigsManager::SplitContig( size_t cg_id,
					 size_t pos1,
					 size_t pos2 )
{
  cout << "SplitContig double" << endl;
  PRINT3( cg_id, pos1, pos2 );
  vec<size_t> new_ids = this->SplitContig( cg_id, pos1 );
  new_ids.push_back( this->SplitContig( new_ids.back( ), pos2 - pos1 )[1] );
  
  return new_ids;
}

/**
 * ContigsManager
 * SplitContig
 */
vec<size_t> ContigsManager::SplitContig( size_t cg_id, size_t pos )
{
  cout << "SplitContig" << endl;
  vec<size_t> new_ids;
  
  // Chunks from original contig.
  const fastavector &original = contigs_[cg_id];
  const size_t original_size = original.size( );
  ForceAssertLt( pos, original_size );
  PRINT3( cg_id, original_size, pos );

  fastavector tg0; 
  fastavector tg1; 
  tg0.SetToSubOf( original, 0, pos );
  tg1.SetToSubOf( original, pos, original_size - pos );
  int tg0len = tg0.size( );
  int tg1len = tg1.size( );
  
  // Update aligns_, index_, and cg2_seqs_.
  
  size_t new_cg_id = contigs_.size( );
  

  cg2seqs_.resize( 1 + contigs_.size( ) );
  // Make a local copy of the read set for the contig, and clear current
  vec<int> seq_set = cg2seqs_[cg_id];
  // Rebuild cg2seqs_ for cg_id and for the last contig (the new one).
  cg2seqs_[cg_id].clear( );
  for (int ii=0; ii<seq_set.isize( ); ii++) {
    int seq_id = seq_set[ii];
    if ( aligns0_index_[seq_id] < 0 ) continue;
    alignlet& al = aligns0_[ aligns0_index_[seq_id] ];
    
    // Read on the left half.
    if ( al.Pos2( ) < (int)pos ) {
      cg2seqs_[cg_id].push_back( seq_id );
      alignlet new_al( al.pos2( ), al.Pos2( ), cg_id, tg0len, al.Fw1( ) );
      aligns0_[ aligns0_index_[seq_id] ] = new_al;
      continue;
    }
    // Read on the right half.
    else if ( al.pos2( ) >= (int)pos ) {
      cg2seqs_[new_cg_id].push_back( seq_id );
      int pos2 = al.pos2( ) - pos;
      int Pos2 = al.Pos2( ) - pos;
      alignlet new_al( pos2, Pos2, new_cg_id, tg1len, al.Fw1( ) );
      aligns0_[ aligns0_index_[seq_id] ] = new_al;
      continue;
    }else{
      // A straddler, remove read (by resetting index_ to -1).
      aligns0_index_[seq_id] = -1;
    }
  }

  
  cg2useqs_.resize( 1 + contigs_.size( ) );
  // Make a local copy of the read set for the contig, and clear current
  vec<int> useq_set = cg2useqs_[cg_id];
  // Rebuild cg2seqs_ for cg_id and for the last contig (the new one).
  cg2useqs_[cg_id].clear( );
  for (int ii=0; ii<useq_set.isize( ); ii++) {
    int useq_id = useq_set[ii];
    for ( size_t ai = 0; ai < ualigns0_index_[useq_id].size(); ai++ ){
      if ( ualigns0_index_[useq_id][ai] < 0 ) continue;
      alignlet& al = ualigns0_[ ualigns0_index_[useq_id][ai] ];
      if ( al.TargetId() != (int)cg_id ) continue;

      // Read on the left half.
      if ( al.Pos2( ) < (int)pos ) {
	cg2useqs_[cg_id].push_back( useq_id );
	alignlet new_al( al.pos2( ), al.Pos2( ), cg_id, tg0len, al.Fw1( ) );
	ualigns0_[ ualigns0_index_[useq_id][ai] ] = new_al;
	continue;
      }
      // Read on the right half.
      else if ( al.pos2( ) >= (int)pos ) {
	cg2useqs_[new_cg_id].push_back( useq_id );
	int pos2 = al.pos2( ) - pos;
	int Pos2 = al.Pos2( ) - pos;
	alignlet new_al( pos2, Pos2, new_cg_id, tg1len, al.Fw1( ) );
	ualigns0_[ ualigns0_index_[useq_id][ai] ] = new_al;
	continue;
      }else{
	// A straddler, remove read (by resetting index_ to -1).
	ualigns0_index_[useq_id][ai] = -1;
      }
    }
  }



  // Update contigs_.
  contigs_[cg_id] = tg0;
  new_ids.push_back( cg_id );
  
  contigs_.push_back( tg1 );
  new_ids.push_back( new_cg_id );
  
  return new_ids;
}

/**
 * ContigsManager
 * CutHead
 */
void ContigsManager::CutHead( size_t cg_id, size_t pos )
{
  cout << "CutHead" << endl;
  // Chunks from original contig.
  const fastavector &original = contigs_[cg_id];
  const size_t original_size = original.size( );
  ForceAssertLt( pos, original_size );
  PRINT3( cg_id, original_size, pos );

  fastavector tg0;
  tg0.SetToSubOf( original, pos, original_size - pos );
  int tg0len = tg0.size( );
  
  // Update aligns_, index_, and cg2_seqs_.
  vec<int> seq_set = cg2seqs_[cg_id];
  cg2seqs_[cg_id].clear( );
  for (int ii=0; ii<seq_set.isize( ); ii++) {
    int seq_id = seq_set[ii];
    if ( aligns0_index_[seq_id] < 0 ) continue;
    alignlet& al = aligns0_[ aligns0_index_[seq_id] ];
    
    // Read on the left half.
    if ( al.pos2( ) < (int)pos ) {
      aligns0_index_[seq_id] = -1;
      continue;
    }
    
    // Read on the right half.
    else {
      cg2seqs_[cg_id].push_back( seq_id );
      int pos2 = al.pos2( ) - pos;
      int Pos2 = al.Pos2( ) - pos;
      alignlet new_al( pos2, Pos2, cg_id, tg0len, al.Fw1( ) );
      aligns0_[ aligns0_index_[seq_id] ] = new_al;
      continue;
    }  
  }
  
  // Update aligns_, index_, and cg2_seqs_.
  vec<int> useq_set = cg2useqs_[cg_id];
  cg2useqs_[cg_id].clear( );
  for (int ii=0; ii<useq_set.isize( ); ii++) {
    int useq_id = useq_set[ii];
    for ( size_t ai = 0; ai < ualigns0_index_[useq_id].size(); ai++ ){
      if ( ualigns0_index_[useq_id][ai] < 0 ) continue;
      alignlet& al = ualigns0_[ ualigns0_index_[useq_id][ai] ];
      if ( al.TargetId() != (int)cg_id ) continue; 

      // Read on the left half.
      if ( al.pos2( ) < (int)pos ) {
	ualigns0_index_[useq_id][ai] = -1;
	continue;
      }

      // Read on the right half.
      else {
	cg2useqs_[cg_id].push_back( useq_id );
	int pos2 = al.pos2( ) - pos;
	int Pos2 = al.Pos2( ) - pos;
	alignlet new_al( pos2, Pos2, cg_id, tg0len, al.Fw1( ) );
	ualigns0_[ ualigns0_index_[useq_id][ai] ] = new_al;
	continue;
      }  
    }
  }
  // Update contigs_.
  contigs_[cg_id] = tg0;
  
}

/**
 * ContigsManager
 * CutTail
 */
void ContigsManager::CutTail( size_t cg_id, size_t pos )
{
  cout << "CutTail" << endl;
  // Chunks from original contig.
  const fastavector &original = contigs_[cg_id];
  const size_t original_size = original.size( );
  ForceAssertLt( pos, original_size );
  PRINT3( cg_id, original_size, pos );

  fastavector tg0;
  tg0.SetToSubOf( original, 0, pos );
  int tg0len = tg0.size( );
  
  // Update aligns_, index_, and cg2seqs_.
  {
    vec<int> seq_set = cg2seqs_[cg_id];
    cg2seqs_[cg_id].clear( );
    
    for (int ii=0; ii<seq_set.isize( ); ii++) {
      int seq_id = seq_set[ii];
      if ( aligns0_index_[seq_id] < 0 ) continue;
      alignlet& al = aligns0_[ aligns0_index_[seq_id] ];
      
      // Seq on the left half.
      if ( al.Pos2( ) <= (int)pos ) {
	cg2seqs_[cg_id].push_back( seq_id );
	alignlet new_al( al.pos2( ), al.Pos2( ), cg_id, tg0len, al.Fw1( ) );
	aligns0_[ aligns0_index_[seq_id] ] = new_al;
	continue;
      }
      
      // Seq on the right half.
      else {
	aligns0_index_[seq_id] = -1;
	continue;
      }
    }
  }

  
  // Update aligns_, index_, and cg2useqs_.
  {
    vec<int> useq_set = cg2useqs_[cg_id];
    cg2useqs_[cg_id].clear( );
    
    for (int ii=0; ii<useq_set.isize( ); ii++) {
      int useq_id = useq_set[ii];
      for ( size_t ai = 0; ai < ualigns0_index_[useq_id].size(); ai++ ){
	if ( ualigns0_index_[useq_id][ai] < 0 ) continue;
	alignlet& al = ualigns0_[ ualigns0_index_[useq_id][ai] ];
	if ( al.TargetId() != (int)cg_id ) continue;
	// Seq on the left half.
	if ( al.Pos2( ) <= (int)pos ) {
	  cg2useqs_[cg_id].push_back( useq_id );
	  alignlet new_al( al.pos2( ), al.Pos2( ), cg_id, tg0len, al.Fw1( ) );
	  ualigns0_[ ualigns0_index_[useq_id][ai] ] = new_al;
	  continue;
	}
	
	// Seq on the right half.
	else {
	  ualigns0_index_[useq_id][ai] = -1;
	  continue;
	}
      }
    }
  }
  
  // Update contigs_.
  contigs_[cg_id] = tg0;
  
}



/**
 * ContigsManager
 * SlideSplitContig (duppplicating part between pos1 and pos2 )
 input:           pos1             pos2
                   |                |
          |--------------------------------------|
 output:  |-------------------------|
 |-----------------------------|
*/
vec<size_t> ContigsManager::SlideSplitContig( size_t cg_id, const size_t pos1, const size_t pos2 ){
  vec<size_t> new_ids;

  cout << "SlideSplitContig" << endl;
  // Chunks from original contig.
  const fastavector &original = contigs_[cg_id];
  const size_t original_size = original.size( );
  PRINT4( cg_id, original_size, pos1, pos2 );
  ForceAssertGe( pos1, 0u );
  ForceAssertLe( pos1, pos2 );
  ForceAssertLe( pos2, original_size );

  fastavector tg0; 
  fastavector tg1; 
  tg0.SetToSubOf( original, 0, pos2 );
  tg1.SetToSubOf( original, pos1, original_size - pos1 );
  int tg0len = tg0.size( );
  int tg1len = tg1.size( );
  
  // Update aligns_, index_, and cg2_seqs_.
  size_t new_cg_id = contigs_.size( );
  cg2seqs_.resize( 1 + contigs_.size( ) );
  
  // Make a local copy of the seq set for the contig, and clear current
  vec<int> seq_set = cg2seqs_[cg_id];
  
  // Rebuild cg2seqs for cg_id and for the last contig (the new one).
  cg2seqs_[cg_id].clear( );
  for (int ii=0; ii<seq_set.isize( ); ii++) {
    int seq_id = seq_set[ii];
    alignlet& al = aligns0_[ aligns0_index_[seq_id] ];
    if ( aligns0_index_[seq_id] < 0 ) continue;
    // Seq on the left half.
    if ( al.Pos2( ) < (int)pos1 ) {
      cg2seqs_[cg_id].push_back( seq_id );
      alignlet new_al( al.pos2( ), al.Pos2( ), cg_id, tg0len, al.Fw1( ) );
      aligns0_[ aligns0_index_[seq_id] ] = new_al;
    }else if ( al.pos2( ) >= (int)pos2 ) {
    // Seq on the right half.
      cg2seqs_[new_cg_id].push_back( seq_id );
      int pos2 = al.pos2( ) - pos1;
      int Pos2 = al.Pos2( ) - pos1;
      alignlet new_al( pos2, Pos2, new_cg_id, tg1len, al.Fw1( ) );
      aligns0_[ aligns0_index_[seq_id] ] = new_al;
    }else{
      // A straddler, remove seq (by resetting index to -1) (sliding regions is duplicated, only one
      // alignment per read is allowed in the read index
      aligns0_index_[seq_id] = -1;
    }
  }

  cg2useqs_.resize( 1 + contigs_.size( ) );
  vec<int> useq_set = cg2useqs_[cg_id];
  
  // Rebuild cg2seqs for cg_id and for the last contig (the new one).
  cg2useqs_[cg_id].clear( );
  for (int ii=0; ii<useq_set.isize( ); ii++) {
    int useq_id = useq_set[ii];
    for ( size_t ai = 0; ai < ualigns0_index_[useq_id].size(); ai++ ){
      if ( ualigns0_index_[useq_id][ai] < 0 ) continue;
      alignlet& al = ualigns0_[ ualigns0_index_[useq_id][ai] ];
      if ( al.TargetId() != (int)cg_id ) continue;
      
      // Seq on the left segment.
      if ( al.Pos2( ) < (int)pos2 ) {
	cg2useqs_[cg_id].push_back( useq_id );
	alignlet new_al( al.pos2( ), al.Pos2( ), cg_id, tg0len, al.Fw1( ) );
	ualigns0_[ ualigns0_index_[useq_id][ai] ] = new_al;
	ForceAssert( new_al.pos2() >=0 && new_al.Pos2() <= tg0len );
      }else if ( al.pos2( ) >= (int)pos1 ) {
	// Seq on the right segment half.
	cg2useqs_[new_cg_id].push_back( useq_id );
	int pos2 = al.pos2( ) - pos1;
	int Pos2 = al.Pos2( ) - pos1;
	alignlet new_al( pos2, Pos2, new_cg_id, tg1len, al.Fw1( ) );
	ualigns0_[ ualigns0_index_[useq_id][ai] ] = new_al;
	ForceAssertGe( new_al.pos2(), 0); ForceAssertLe( new_al.Pos2(), tg1len );
      }else{	
	// A straddler, remove seq (by resetting index to -1).
	ualigns0_index_[useq_id][ai] = -1;
      }
    }
  }

  cout << "Updating contigs_" << endl;
  // Update contigs_.
  contigs_[cg_id] = tg0;
  new_ids.push_back( cg_id );
  
  contigs_.push_back( tg1 );
  new_ids.push_back( new_cg_id );
  
  return new_ids;
}


/**
 * ContigsManager
 * MergeContigds
 *
 * Contig cg2 will be resized to 0, and the merged contig will be
 * saved in cg1. Both cg1 and cg2 may be reverse-complemented.  The
 * actual merge is performed triviallyI if cg2 is embedded in cg1,
 * then cg2 is simply removed. Otherwise, extend cg1 with the portion
 * of cg2 specified by the alignment.
 */
void ContigsManager::MergeContigs( int sig_cg1,
				   int sig_cg2,
				   const alignment &al )
{
  cout << "MergeContgis" << endl;
  int cg1 = ( sig_cg1 < 0 ? - sig_cg1 - 1 : sig_cg1 );
  int cg2 = ( sig_cg2 < 0 ? - sig_cg2 - 1 : sig_cg2 );

  // Contig cg2 is embedded in cg1.
  ForceAssertEq( al.pos2( ), 0 );
  if ( al.Pos2( ) == (int)contigs_[cg2].size( ) ) {
    this->RemoveAligns( cg2 );
    contigs_[cg2].clear( );
    return;
  }

  // If needed, flip contigs (for easier bookkeping later on).
  if ( sig_cg1 < 0 ) this->ReverseComplement( (size_t)cg1 );
  if ( sig_cg2 < 0 ) this->ReverseComplement( (size_t)cg2 );
  
  // Remove aligns on cg2 between 0 and Pos2.
  ho_interval alwin( 0, al.Pos2( ) );
  this->RemoveAligns( (size_t)cg2, &alwin );

  // Adjust contigs.
  int len1orig = (int)contigs_[cg1].size( );
  int len2orig = (int)contigs_[cg2].size( );
  fastavector chunk_contig;
  chunk_contig.SetToSubOf( contigs_[cg2], al.Pos2( ), len2orig - al.Pos2( ) );
  contigs_[cg1].Append( chunk_contig );
  contigs_[cg2].clear( );

  int merged_len = len1orig + ( len2orig - al.Pos2( ) );
  ForceAssertEq( merged_len, (int)contigs_[cg1].size( ) );
  
  // Adjust alignments (first on cg1, then on the remaining aligns of cg2).
  {
    const vec<int> &ids1 = cg2seqs_[cg1];
    for (int ii=0; ii<ids1.isize( ); ii++) {
      int idx = aligns0_index_[ ids1[ii] ];
      if ( idx < 0 ) continue;
      aligns0_[idx].SetTargetLength( merged_len, true );
    }
    
    const vec<int> &ids2 = cg2seqs_[cg2];
    for (int ii=0; ii<ids2.isize( ); ii++) {
      int idx = aligns0_index_[ ids2[ii] ];
      if ( idx < 0 ) continue;
      aligns0_[idx].SetTargetId( cg1 );
      aligns0_[idx].SetTargetLength( merged_len, true );
      aligns0_[idx].Shift( len1orig - al.Pos2( ) );
    }
  }

  {
    const vec<int> &ids1 = cg2useqs_[cg1];
    for (int ii=0; ii<ids1.isize( ); ii++) {
      for ( size_t ai = 0; ai < ualigns0_index_[ ids1[ii] ].size(); ai++ ){
	int idx = ualigns0_index_[ ids1[ii] ][ai];
	if ( idx < 0 ) continue;
	if ( ualigns0_[idx].TargetId() != (int)cg1 ) continue;
	ualigns0_[idx].SetTargetLength( merged_len, true );
      }
    }
    
    const vec<int> &ids2 = cg2seqs_[cg2];
    for (int ii=0; ii<ids2.isize( ); ii++) {
      for ( size_t ai = 0; ai < ualigns0_index_[ ids2[ii] ].size(); ai++ ){
	int idx = ualigns0_index_[ ids2[ii] ][ai];
	if ( idx < 0 ) continue;
	if ( ualigns0_[idx].TargetId() != (int)cg2 ) continue;
	ualigns0_[idx].SetTargetId( cg1 );
	ualigns0_[idx].SetTargetLength( merged_len, true );
	ualigns0_[idx].Shift( len1orig - al.Pos2( ) );
      }
    }
  }
}

/**
 * ContigsManager
 * IsolateContig
 */
void ContigsManager::IsolateContig( size_t cg_id )
{
  cout << "IsolateContig" << endl;
  this->RemoveAligns( cg_id );
}



/**
 * ContigsManager
 * ReverseComplement
 * private
 */
void ContigsManager::ReverseComplement( size_t cg_id )
{
  cout << "ReverseComplement contig" << endl;
  contigs_[cg_id].ReverseComplement( );
  {
    const vec<int> &seqids = cg2seqs_[cg_id];
    for (int ii=0; ii<seqids.isize( ); ii++) {
      int seqidx = aligns0_index_[ seqids[ii] ];
      if ( seqidx < 0 ) continue;
      aligns0_[seqidx].Reverse( );
    }
  }
  
  {
    const vec<int> &useqids = cg2useqs_[cg_id];
    for (int ii=0; ii<useqids.isize( ); ii++) {
      for ( size_t ai = 0; ai < ualigns0_index_[ useqids[ii] ].size(); ai++ ){
	int useqidx = ualigns0_index_[ useqids[ii] ][ai];
	if ( useqidx < 0 ) continue;
	ualigns0_[useqidx].Reverse( );
      }
    }
  }
}

/**
 * ContigsManager
 * RemoveAligns
 * private
 *
 * Remove (ie reset index to -1) all seqs aligning cg_id on win. If
 * win is null, it removes all seqs on cg_id.
 */
void ContigsManager::RemoveAligns( size_t cg_id, const ho_interval *win )
{
  cout << "RemoveAligns in contig" << endl;
  // Make a local copy of the seq set for the contig.
  {
    vec<int> ids = cg2seqs_[cg_id];  
    cg2seqs_[cg_id].clear( );
    
    // Remove aligns (by resetting index to -1).
    for (int ii=0; ii<ids.isize( ); ii++) {
      int idx = aligns0_index_[ ids[ii] ];
      if ( idx < 0 ) continue;
      const alignlet &al = aligns0_[idx];
      ho_interval alwin( al.pos2( ), al.Pos2( ) );
      if ( ( !win ) || Overlap( *win, alwin ) > 0 ) {
	aligns0_index_[ ids[ii] ] = -1;
	continue;
      }
      cg2seqs_[cg_id].push_back( ids[ii] );
    }
  }

  {  
    vec<int> ids = cg2useqs_[cg_id];  
    cg2useqs_[cg_id].clear( );
    
    // Remove aligns (by resetting index to -1).
    for (int ii=0; ii<ids.isize( ); ii++) {
      for ( size_t ai = 0; ai < ualigns0_index_[ ids[ii] ].size(); ai++ ){
	int idx = ualigns0_index_[ ids[ii] ][ai];
	if ( idx < 0 ) continue;
	if ( ualigns0_[idx].TargetId() != (int)cg_id ) continue;
	const alignlet &al = ualigns0_[idx];
	ho_interval alwin( al.pos2( ), al.Pos2( ) );
	if ( ( !win ) || Overlap( *win, alwin ) > 0 ) {
	  ualigns0_index_[ ids[ii] ][ai] = -1;
	  continue;
	}
	cg2useqs_[cg_id].push_back( ids[ii] );
      }
    }
    
  }
}



/**
 * ContigsManager
 * Setup
 * private
 */
void ContigsManager::Setup( )
{
  cg2seqs_.resize( contigs_.size( ) );
  for (size_t seq_id=0; seq_id < aligns0_index_.size(); seq_id++) {
    if ( aligns0_index_[seq_id] < 0 ) continue;
    cg2seqs_[ aligns0_[ aligns0_index_[seq_id] ].TargetId( ) ].push_back( seq_id );
  }
  
  cg2useqs_.resize( contigs_.size( ) );
  for (size_t useq_id=0; useq_id < ualigns0_index_.size(); useq_id++) {
    for (size_t ai=0; ai < ualigns0_index_[useq_id].size(); ai++) {
      if ( ualigns0_index_[useq_id][ai] < 0 ) continue;
      cg2useqs_[ ualigns0_[ ualigns0_index_[useq_id][ai] ].TargetId( ) ].push_back( useq_id );
    }
  }
}
