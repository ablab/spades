///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "CoreTools.h"
#include "IndexedAlignmentPlusVector.h"
#include "PackAlign.h"
#include "STLExtensions.h"
#include "VecAlignmentPlus.h"



/*
 * vec_alignment_plus
 * constructor
 */
vec_alignment_plus::vec_alignment_plus( const String &alignments_file,
					const vec<int> &read_lengths ) :
  read_lengths_ ( &read_lengths )
{
  ForceAssert( IsRegularFile( alignments_file ) );
  this->Load( alignments_file );
}



/*
 * vec_alignment_plus
 * constructor
 */
vec_alignment_plus::vec_alignment_plus( const vec<alignment_plus> &orig_aligns,
					const vec<int> &read_lengths ) :
  read_lengths_ ( &read_lengths )
{
  ForceAssert( is_sorted( orig_aligns.begin( ), orig_aligns.end( ) ) );

  // Reserve memory.
  int n_alignments = orig_aligns.size( );
  
  alignments_.reserve( n_alignments / 2 );
  all_aligns_ids_.reserve ( n_alignments );
  all_aligns_flip_.reserve ( n_alignments );

  // first_align_id[ii] is the location of the first alignment in the
  //  vector alignments_ which contains read ii.
  vec<int> first_align_id;
  first_align_id.resize( read_lengths_->size(), -1 );

  // Fill data structures.
  for (int ii=0; ii<n_alignments; ii++) {
    const alignment_plus *al_plus = &( orig_aligns[ii] );
    
    if ( al_plus->Id1() < al_plus->Id2() ) {
      alignments_.push_back( *al_plus );
      all_aligns_flip_.push_back ( False );
      all_aligns_ids_.push_back( alignments_.size() - 1 );
      
      if ( first_align_id[ al_plus->Id1() ] < 0 )
	first_align_id[ al_plus->Id1() ] = alignments_.size() - 1;
    }
    else {
      int al_id = first_align_id[ al_plus->Id2() ];
      while ( alignments_[al_id].Id2() != al_plus->Id1() ) {
	ForceAssert( alignments_[al_id].Id1() == al_plus->Id2() );
	++al_id;
      }
      
      all_aligns_flip_.push_back ( True );
      all_aligns_ids_.push_back( al_id );
    }
  }
  
  // Check sizes.
  ForceAssert ( all_aligns_ids_.size() == all_aligns_flip_.size() );
}



/*
 * vec_alignment_plus
 * GetAlignment
 */
void vec_alignment_plus::GetAlignment( alignment_plus &al_plus, int align_id ) const
{
  if ( all_aligns_flip_[ align_id ] ) {
    const alignment_plus &local_al = alignments_[ all_aligns_ids_[align_id] ];
    int rd1_length = (*read_lengths_)[ local_al.Id1() ];
    int rd2_length = (*read_lengths_)[ local_al.Id2() ];
    
    al_plus.SetToSwapOf(local_al, rd1_length, rd2_length);
  }
  else {
    int al_id = all_aligns_ids_[align_id];
    al_plus = alignments_[ al_id ];
  }  
  
}



/*
 * vec_alignment_plus
 * GetAlignsIndex
 */
int vec_alignment_plus::GetAlignsIndex( int read_id ) const
{
  // Fill all_aligns_index_ the first time GetAlignsIndex is called.
  if ( all_aligns_index_.empty() ) {
    all_aligns_index_.resize( read_lengths_->size(), -1 );

    // Fill vector.
    for ( int ii = (int)all_aligns_ids_.size()-1; ii>=0; ii-- ) {
      const alignment_plus &al_plus = alignments_[ all_aligns_ids_[ii] ];      
      int id = all_aligns_flip_[ii] ? al_plus.Id2() : al_plus.Id1();

      all_aligns_index_[ id ] = ii;
    }
  }

  // Return.
  return all_aligns_index_[ read_id ];

}



/*
 * vec_alignment_plus
 * GetAlignmentId1
 */
int vec_alignment_plus::GetAlignmentId1( int align_id ) const
{
  int al_pos = all_aligns_ids_[ align_id ];
  bool al_flip = all_aligns_flip_[ align_id ];

  if ( al_flip )
    return alignments_[ al_pos ].Id2();
  else
    return alignments_[ al_pos ].Id1();

}



/*
 * vec_alignment_plus
 * GetAlignmentId2
 */
int vec_alignment_plus::GetAlignmentId2( int align_id ) const
{
  int al_pos = all_aligns_ids_[ align_id ];
  bool al_flip = all_aligns_flip_[ align_id ];

  if ( al_flip )
    return alignments_[ al_pos ].Id1();
  else
    return alignments_[ al_pos ].Id2();

}



/*
 * vec_alignment_plus
 * GetAlignmentScore
 */
float vec_alignment_plus::GetAlignmentScore( int align_id ) const
{
  int al_pos = all_aligns_ids_[ align_id ];
  return alignments_[ al_pos ].score;
}



/*
 * vec_alignment_plus
 * GetAlignmentLength
 *
 * Alignment length is defined either as Pos1-pos1, or as Pos2-pos2 (set
 * use_pos2=true for the latter).
 */
int vec_alignment_plus::GetAlignmentLength( int align_id,
					    bool use_pos2 ) const
{
  int al_pos = all_aligns_ids_[ align_id ];

  if ( use_pos2 )
    return alignments_[ al_pos ].a.Pos2( ) - alignments_[ al_pos ].a.pos2( );

  return alignments_[ al_pos ].a.Pos1( ) - alignments_[ al_pos ].a.pos1( );
}



/*
 * vec_alignment_plus
 * GetSymmetricAlignmentLength
 *
 * It returns the average between the two possible alignment lengths.
 */
int vec_alignment_plus::GetSymmetricAlignmentLength( int align_id ) const
{
  int al_pos = all_aligns_ids_[align_id];
  int len1 = alignments_[al_pos].a.Pos1( ) - alignments_[al_pos].a.pos1( );
  int len2 = alignments_[al_pos].a.Pos2( ) - alignments_[al_pos].a.pos2( );

  return ( ( len1 + len2 ) / 2 );
}



/*
 * vec_alignment_plus
 * Rc2
 */
Bool vec_alignment_plus::Rc2( int align_id ) const
{
  int al_pos = all_aligns_ids_[ align_id ];
  return alignments_[ al_pos ].Rc2( );
}



/*
 * vec_alignment_plus
 * SaveAlignments
 */
void vec_alignment_plus::SaveAlignments( const String &filename )
{
  VecAlignmentPlusWriter( filename ).Write( *this );
}



/*
 * vec_alignment_plus
 * SetPlainAlignment
 */
void vec_alignment_plus::SetPlainAlignment( int align_id,
					    int pos1,
					    int pos2,
					    int errors,
					    const avector<int>& gaps,
					    const avector<int>& lengths,
					    int nblocks )
{
  // Position in the alignments_ vector.
  int al_pos = all_aligns_ids_[ align_id ];

  if ( all_aligns_flip_[ align_id ] ) {
    alignment_plus &the_align = alignments_[ al_pos ];
    int rd1_length = (*read_lengths_)[ the_align.Id1() ];
    int rd2_length = (*read_lengths_)[ the_align.Id2() ];

    // Flip.
    alignment_plus local_al;
    local_al.SetToSwapOf( the_align, rd1_length, rd2_length );

    // Set.
    local_al.a.Set( pos1, pos2, errors, gaps, lengths, nblocks );
    
    // Flip back.
    the_align.SetToSwapOf( local_al, rd2_length, rd1_length );
  }
  else
    alignments_[ al_pos ].a.Set( pos1, pos2, errors, gaps, lengths, nblocks );
  
}



/*
 * vec_alignment_plus
 * SetPlainAlignment
 */
void vec_alignment_plus::SetPlainAlignment( int align_id,
					    alignment &plain_al )
{
  // Position in the alignments_ vector.
  int al_pos = all_aligns_ids_[ align_id ];
  
  if ( all_aligns_flip_[ align_id ] ) {
    alignment_plus &the_align = alignments_[ al_pos ];
    int rd1_length = (*read_lengths_)[ the_align.Id1() ];
    int rd2_length = (*read_lengths_)[ the_align.Id2() ];
    
    // Flip.
    alignment_plus local_al;
    local_al.SetToSwapOf( the_align, rd1_length, rd2_length );
    
    // Set.
    local_al.a = plain_al;
    
    // Flip back.
    the_align.SetToSwapOf( local_al, rd2_length, rd1_length );
  }
  else
    alignments_[ al_pos ].a = plain_al;
  
}



/*
 * vec_alignment_plus
 * SetPlainAlign
 */
void vec_alignment_plus::SetPlainAlign( int align_id,
					align &plain_al )
{
  // Position in the alignments_ vector.
  int al_pos = all_aligns_ids_[ align_id ];
  
  if ( all_aligns_flip_[ align_id ] ) {
    alignment_plus &the_align = alignments_[ al_pos ];
    int rd1_length = (*read_lengths_)[ the_align.Id1() ];
    int rd2_length = (*read_lengths_)[ the_align.Id2() ];
    
    // Flip.
    alignment_plus local_al;
    local_al.SetToSwapOf( the_align, rd1_length, rd2_length );
    
    // Set.
    local_al.a.Set( plain_al );
    
    // Flip back.
    the_align.SetToSwapOf( local_al, rd2_length, rd1_length );
  }
  else
    alignments_[ al_pos ].a.Set( plain_al );
  
}



/*
 * vec_alignment_plus
 * KillImproperAligns
 */
int vec_alignment_plus::KillImproperAligns( )
{
  ForceAssert( read_lengths_ );

  int tot_killed = 0;
  for (int ii=0; ii<(int)alignments_.size( ); ii++) {
    alignment_plus &ap = alignments_[ii];
    if ( ap.a.Dead( ) ) continue;
    if ( RequireProper( ap, *read_lengths_, 666, False ) ) continue;
    ap.a.Kill( );
    tot_killed++;
  }

  return 2 * tot_killed;
}



/*
 * vec_alignment_plus
 * SetAlignmentScore
 */
void vec_alignment_plus::SetAlignmentScore( int align_id, float score )
{
  // Position in the alignments_ vector.
  int al_pos = all_aligns_ids_[ align_id ];
  
  // Set score.
  alignments_[ al_pos ].score = score;
  
}



/*
 * vec_alignment_plus
 * Load
 */
void vec_alignment_plus::Load( const String &alignments_file )
{
  VecAlignmentPlusReader( alignments_file ).ReadHalf( alignments_, 
                                                      all_aligns_flip_, 
                                                      all_aligns_ids_ );
}

// =================================================================================
//
// alignment index code
//
// =================================================================================

int rai_fd(-1);

void BuildAlignsIndexFin( int n_reads, const vec<int>& aligns_to,
     vec<int>& aligns_to_index, const String& run_dir )
{    for ( int i = 1; i <= n_reads; i++ )
          if ( aligns_to_index[i] < 0 ) aligns_to_index[i] = aligns_to_index[i-1];
     int fd = Open( run_dir + "/aligns.index", O_WRONLY | O_CREAT );
     WriteBytes( fd, &aligns_to_index[0], (n_reads + 1) * sizeof(int) );
     WriteBytes( fd, &aligns_to[0],
          (longlong) aligns_to.size( ) * (longlong) sizeof(int) );
     close(fd);
     if ( rai_fd >= 0 )
     {    close(rai_fd);
          rai_fd = -1;    }    }

void BuildAlignsIndex( const String& run_dir, const vec_alignment_plus& all_aligns,
     int n_reads )
{    int n_aligns = all_aligns.GetNumberAlignments( );
     vec<int> aligns_to(n_aligns), aligns_to_index( n_reads + 1, -1 );
     aligns_to_index[0] = 0;
     for ( int i = 0; i < n_aligns; i++ )
     {    int j, id1 = all_aligns.GetAlignmentId1(i);
          for ( j = i + 1; j < n_aligns; j++ )
               if ( id1 != all_aligns.GetAlignmentId1(j) ) break;
          for ( int k = i; k < j; k++ )
               aligns_to[k] = all_aligns.GetAlignmentId2(k);
          sort( aligns_to.begin( ) + i, aligns_to.begin( ) + j );
          aligns_to_index[id1] = i, aligns_to_index[ id1 + 1 ] = j;
          i = j - 1;    }
     BuildAlignsIndexFin( n_reads, aligns_to, aligns_to_index, run_dir );    }

void BuildAlignsIndex( const String& run_dir, const vec<alignment_plus>& all_aligns,
     int n_reads )
{    int n_aligns = all_aligns.size( );
     vec<int> aligns_to(n_aligns), aligns_to_index( n_reads + 1, -1 );
     aligns_to_index[0] = 0;
     for ( int i = 0; i < n_aligns; i++ )
     {    int j, id1 = all_aligns[i].Id1( );
          for ( j = i + 1; j < n_aligns; j++ )
               if ( id1 != all_aligns[j].Id1( ) ) break;
          for ( int k = i; k < j; k++ )
               aligns_to[k] = all_aligns[k].Id2( );
          sort( aligns_to.begin( ) + i, aligns_to.begin( ) + j );
          aligns_to_index[id1] = i, aligns_to_index[ id1 + 1 ] = j;
          i = j - 1;    }
     BuildAlignsIndexFin( n_reads, aligns_to, aligns_to_index, run_dir );    }

void AlignsIndexReader::readIndex( int id1, vec<int>& to ) const
{
    int startstop[2];
    mFR.seek(id1 * sizeof(int));
    mFR.read(startstop, sizeof(startstop));
    to.resize(startstop[1] - startstop[0]);
    if ( !to.empty() )
    {
        mFR.seek( sizeof(int) * (mNReads + 1ul + startstop[0]) );
        mFR.read(&to[0], to.size() * sizeof(int));
    }
}
