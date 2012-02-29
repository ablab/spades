// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Basevector.h"
#include "SuperMap.h"
#include "Vec.h"



/*
 * super_map
 * Constructor
 */
super_map::super_map( const super &new_super,
		      const vec<int> &contig_lengths )
{
  const vec<int> *the_lengths = &contig_lengths;
  vecbasevector *the_contigs = 0;

  this->Setup( new_super, the_lengths, the_contigs );
}



/*
 * super_map
 * Constructor
 */
super_map::super_map( const super &new_super,
		      const vecbasevector &contigs )
{
  vec<int> *the_lengths = 0;
  const vecbasevector *the_contigs = &contigs;
  
  this->Setup( new_super, the_lengths, the_contigs );
}



/*
 * super_map
 * Setup
 */
void super_map::Setup( const super &new_super,
		       const vec<int> &contig_lengths )
{
  const vec<int> *the_lengths = &contig_lengths;
  vecbasevector *the_contigs = 0;

  this->Setup( new_super, the_lengths, the_contigs );
}



/*
 * super_map
 * Setup
 */
void super_map::Setup( const super &new_super,
		       const vecbasevector &contigs )
{
  vec<int> *the_lengths = 0;
  const vecbasevector *the_contigs = &contigs;
  
  this->Setup( new_super, the_lengths, the_contigs );
}



/*
 * super_map
 * InSupercontig
 */
bool super_map::InSupercontig( int contig_id ) const
{
  map<int, int>::const_iterator iter = contig_to_pos_.find( contig_id );
  return ( iter != contig_to_pos_.end( ) );
}



/*
 * super_map
 * GetNumberContigs
 */
int super_map::GetNumberContigs( ) const
{
  return (int)ids_.size();
}



/*
 * super_map
 * GetSupercontigLength
 */
int super_map::GetSupercontigLength( ) const
{
  return *max_element( stop_.begin( ), stop_.end( ) ) -
    *min_element( start_.begin(), start_.end() ) + 1;
}



/*
 * super_map
 * GetUngappedSupercontigLength
 */
int super_map::GetUngappedSupercontigLength( ) const
{
  int tot_length = 0;
  for (int jj=0; jj<(int)ids_.size(); jj++)
    tot_length += 1 + stop_[jj] - start_[jj];
  
  return tot_length;
}



/*
 * super_map
 * GetContigLength
 */
int super_map::GetContigLength( int contig_id ) const
{
  int pos = this->FindPos( contig_id );
  return 1 + stop_[pos] - start_[pos];
}



/*
 * super_map
 * GetContigLengthPos
 */
int super_map::GetContigLengthPos( int contig_pos ) const
{
  return 1 + stop_[contig_pos] - start_[contig_pos];
}



/*
 * super_map
 * GetContigLengths
 */
vec<int> super_map::GetContigLengths( ) const
{
  vec<int> lengths( ids_.size() );
  for (int jj=0; jj<(int)ids_.size(); jj++)
    lengths[jj] = 1 + stop_[jj] - start_[jj];
  return lengths;
}
  
 

/*
 * super_map
 * GetPositionInSupercontig
 */
int super_map::GetPositionInSupercontig( int contig_id ) const
{
  return this->FindPos( contig_id );
}



/*
 * super_map
 * GetContigId
 */
int super_map::GetContigId( int pos_in_supercontig ) const
{
  return ids_[pos_in_supercontig];
}



/*
 * super_map
 * GetContigIds
 */
vec<int> super_map::GetContigIds( ) const
{
  return ids_;
}



/*
 * super_map
 * GetGap
 */
int super_map::GetGap( int after_pos_in_supercontig ) const
{
  int pos = after_pos_in_supercontig;
  ForceAssert( pos < (int)ids_.size( ) - 1 );
  return ( start_[pos + 1] - stop_[pos] - 1 );
}



/*
 * super_map
 * GetGaps
 */
vec<int> super_map::GetGaps( ) const
{
  vec<int> gaps( ids_.size() - 1 );
  for ( unsigned int ii = 0; ii < gaps.size(); ++ii )
    gaps[ii] = start_[ii+1] - stop_[ii] - 1;
  
  return gaps;
}



/*
 * super_map
 * GetStartOnSupercontig
 */
int super_map::GetStartOnSupercontig( int contig_id ) const
{
  return start_[ this->FindPos( contig_id ) ];
}



/*
 * super_map
 * GetStartOnSupercontigPos
 */
int super_map::GetStartOnSupercontigPos( int contig_pos ) const
{
  return start_[ contig_pos ];
}



/*
 * super_map
 * GetUngappedStartOnSupercontig
 */
int super_map::GetUngappedStartOnSupercontig( int contig_id ) const
{
  int position_on_super = this->FindPos( contig_id );
  
  int start = 0;
  for ( int ii = 0; ii < position_on_super; ++ii )
    start += ( stop_[ii] - start_[ii] + 1 );
  
  return start;
}



/*
 * super_map
 * GetUngappedStartOnSupercontigPos
 */
int super_map::GetUngappedStartOnSupercontigPos( int contig_pos ) const
{
  int start = 0;
  for ( int ii = 0; ii < contig_pos; ++ii )
    start += ( stop_[ii] - start_[ii] + 1 );
  
  return start;
}



/*
 * super_map
 * GetStopOnSupercontig
 */
int super_map::GetStopOnSupercontig( int contig_id ) const
{
  return stop_[ this->FindPos( contig_id ) ];
}



/*
 * super_map
 * GetStopOnSupercontigPos
 */
int super_map::GetStopOnSupercontigPos( int contig_pos ) const
{
  return stop_[ contig_pos ];
}



/*
 * super_map
 * GetUngappedStopOnSupercontig
 */
int super_map::GetUngappedStopOnSupercontig( int contig_id ) const
{
  int position_on_super = this->FindPos( contig_id );

  int stop = 0;
  for ( int ii = 0; ii < position_on_super+1; ++ii )
    stop += ( stop_[ii] - start_[ii] + 1 );
  
  return stop - 1;
}



/*
 * super_map
 * GetUngappedStopOnSupercontigPos
 */
int super_map::GetUngappedStopOnSupercontigPos( int contig_pos ) const
{
  int stop = 0;
  for ( int ii = 0; ii < contig_pos+1; ++ii )
    stop += ( stop_[ii] - start_[ii] + 1 );
  
  return stop - 1;
}



/*
 * super_map
 * SortContigsByStart
 */
void super_map::SortContigsByStart( vec<int> &contig_ids ) const
{
  contig_ids.clear( );
  
  int n_contigs = this->GetNumberContigs( );
  contig_ids.reserve( n_contigs );
  for (int jj=0; jj<n_contigs; jj++)
    contig_ids.push_back( jj );
  
  order_contig_pos_Start sorter( this );
  sort( contig_ids.begin( ), contig_ids.end( ), sorter );

  for (int jj=0; jj<n_contigs; jj++)
    contig_ids[jj] = ids_[ contig_ids[jj] ];
}



/*
 * super_map
 * SortContigsByMaxGap
 */
void super_map::SortContigsByMaxGap( vec<int> &contig_ids ) const
{
  contig_ids.clear( );
  
  int n_contigs = this->GetNumberContigs( );
  contig_ids.reserve( n_contigs );
  for (int jj=0; jj<n_contigs; jj++)
    contig_ids.push_back( jj );
  
  order_contig_pos_MaxGap sorter( this );
  sort( contig_ids.begin( ), contig_ids.end( ), sorter );
  
  for (int jj=0; jj<n_contigs; jj++)
    contig_ids[jj] = ids_[ contig_ids[jj] ];
}



/*
 * super_map
 * Setup
 *
 * To avoid code duplication, all Setups converge here. If contig_lengths
 * is null, then contis is used to determine contig lengths (vice versa,
 * if contigs is null, then contig_lengts is used).
 */
void super_map::Setup( const super &new_super,
		       const vec<int> *contig_lengths,
		       const vecbasevector *contigs )
{
  int n_contigs = (int)new_super.mtig.size( );
  
  this->Clear( );
  this->Reserve( n_contigs );
  
  int pos_on_scg = -1;
  for (int jj=0; jj<n_contigs; jj++) {
    int id = new_super.mtig[jj];
    int length = (contigs) ? (*contigs)[id].size( ) : (*contig_lengths)[id];
    int start = pos_on_scg + 1;
    int stop = pos_on_scg + length;
    
    this->AddContig( id, start, stop );
      
    pos_on_scg = stop;
    if ( jj<n_contigs - 1 )
      pos_on_scg += new_super.gap[jj];
  }
}



/*
 * super_map
 * FindPos
 */
int super_map::FindPos( int contig_id ) const
{
  map<int, int>::const_iterator iter = contig_to_pos_.find( contig_id );
  ForceAssert ( iter != contig_to_pos_.end( ) );
  return iter->second;
}



/*
 * super_map
 * Clear
 */
void super_map::Clear( )
{
  ids_.clear( );
  start_.clear( );
  stop_.clear( );
  contig_to_pos_.clear( );
}



/*
 * super_map
 * Reserve
 */
void super_map::Reserve( unsigned int size )
{
  ids_.reserve( size );
  start_.reserve( size );
  stop_.reserve( size );
}



/*
 * super_map
 * AddContig
 */
void super_map::AddContig( int contig_id, int start, int stop )
{
  ids_.push_back( contig_id );
  start_.push_back( start );
  stop_.push_back( stop );

  ForceAssert( ids_.size( ) == start_.size( ) &&
	       start_.size( ) == stop_.size( ) );

  contig_to_pos_[ contig_id ] = (int)ids_.size( ) - 1;
}



