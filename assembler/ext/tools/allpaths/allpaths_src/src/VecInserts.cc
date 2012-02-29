// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//

#include "InsertEnds.h"
#include "ReadLocation.h"
#include "ReadPairing.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"
#include "VecInserts.h"



/*
 * vec_inserts
 * Constructor
 */
vec_inserts::vec_inserts( const assembly *the_assembly ) :
  the_assembly_ ( the_assembly )
{ }



/*
 * vec_inserts
 * Clear
 */
void vec_inserts::Clear( )
{
  contig_ids_.clear( );
  inserts_.clear( );
}



/*
 * vec_inserts
 * AddAllContigs
 */
void vec_inserts::AddAllContigs( )
{
  // Wipe out old existing data.
  contig_ids_.clear( );
  inserts_.clear( );

  // Call AddContigs.
  int n_contigs = the_assembly_->reads_orig_index.size( );
  vec<int> contig_ids;
  contig_ids.reserve( n_contigs );
  for (int ii=0; ii<n_contigs; ii++)
    contig_ids.push_back( ii );

  this->AddContigs( contig_ids );
}



/*
 * vec_inserts
 * AddSupercontigs
 */
void vec_inserts::AddSupercontigs( const vec<int> &super_ids )
{
  int n_contigs = 0;
  for (int ii=0; ii<(int)super_ids.size( ); ii++)
    n_contigs += (int)the_assembly_->supers[ super_ids[ii] ].mtig.size( );
  
  vec<int> contig_ids;
  contig_ids.reserve( n_contigs );

  for (int ii=0; ii<(int)super_ids.size( ); ii++) {
    const vec<int> &mtigs = the_assembly_->supers[ super_ids[ii] ].mtig;
    copy( mtigs.begin( ), mtigs.end( ), back_inserter( contig_ids ) );
  }

  this->AddContigs( contig_ids );
}



/*
 * vec_inserts
 * AddContigs
 */
void vec_inserts::AddContigs( const vec<int> &contig_ids )
{
  const vec<int> &old_contigs = contig_ids_;
  const vec<int> &new_contigs = contig_ids;

  const vec< read_pairing > &pairs = the_assembly_->pairs;
  const vec< read_location > &locs = the_assembly_->reads_orig;
  const vec< vec<int> > &reads_index = the_assembly_->reads_orig_index;
  const vec< int > &pairs_index = the_assembly_->pairs_index;

  // Create (and sort) the id_locs map (pairs of <read_id, location_id>).
  vec< const vec<int>* > all_contigs;
  all_contigs.push_back( &old_contigs );
  all_contigs.push_back( &new_contigs );
  
  int pairs_count = 0;
  for (int ii=0; ii<2; ii++) {
    const vec<int> &cgs = *all_contigs[ii];
    for (int jj=0; jj<(int)cgs.size( ); jj++) 
      pairs_count += reads_index[ cgs[jj] ].size( );
  }

  vec< pair<int, int> > id_locs;
  id_locs.reserve( pairs_count );

  for (int ii=0; ii<2; ii++) {
    const vec<int> &cgs = *all_contigs[ii];
    for (int jj=0; jj<(int)cgs.size( ); jj++) {
      const vec<int> &loc_ids = reads_index[ cgs[jj] ];
      for (int kk=0; kk<(int)loc_ids.size( ); kk++) {
	int loc_id = loc_ids[kk];
	int read_id = locs[ loc_ids[kk] ].ReadId( );
	id_locs.push_back( pair<int, int>( read_id, loc_id ) );
      }
    }
  }

  sort( id_locs.begin( ), id_locs.end( ) );

  // Estimate number of links and reserve memory (taking into account
  //  the fact there may be multiply placed reads).
  int links_capacity = inserts_.capacity( );
  if ( links_capacity < (int)( 1.3 * float ( id_locs.size( ) ) ) ) {
    int estimated_links_count = (int)( 1.5 * float( id_locs.size( ) ) );
    inserts_.reserve( estimated_links_count );
  }
    
  // Add new inserts. To avoid duplications, keep inserts iff they are in
  //  the same read_id order as in the read_pairing they belong to.
  vec< pair<int,int> >::iterator iter;
  
  for (int ii=0; ii<(int)new_contigs.size( ); ii++) {
    const vec<int> &loc_ids = reads_index[ new_contigs[ii] ];
    for (int jj=0; jj<(int)loc_ids.size( ); jj++) {
      int loc_id = loc_ids[jj];
      int read_id = locs[ loc_ids[jj] ].ReadId( );
      int pair_id = pairs_index[read_id];
      
      // Unpaired read.
      if ( pair_id < 0 )
	continue;
      
      // Partner read and its locations.
      vec<int> p_loc_ids;

      const read_pairing &the_pair = pairs[ pair_id ];
      int p_read_id = (the_pair.id1 == read_id) ? the_pair.id2 : the_pair.id1;
      pair<int, int> test_pair( p_read_id, -1 );
      iter = lower_bound( id_locs.begin( ), id_locs.end( ), test_pair );
      while( iter != id_locs.end( ) && iter->first == p_read_id ) {
	p_loc_ids.push_back( iter->second );
	iter++;
      }

      // No read_id's partners found in contig_ids_.
      if ( p_loc_ids.size( ) < 1 )
	continue;

      // Analyze each pair of locations.
      for (int kk=0; kk<(int)p_loc_ids.size( ); kk++) {
	pair<vec<int>::const_iterator, vec<int>::const_iterator> result;
	int p_loc_id = p_loc_ids[kk];
	int p_cg = locs[p_loc_id].Contig( );
	result = equal_range( old_contigs.begin( ), old_contigs.end( ), p_cg );

	int id1 = the_pair.id1;
	int id2 = the_pair.id2;
	insert_ends new_ins;
	
	// Partner location belongs to the old contigs.
	if ( result.first < old_contigs.end( ) ) {
	  if ( locs[loc_id].ReadId( ) == id1 )
	    new_ins.Set( loc_id, p_loc_id, pair_id );
	  else
	    new_ins.Set( p_loc_id, loc_id, pair_id );
	  inserts_.push_back( new_ins );
	  continue;
	}
	
	// Partner location belongs to the new contigs;
	if ( locs[loc_id].ReadId( ) == id1 ) {
	  new_ins.Set( loc_id, p_loc_id, pair_id );
	  inserts_.push_back( new_ins );
	}
      }
    }
  }

  // Copy contig_ids at the end of contig_ids_ and sort.
  copy( new_contigs.begin(), new_contigs.end(), back_inserter( contig_ids_ ) );
  sort( contig_ids_.begin(), contig_ids_.end() );
  
  // Sort inserts_.
  order_InsertEnds_Contigs sorter( locs );
  sort( inserts_.begin( ), inserts_.end( ), sorter );
}



/*
 * vec_inserts
 * FindGapLinks
 *
 * It finds all the inserts spanning gap after contig pos in super_id.
 * All means really all: illogical, invalid and valid.  Output is
 * insert_ids (ids in *this of the inserts spanning gap).
 */
void vec_inserts::FindGapLinks( int super_id,
				int pos,
				vec<int> &insert_ids ) const
{
  insert_ids.clear( );
  
  for (int ins_id=0; ins_id<(int)inserts_.size( ); ins_id++) {
    const insert_ends &the_insert = inserts_[ins_id];
    const vec<read_location> &locs = the_assembly_->reads_orig;
    int cg1 = locs[the_insert.Loc1( )].Contig( );
    int cg2 = locs[the_insert.Loc2( )].Contig( );
    map<int,int>::const_iterator iter;
    iter = the_assembly_->mtigs_to_supers.find( cg1 );
    if ( iter == the_assembly_->mtigs_to_supers.end( ) )
      continue;
    int sup1 = iter->second;
    if ( sup1 != super_id )
      continue;
    iter = the_assembly_->mtigs_to_supers.find( cg2 );
    if ( iter == the_assembly_->mtigs_to_supers.end( ) )
      continue;
    int sup2 = iter->second;
    if ( sup2 != super_id )
      continue;
    iter = the_assembly_->mtigs_to_super_pos.find( cg1 );
    ForceAssert( iter != the_assembly_->mtigs_to_super_pos.end( ) );
    int pos1 = iter->second;
    iter = the_assembly_->mtigs_to_super_pos.find( cg2 );
    ForceAssert( iter != the_assembly_->mtigs_to_super_pos.end( ) );
    int pos2 = iter->second;
    bool good_12 = ( pos1 <= pos && pos2 > pos );
    bool good_21 = ( pos2 <= pos && pos1 > pos );
    if ( ! ( good_12 || good_21 ) )
      continue;

    insert_ids.push_back( ins_id );
  }
}



/*
 * vec_inserts
 * FindSpotLinks
 *
 * Similar to FindGapLinks, the difference being that instead of finding
 * links over a gap, it finds links spanning a spot (pivot) on a given
 * contig_id. Output is stored in insert_ids.
 *
 * Definition of spanning link: call pos the position in super of contig_id,
 * and leg_left, leg_right the two legs of the insert. A link is spanning
 * if leg_left is placed either on [0, pos), or it is placed on pos but it
 * starts before pivot; and if leg_right is either placed > pos, or it is
 * placed at pos but it stops after pivot.
 */
void vec_inserts::FindSpotLinks( int contig_id,
				 int pivot,
				 vec<int> &insert_ids ) const
{
  insert_ids.clear( );
  
  map<int,int>::const_iterator iter;
  iter = the_assembly_->mtigs_to_supers.find( contig_id );
  ForceAssert( iter != the_assembly_->mtigs_to_supers.end( ) );
  int super_id = iter->second;
  iter = the_assembly_->mtigs_to_super_pos.find( contig_id );
  ForceAssert( iter != the_assembly_->mtigs_to_super_pos.end( ) );
  int pos = iter->second;

  for (int ins_id=0; ins_id<(int)inserts_.size( ); ins_id++) {
    const insert_ends &the_insert = inserts_[ins_id];
    const vec<read_location> &locs = the_assembly_->reads_orig;
    int loc1_id = the_insert.Loc1( );
    int loc2_id = the_insert.Loc2( );
    int cg1 = locs[loc1_id].Contig( );
    int cg2 = locs[loc2_id].Contig( );
    iter = the_assembly_->mtigs_to_supers.find( cg1 );
    if ( iter == the_assembly_->mtigs_to_supers.end( ) )
      continue;
    int sup1 = iter->second;
    if ( sup1 != super_id )
      continue;
    iter = the_assembly_->mtigs_to_supers.find( cg2 );
    if ( iter == the_assembly_->mtigs_to_supers.end( ) )
      continue;
    int sup2 = iter->second;
    if ( sup2 != super_id )
      continue;
    iter = the_assembly_->mtigs_to_super_pos.find( cg1 );
    ForceAssert( iter != the_assembly_->mtigs_to_super_pos.end( ) );
    int pos1 = iter->second;
    iter = the_assembly_->mtigs_to_super_pos.find( cg2 );
    ForceAssert( iter != the_assembly_->mtigs_to_super_pos.end( ) );
    int pos2 = iter->second;
    int start1 = locs[loc1_id].StartOnContig( );
    int start2 = locs[loc2_id].StartOnContig( );
    int stop1 = locs[loc1_id].StopOnContig( );
    int stop2 = locs[loc2_id].StopOnContig( );
    if ( pos2 < pos1 || ( pos1 == pos2 && start2 < start1 ) ) {
      swap( pos1, pos2 );
      swap( start1, start2 );
      swap( stop1, stop2 );
    }
    if ( pos1 > pos || ( pos1 == pos && start1 > pivot ) )
      continue;
    if ( pos2 < pos || ( pos2 == pos && stop2 < pivot ) )
      continue;

    insert_ids.push_back( ins_id );
  }
}



