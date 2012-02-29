// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include <fstream>

#include "Loader.h"
#include "STLExtensions.h"


/*
 * LoadNewSuper
 */
void LoadNewSupers( const String &asupers_file,
		    vec<super> &new_supers )
{
  READ( asupers_file, vec<annotated_supercontig>, asupers );
  LoadNewSupers( asupers, new_supers );
}



/*
 * LoadNewSuper
 */
void LoadNewSupers( vec<annotated_supercontig> &asupers,
		    vec<super> &new_supers )
{
  new_supers.clear( );
  new_supers.resize( asupers.size( ) );

  for (int ii=0; ii<(int)new_supers.size( ); ii++) {
    annotated_supercontig& as = asupers[ii];
    
    for (int jj=0; jj<as.NumContigs( ); jj++) {
      new_supers[ii].mtig.push_back( as.Contig(jj).ID( ) );
      if ( jj < as.NumContigs( ) - 1 )
	new_supers[ii].gap.push_back( as.Gap(jj) );
    }
  }
}



/*
 * LoadContigLengths
 */
void LoadContigLengths( const String &contigs_file,
			vec<int> &lengths )
{
  vecbasevector contigs;
  contigs.ReadAll( contigs_file );
  LoadContigLengths( contigs, lengths );
}



/*
 * LoadContigLengths
 */
void LoadContigLengths( const vecbasevector &contigs,
			vec<int> &lengths )
{
  lengths.resize( contigs.size( ) );
  for (int ii=0; ii<(int)contigs.size( ); ii++)
    lengths[ii] = contigs[ii].size( );
}



/*
 * LoadFirstLocs
 */
void LoadFirstLocs( const vec<read_location> &locs,
		    vec<int> &first_locs )
{
  first_locs.clear( );
  
  ForceAssert( is_sorted( locs.begin(), locs.end() ) );
  int max_cg_id = locs[locs.size()-1].Contig();
  first_locs.resize( 1 + max_cg_id , -1 );
  for (int ii=0; ii<(int)locs.size(); ii++) {
    int contig = locs[ii].Contig();
    if ( contig > -1 && first_locs[contig] < 0 )
      first_locs[contig] = ii;
  }
}



/*
 * LoadContigCoverages
 */
void LoadContigCoverages( const vec<read_location> &locs,
			  const vec<int> &first_locs,
			  vec<float> &coverages )
{
  coverages.clear( );
  coverages.resize( first_locs.size( ), 0.0 );

  for (int contig_id=0; contig_id<(int)first_locs.size( ); contig_id++) {
    int first_loc = first_locs[contig_id];
    if ( first_loc < 0 )
      continue;
    
    longlong contig_len = locs[first_loc].LengthOfContig( );
    if ( contig_len < 1 )
      continue;
    
    longlong tot_len = 0;
    for (int ii=first_loc; ii<(int)locs.size( ); ii++) {
      if ( locs[ii].Contig( ) != contig_id )
	break;
      tot_len += locs[ii].LengthOfRead( );
    }
    
    coverages[contig_id] = SafeQuotient( tot_len, contig_len );
  }
}



/*
 * LoadLocs
 */
void LoadLocs( const String &locs_file,
	       vec<read_location> &locs,
	       const vec<Bool> *is_transposon,
	       bool check_sorted )
{
  locs.clear( );
  READX( locs_file, locs );

  if ( check_sorted )
    if ( ! is_sorted( locs.begin( ), locs.end( ) ) )
      sort( locs.begin( ), locs.end( ) );
  
  if ( is_transposon ) {
    vec<read_location>::iterator rl_iter;
    for ( rl_iter = locs.begin(); rl_iter != locs.end(); ++rl_iter ) {
      if ( (*is_transposon)[ rl_iter->ReadId() ] ) {
	rl_iter->SetOrientationOnContig( ( rl_iter->OrientationOnContig()
					   == ForwardOr
					   ? ReverseOr
					   : ForwardOr ) );
      }
    }
  }
}



/*
 * LoadLocsContigPairs
 */
void LoadContigPairs( const vec_inserts &ins,
		      const assembly &the_assembly,
		      vec< contig_pair > &cg_pairs )
{
  const vec<read_location> &all_locs = the_assembly.reads_orig;
  const map<int, int> &mts = the_assembly.mtigs_to_supers;
  map< contig_pair, int> pair_to_pos;
  map< contig_pair, int>::iterator cp_iter;
  map<int, int>::const_iterator mts_iter;
  
  // Loop over all the insert.
  for (int ii=0; ii<(int)ins.size( ); ii++) {
    const insert_ends &the_ins = ins[ii];
    
    int loc1 = the_ins.Loc1( );
    int loc2 = the_ins.Loc2( );
    int cg1 = all_locs[loc1].Contig( );
    int cg2 = all_locs[loc2].Contig( );
    
    mts_iter = mts.find( cg1 );
    ForceAssert ( mts_iter != mts.end( ) );
    int super1 = mts_iter->second;
    
    mts_iter = mts.find( cg2 );
    ForceAssert ( mts_iter != mts.end( ) );
    int super2 = mts_iter->second;
    
    // Links between different supercontigs.
    if ( super1 != super2 )
      continue;

    // Find (or create) contig_pair for this link.
    contig_pair the_pair( super1, cg1, cg2 );
    
    int pair_id = -1;
    cp_iter = pair_to_pos.find( the_pair );
    if ( cp_iter != pair_to_pos.end( ) )
      pair_id = cp_iter->second;
    else {
      cg_pairs.push_back( the_pair );
      pair_id = (int)cg_pairs.size( ) - 1;
      pair_to_pos[ the_pair ] = pair_id;
    }
    
    // Add link to contig_pair;
    cg_pairs[pair_id].AddInsert( &the_ins );
  }

  // Sort.
  sort( cg_pairs.begin( ), cg_pairs.end( ) );
}





