// Copyright (c) 2003 Whitehead Institute for Biomedical Research
//


// Copyright (c) 2003 Whitehead Institute for Biomedical Research
#include <map>

#include "ReadLocationUtil.h"
#include "STLExtensions.h"
#include "VecAlignmentPlus.h"



int AlignsCount( const vec<read_location> &locs,
		 int loc_id,
                 int *n_aligns_found,
                 AlignsIndexReader const& alignsIndexRdr,
		 int min_align,
		 int max_skip )
{
  int aligns_count = 0;
  int contig_id = locs[loc_id].Contig( );
  
  // Parse aligns.index to determine the other reads loc_id aligns to.
  vec<int> to;
  alignsIndexRdr.readIndex(locs[loc_id].ReadId(),to);
  *n_aligns_found = 0;
  
  // Left of loc_id.
  int skip_count = 0;
  for (int jj=loc_id-1; jj>=0; jj--) {
    if ( locs[jj].Contig( ) != contig_id )
      break;
    if ( skip_count >= max_skip )
      break;
    int locs_overlap = LocsOverlap( locs, jj, loc_id );
    if ( locs_overlap < min_align ) {
      skip_count++;
      continue;
    }
    aligns_count++;
    if ( binary_search( to.begin( ), to.end( ), locs[jj].ReadId( ) ) )
      *n_aligns_found += 1;
  }
  
  // Right of loc_id.
  skip_count = 0;
  for (int jj=loc_id+1; jj<(int)locs.size( ); jj++) {
    if ( locs[jj].Contig( ) != contig_id )
      break;
    if ( skip_count >= max_skip )
      break;
    int locs_overlap = LocsOverlap( locs, loc_id, jj);
    if ( locs_overlap < min_align ) {
      skip_count++;
      continue;
    }
    aligns_count++;
    if ( binary_search( to.begin( ), to.end( ), locs[jj].ReadId( ) ) )
      *n_aligns_found += 1;
  }

  // Return.
  return aligns_count;
}



void ReadsInContig( const vec<read_location> &locs, vec<int> &reads_in_contig )
{
  reads_in_contig.clear( );

  // This is only a guess, since locs may not be sorted.
  int n_contigs = 1 + locs.back( ).Contig( );

  reads_in_contig.resize( n_contigs, 0 );
  for (int ii=0; ii<(int)locs.size( ); ii++) {
    if ( (int)reads_in_contig.size( ) <= locs[ii].Contig( ) )
      reads_in_contig.resize( 1 + locs[ii].Contig( ), 0 );
    reads_in_contig[locs[ii].Contig( )] += 1;
  }
}



int LocsFileContigsCount( const String &locs_file )
{
  int largest_contig_id = -1;

  ifstream in( locs_file.c_str( ) );
  int n_locs;
  in >> n_locs;
  char c;
  in.get( c );
  for (int ii=0; ii<n_locs; ii++) {
    read_location loc;
    in.read( (char*) &loc, sizeof( read_location ) );
    largest_contig_id = Max( largest_contig_id, loc.Contig( ) );
  }
  in.close( );

  return largest_contig_id + 1;
}



int LocsOverlap( const vec<read_location> &locs, int id1, int id2 )
{
  const read_location &loc1 = locs[id1];
  const read_location &loc2 = locs[id2];
  if ( loc1.Contig( ) != loc2.Contig( ) )
    return 0;
  int left = Max( loc1.StartOnContig( ), loc2.StartOnContig( ) );
  int right = 1 + Min( loc1.StopOnContig( ), loc2.StopOnContig( ) );
  
  return Max( 0 , right - left );
}



int SignedLocsOverlap( const vec<read_location> &locs, int id1, int id2 )
{
  const read_location &loc1 = locs[id1];
  const read_location &loc2 = locs[id2];
  ForceAssert( loc1.Contig( ) == loc2.Contig( ) );
  int left = Max( loc1.StartOnContig( ), loc2.StartOnContig( ) );
  int right = 1 + Min( loc1.StopOnContig( ), loc2.StopOnContig( ) );
  
  return ( right - left );
}



int LocsOffset( const vec<read_location> &locs, int id1, int id2 )
{ 
  const read_location &loc1 = locs[id1];
  const read_location &loc2 = locs[id2];
  ForceAssert( loc1.Contig( ) == loc2.Contig( ) );
  
  return loc2.StartOnContig( ) - loc1.StartOnContig( );
}



void MergeLocations( const vec<read_location> &unmerged_locs, 
		     vec<read_location> &mapping_locs,
		     vec<read_location> &merged_locs )
{

  if ( mapping_locs.empty() )
  {
    merged_locs = unmerged_locs;
    return;
  }

  sort( mapping_locs.begin(), 
        mapping_locs.end(), 
        order_read_locations_by_readid() );
 
  bool is_sorted_by_contig = is_sorted( unmerged_locs.begin(), unmerged_locs.end(),
                                        order_read_locations_by_contig() );

  pair< vec<read_location>::iterator, vec<read_location>::iterator > range;
  int num_locs( unmerged_locs.size() );
  for ( int i=0; i<num_locs; ++i )
  {
      const read_location &unmerged_loc = unmerged_locs[i];
      
      range = equal_range( mapping_locs.begin(), 
                           mapping_locs.end(),
                           unmerged_loc,
                           order_read_locations_by_readid() );
      
      //  if read at this location did not take part in a merger, add the location to the list
      if ( range.first == range.second )
          merged_locs.push_back( unmerged_loc );

      //  else calculate the merged read's location from the unmerged
      //  read's location and the mapping information.
      else
      {
          ForceAssertEq( distance( range.first, range.second ), 1 );

          // we want a copy of the mapping location
          read_location mapping_loc = *range.first;

          read_location merged_loc( unmerged_loc );
          merged_loc.SetReadId( mapping_loc.Contig() );
          merged_loc.SetLengthOfRead( mapping_loc.LengthOfContig() );

          // If the unmerged read's orientation agrees with its
          // orientation in the map, then the merged read is forward
          // on the contig.
          if ( unmerged_loc.OrientationOnContig() == mapping_loc.OrientationOnContig() )
          {
              merged_loc.SetOrientationOnContig( ForwardOr );
          }

          // If they disagree, the merged read is reversed on the
          // contig and we flip the copy of the mapping location to
          // make our calculations easier later on.
          else
          {
              mapping_loc.Reverse();
              // read_location::Reverse() has a bug
              mapping_loc.SetStartOnContig( mapping_loc.StartOnContig() - 1 );
              merged_loc.SetOrientationOnContig( ReverseOr );
          }

          merged_loc.SetStartOnContig( unmerged_loc.StartOnContig() - 
                                       mapping_loc.StartOnContig() );

          // We now go hunting for this read in the already merged
          // locations.  Three possibilities exist:

          // 1. A placement of this merged read exists but does not
          // agree exactly with this one: we barf.
          
          // 2. A placement of this merged read exists and agrees
          // exactly: we do nothing.

          // 3. No placement of this merged read exists: we add this
          // placement to the list of merged locs.
          
          const int tolerance = 0;

          bool found = false;
          
          // If we know that the unmerged_locs are sorted by contig,
          // we can stop searching when we hit a different contig.  If
          // we don't know that, we have to go all the way through the
          // whole vector.

          vec<read_location>::reverse_iterator merged_iter = merged_locs.rbegin();

          for ( ; merged_iter != merged_locs.rend() && ! found; ++merged_iter )
          {
              if ( is_sorted_by_contig && merged_iter->Contig() != merged_loc.Contig() )
                  break;

              if ( merged_iter->ReadId() == merged_loc.ReadId() &&
                   merged_iter->Contig() == merged_loc.Contig() )
              {
                  if ( merged_iter->OrientationOnContig() != merged_loc.OrientationOnContig() ||
                       Abs( merged_iter->StartOnContig() - 
                            merged_loc.StartOnContig() ) > tolerance ||
                       Abs( merged_iter->StopOnContig() -
                            merged_loc.StopOnContig() ) > tolerance )
                  {
                      cout << endl;
                      cout << "Unable to merge reads." << endl;

                      cout << "Unmerged location: " << unmerged_loc;
                      cout << "Mapping location: " << mapping_loc;
                      cout << "Merged location: " << merged_loc;

                      cout << "Existing merged location: " << *merged_iter;
                      exit(-1);
                  }
                  else
                      found = true;
              }
          }

          if ( ! found )
              merged_locs.push_back( merged_loc );
      }
  }

  sort(merged_locs.begin(), merged_locs.end());

}

void UnmergeLocations( const vec<read_location> &merged_locs,
		       vec<read_location> &mapping_locs,
		       vec<read_location> &unmerged_locs )
{
  
  //  contig id in map corresponds to the read id of one of the merged/rc-ed reads
  sort( mapping_locs.begin(), 
	mapping_locs.end(), 
	order_read_locations_by_contig() );
  
  pair< vec<read_location>::iterator, vec<read_location>::iterator > range;

  int num_locs = merged_locs.size();
  for ( int i=0; i<num_locs; ++i )
  {
      const read_location &merged_loc = merged_locs[i];
      
      read_location target;
      target.SetContig( merged_loc.ReadId() );

      range = equal_range( mapping_locs.begin(), 
                           mapping_locs.end(),
                           target,
                           order_read_locations_by_contig() );

      if ( range.first == range.second )
          unmerged_locs.push_back( merged_loc );
      else
      {
          while ( range.first != range.second )
          {
              read_location rl(*range.first);

              if ( merged_locs[i].OrientationOnContig() == ReverseOr )
              {
                  rl.Reverse();
                  // read_location::Reverse() has a bug
                  rl.SetStartOnContig( rl.StartOnContig() - 1 );
              }

              rl.SetContig( merged_locs[i].Contig() );
              rl.SetLengthOfContig( merged_locs[i].LengthOfContig() );

              rl.SetStartOnContig( merged_locs[i].StartOnContig() + rl.StartOnContig() );

              unmerged_locs.push_back( rl );
              ++range.first;
          }
      }
  }

  sort( unmerged_locs.begin(), unmerged_locs.end() );
  
}

bool Stretch( const read_location &locA,
	      const read_location &locB,
	      const read_pairing &pair,
	      float &stretch )
{
  bool is_logic = false;
  stretch = 0.0;
  
  // Check pair.
  int id1 = locA.ReadId( );
  int id2 = locB.ReadId( );
  if ( id1 != pair.id1 )
    swap( id1, id2 );
  ForceAssert( id1 == pair.id1 && id2 == pair.id2 );

  // The two locs belong to different contigs.
  if ( locA.Contig( ) != locB.Contig( ) )
    return is_logic;

  // Same orientation (really illogical).
  if ( locA.OrientationOnContig( ) == locB.OrientationOnContig( ) )
    return is_logic;

  is_logic = true;

  // Which is fw, which is rc.
  orientation orA = locA.OrientationOnContig( );
  orientation orB = locB.OrientationOnContig( );
  const read_location &locFw = orA == ForwardOr ? locA : locB;
  const read_location &locRc = orA == ForwardOr ? locB : locA;

  // Find stretch and return.
  int rc_start = locRc.StartOnContig( );
  int fw_end = locFw.StopOnContig( ) + 1;
  int sep_observed = rc_start - fw_end;
  int sep_given = pair.sep;
  int stdev = pair.sd;
  ForceAssert( stdev > 0 );

  stretch = ( (float)( sep_observed - sep_given ) / (float)stdev );
  
  return is_logic;
}

void AllStretches( int contig_id,
		   const vec<read_location> &locs,
		   const vec<int> &flocs,
		   const vec<int> &to_pair,
		   const vec<read_pairing> &pairs,
		   vec< pair<int,float> > &pair2stretch,
		   bool append )
{
  if ( ! append )
    pair2stretch.clear( );

  int floc = flocs[contig_id];
  if ( floc < 0 )
    return;

  map<int,int> to_loc;
  for (int ii=floc; ii<(int)locs.size( ); ii++) {
    if ( locs[ii].Contig( ) != contig_id )
      break;
    to_loc[locs[ii].ReadId( )] = ii;
  }
  
  for (int ii=floc; ii<(int)locs.size( ); ii++) {
    if ( locs[ii].Contig( ) != contig_id )
      break;

    int read_id = locs[ii].ReadId( );
    int pair_id = to_pair[read_id];
    if ( pair_id < 0 )
      continue;
    const read_pairing &pair = pairs[pair_id];
    int partner_id = pair.id1 == read_id ? pair.id2 : pair.id1;
    map<int,int>::iterator it = to_loc.find( partner_id );
    if ( it == to_loc.end( ) )
      continue;
    int jj = it->second;
    if ( jj < ii )
      continue;
    
    float stretch = 0.0;
    const read_location &locA = locs[ii];
    const read_location &locB = locs[jj];
    if ( !Stretch( locA, locB, pair, stretch ) )
      continue;
    
    pair2stretch.push_back( make_pair( pair_id, stretch ) );
  }
}

void AllStretches( const vec<read_location> &locs,
		   const vec<int> &flocs,
		   const vec<int> &to_pair,
		   const vec<read_pairing> &pairs,
		   vec< pair<int,float> > &pair2stretch )
{
  pair2stretch.clear( );
  pair2stretch.reserve( pairs.size( ) );

  ForceAssert( is_sorted( locs.begin( ), locs.end( ) ) );
  
  for (int ii=0; ii<(int)flocs.size( ); ii++)
    AllStretches( ii, locs, flocs, to_pair, pairs, pair2stretch, true );

  sort( pair2stretch.begin( ), pair2stretch.end( ) );
}

