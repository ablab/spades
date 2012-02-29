///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PrettyPrintTable.h"
#include "ReadLocationLG.h"
#include "LocsHandlerLG.h"
#include "STLExtensions.h"
#include "String.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"

/**
 * LocsHandlerLG
 * Constructor
 */
LocsHandlerLG::LocsHandlerLG( ) :
  rlens_ ( 0 ),
  clens_ ( 0 ),
  to_rc_ ( 0 ),
  locs_ ( 0 )
{ }

/**
 * LocsHandlerLG
 * Constructor
 */
LocsHandlerLG::LocsHandlerLG( const int K, const String &reads_head )
{
  this->LoadFromKmers( K, reads_head );
  this->SetIndices( );
}

/**
 * LocsHandlerLG
 * Constructor
 */
LocsHandlerLG::LocsHandlerLG( const vec<int> *rlens,
		      const vec<int> *clens,
		      const vec<int> *to_rc,
		      const vec<ReadLocationLG> *locs ) :
  rlens_ ( rlens ),
  clens_ ( clens ),
  to_rc_ ( to_rc ),
  locs_ ( locs )
{
  this->SetIndices( );
}

/**
 * LocsHandlerLG
 * LoadFromKmers
 */
void LocsHandlerLG::LoadFromKmers( const int K, const String &reads_head )
{
  // File names.
  String strK = ToString( K );
  String paths_file = reads_head + ".paths.k" + strK;
  String unipaths_file = reads_head + ".unipaths.k" + strK;
  String unipathsdb_file = reads_head + ".unipathsdb.k" + strK;
  String unilocs_file = reads_head + ".unilocs." + strK + ".10.1";

  // Mapped data.
  vecKmerPath paths(paths_file);
  vecKmerPath unipaths(unipaths_file);
  vec<tagged_rpint> unipathsdb;
  BinaryRead3(unipathsdb_file, unipathsdb);

  // Load core_rlens_, core_clens_, and core_to_rc_.
  core_rlens_.resize( paths.size( ) );
  for (size_t ii=0; ii<paths.size( ); ii++)
    core_rlens_[ii] = paths[ii].KmerCount( );
  
  core_clens_.resize( unipaths.size( ) );
  for (size_t ii=0; ii<unipaths.size( ); ii++)
    core_clens_[ii] = unipaths[ii].KmerCount( );
  
  UnipathInvolution( unipaths, unipathsdb, core_to_rc_ );

  // Load and sort locs.
  BinaryRead2( unilocs_file, core_locs_ );
  if ( ! is_sorted( core_locs_.begin( ), core_locs_.end( ) ) )
    sort( core_locs_.begin( ), core_locs_.end( ) );
  
  // Assign pointers.
  rlens_ = &core_rlens_;
  clens_ = &core_clens_;
  to_rc_ = &core_to_rc_;
  locs_ = &core_locs_;
}

/**
 * LocsHandlerLG
 * BracketOnRead
 *
 * Find bracket on read of overlap between read and contig. This is
 * returned as [begin,end), where end is in STL style (not stop).
 */
pair<int,int> LocsHandlerLG::BracketOnRead( longlong loc_id ) const
{
  const ReadLocationLG &rloc = (*locs_)[loc_id];
  int start = rloc.Start( );
  int contig_len = (*clens_)[rloc.Contig( )];
  int read_len = (*rlens_)[rloc.ReadId( )];
  
  int begin = Max( - start, 0 );
  int end = Min( contig_len - start, read_len );

  return make_pair( begin, end );
}

/**
 * LocsHandlerLG
 * BracketOnContig
 *
 * Find bracket on contig of overlap between read and contig. This is
 * returned as [begin,end), where end is in STL style (not stop).
 */
pair<int,int> LocsHandlerLG::BracketOnContig( longlong loc_id ) const
{
  const ReadLocationLG &rloc = (*locs_)[loc_id];
  int start = rloc.Start( );
  int contig_len = (*clens_)[rloc.Contig( )];
  int read_len = (*rlens_)[rloc.ReadId( )];
  
  int begin = Max( 0, start );
  int end = Min( contig_len, read_len + start );
  
  return make_pair( begin, end );
}

/**
 * GetFwChain
 *
 * A chain of locs for a given read is a (non overlapping and)
 * consistent partitioning of the read in chunks, where "consistent"
 * is defined as follows. Say that the read is partitioned in locs as:
 *
 *         ...[__)[____)...[________)       read
 *
 * where heach interval is a loc in the chain (the dots represent
 * portions of the read for which there is no loc, ie "holes"). Then
 * the contigs spcified by the locs either extend on the left (or on
 * the right) the first (or the last) loc, or match exactly the loc on
 * the read, as in:
 *
 *         ...[__)[____)...[________)       read
 *         ...[__)[____)...[___________)    contig
 *
 * Returns false on error (if, for example, the chain is not consistent).
 */
bool LocsHandlerLG::GetFwChain( longlong read_id, vec<longlong> &locpos ) const
{
  // The set of locs for read_id.
  locpos = this->GetFwPlacements( read_id );

  // No locs found (not an error).
  if ( locpos.size( ) < 1 )
    return true;
  
  // Match each loc in the chain with the bracket on the read (sort locpos).
  {
    vec< pair< pair<int,int>, longlong > > brackets2locpos;
    brackets2locpos.reserve( locpos.size( ) );
    for (longlong ii=0; ii<(longlong)locpos.size( ); ii++) {
      pair<int,int> brack = BracketOnRead( locpos[ii] );
      brackets2locpos.push_back( make_pair( brack, locpos[ii] ) );
    }
    sort( brackets2locpos.begin( ), brackets2locpos.end( ) );
    for (longlong ii=0; ii<(longlong)locpos.size( ); ii++)
      locpos[ii] = brackets2locpos[ii].second;
  }
  
  // Get brackets on reads and contigs.
  vec< pair<int,int> > rBracks;
  vec< pair<int,int> > cBracks;
  rBracks.reserve( locpos.size( ) );
  cBracks.reserve( locpos.size( ) );

  for (longlong ii=0; ii<(longlong)locpos.size( ); ii++) {
    const ReadLocationLG &rl = (*locs_)[ locpos[ii] ];
    rBracks.push_back( this->BracketOnRead( locpos[ii] ) );
    cBracks.push_back( this->BracketOnContig( locpos[ii] ) );
  }

  // Test non-overlapping.
  bool is_valid = true;
  for (longlong ii=0; ii<(longlong)rBracks.size( )-1; ii++) {
    if ( rBracks[ii].second > rBracks[ii+1].first ) {
      is_valid = false;
      break;
    }
  }

  // Test consistency: 1. One bracket only, but there are holes.
  if ( is_valid && cBracks.size( ) == 1 ) {
    if ( rBracks[0].first > 0 ) {
      if ( cBracks[0].first > 0 )
	is_valid = false;
    }
    else if ( rBracks[0].second < (*rlens_)[read_id] ) {
      int contig_id = (*locs_)[ locpos[0] ].Contig( );
      if ( cBracks[0].second < (*clens_)[contig_id] )
	is_valid = false;
    }
  }
  
  // 2. Two brackets or more.
  if ( is_valid ) {
    for (longlong ii=0; ii<cBracks.isize( ); ii++) {
      int contig_id = (*locs_)[ locpos[ii] ].Contig( );

      bool read_begin = ( ii == 0 && rBracks[0].first  == 0 );
      if ( ( ! read_begin ) && cBracks[ii].first != 0 ) {
	is_valid = false;
	break;
      }
      
      int brack_end = rBracks[ rBracks.size( )-1 ].second;      
      bool is_last = ( ii == cBracks.isize( ) - 1 );
      bool read_end = ( is_last && brack_end == (*rlens_)[read_id] );
      if ( ( ! read_end ) && cBracks[ii].second != (*clens_)[contig_id] ) {
	is_valid = false;
	break;
      }
    }
  }
  
  // Return.
  return is_valid;

}

/**
 * LocsHandlerLG
 * PrintChain
 *
 * It does not check if the chain is valid.
 */
void LocsHandlerLG::PrintChain( vec<longlong> &locpos, ostream &out ) const
{
  if ( locpos.size( ) < 1 ) {
    out << "Empty chain: skip printing.\n";
    return;
  }
  
  vec< vec<String> > table;
  vec<String> line;
  longlong read_id = (*locs_)[ locpos[0] ].ReadId( );
  
  vec< pair<int,int> > rBracks;
  vec< pair<int,int> > cBracks;
  rBracks.reserve( locpos.size( ) );
  cBracks.reserve( locpos.size( ) );

  for (longlong ii=0; ii<(longlong)locpos.size( ); ii++) {
    const ReadLocationLG &rl = (*locs_)[ locpos[ii] ];
    rBracks.push_back( this->BracketOnRead( locpos[ii] ) );
    cBracks.push_back( this->BracketOnContig( locpos[ii] ) );
  }

  for (int level=0; level<3; level++) {
    line.clear( );
    for (longlong ii=0; ii<(longlong)locpos.size( ); ii++) {
      const ReadLocationLG &rl = (*locs_)[ locpos[ii] ];
      const pair<int,int> &bra = level == 0 ? rBracks[ii] : cBracks[ii];
      longlong read_id = rl.ReadId( );
      String str1 = ToString( bra.first );
      String str2 = ToString( bra.second );
      String entry = "[" + str1 + "," + str2 + ")";
      if ( level == 0 ) {
	if ( ii == (longlong)locpos.size( )-1 )
	  entry += "_" + ToString( (*rlens_)[read_id] );
      }
      else if ( level == 1 )
	entry += "_" + ToString( (*clens_)[rl.Contig( )] );
      else
	entry = "c" + ToString( rl.Contig( ) );
      line.push_back( entry );
    }
    table.push_back( line );
  }
  
  out << "read_" << read_id << "\n";
  BeautifyAndPrintTable( table, out );
}

/**
 * LocsHandlerLG
 * GetAllPlacements
 *
 * Notice this makes a copy of the vector of locs of read_id.
 */
vec<longlong> LocsHandlerLG::GetAllPlacements(const longlong read_id) const
{
  vec<longlong> all_locs;
  const Int64Vec & iv = to_loc_[read_id];
  all_locs.reserve(iv.size());
  for (Int64Vec::const_iterator iv_it = iv.begin(); iv_it != iv.end(); iv_it++)
    all_locs.push_back(*iv_it);

  return all_locs;
}

/**
 * LocsHandlerLG
 * GetFwPlacements
 */
vec<longlong> LocsHandlerLG::GetFwPlacements(const longlong read_id) const
{
  vec<longlong> fw_locs;
  const Int64Vec & iv = to_loc_[read_id];
  fw_locs.reserve(iv.size() / 2);
  for (Int64Vec::const_iterator iv_it = iv.begin(); iv_it != iv.end(); iv_it++)
    if ((*locs_)[*iv_it].Fw())
      fw_locs.push_back(*iv_it);
			
  return fw_locs;
}

/**
 * LocsHandlerLG
 * GetRcPlacements
 */
vec<longlong> LocsHandlerLG::GetRcPlacements(const longlong read_id) const
{
  vec<longlong> rc_locs;
  const Int64Vec & iv = to_loc_[read_id];
  rc_locs.reserve(iv.size() / 2);
  for (Int64Vec::const_iterator iv_it = iv.begin(); iv_it != iv.end(); iv_it++)
    if ((*locs_)[*iv_it].Rc())
      rc_locs.push_back(*iv_it);
			
  return rc_locs;
}

/**
 * LocsHandlerLG
 * SetIndices
 */
void LocsHandlerLG::SetIndices( )
{
  // Generate to_loc_.
  to_loc_.resize( rlens_->size() );
  for ( size_t i = 0; i < locs_->size(); i++ )
    // Good memory management: Increment SerfVec capacity by 2 - since the
    // number of locs per read is almost always even and is usually 2.
    to_loc_[ (*locs_)[i].ReadId()  ].push_back( i, 1, 2 );
}

