///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "Superb.h"
#include "STLExtensions.h"
#include "VecTemplate.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/reporting/COffset.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * class CSuperLinks
 * Constructor
 */
CSuperLinks::CSuperLinks( const PairsManager *pairs,
			  const vec<superb> *supers,
			  const vec<alignlet> *aligns,
			  const vec<int> *index ) :
  pairs_ ( pairs ),
  supers_ ( supers ),
  aligns_ ( aligns ),
  index_ ( index )
{
  this->GenerateMaps( );
}

/**
 * class CSuperLinks
 * SetPointers
 */
void CSuperLinks::SetPointers( const PairsManager *pairs,
			       const vec<superb> *supers,
			       const vec<alignlet> *aligns,
			       const vec<int> *index )
{
  pairs_ = pairs;
  supers_ = supers;
  aligns_ = aligns;
  index_ = index;
  
  this->GenerateMaps( );
}

/**
 * class CSuperLinks
 * AllLinks
 */
COffset CSuperLinks::AllLinks( int super1, int super2, bool rc2 ) const
{
  COffset empty_offset;
  vec<COffset> offsets;
  this->AllLinks( super1, offsets );
  for (int ii=0; ii<offsets.isize( ); ii++) {
    if ( offsets[ii].Super2( ) != super2 ) continue;
    if ( offsets[ii].Rc2( ) != rc2 ) continue;
    return offsets[ii];   // Found valid links to super2 (rc2): return.
  }

  // No valid link found.
  return empty_offset;
}

/**
 * class CSuperLinks
 * AllLinks
 *
 * slop: if >0, only consider a link if the ends of the two reads are
 *       closer to the edge of the super than ( sep * slop )
 * stretch: if not NULL, pass this arg to ClusterLinks( )
 */
void CSuperLinks::AllLinks( int super_id,
			    vec<COffset>& out,
			    float slop,
			    float *stretch ) const
{
  COffset tmp(super_id,0,false,true_end_[super_id]-true_begin_[super_id],0);
  set<COffset> offSet;
  superb const& sup = (*supers_)[super_id];
  for ( int cgpos = 0; cgpos < sup.Ntigs(); cgpos++ ) {
    typedef IntVec::const_iterator Itr;
    IntVec const& ids = read_ids_[sup.Tig(cgpos)];
    for ( Itr itr(ids.begin()), end(ids.end()); itr != end; ++itr ) {
      int read_id = *itr;
      int p_read_id = pairs_->getPartnerID(read_id);
      int idx = (*index_)[p_read_id];
      if ( idx < 0 )
	continue;
      
      alignlet const& alignlet2 = (*aligns_)[idx];
      int p_super_id = super_id_[alignlet2.TargetId()];
      if ( p_super_id < 0 || super_id == p_super_id )
	continue;
      
      alignlet const& alignlet1 = (*aligns_)[(*index_)[read_id]];
      tmp.SetSupers(tmp.Super1(),
		    p_super_id,
		    alignlet1.Fw1()==alignlet2.Fw1(),
		    tmp.Slen1(),
		    true_end_[p_super_id]-true_begin_[p_super_id]);
      
      // const_cast is legitimate because adding links doesn't alter order
      COffset& cOffset = const_cast<COffset&>(*offSet.insert(tmp).first);
      
      int offset = super_begin_[alignlet1.TargetId()];
      ho_interval w1(offset+alignlet1.pos2(), offset+alignlet1.Pos2());
      offset = super_begin_[alignlet2.TargetId()];
      ho_interval w2(offset+alignlet2.pos2(), offset+alignlet2.Pos2());
      
      longlong pair_id = pairs_->getPairID(read_id);
      int sep = pairs_->sep( pair_id );
      int sd = pairs_->sd( pair_id );
      if (slop >= 0) {
	int rc1 = alignlet1.Fw1() ? 0 : 1;
	int rc2 = alignlet2.Fw1() ? 0 : 1;
	int dist_to_end1 = rc1 ? w1.Start() : (tmp.Slen1() - w1.Stop());
	int dist_to_end2 = rc2 ? w2.Start() : (tmp.Slen2() - w2.Stop());
	int max_dist = sep + sd*slop;
	if (dist_to_end1 > max_dist) continue;
	if (dist_to_end2 > max_dist) continue;
      }
      int raw_offset
	= alignlet1.Fw1()
	? w1.Stop()+sep-(tmp.Rc2()?(tmp.Slen2()-w2.Stop()):w2.Start())
	: w1.Start()-sep-(tmp.Rc2()?(tmp.Slen2()-w2.Start()):w2.Stop());
      normal_distribution nd( raw_offset, sd );
      cOffset.AddLink( SLink( nd, w1, w2, pair_id ) );
    }
  }
  
  typedef set<COffset>::iterator SItr;
  for ( SItr itr(offSet.begin()), end(offSet.end()); itr != end; ++itr )
    stretch ? itr->ClusterLinks( *stretch ) : itr->ClusterLinks( );
  
  out.assign(offSet.begin(),offSet.end());
}

/**
 * class CSuperLinks
 * PrintAllLinks
 */
void CSuperLinks::PrintAllLinks( ostream &out, int super_id, bool full ) const
{
  vec<COffset> offsets;
  this->AllLinks( super_id, offsets );
  const PairsManager *ppairs = full ? pairs_ : NULL;
  for (size_t ii=0; ii<offsets.size( ); ii++)
    offsets[ii].Print( out, ppairs );
}

/**
 * class CSuperLinks
 * AllPairs
 */
vec<int> CSuperLinks::AllPairs( int super_id ) const
{
  vec<int> pids;
  
  const superb &sup = (*supers_)[super_id];
  for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
    int edge_id = sup.Tig( cgpos );
    for (int ii=0; ii<(int)read_ids_[edge_id].size( ); ii++) {
      int read_id = read_ids_[edge_id][ii];
      int partner_id = pairs_->getPartnerID( read_id );
      if ( (*index_)[partner_id] < 0 ) continue;
      int p_edge_id = (*aligns_)[ (*index_)[partner_id] ].TargetId( );
      int p_super_id = super_id_[p_edge_id];
      if ( p_super_id < 0 ) continue;
      if ( super_id != p_super_id )
	pids.push_back( pairs_->getPairID( read_id ) );
    }
  }
  
  return pids;
}

/**
 * class CSuperLinks
 * AddLink
 *
 * It asserts if pair_id does not link the supers in offset.
 */
void CSuperLinks::AddLink( int pair_id, COffset &offset ) const
{
  int id1 = pairs_->ID1( pair_id );
  int id2 = pairs_->ID2( pair_id );
  int s1 = this->SuperId( id1 );
  int s2 = this->SuperId( id2 );
  if ( s1 != offset.Super1( ) ) {
    swap( id1, id2 );
    swap( s1, s2 );
  }
  ForceAssertEq( s1, offset.Super1( ) );
  ForceAssertEq( s2, offset.Super2( ) );

  // Windows of reads on their supers.
  pair<int,int> win1 = this->WinOnSuper( id1 );
  pair<int,int> win2 = this->WinOnSuper( id2 );
  int truelen1 = true_end_[s1] - true_begin_[s1];
  int truelen2 = true_end_[s2] - true_begin_[s2];
  int beg1 = win1.first;
  int end1 = win1.second;
  int beg2 = win2.first;
  int end2 = win2.second;

  // Orientation of s2 with respect to s1.
  bool rc2 = ( this->FwOnSuper( id1 ) == this->FwOnSuper( id2 ) );

  // Offset and stdev.
  int sep = pairs_->sep( pair_id );
  int stdev = pairs_->sd( pair_id );
  int raw_offset
    = this->FwOnSuper( id1 )
    ? end1 + sep - ( rc2 ? ( truelen2 - end2 ) : beg2 )
    : beg1 - sep - ( rc2 ? ( truelen2 - beg2 ) : end2 );
  normal_distribution nd( raw_offset, stdev );
  
  // Windows.
  ho_interval w1( beg1, end1 );
  ho_interval w2( beg2, end2 );
  
  // Add link and return.
  offset.AddLink( SLink( nd, w1, w2, pair_id ) );
  
}

/**
 * class CSuperLinks
 * GenerateMaps
 * private
 */
void CSuperLinks::GenerateMaps( )
{
  super_id_.clear( );
  super_pos_.clear( );
  super_begin_.clear( );
  super_end_.clear( );
  read_ids_.clear( );

  true_begin_.clear( );
  true_end_.clear( );

  int n_edges = 0;
  for (int ii=0; ii<supers_->isize( ); ii++)
    for (int jj=0; jj<(*supers_)[ii].Ntigs( ); jj++)
      n_edges = Max( n_edges, (*supers_)[ii].Tig( jj ) );
  for (int read_id=0; read_id<(int)index_->size( ); read_id++) {
    if ( (*index_)[read_id] < 0 ) continue;
    n_edges = Max( n_edges, (*aligns_)[ (*index_)[read_id] ].TargetId( ) );
  }
  n_edges++;
  
  super_id_.resize( n_edges, -1 );
  super_pos_.resize( n_edges, -1 );
  for (int ii=0; ii<supers_->isize( ); ii++) {
    for (int jj=0; jj<(*supers_)[ii].Ntigs( ); jj++) {
      int edge_id = (*supers_)[ii].Tig( jj );
      super_id_[edge_id] = ii;
      super_pos_[edge_id] = jj;
    }
  }
  
  super_begin_.resize( n_edges, 0 );
  super_end_.resize( n_edges, 0 );
  for (int ii=0; ii<(int)supers_->size( ); ii++) {
    int pos = 0;
    for (int jj=0; jj<(int)(*supers_)[ii].Ntigs( ); jj++) {
      super_begin_[ (*supers_)[ii].Tig( jj ) ] = pos;
      pos += (*supers_)[ii].Len( jj );
      super_end_[ (*supers_)[ii].Tig( jj ) ] = pos;
      if ( jj < (int)(*supers_)[ii].Ntigs( ) - 1 )
	pos += (*supers_)[ii].Gap( jj );
    }
  }
  
  read_ids_.resize( n_edges );
  for (int read_id=0; read_id<(int)index_->size( ); read_id++) {
    if ( (*index_)[read_id] < 0 ) continue;
    const alignlet &al = (*aligns_)[ (*index_)[read_id] ];
    read_ids_[ al.TargetId( ) ].push_back( read_id );
  }
  
  true_begin_.resize( supers_->size( ), 0 );
  true_end_.resize( supers_->size( ), 0 );
  for (int ii=0; ii<(int)supers_->size( ); ii++) {
    true_begin_[ii] = (*supers_)[ii].TrueBegin( );
    true_end_[ii] = (*supers_)[ii].TrueEnd( );
  }

}

/**
 * class CSuperLinks
 * SuperId
 * private
 */
int CSuperLinks::SuperId( int read_id ) const
{
  const alignlet &al = (*aligns_)[ (*index_)[read_id] ];
  return super_id_[ al.TargetId( ) ];
}

/**
 * class CSuperLinks
 * FwOnSuper
 * private
 */
bool CSuperLinks::FwOnSuper( int read_id ) const
{
  const alignlet &al = (*aligns_)[ (*index_)[read_id] ];
  return al.Fw1( );
}

/**
 * class CSuperLinks
 * WinOnSuper
 * private
 */
pair<int,int> CSuperLinks::WinOnSuper( int read_id ) const
{
  const alignlet &al = (*aligns_)[ (*index_)[read_id] ];
  int offset = super_begin_[ al.TargetId( ) ];
  return make_pair( offset + al.pos2( ), offset + al.Pos2( ) );
}

