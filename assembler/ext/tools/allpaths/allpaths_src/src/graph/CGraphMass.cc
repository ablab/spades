/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "Vec.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "graph/CGraphMass.h"

/**
 * class CGraphMass
 * Constructor
 */
CGraphMass::CGraphMass( const HyperKmerPath *hkp, int cap ) :
  hkp_ ( hkp ),
  cap_ ( cap )
{
  if ( hkp_ ) this->BuildIndexes( );
}

/**
 * class CGraphMass
 * SetPointers
 */
void CGraphMass::SetPointers( const HyperKmerPath *hkp )
{
  hkp_ = hkp;
  this->BuildIndexes( );
}

/**
 * class CGraphMass
 * DirectedKMass
 */
int CGraphMass::DirectedKMass( int edge_id, bool fw ) const
{
  // Add mass of edge (initialize recursion).
  int mass = 0 ;
  vec<int> seen;
  if ( this->AddEdge( edge_id, mass, seen ) ) return mass;
  
  // Add edges "parallel" to edge_id to the seen list.
  int vertex_id = fw ? toright_idx_[edge_id] : toleft_idx_[edge_id];
  vec<int> from_edge_ids = hkp_->FromEdgeObj( vertex_id );
  vec<int> to_edge_ids = hkp_->ToEdgeObj( vertex_id );
  const vec<int> &par_edges = fw ? to_edge_ids : from_edge_ids;
  copy( par_edges.begin( ), par_edges.end( ), back_inserter( seen ) );

  // Add mass of edges connected to edge_id.
  const vec<int> &edge_ids = fw ? from_edge_ids : to_edge_ids;
  for (int ii=0; ii<edge_ids.isize( ); ii++)
    if ( this->KMass( edge_ids[ii], mass, seen ) ) return mass;
  
  // Done, return.
  return mass;
}

/**
 * class CGraphMass
 * BuildIndexes
 */
void CGraphMass::BuildIndexes( )
{
  hkp_->ToLeft( toleft_idx_ );
  hkp_->ToRight( toright_idx_ );
}

/**
 * class CGraphMass
 * KMass
 *
 * This can rapidly become very inefficient (if there are many short
 * unipaths, or cap_ is too large).
 */
bool CGraphMass::KMass( int edge_id, int &mass, vec<int> &seen ) const
{
  // Leave now if mass is already over cap_.
  if ( mass >= cap_ ) return true;
  
  // This edge has already been added to the current mass.
  vec<int>::iterator it = find( seen.begin( ), seen.end( ), edge_id );
  bool skip = ( it != seen.end( ) );
  
  // All edges connected with this edge (both fw and bw).
  vec<int> fw_edges = FindAllEdges( toright_idx_[edge_id] );
  vec<int> bw_edges = FindAllEdges( toleft_idx_[edge_id] );

  // Select edges to be added.
  vec<int> edges;
  edges.reserve( fw_edges.size( ) + bw_edges.size( ) );
  for (int direction=0; direction<2; direction++) {
    const vec<int> &current = ( direction == 0 ) ? fw_edges : bw_edges;
    for (int ii=0; ii<current.isize( ); ii++) {
      vec<int>::iterator it = find( seen.begin( ), seen.end( ), current[ii] );
      if ( it == seen.end( ) )
	edges.push_back( current[ii] );
    }
  }
  
  // Add mass of each of the connecting edges (recursive).
  for (int ii=0; ii<edges.isize( ); ii++) {
    if ( mass >= cap_ ) return true;
    if ( this->AddEdge( edges[ii], mass, seen ) ) return true;
    if ( this->KMass( edges[ii], mass, seen ) ) return true;
  }
  
  // Return.
  return ( mass >= cap_ );
}

/**
 * class CGraphMass
 * AddEdge
 */
bool CGraphMass::AddEdge( int edge_id, int &mass, vec<int> &seen ) const
{
  // This edge has already been added to the current mass.
  vec<int>::iterator it = find( seen.begin( ), seen.end( ), edge_id );
  if ( it != seen.end( ) ) return ( mass >= cap_ );

  // Add edge and return.
  int edge_kmass = hkp_->Edges( )[edge_id].KmerCount( );
  mass = Min( cap_, mass + edge_kmass );
  seen.push_back( edge_id );
  return ( mass >= cap_ );
}

/**
 * class CGraphMass
 * FindAllEdges
 *
 * Find all edges From and To the given vertex_id.
 */
vec<int> CGraphMass::FindAllEdges( int vertex_id ) const
{
  vec<int> edges;

  const vec<int> from_ids = hkp_->FromEdgeObj( vertex_id );
  const vec<int> to_ids = hkp_->ToEdgeObj( vertex_id );
  edges.reserve( from_ids.size( ) + to_ids.size( ) );

  for (int direction=0; direction<2; direction++) {
    const vec<int> &ids = ( direction == 0 ) ? from_ids : to_ids;
    for (int ii=0; ii<ids.isize( ); ii++)
      edges.push_back( ids[ii] );
  }

  sort( edges.begin( ), edges.end( ) );
  edges.erase( unique( edges.begin( ), edges.end( ) ), edges.end( ) );
  
  return edges;
}

