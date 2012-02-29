// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/MuxGraph.h"

#include "STLExtensions.h"

#include <set>


// Functions to translate pathIds to nodeIds and back.

OrientedKmerPathId  PathIdFromNodeId( int nodeId ) 
{
  return OrientedKmerPathId( nodeId/2, nodeId&1 ); 
}
  
int  NodeIdFromPathId( const OrientedKmerPathId& pathId )
{
  return pathId.GetId()*2 + pathId.IsRc(); 
}


int
MuxGraph::size() const
{
  return m_singleEdgeNodes.size() / 2;
}


// Accomodate the forward and reverse versions of numReads paths.
void 
MuxGraph::resize( const int numReads )
{
  m_singleEdgeNodes.resize( 2*numReads, Mux( OrientedKmerPathId(), s_empty, 0 ) );
}


// Methods to get/set muxs.

void
MuxGraph::GetMuxesOf( const OrientedKmerPathId& pathId, vec<Mux>& muxes ) const
{
  muxes.clear();

  const int nodeId = NodeIdFromPathId( pathId );
  const Mux& aMux = m_singleEdgeNodes[ nodeId ];

  if ( ! IsSpecial( aMux ) )
    muxes.resize( 1, aMux );

  else if ( ! IsEmpty( aMux ) )
  {
    int multiStart, numMuxes;
    GetFromMulti( aMux, multiStart, numMuxes );

    ForceAssertGt( numMuxes, 1 );
    muxes.reserve( numMuxes );
    muxes.insert( muxes.end(),
                  m_multiEdgeNodes.begin() + multiStart,
                  m_multiEdgeNodes.begin() + multiStart + numMuxes );
  }
}

int 
MuxGraph::NumMuxesOf( const OrientedKmerPathId& pathId ) const
{
  
  const int nodeId = NodeIdFromPathId( pathId );
  const Mux& aMux = m_singleEdgeNodes[ nodeId ];

  if ( IsSingle( aMux ) )
    return 1;
  else
    return GetNumMuxesFromSpecial( aMux );
}

void
MuxGraph::SetMuxOf( const OrientedKmerPathId& pathId, const Mux& newMux )
{
  Assert( ! IsSpecial( newMux ) );
  
  m_singleEdgeNodes[ NodeIdFromPathId( pathId ) ] = newMux;
}

void 
MuxGraph::SetMuxesOf( const OrientedKmerPathId& pathId, const vec<Mux>& newMuxes )
{
  // No effort is made to remove entries from the multiEdgeNodes table
  // if they exist, as it would not only be grossly inefficient
  // itself, it would also require updating the start position
  // information for all the multi nodes that point into later
  // positions in the multiEdgeNodes table.

  if ( newMuxes.size() == 1 )
  {
    this->SetMuxOf( pathId, newMuxes.front() );
    return;
  }

  const int nodeId = NodeIdFromPathId( pathId );
  Mux& aMux = m_singleEdgeNodes[ nodeId ];

  // Add new edges, if any.
  if ( newMuxes.empty() )
  {
    SetToEmpty( aMux );
  }
  else
  {
    SetToMulti( aMux, m_multiEdgeNodes.size(), newMuxes.size() );

    m_multiEdgeNodes.insert( m_multiEdgeNodes.end(),
                             newMuxes.begin(), newMuxes.end() );
  }
}


void
MuxGraph::Write( const String& filename ) const
{
  BinaryWrite2( filename + ".single", m_singleEdgeNodes );
  BinaryWrite2( filename + ".multi", m_multiEdgeNodes );
}



void
MuxGraph::Read( const String& filename )
{
  m_singleEdgeNodes.clear();
  m_multiEdgeNodes.clear();

  BinaryRead2( filename + ".single", m_singleEdgeNodes );
  BinaryRead2( filename + ".multi", m_multiEdgeNodes );
}

bool
MuxGraph::FilesExist( const String& filename ) const
{
  return( IsRegularFile( filename + ".single" ) &&
	  IsRegularFile( filename + ".multi" ) );
}



bool
MuxGraph::VerifySameAs( const MuxGraph& other ) const
{
  if ( m_singleEdgeNodes.size() != other.m_singleEdgeNodes.size() )
  {
    PRINT2( m_singleEdgeNodes.size(), other.m_singleEdgeNodes.size() );
    return false;
  }

  vec<Mux> thisMuxes, otherMuxes;
  int numNodes = m_singleEdgeNodes.size();
  for ( int i = 0; i < numNodes; ++i )
  {
    OrientedKmerPathId pathId = PathIdFromNodeId( i );
    this->GetMuxesOf( pathId, thisMuxes );
    other.GetMuxesOf( pathId, otherMuxes );

    sort( thisMuxes.begin(), thisMuxes.end() );
    sort( otherMuxes.begin(), otherMuxes.end() );

    if ( thisMuxes != otherMuxes )
    {
      PRINT2( i, pathId );
      PRINT( thisMuxes );
      PRINT( otherMuxes );

      return false;
    }
  }
  
  return true;
}


void PrintNode( ostream& out, const OrientedKmerPathId& okpid, 
                set<OrientedKmerPathId>& seen, const int partition ) {
  if ( ! seen.count( okpid ) ) {
    out << "  " << NodeIdFromPathId( okpid ) 
        << " ["
        << "label=\"" << okpid << "\"";
    if ( okpid.GetId() < partition )
      out << "fill=gray";
    out << "];" << "\n";
    seen.insert( okpid );
  }
}
  

void
MuxGraph::PrintDot( const String& filename, const int partition ) const 
{
  ofstream out( filename.c_str() );

  out << "digraph MuxGraph {" << "\n";
  out << "  rankdir=LR;" << "\n";
  out << "  edge [dir=back];" << "\n";
  
  set<OrientedKmerPathId> seen;

  vec<Mux> muxes;
  for ( int i = 0; i < size(); ++i ) {
    for ( int rc = 0; rc < 2; ++rc ) {
      OrientedKmerPathId thisOkpid( i, rc );
      GetMuxesOf( thisOkpid, muxes );
      if ( ! muxes.empty() )
        PrintNode( out, thisOkpid, seen, partition );
      for ( int m = 0; m < muxes.isize(); ++m ) {
        OrientedKmerPathId muxOkpid( muxes[m].GetPathId() );
        PrintNode( out, muxOkpid, seen, partition );
        out << "  " << NodeIdFromPathId( muxOkpid ) 
            << " -> " << NodeIdFromPathId( thisOkpid )
            << " [label=\"" << muxes[m].GetNumKmers() << "\"];" << "\n";
      }
    }
  }
  
  out << "}" << "\n";
}
