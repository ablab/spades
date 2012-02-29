// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_MUX_H
#define PATHS_MUX_H

#include "paths/OrientedKmerPathId.h"

/// A Mux represents one (identified) path extending past the end of
/// another (anonymous) one.  (That is, a Mux can tell you who is
/// doing the extending, but not who is being extended.)  It knows how
/// many kmers and how many segments out the extension reaches.
///
/// For an overview of Mux-based insert walking, see paths/doc/mux.pdf

class Mux {

public:
  Mux()
    : m_segment( -1 )
  {}

  Mux( const OrientedKmerPathId& okpi, int segment, int numKmers )
    : m_pathId   ( okpi ),
      m_segment  ( segment ),
      m_numKmers ( numKmers )
  {}

  OrientedKmerPathId  GetPathId() const { return m_pathId; }
  void SetPathId( const OrientedKmerPathId& id ) { m_pathId = id; }
  
  int  GetSegment() const  { return m_segment; }
  void SetSegment( const int segment )  { m_segment = segment; }
  
  int  GetNumKmers() const { return m_numKmers; }
  void SetNumKmers( const int numKmers ) { m_numKmers = numKmers; }

  bool operator< ( const Mux& other ) const
  {
    return ( m_pathId < other.m_pathId ||
             m_pathId == other.m_pathId && m_segment < other.m_segment );
  }
  
  bool operator== ( const Mux& other ) const
  {
    return ( m_pathId == other.m_pathId && m_segment == other.m_segment ); 
  }

  friend struct OrderBySegment;
  struct OrderBySegment : public binary_function<Mux,Mux,bool>
  {
    bool operator() ( const Mux& lhs, const Mux& rhs ) const
    {
      return ( lhs.m_segment < rhs.m_segment ||
               lhs.m_segment == rhs.m_segment && lhs.m_pathId < rhs.m_pathId );
    }
  };

  friend ostream& operator<< ( ostream& out, const Mux& aMux )
  {
    return out << aMux.m_pathId << "(s" << aMux.m_segment << ",k" << aMux.m_numKmers << ")";
  }

  void Write( ostream& out ) const;
  void Read( istream& in );

private:
  OrientedKmerPathId m_pathId;
  int m_segment;
  int m_numKmers;

};

inline
void
Mux::Write( ostream& out ) const
{
  m_pathId.Write( out );
  BinWrite( out, m_segment );
  BinWrite( out, m_numKmers );
}

inline
void
Mux::Read( istream& in )
{
  m_pathId.Read( in );
  BinRead( in, m_segment );
  BinRead( in, m_numKmers );
}



#endif
