// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_ORIENTEDKMERPATHID_H
#define PATHS_ORIENTEDKMERPATHID_H

#include "paths/KmerPath.h"

// A tuple of path id and orientation.

class OrientedKmerPathId
{
 public:
  OrientedKmerPathId()
    : m_idRc( INT_MAX ) {}

  explicit
  OrientedKmerPathId( int idRc )
    : m_idRc( idRc ) {}

  OrientedKmerPathId( int id, bool rc )
    : m_idRc( rc ? -id-1 : id ) {}

  bool IsValid() const { return ( m_idRc != INT_MAX ); }

  bool IsFw() const { return m_idRc >= 0; }
  bool IsRc() const { return m_idRc < 0; }
  int  GetId() const { return ( IsRc() ? -m_idRc-1 : m_idRc ); }
  int  GetIdRc() const { return m_idRc; }

  template <class T, class vecT>
  const T* GetPtr( const vecT& fw, const vecT& rc ) const
  {
    return &( IsRc() ? rc[ -m_idRc-1 ] : fw[ m_idRc ] );
  }

  const KmerPath * GetPathPtr( const vecKmerPath &paths, 
                               const vecKmerPath &paths_rc ) const 
  {
    return GetPtr<KmerPath,vecKmerPath>( paths, paths_rc );
  }

  bool operator< ( const OrientedKmerPathId &other ) const
  {
    return ( this->m_idRc < other.m_idRc );
  }

  bool operator== ( const OrientedKmerPathId &other ) const
  {
    return ( this->m_idRc == other.m_idRc );
  }

  void Write( ostream& out ) const
  {
    out.write( (char*)&m_idRc, sizeof(m_idRc) );
  }

  void Read( istream& in )
  {
    in.read( (char*)&m_idRc, sizeof(m_idRc) );
  }
  
 private:
  int m_idRc; // TODO: potentially dangerous truncation of index
};

inline
ostream& operator<< ( ostream& out, const OrientedKmerPathId& okpd )
{
  return out << okpd.GetId() << ( okpd.IsRc() ? "rc" : "fw" );
}


// Define hash<OKPI>, so we can use it as a key in hash maps and sets
#include <sys/types.h>
#if __GNUC__ > 2
namespace __gnu_cxx {
#endif

template <class T> struct hash;

template<>
struct hash<OrientedKmerPathId> {
  size_t 
  operator() ( const OrientedKmerPathId& okpid ) const
  { return okpid.GetIdRc(); }
};

#if __GNUC__ > 2
}
#endif




#endif
