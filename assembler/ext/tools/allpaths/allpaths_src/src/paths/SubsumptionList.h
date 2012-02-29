// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_SUBSUMPTIONLIST_H
#define PATHS_SUBSUMPTIONLIST_H

#include "paths/OrientedKmerPathId.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

class BriefSubsumptionRecord
{
 public:
  BriefSubsumptionRecord()
    : m_leftOverhang( -1 )
  {}
  
  BriefSubsumptionRecord( const OrientedKmerPathId &superPath,
                          const int leftOverhang )
    : m_superPath( superPath ),
      m_leftOverhang( leftOverhang )
  {}
  
  OrientedKmerPathId GetSuperPathId() const { return m_superPath; }
  int GetLeftOverhang() const { return m_leftOverhang; }
  
  void Write( ostream& out ) const;
  void Read( istream& in );
  
  bool operator== ( const BriefSubsumptionRecord& other ) const
  {
    return ( m_superPath == other.m_superPath && m_leftOverhang == other.m_leftOverhang );
  }

  bool operator!= ( const BriefSubsumptionRecord& other ) const
  {
    return ! ( *this == other );
  }
  
 private:
  OrientedKmerPathId m_superPath;
  int m_leftOverhang;
};

TRIVIALLY_SERIALIZABLE(BriefSubsumptionRecord);
typedef SerfVec<BriefSubsumptionRecord> BriefSubsumptionRecordVec;
typedef MasterVec<BriefSubsumptionRecordVec> VecBriefSubsumptionRecordVec;

class SubsumptionRecord
{
 public:
  SubsumptionRecord()
    : m_leftOverhang( -1 )
  {}

  SubsumptionRecord( const OrientedKmerPathId &subPath,
                     const OrientedKmerPathId &superPath,
                     const int leftOverhang )
    : m_subPath( subPath ),
      m_superPath( superPath ),
      m_leftOverhang( leftOverhang )
  {}

  SubsumptionRecord( const OrientedKmerPathId &subPath,
                     const BriefSubsumptionRecord &record )
    : m_subPath( subPath ),
      m_superPath( record.GetSuperPathId() ),
      m_leftOverhang( record.GetLeftOverhang() )
  {}

  OrientedKmerPathId GetSubPathId() const { return m_subPath; }
  OrientedKmerPathId GetSuperPathId() const { return m_superPath; }
  int GetLeftOverhang() const { return m_leftOverhang; }

  void Write( ostream& out ) const;
  void Read( istream& in );

  bool operator== ( const SubsumptionRecord& other ) const
  {
    return ( m_subPath == other.m_subPath &&
             m_superPath == other.m_superPath && 
             m_leftOverhang == other.m_leftOverhang );
  }

  bool operator< ( const SubsumptionRecord& other ) const
  {
    if ( m_subPath < other.m_subPath ) return true;
    if ( other.m_subPath < m_subPath ) return false;
    if ( m_superPath < other.m_superPath ) return true;
    if ( other.m_superPath < m_superPath ) return false;
    return ( m_leftOverhang < other.m_leftOverhang );
  }

 private:
  OrientedKmerPathId m_subPath;
  OrientedKmerPathId m_superPath;
  int m_leftOverhang;
};


class SubsumptionList {
 public:
  SubsumptionList() {}
  
  SubsumptionList( const int nReads )
    : m_list( nReads*2 )
  {}

  void resize( const int nReads ) { m_list.resize( nReads*2 ); }

  void
  SetBriefRecordsFor( const OrientedKmerPathId& subsumedRead,
                      const vec<BriefSubsumptionRecord>& records );

  void
  GetFullRecordsFor( const OrientedKmerPathId& subsumedRead,
                     vec<SubsumptionRecord>& records ) const;

  void Write( const String& filename ) const;
  void Read( const String& filename );

  bool operator== ( const SubsumptionList& other ) const 
  {
    return ( m_list == other.m_list );
  }

 private:
  VecBriefSubsumptionRecordVec m_list;
};

#endif
