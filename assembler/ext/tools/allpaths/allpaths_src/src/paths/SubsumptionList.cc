// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/SubsumptionList.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

template class SmallVec< BriefSubsumptionRecord, MempoolAllocator<BriefSubsumptionRecord> >;
template class OuterVec<BriefSubsumptionRecordVec>;

// Functions to translate pathIds to indexes and back.

OrientedKmerPathId  PathIdFromIndex( int index ) 
{
  return OrientedKmerPathId( index/2, index&1 ); 
}
  
int IndexFromPathId( const OrientedKmerPathId& pathId )
{
  return pathId.GetId()*2 + pathId.IsRc(); 
}

void
BriefSubsumptionRecord::Write( ostream& out ) const
{
  m_superPath.Write( out );
  BinWrite( out, m_leftOverhang );
}


void
BriefSubsumptionRecord::Read( istream& in ) 
{
  m_superPath.Read( in );
  BinRead( in, m_leftOverhang );
}


void
SubsumptionList::SetBriefRecordsFor( const OrientedKmerPathId& subsumedRead,
                                     const vec<BriefSubsumptionRecord>& records )
{
  BriefSubsumptionRecordVec& ownedVec = m_list[ IndexFromPathId( subsumedRead ) ];

  ownedVec.clear();
  ownedVec.reserve( records.size() );
  copy( records.begin(), records.end(),
        back_inserter( ownedVec ) );
}


void
SubsumptionList::GetFullRecordsFor( const OrientedKmerPathId& subsumedRead,
                                    vec<SubsumptionRecord>& records ) const
{
  const BriefSubsumptionRecordVec& briefRecords = m_list[ IndexFromPathId( subsumedRead ) ];
  records.clear();
  records.reserve( briefRecords.size() );
  for ( BriefSubsumptionRecordVec::const_iterator briefIter = briefRecords.begin();
        briefIter != briefRecords.end(); ++briefIter )
    records.push_back( SubsumptionRecord( subsumedRead, *briefIter ) );
}


void
SubsumptionList::Write( const String& filename ) const
{
  m_list.WriteAll( filename );
}


void
SubsumptionList::Read( const String& filename )
{
  m_list.ReadAll( filename );
}
