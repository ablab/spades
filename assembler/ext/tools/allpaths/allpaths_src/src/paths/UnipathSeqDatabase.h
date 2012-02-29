/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_UNIPATHSEQDATABASE_H
#define PATHS_UNIPATHSEQDATABASE_H

// This class allows one to quickly find all the occurences of a given
// unipath in a set of UnipathSeqs.

class UnipathSeqDatabase {
 public:
  UnipathSeqDatabase( const vecUnipathSeq& unipathSeqs );
  //TODO: potentially dangerous truncation of ids
  struct Record {
    int unipath;
    int seqId;
    int index;
    
    Record() {}

    Record( int u, int s, int i )
      : unipath(u), seqId(s), index(i) {}
    
    bool operator< ( const Record& other ) const {
      return ( this->unipath < other.unipath );
    }

    void Print( ostream& out, const vecUnipathSeq& unipathSeqs ) const {
      for ( UnipathSeq::size_type i = 0; i < unipathSeqs[seqId].size(); ++i ) {
        if ( i > 0 ) out << ".";
        String unipathName = BaseAlpha( unipathSeqs[seqId][i] );
        if ( i != static_cast<UnipathSeq::size_type>(index) )
          unipathName.ToLower();
        out << unipathName;
      }
    }
  };

  void Find( const int unipath, vec<Record>& records ) const;
  
 private:
  vec<Record> m_records;
};


inline
UnipathSeqDatabase::UnipathSeqDatabase( const vecUnipathSeq& unipathSeqs ) {
  int numRecords = 0;
  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqs.size(); ++i )
    numRecords += unipathSeqs[i].size();

  m_records.reserve( numRecords );
  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqs.size(); ++i )
    for ( UnipathSeq::size_type j = 0; j < unipathSeqs[i].size(); ++j )
      m_records.push_back( Record( unipathSeqs[i][j], i, j ) );

  sort( m_records.begin(), m_records.end() );
}


inline
void UnipathSeqDatabase::Find( const int unipath, vec<UnipathSeqDatabase::Record>& records ) const
{
  pair<vec<Record>::const_iterator,vec<Record>::const_iterator> range;
  
  range = equal_range( m_records.begin(), m_records.end(),
                       Record( unipath, 0, 0 ) );

  records.resize( distance( range.first, range.second ) );
  copy( range.first, range.second, records.begin() );
}

#endif
