/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/ExtendUnipathSeqs.h"

#include "paths/UnipathSeqDatabase.h"

#include <set>

void ExtendUnipathSeqs( const vecKmerPath& unipaths,
                        const vecUnipathSeq& unipathSeqs,
                        vecUnipathSeq& extendedUnipathSeqs,
                        vec<Mux>& muxes ) 
{
  ExtendUnipathSeqs( unipaths, unipathSeqs, 
                     unipathSeqs, extendedUnipathSeqs, muxes );
}


void ExtendUnipathSeqs( const vecKmerPath& unipaths,
                        const vecUnipathSeq& unipathSeqs,
                        const vecUnipathSeq& unipathSeqsToExtend,
                        vecUnipathSeq& extendedUnipathSeqs,
                        vec<Mux>& muxes ) 
{
  // cout << "Extending unipath seqs... " << endl;

  UnipathSeqDatabase unipathSeqDb( unipathSeqs );
  
  int numSeqs = unipathSeqsToExtend.size();
  
  int seqsPerDot = 1;
  while ( numSeqs / seqsPerDot > 100 )
    seqsPerDot *= 10;
  
  int numDots = numSeqs / seqsPerDot + 1;

  vecUnipathSeq newUnipathSeqs;
  muxes.clear();

  // cout << "Processing seqs in " << numDots << " passes." << endl;
  
  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqsToExtend.size(); ++i ) {
    // if ( i % seqsPerDot == 0 ) Dot( cout, i / seqsPerDot );
    
    // The left extension will be constructed in reverse order.
    vec<int> leftExtension;
    vec<int> rightExtension;
    
    for ( int pass = 0; pass < 2; ++pass ) {
      if ( unipathSeqsToExtend[i].empty() ) continue;

      set<int> unipathsInExtension;

      vec<int>& extension = ( pass == 0 ? leftExtension : rightExtension );
      int offset = ( pass == 0 ? -1 : 1 );
      
      int unipath = ( pass == 0 ? unipathSeqsToExtend[i].front() : unipathSeqsToExtend[i].back() );
      
      const int k_noUnipath = -1;
      const int k_multipleUnipaths = -2;
      
      while ( unipath >= 0 ) {
        int newUnipath = k_noUnipath;
        
        vec<UnipathSeqDatabase::Record> records;
        unipathSeqDb.Find( unipath, records );
            
        for ( unsigned int j = 0; j < records.size(); ++j ) {
          int seqId = records[j].seqId;
          int index = records[j].index+offset;
          if ( index >= 0 &&
               static_cast<unsigned>(index) < unipathSeqs[ seqId ].size() ) {
            int otherUnipath = unipathSeqs[ seqId ][ index ];
            if ( newUnipath == k_noUnipath )
              newUnipath = otherUnipath;
            else if ( newUnipath != otherUnipath ) {
              newUnipath = k_multipleUnipaths;
              break;
            }
          }
          
          if ( newUnipath == k_multipleUnipaths )
            break;
        }

        if ( newUnipath >= 0 ) {

          // If the only neighboring unipath is itself, then we've
          // entered a terminal loop in the unipath graph and we should
          // break.
          if ( newUnipath == unipath )
            break;

          // If the neighboring unipath has already been encountered,
          // then we've entered a loop and should terminate.  This
          // will only happen if the unipath sequences we're extending
          // with are a proper subset of the sequences used to derive
          // the unipaths.
          pair<set<int>::iterator,bool> insertResult = unipathsInExtension.insert( newUnipath );
          if ( insertResult.second ) 
            extension.push_back( newUnipath );
          else
            break;
        }
        unipath = newUnipath;
      }
    }

    if ( leftExtension.empty() && rightExtension.empty() ) {
      newUnipathSeqs.push_back_reserve( unipathSeqsToExtend[i], 0, 2.0 );
      Mux theMux;
      theMux.SetSegment( 0 );
      theMux.SetNumKmers( 0 );
      muxes.push_back( theMux );
    }
    else {
      UnipathSeq newUnipathSeq;
      int newSize = leftExtension.size() + unipathSeqsToExtend[i].size() + rightExtension.size();
      newUnipathSeq.reserve( newSize );
      
      // Note rbegin and rend here, since leftExtension was constructed backwards.
      copy( leftExtension.rbegin(), leftExtension.rend(),
            back_inserter( newUnipathSeq ) );
      copy( unipathSeqsToExtend[i].begin(), unipathSeqsToExtend[i].end(),
            back_inserter( newUnipathSeq ) );
      copy( rightExtension.begin(), rightExtension.end(),
            back_inserter( newUnipathSeq ) );
      
      newUnipathSeqs.push_back_reserve( newUnipathSeq, 0, 2.0 );
      
      Mux theMux;
      theMux.SetSegment( leftExtension.size() );
      int numKmers = 0;
      for ( unsigned int j = 0; j < leftExtension.size(); ++j )
        numKmers += unipaths[ leftExtension[j] ].KmerCount();
      theMux.SetNumKmers( numKmers );
      muxes.push_back( theMux );
    }
  }

  newUnipathSeqs.swap( extendedUnipathSeqs );

  // cout << endl;
  // cout << "done." << endl;
}
