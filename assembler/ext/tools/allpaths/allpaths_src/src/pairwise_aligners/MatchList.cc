/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "pairwise_aligners/MatchList.h"

bool MatchList::FindMatchInSorted( const int id1, const Match& newMatch ) const 
{
  const vec<Match>& sortedMatches = m_sortedMatches[id1];

  pair<vec<Match>::const_iterator,vec<Match>::const_iterator> range;
  range = equal_range( sortedMatches.begin(), sortedMatches.end(), newMatch );
  
  // We can use QuickContains here because the equal_range ensures that all the
  // elements in the range have the same id2, rc, and offset.
  for ( ; range.first != range.second; ++range.first )
    if ( range.first->QuickContains( newMatch ) )
      return true;

  return false;
}

bool MatchList::FindMatchInUnsorted( const int id1, const Match& newMatch ) const 
{
  const vec<Match>& unsortedMatches = m_unsortedMatches[id1];

  for ( vec<Match>::const_iterator matchIter = unsortedMatches.begin();
        matchIter != unsortedMatches.end(); ++matchIter )
    if ( matchIter->Contains( newMatch ) )
      return true;

  return false;
}

void MatchList::UnsortedToSorted( const int id1 ) 
{
  vec<Match>& sortedMatches = m_sortedMatches[id1];
  vec<Match>& unsortedMatches = m_unsortedMatches[id1];

  if ( unsortedMatches.empty() )
    return;

  unsigned int sortedSize = sortedMatches.size();
  unsigned int unsortedSize = unsortedMatches.size();

  sort( unsortedMatches.begin(), unsortedMatches.end() );

  if ( sortedMatches.empty() )
    // This has the side effect of clearing the unsorted vector.
    sortedMatches.swap( unsortedMatches );
  
  else {
    // We do a merge sort of the two sorted vectors.

    if ( sortedMatches.capacity() < sortedMatches.size() + unsortedMatches.size() ) {
      int newCapacity = 2 * max( sortedMatches.size(), unsortedMatches.size() );
      sortedMatches.reserve( newCapacity );
    }

    int oldSize = sortedMatches.size();
    int newSize = oldSize + unsortedMatches.size();

    sortedMatches.resize( newSize );
    
    int newIndex = newSize - 1;
    int oldSortedIndex = oldSize - 1;
    int oldUnsortedIndex = unsortedMatches.size() - 1;
    
    while ( newIndex >= 0 ) {
      if ( oldSortedIndex >= 0 ) 
        if ( oldUnsortedIndex >= 0 ) 
          // In this case, we need to compare the two matches and copy the greater one
          // first (since we're going backwards through the vectors).
          if ( sortedMatches[oldSortedIndex] > unsortedMatches[oldUnsortedIndex] ) 
            sortedMatches[newIndex--] = sortedMatches[oldSortedIndex--];
          else
            sortedMatches[newIndex--] = unsortedMatches[oldUnsortedIndex--];
        else { // oldUnsortedIndex < 0
          // In this case, we've finished with the new "unsorted"
          // matches, so we're just copying the pre-existing matches.
          // We can stop if we'd just be copying values onto
          // themselves.
          if ( newIndex == oldSortedIndex )
            break;
          sortedMatches[newIndex--] = sortedMatches[oldSortedIndex--];
        }
      else // oldSortedIndex < 0 
        // Lastly, if we've finished copying all the pre-existing matches and there are
        // still some new "unsorted" matches left, copy them.
        sortedMatches[newIndex--] = unsortedMatches[oldUnsortedIndex--];
    }

    // All the Matches in the unsorted vector are now in the sorted vector.
    unsortedMatches.clear();
  }
}  

void MatchList::AddMatch( const int id1, const Match& newMatch ) 
{
  vec<Match>& unsortedMatches = m_unsortedMatches[id1];

  unsortedMatches.push_back( newMatch );

  if ( unsortedMatches.size() == m_sortedBatchSize )
    this->UnsortedToSorted( id1 );
}

void MatchList::ProcessMatchingKmers( int pos1, int pos2, int len,
                                      int id1, int id2,
                                      const basevector* pSeq1, const basevector* pSeq2 )
{
  // Is the match rc?
  Bool rc = (pos1 < 0) ^ (pos2 < 0);

  // Decode the pos values.
  if ( pos1 < 0 ) pos1 = -pos1;
  if ( pos2 < 0 ) pos2 = -pos2;
  --pos1;
  --pos2;

  // Where should the match be stored?  The idea here is that we
  // should store matches under the shorter sequence, which is likely
  // to have fewer matches overall.  In the event that both sequences
  // are the same length, then if the difference between the ids is
  // odd, we store it under whichever id is greater; even, under
  // whichever is less.

  /*
   * bool swapSeqs = ( pSeq1->size() > pSeq2->size() );
   *
   * if ( ! swapSeqs && pSeq1->size() == pSeq2->size() )
   */
  bool swapSeqs = false;
  
    // We have to be careful here, since modulo on negative numbers is weird.
    if ( id2 > id1 )
      swapSeqs = ( (id2-id1)%2 == 1 );
    else if ( id1 > id2 )
      swapSeqs = ( (id1-id2)%2 == 0 );
    else // id1 == id2
      swapSeqs = ( pos1 > pos2 );

  if ( swapSeqs ) {
    swap( pos1, pos2 );
    swap( id1, id2 );
    swap( pSeq1, pSeq2 );
  }

  if (rc) pos2 = pSeq2->size( ) - len - pos2;

  Match newMatch( id2, rc, pos1, pos2, len );

  // Check sorted first [O(log n)].
  if ( this->FindMatchInSorted( id1, newMatch ) )
    return;

  // Check unsorted last [O(n)].
  if ( this->FindMatchInUnsorted( id1, newMatch ) )
    return;

  // If it's not in either, then compute its maximal extension and store that
  // match.
  int errors;
  if (!rc) 
    MaxMutmerFromMer( pos1, pos2, len, errors, *pSeq1, *pSeq2, True );
  else
    MaxMutmerFromMerRev( pos1, pos2, len, errors, *pSeq1, *pSeq2, True );

  this->AddMatch( id1, Match( id2, rc, pos1, pos2, len ) );
}

void MatchList::GetMatches( const int id1, vec<Match>& matches ) const 
{
  const vec<Match>& sortedMatches = m_sortedMatches[id1];
  const vec<Match>& unsortedMatches = m_unsortedMatches[id1];

  matches.clear();
  matches.reserve( sortedMatches.size() + unsortedMatches.size() );
  copy( sortedMatches.begin(), sortedMatches.end(),
        back_inserter( matches ) );
  copy( unsortedMatches.begin(), unsortedMatches.end(),
        back_inserter( matches ) );
}

void MatchList::GetMatchCounts( vec<int>& matchCounts ) const 
{
  matchCounts.clear();
  matchCounts.resize( m_sortedMatches.size() );
  
  for ( unsigned int i = 0; i < m_sortedMatches.size(); ++i )
    matchCounts[i] = m_sortedMatches[i].size() + m_unsortedMatches[i].size();
}

