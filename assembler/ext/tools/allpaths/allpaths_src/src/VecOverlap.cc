/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "VecOverlap.h"

#include "STLExtensions.h"
#include <numeric>

// Represents the suffix of the index'th word starting at the pos'th
// character.
template <class T>
struct Suffix {
  Suffix()
    : index(-1), pos(-1) {}
  Suffix( int i, int p )
    : index(i), pos(p) {}

  int index;
  int pos;

  int GetLength( const vec< vec<T> >& words ) const
  {
    return words[index].size() - pos;
  }

  T GetValue( const vec< vec<T> >& words, int offset = 0 ) const
  {
    if ( pos+offset >= words[index].isize() ) 
      return -1;
    else
      return words[index][pos+offset];
  }
  
  int GetLongestCommonPrefix( const vec< vec<T> >& words, const Suffix<T>& other ) const
  {
    int lcp = 0;
    while(1) {
      T thisValue = this->GetValue( words, lcp );
      T prevValue = other.GetValue( words, lcp );
      if ( thisValue >= 0 && thisValue == prevValue ) // implies prevValue >= 0, too
        ++lcp;
      else
        break;
    }
    return lcp;
  }

  void Print( ostream& out, const vec< vec<T> >& words ) const {
    copy( words[index].begin() + pos, words[index].end(),
          ostream_iterator<T>( out, "." ) );
    out << "\n";
  }
};


// Compare two Suffix objects that refer to suffixes in "words",
// starting at the start'th character (of each *suffix*) and comparing
// no more than length characters.
template <class T>
struct SuffixCompare {
  SuffixCompare( const vec< vec<T> >& words, int start, int length )
    : m_words( words ), m_start( start ), m_length( length ) {}
  
  bool operator() ( const Suffix<T>& lhs, const Suffix<T>& rhs ) {
    for ( int i = 0; i < m_length; ++i ) {
      if ( rhs.GetValue( m_words, m_start+i ) < 0 )
        return false;
      else if ( lhs.GetValue( m_words, m_start+i ) < 0 )
        return true;
      else if ( lhs.GetValue( m_words, m_start+i ) < rhs.GetValue( m_words, m_start+i ) )
        return true;
      else if ( lhs.GetValue( m_words, m_start+i ) > rhs.GetValue( m_words, m_start+i ) )
        return false;
    }
    return false;
  }
  
 private:
  const vec< vec<T> >& m_words;
  const int m_start;
  const int m_length;
};


template <class T> 
class vec_overlap<T>::imp {
 public:
  imp( const vec< vec<T> >& x );

  void GetOverlaps( const vec<T>& y, 
                    vec< pair<int,int> >& overlaps,
                    Bool allow_subset_x = True,
                    Bool recurse = True );
  
 private:
  vec< vec<T> > m_words;
  vec< Suffix<T> > m_suffixArray;
  vec< int > m_lcp; // longest common prefix
  int m_maxWordLength;
};

template <class T>
void RadixSort( typename vec< Suffix<T> >::iterator begin, 
                typename vec< Suffix<T> >::iterator end,
                const vec< vec<T> >& words, const int offset, 
                const int maxValue, const int maxLength ) 
{
  if ( distance( begin, end ) < 30 ) {
    sort( begin, end,
          SuffixCompare<T>( words, offset, maxLength+1 ) );
    return;
  }

  const int numBuckets = maxValue+2;

  vec<int> bucketUpperBound( numBuckets, 0 );

  // Fill bucketUpperBound.
  for ( typename vec< Suffix<T> >::iterator i = begin; i != end; ++i )
    ++bucketUpperBound[ i->GetValue(words,offset) + 1 ];

  int minBucket = numBuckets;
  for ( int bucket = 0; bucket < numBuckets; ++bucket )
    if ( bucketUpperBound[bucket] > 0 ) {
      minBucket = bucket;
      break;
    }

  int maxBucket = -1;
  for ( int bucket = maxValue+1; bucket >= 0; --bucket )
    if ( bucketUpperBound[bucket] > 0 ) {
      maxBucket = bucket;
      break;
    }

  partial_sum( bucketUpperBound.begin(), bucketUpperBound.end(),
               bucketUpperBound.begin() );

  vec<int> bucketLowerBound( numBuckets );
  bucketLowerBound[0] = 0;
  copy( bucketUpperBound.begin(), bucketUpperBound.end() - 1,
        bucketLowerBound.begin() + 1 );

  for ( int bucket = minBucket; bucket <= maxBucket; ++bucket )
  {
    // Now we go through all the positions where the elements in
    // this bucket should go.  If that position contains an element
    // from another bucket, swap that element into the first
    // available position in the bucket it belongs to.  Check the
    // newly swapped in element and repeat, until the element
    // swapped in belongs in the current bucket.  At that point,
    // move onto the next position, continuing until all the
    // elements in the bucket have been correctly placed.
    
    // This strategy has the advantage that every move of an element
    // puts it into the correct place, so every element is moved at
    // most once.
    while ( bucketLowerBound[bucket] != bucketUpperBound[bucket] ) {
      const int realBucket = (begin+bucketLowerBound[bucket])->GetValue(words,offset) + 1;
      if ( realBucket == bucket )
        ++bucketLowerBound[bucket];
      else {
        swap( *(begin+bucketLowerBound[bucket]), *(begin+bucketLowerBound[realBucket]) );
        ++bucketLowerBound[realBucket];
      }
    }
  }
   
  // We recurse through the buckets, skipping bucket 0, which is the
  // bucket of sequences that have already ended.
  for ( int bucket = max(1,minBucket); bucket <= maxBucket; ++bucket )
  {
    // We have to recover the lowerBound from the bucketUpperBound of
    // the previous bucket, since we clobbered bucketLowerBound in the
    // swapping process above.
    const int lowerBound = bucketUpperBound[ bucket-1 ];
    const int upperBound = bucketUpperBound[ bucket ];
    const int numElements = upperBound - lowerBound;

    if ( numElements > 1 )
      RadixSort( begin + lowerBound, begin + upperBound, words, offset+1, maxValue, maxLength );
  }
}


template <class T> 
vec_overlap<T>::imp::imp( const vec< vec<T> >& words ) 
  : m_words( words )
{
  m_maxWordLength = 0;
  int numSuffixes = 0;
  for ( unsigned int i = 0; i < m_words.size(); ++i ) {
    numSuffixes += m_words[i].size();
    m_maxWordLength = max<int>( m_maxWordLength, m_words[i].size() );
  }
  m_suffixArray.reserve( numSuffixes );

  T maxValue = -1;
  for ( unsigned int i = 0; i < m_words.size(); ++i )
    for ( unsigned int j = 0; j < m_words[i].size(); ++j ) {
      m_suffixArray.push_back( Suffix<T>( i, j ) );
      maxValue = max<T>( maxValue, m_words[i][j] );
    }

  RadixSort<T>( m_suffixArray.begin(), m_suffixArray.end(), m_words, 
                0, maxValue, m_maxWordLength );

  m_lcp.resize( m_suffixArray.size(), 0 );
  for ( unsigned int i = 1; i < m_suffixArray.size(); ++i )
    m_lcp[i] = m_suffixArray[i].GetLongestCommonPrefix( m_words, m_suffixArray[i-1] );
}


template <class T>
void vec_overlap<T>::imp::GetOverlaps( const vec<T>& y, 
                                       vec< pair<int,int> >& overlaps,
                                       Bool allow_subset_x,
                                       Bool recurse ) 
{
  overlaps.clear();

  int maxWordLength = max<int>( m_maxWordLength, y.size() );
  
  m_words.push_back( y );
  
  Suffix<T> ySuffix( m_words.size()-1, 0 );

  typename vec< Suffix<T> >::iterator insertionPoint =
    upper_bound( m_suffixArray.begin(), m_suffixArray.end(),
                 ySuffix,
                 SuffixCompare<T>( m_words, 0, maxWordLength + 1 ) );
  
  int index = distance( m_suffixArray.begin(), insertionPoint );

  int backIndex = index - 1;
  if ( backIndex >= 0 ) {
    int lcp = ySuffix.GetLongestCommonPrefix( m_words, m_suffixArray[backIndex] );
    
    while ( lcp > 0 ) {
      if ( lcp == m_suffixArray[backIndex].GetLength( m_words ) &&
           ( allow_subset_x || 
             m_suffixArray[backIndex].pos > 0 ||
             m_suffixArray[backIndex].GetLength( m_words ) > y.isize() ) )
        overlaps.push_back( make_pair( m_suffixArray[backIndex].index,
                                       - m_suffixArray[backIndex].pos ) );
      
      if ( backIndex == 0 )
        break;
      else {
        lcp = min( m_lcp[backIndex], lcp );
        backIndex--;
      }
    }
  }

  unsigned int forwIndex = index;
  if ( forwIndex < m_suffixArray.size() ) {
    unsigned int lcp = ySuffix.GetLongestCommonPrefix( m_words, m_suffixArray[forwIndex] );
    
    while ( lcp >= y.size() ) {
      if ( allow_subset_x ||
           m_suffixArray[forwIndex].pos > 0 ||
           m_suffixArray[forwIndex].GetLength( m_words ) > y.isize() )
        overlaps.push_back( make_pair( m_suffixArray[forwIndex].index,
                                     - m_suffixArray[forwIndex].pos ) );

      ++forwIndex;
      if ( forwIndex == m_suffixArray.size() ) 
        break;
      else
        lcp = m_lcp[forwIndex];
    }
  }

  m_words.pop_back();

  if ( recurse ) {
    typename vec_overlap<T>::imp yOverlapper( vec< vec<T> >( 1, y ) );
  
    vec< pair<int, int> > inverseOverlaps;
    for ( unsigned int x = 0; x < m_words.size(); ++x ) {
      yOverlapper.GetOverlaps( m_words[x], inverseOverlaps, allow_subset_x, False );
      for ( unsigned int o = 0; o < inverseOverlaps.size(); ++o )
        if ( inverseOverlaps[o].second != 0 &&
             ( allow_subset_x || - inverseOverlaps[o].second + m_words[x].isize() > y.isize() ) )
          overlaps.push_back( make_pair( x, - inverseOverlaps[o].second ) );
    }
  }
}


template <class T>
vec_overlap<T>::vec_overlap( const vec< vec<T> >& words )
  : m_pImp( new imp( words ) )
{}


template <class T>
vec_overlap<T>::~vec_overlap()
{
  delete m_pImp;
}


template <class T>
void vec_overlap<T>::GetOverlaps( const vec<T>& y, 
                                   vec< pair<int,int> >& overlaps,
                                   Bool allow_subset_x ) const
{
  m_pImp->GetOverlaps( y, overlaps, allow_subset_x, True );
}


#define INSTANTIATE(X) \
template vec_overlap<X>::vec_overlap( const vec< vec<X> >& ); \
template vec_overlap<X>::~vec_overlap(); \
template void vec_overlap<X>::GetOverlaps( const vec<X>&, vec< pair<int,int> >&, Bool ) const;

INSTANTIATE(int)
