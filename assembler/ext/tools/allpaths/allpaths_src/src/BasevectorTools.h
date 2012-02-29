/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// File: BasevectorTools.h
//
// This file defines some handy and generic enough methods to be used with/
// applied to basevectors.

#pragma once
#ifndef __BASEVECTOR_TOOLS_H
#define __BASEVECTOR_TOOLS_H

#include "math/Functions.h"
#include "Qualvector.h"
#include "Basevector.h"

/** \brief Implements trivial base comparator: two bases "match" only iff they are
 * the same. @see MismatchCount .
 *
 * This is a "functor" class; it has no data and only one method - operator() -
 * defined.
 */
struct ExactBaseMatch {
  bool operator() (base_t q, base_t t ) const {return q==t; }
};


/** \brief Implements generalized base comparator: works on \em letters (ACTG...),
 * \em not on numeric representation (hence the tag _Char)!!
 * This comparator is now symmetric:  either q or t or both can be generalized
 * bases in upper or lower case.
 * All unusual characters get mapped to "not a base", and don't match anything.
 * @see MismatchCount .
 *
 * This is a "functor" class; it has no data and only one method - operator() -
 * defined.
 */
struct GeneralizedBaseMatch_Char {
  bool operator() (char q, char t ) const {
      return GeneralizedBase::isGeneralizedBase(q) &&
             GeneralizedBase::isGeneralizedBase(t) &&
             GeneralizedBase::fromChar(q).matches(GeneralizedBase::fromChar(t));
  }
};




/** \brief Counts (the effective number of) mismatches in
 * the given alignment (as specified by start positions and length, no indels)
 * between the query sequence \c query and the target reference \c target.
 *
 * An alignment between the query sequence \c query and the target
 * reference \c target is defined by its \c length and start positions on
 * the query (\c q_start) and the target (\c t_start) sequences. For an
 * alignment so defined, this method calculates the "effective number of
 * mismatches":
 *    - If the qualities of the query bases are not specified, then
 *      the score is simply the number of generalized mismatches, where the
 *      "match" vs. "mismatch" decision for a pair of bases at the same
 *      offsets from the respective alignment starts on the query and target
 *      sequences is made by the client-specified comparator functor.
 *    - If the quality scores (non-null \c quals) are passed, then each
 *      generalized mismatch (as defined by the \c BaseComparator functor)
 *      is multiplied by max(quality at that base, \c minQual), such weighted
 *      mismatches are summed up (which gives, approximately, negative
 *      log-likelihood of the alignment), and the total error score is then
 *      divided by \c expectedQual and rounded. The output thus shows how many
 *      mismatches with the fixed quality \c expectedQual would result in the
 *      same probability of error.
 *
 * The template parameter (comparator functor) must have binary
 * bool operator()(base_t q, base_t t) defined; it's contract is to return
 * \c true when the base \c q on the query sequence is considered to be
 * matching the base \c t on the target sequence. Functors less trivial than
 * ( q==t ) can be defined for specific analysis needs; note also that the
 * functor is \e not assumed to be or treated as a symmetric one (the order of
 * the arguments matters): for instance, the functor is allowed to define that
 * T on a query sequence matches both T and C on the target, but C on the
 * query matches only C on the target. DO NOT forget to inline the operator()
 * in custom functors or performance hit can be expected!
 *
 *  The return value of this function is:
 *  if no quality scores are specified, then
 *  - the total number of mismatches if it does not exceed max. count cutoff
 *  - \c max_mismatches + 1 is returned as soon as the number of actually
 *    discovered mismatches exceeds the specified cutoff
 *    \c max_mismatches (this can save time)
 *  if quality scores are provided:
 *  - the alignment score ( sum of mismatch*
 *    max(\c minQual, \c quals[mismatch_position]) over all mismatches),
 *    divided by \c expectedQual; if the score does not exceed
 *    \c max_mismatches * \c expectedQual. [ Note that the score (before
 *    normalization) is approx. -10*(log-likelihood of the alignment)].
 *  - \c max_mismatches+1 is immediately returned as soon as the running
 *    cumulative score exceeds \c max_mismatches * \c expectedQual cutoff
 *    (can save time as the rest of the sequences is not scanned in that case)
 *
 *  NOTE: this method is \e unchecked. If the alignment as defined by start
 *  positions and the length run out of the sequence or qaulity vector boundaries,
 *  the behavior is undefined.
 *
 *  @param query query dna sequence
 *  @param target target dna sequence
 *  @param q_start start of the alignment on the \c query sequence
 *  @param t_start start of the alignment on the \c target sequence
 *  @param length length of the alignment
 *  @param max_mismatches (default 1,000,000) the moment the number of
 *         discovered mismatches exceeds this limit,
 *         the method immediately returns with \c max_mismatches+1
 *  @param[in] quals points to vector of base quality scores for the query
 *            sequence; can be \c null (default)
 *  @param expectQual used to normalize alignment score only when quality
 *         scores are passed (<tt> qual!=null</tt>); expected  quality (each
 *         mismatch is added scaled by a factor (this base quality) / \c expectedQual).
 *  @param minQual used only when quality scores are used; if quality of a base
 *          is worse than \c minQual, then \c minQual is used instead.
*/

template <class BaseComparator>
inline int MismatchCount(const basevector &query,
		      const basevector &target,
		      unsigned int q_start,
		      unsigned int t_start,
		      unsigned int length,
		      int max_mismatches = 1000 * 1000,
		      const qualvector * quals = 0,
		      double expectedQual = 30,
		      unsigned char minQual = 10
)
 {

   return MismatchCount(BaseComparator(),
			query, target, q_start, t_start, length,
			max_mismatches, quals, expectedQual, minQual);
 }


/// This is exaclty the same mismatch counter as the other
/// version (@see MismatchCount(basevector &, basevector &,...)), except
/// that it does not have to be invoked with explicit template argument
/// but rather takes an instance of base comparator as its first argument
/// (which is more STL-style and in line with STL algorithm implementations).
/// These two flavors will probably coexist until a decisive design choice
/// will be made (if ever). See the other implementation of MismatchCount
/// for the detailed description of the arguments.

template <class BaseComparator>
int MismatchCount(    const BaseComparator & equal,
		      const basevector &query,
		      const basevector &target,
		      unsigned int q_start,
		      unsigned int t_start,
		      unsigned int length,
		      int max_mismatches = 1000 * 1000,
		      const qualvector * quals = 0,
		      double expectedQual = 30,
		      unsigned char minQual = 10
)
 {
   basevector::iterator query_iter = query.Begin(q_start);
   basevector::iterator target_iter = target.Begin(t_start);
   const basevector::iterator query_end = query.End();
  if ( ! quals ) {
    return MismatchCount(equal,query_iter,query_end,target_iter,max_mismatches);
  } else {
    return MismatchScore(equal,query_iter,query_end,target_iter,
			 quals->begin()+q_start,max_mismatches*expectedQual,expectedQual,minQual);
  }

}


template < class BaseComparator, class InputIterator1, class InputIterator2 >
int MismatchCount( BaseComparator equal, InputIterator1 query_iter,
		      const InputIterator1 & query_end, InputIterator2 target_iter,
		      int max_mismatches = 1000 * 1000)
{
    int mismatches = 0; // mismatches discovered so far
    for ( ; query_iter != query_end; ++query_iter, ++target_iter )
    {
        // if comparator thinks that the bases "match"
        // go fetch next bases
        if (equal(*query_iter, *target_iter))
            continue;
        ++mismatches;    
        if (mismatches > max_mismatches) 
            break;
    }
   return mismatches;
}


template <class BaseComparator, class QueryIterator>
int MismatchScore(BaseComparator equal,
          QueryIterator query_iter,
		  QueryIterator query_end,
		  basevector::const_iterator target_iter,
		  qualvector::const_iterator qual_iter,
		  double maxBad,
		  double expectedQual = 30,
		  unsigned char minQual = 10)
{
  float bad = 0;
  for ( ; query_iter != query_end ;  ++qual_iter, ++query_iter, ++target_iter )  {
      if (  equal( *query_iter, *target_iter)  ) continue;
      bad += Max(minQual,(*qual_iter));
      if ( bad > maxBad ) {
	bad+=expectedQual;
	break;
      }
  }
  return int( round(bad/ expectedQual));
}



/// Scans sinchronously along the two iterators up until \c length steps are made
/// and returns a string of length \c length. In the returned string, the positions
/// where the same value occured in both sequences (represented by their iterators
/// passed to this method) are set to ' ' (space), and all the positions where
/// mismatch occured are set to \c mismatch_symbol (default is '*')
/// [the string thus represents a \em gapless
/// alignment]. NOTE: the iterators must implement operator++() and operator*()
/// (basevector iterators conform to this contract).

template <class queryiterator, class targetiterator>
String MismatchString( queryiterator qi, targetiterator ti, int length, char mismatch_symbol='*' ) {
  String S(length);
  for ( int i = 0 ; i < length ; i++, ++ti, ++qi) {
    if ( *ti == *qi ) S[i] = ' ';
    else S[i] = mismatch_symbol;
  }
  return S;
}

/// Same as MismatchString(queryiterator, targetiterator, int), but prints into
/// the passed String instance (and returns the reference to the same String
/// instance it is supplied with). The string \s will be resized to \c length.
/// Additional parameter \c overwrite controls whether the string will be completely
/// overwritten at every position (\c overwrite = \c true, default) or only mismatch
/// positions will be marked and matching positions will be left untouched (note that
/// positions in \c s where \c *qi and \c *ti do not match will be \em always overwritten)
template <class queryiterator, class targetiterator>
String & MismatchString( queryiterator qi, targetiterator ti, int length,
			   String &s, bool overwrite=true, char mismatch_symbol='*') {
  s.resize(length);
  for ( int i = 0 ; i < length ; i++, ++ti, ++qi ) {
    if ( *ti == *qi ) {
      if (overwrite) s[i] = ' ';
    }
    else s[i] = mismatch_symbol;
  }
  return s;
}

/// Same as MismatchString, but prints the mismatches directly into the
/// specified output stream. Prints exactly \c length characters:
/// ' ' (space) if the values at the corresponding position in the sequences
/// (represented by the passed iterators) are the same, or '*' if the values
/// at the corresponding position in the sequences differ (mismatch). NOTE:
/// this method \em does \em not print newline after the mismatch string!

template <class queryiterator, class targetiterator>
void PrintMismatchString( ostream & out, queryiterator qi, targetiterator ti, int length ) {
  for ( int i = 0 ; i < length ; i++, ++ti, ++qi ) {
    if ( *ti == *qi ) out << ' ';
    else out << '*';
  }
}



/// \brief Accumulates statistics of G and C counts observed across windows of
/// specified length \c W  placed uniformly, one at each
/// position on the specified \c reference sequence.
///
/// The GC content of a signle window is the count of G and C bases in that
/// window (thus it is a number in the range [0, W]). The statistics computed
/// is the numbers of windows that exhibited each individual GC count in that
/// range.
///
/// This method \em does \em not clear the output counter \c GC_counter, but
/// \em adds newly observed counts to the previously recorded ones (i.e.
/// multiple calls on, e.g., different reference sequences with the same
/// output counter will result in accumulation of counts).
///
/// Note: simple flat array or vector do conform to the contracts for the
/// counter (see below): \c GC_counter can be an
/// array/vector of W+1 elements. Design/performance/memory
/// considerations, however, may require counters backed by a
/// sparse vector/map etc.
///
/// @param W window size
/// @param reference[in] reference sequence
/// @param GC_counter[out] template argument, an indexable container; upon
/// return it will hold at position \c n the number of windows with
/// GC count equal to n found on the reference. The contract is that
/// 1) CounterOut::operator[unsigned int] is defined and returns (mutable)
/// reference that supports operator ++() (in a typical situation it will
/// simply return 'int &'); 2) CounterOut::operator[] call is valid with any
/// argument in the range [0...W] (i.e. the storage is pre-allocated).
/// @return window GC count statistics on the reference;
/// this method will \em increment all the elements of \c GC_counter at
/// positions [0,W] by the numbers of windows having the corresponding
/// GC count. If the output container have elements
/// at positions W+1 and above legitimately accessible, they will remain
/// accessible and the values at those positions will stay unchanged.
template <class CounterOut>
  void ComputeWindowGCCounts (CounterOut & GC_counter, ///< GC stats (output)
			      int W,  ///< window size
			      const basevector & reference ///< windows are placed on this sequence
			     )
{
      unsigned int contig_size = reference.size();

     // contig too short to accomodate a single window, nothing to do!
     if ( contig_size < (unsigned int)W ) return;

     // Simple (but inefficient) strategy equivalent to what this method
     // actually does is:
     // for each position p on the reference r {
     //    consider window w=r[p,p+W)
     //    local_gc_count = 0;
     //    for each base 0...W-1 in the window w {
     //          if ( base==C or base==G ) local_gc_count++;
     //    }
     //    // we just discovered yet another window with local_gc_count gc's:
     //    GC_counter[local_gc_count]++;
     // }
     //
     // or a similar nested loop for anchored windows, except that
     // the outer loop would be 'for each anchor position pi'
     // and GC[local_gc_count] should be incremented by n(pi)
     // (anchor multiplicity) rather than 1.
     //
     // Efficiency considerations (applies to counting all windows on the ref;
     // does not necessarily apply to counting anchored windows since there
     // maybe only a few of those):
     // 1) it is highly inefficient to re-scan the new window w every time:
     //    consider right-shifting a window along the reference by x
     //    positions (so far, we actually use only x=1 in this method!!):
     //    window [p,p+W) ---> new window [p+x,p+x+W). The gc count for the
     //    new window is equal to the gc count in the old window
     //    minus the number of GCs among bases [p,p+x) (leftmost in the old
     //    window) plus the number of GCs among bases [p+W,p+W+x) (rightmost
     //    in the new window) on the reference. This is an O(x) (compared to
     //    O(W) required to re-read all the bases in the new window from
     //    scratch  -  if W=50 and x=1, it's a lot).
     // 2) switching between accessing old window leftmost bases and new
     //    window rightmost bases directly on the reference requires direct,
     //    non-consequtive indexing of the basevector. The alternative
     //    is to store the required bases locally, then we can make only
     //    one linear pass through the reference, and for that we can use
     //    BasevectorIterator, which is much faster than
     //    basevector::operator[int].
     // 3) We can not just store the
     //    leftmost base (i=0) of the current window: after the window shifts,
     //    the i=1 base of the old window will be the leftmost and we will need
     //    it too, and so on. At any point, every base at position k in the
     //    current window will be the leftmost one after k shifts, so we have
     //    to remember them all for future use once we've read them.
     //    The leftmost base that gets pushed out of the window
     //    will never be needed again, so we have to
     //    remember exactly W bases at each point.
     // 4) It is unfeasible to keep a
     //    simple array/vector of the current window sequence: updating it for
     //    the new window would require left-shifting the sequence and adding
     //    new base(s) from the reference to the right. Shifting the
     //    whole sequence is O(W) and we loose all the benefits of reusing
     //    old window sequence/gc-count. To resolve this problem we use the
     //    *circular* array - an array of W bases, in which the logical start
     //    position can be anywhere inside the array. The sequnce represented
     //    by a circular array A with logical start position s is defined as
     //    A[s,W) followed by A[0,s) (wrapping over). In such representation,
     //    physical shift of the array elemenst is not required: when the window
     //    shifts (by 1), we simply increment the logical start position and
     //    update the base at the old logical start position with the new
     //    base form the reference. It is easy to see that {A[s+1,W),A[0,s+1)}
     //    correctly represents the sequence in the new window, provided A[s]
     //    is updated with the new window's rightmost base. Updating circular
     //    array is thus O[1] (update A[s] and increment s), and we can still
     //    enjoy the benefits of fast sequential access to the reference
     //    basevector.




//       ForceAssertLe(W,100);
       base_t * bases = new base_t[W]; // circular array
       int bases_logical_start = 0; // start of the sequence in the
                                    // circularly wrapped array 'bases'

       // start iterating from the first position on the contig:
       basevector::iterator b_iter = reference.Begin();

       int gc = 0 ; // will keep local gc-count within the current window

       // position of the current window on the contig:
       unsigned int window_position = 0;

       // initialization:
       // pre-compute gc count in the first window on the contig
       // and copy bases from the first window into the 'bases' array
       // (we checked above that contig size >= W, so the loop is safe):
       for ( int pos = 0 ; pos < W ; pos++, ++b_iter ) {
	 base_t base = bases[pos] = *b_iter;
	 if ( base == BASE_C || base == BASE_G ) gc++;
       }

       // don't forget to count in the gc value observed in the first window:
       ++GC_counter[gc];

       // now walk through all the (remaining) windows on the contig:
       for ( window_position++ ; window_position <= contig_size - W ; window_position++, ++b_iter ) {

	 // window has just shifted, get the next base -
	 // the rightmost one in the new window (we already have the others):
         base_t base = *b_iter;

	 // if a G or C base just moved into the window, increase the gc count
         if ( base == BASE_C || base == BASE_G ) gc++;

	 // the first base in the old window is pushed out;
	 // if it is G or C, then the gc count in the new window will decrease:
	 base_t old_base = bases[bases_logical_start];
	 if ( old_base == BASE_C || old_base == BASE_G ) gc--;

	 // update circular array so it is synchronised with the new window:
	 bases[bases_logical_start++] = base;

	 // make sure logical start position does not run away:
	 // wrap around (performing 'if'
	 // is probably faster than doing %=W, but who knows..)
	 if ( bases_logical_start >= W ) bases_logical_start -= W;

	 ++GC_counter[gc]; // count in the occurence of the count on the ref.

       } // end of for ( window_position ) loop
       // all windows on the contig are counted now

       delete [] bases; // free memory!
}



/// \brief Accumulates statistics of G and C counts observed across windows of
/// specified length \c W
/// placed at positions specified in \c anchor_pos_counts, with the specified
/// multiplicities (number of windows per anchor positions).
///
/// The GC content of a signle window is the count of G and C bases in that
/// window (thus it is a number in the range [0, W]). The statistics computed
/// is the numbers of windows that exhibited each individual GC count in that
/// range.
///
/// The array-like (or rather array-like addressable) structure
/// \c anchor_pos_counts should specify multiplicity count n(p_i) at each
/// position p_i where window(s) are anchored, and should return 0 at any
/// other position on the reference. Exactly n(p_i) copies
/// of the window at position p_i will be considered giving rise to n(p_i)
/// observations of the same GC count associated with that window.
///
/// This method \em does \em not clear the output counter \c GC_counter, but
/// \em adds newly observed counts to the previously recorded ones (i.e.
/// multiple calls on, e.g., different reference sequences with the same
/// output counter will result in accumulation of counts).
///
/// Note: simple flat arrays or vectors do conform to the contracts for the
/// in and out counters (see below), e.g. \c anchor_pos_counts can be an
/// array/vector of size(reference) elements with elements at p_i set
/// to n(p_i) at all other elements set to 0. Performance vs memory
/// considerations, however, may require counters backed by a
/// sparse vector/map etc.
///
/// @param W window size
/// @param reference[in] reference sequence
/// @param anchor_pos_counts[in] indexable container that provides counts n(pi)
/// at each anchor position pi. The actual argument (type) must conform to the
/// following contract: 1) it must have ::operator[unsigned int] defined
/// and returning an int (or a type automatically convertible to int),
/// 2) it must be legal to invoke this operator with any argument in the range
/// [ 0, size(reference) ), and 3) ::operator[pi] is expected to return anchor
/// multiplicity n(pi) for all positions pi on the reference contig, where the
/// anchor is present, and it must return 0 for all other positions.
/// @param GC_counter[out] template argument, an indexable container; upon
/// return it will hold at position \c n the number of windows with
/// GC count equal to n found on the reference (with multiplicities counted in)
/// The contract is that
/// 1) CounterOut::operator[unsigned int] is defined and returns (mutable)
/// reference that supports operator +=(int) (in a typical situation it will
/// simply return 'int &'); 2) CounterOut::operator[] call is valid with any
/// argument in the range [0...W] (i.e. the storage is pre-allocated).
/// @return window GC count statistics on the reference computed across
/// anchored windows with multiplicities;
/// this method will \em increment all the elements of \c GC_counter at
/// positions [0,W] by the numbers of windows having the corresponding
/// GC count. If the output container have elements
/// at positions W+1 and above legitimately accessible, they will remain
/// accessible and the values at those positions will stay unchanged.
template <class CounterIn, class CounterOut>
  void ComputeWindowGCCounts (CounterOut & GC_counter, ///< GC stats (output)
			      int W,  ///< window size
			      const basevector & reference, ///< windows are placed on this sequence
			      const CounterIn & anchor_pos_counts) ///< window positions/multiplicities
{
      unsigned int contig_size = reference.size();

     // contig too short to accomodate a single window, nothing to do!
     if ( contig_size < (unsigned int)W ) return;

     //  see the other implementation for the description of the algorithm!

     // count the total number of anchors (without multiplicities!) requested:
     unsigned int total_pos_cnt = 0;

     for ( unsigned int i = 0 ; i <= contig_size - W ; i++ ) {
       if ( anchor_pos_counts[i] > 0 ) total_pos_cnt++;
     }

     // check how dense, on average, the anchors on the reference are.
     // if anchors occur on average more often that every W/2 bases,
     // we will simply scan the whole reference; if the anchors are
     // far apart ( > W/2 on average ), we will independently read
     // all bases for each window from the reference. [the W/2 cutoff
     // is random and not tested for being optimal...]
     if ( total_pos_cnt < 2*(contig_size/W) ) {
       // anchors are sparser than one per W/2 bases:

       for ( unsigned int window_pos = 0 ; window_pos <= contig_size - W ; window_pos++ ) {
	 unsigned int n = (unsigned int) anchor_pos_counts[window_pos];
	 if ( n == 0 ) continue; // no windows anchored at this position

	 unsigned int gc_cnt = 0;

	 // initialize iterator at current position on the reference
	 basevector::iterator b_iter = reference.Begin(window_pos);
	 // read all bases for the current window and count gc's:
	 for ( int pos = 0 ; pos < W ; pos++, ++b_iter ) {
	   base_t base = *b_iter;
	   if ( base == BASE_C || base == BASE_G ) gc_cnt++;
	 }
	 GC_counter[gc_cnt]+=n;
       }
     } else {
       // anchor are dense, let's read through the whole reference in one pass

       base_t * bases = new base_t[W]; // circular array
       int bases_logical_start = 0; // start of the sequence in the
                                    // circularly wrapped array 'bases'



       // start iterating from the first position on the contig:
       basevector::iterator b_iter = reference.Begin();

       int gc_cnt = 0 ; // will keep local gc-count within the current window

       // position of the current window on the contig:
       unsigned int window_position = 0;

       // initialization:
       // pre-compute gc count in the first window on the contig
       // and copy bases from the first window into the 'bases' array
       // (we checked above that contig size >= W, so the loop is safe):
       for ( int pos = 0 ; pos < W ; pos++, ++b_iter ) {
	 base_t base = bases[pos] = *b_iter; // get and store current base
	 if ( base == BASE_C || base == BASE_G ) gc_cnt++; // count Cs and Gs
       }
       // now we got gc count in the first window and window bases are stored;
       // done with initialization!

       // don't forget to count in the gc value observed if there is a window
       // anchored at window_position=0:
       GC_counter[gc_cnt] += anchor_pos_counts[window_position];

       // now walk through all the (remaining) windows on the contig:
       for ( window_position++ ; window_position <= contig_size - W ; window_position++, ++b_iter ) {

	 // window has just shifted, so we need to get the next base -
	 // the rightmost base in the new window (we already have the others):
         base_t base = *b_iter;

	 // if a G or C base just moved into the window, increase the gc count
         if ( base == BASE_C || base == BASE_G ) gc_cnt++;

	 // if G or C was on the left of the old window (now pushed out), decrease count:
	 base_t old_base = bases[bases_logical_start];
	 if ( old_base == BASE_C || old_base == BASE_G ) gc_cnt--;

	 // update circular array so that it is synchronised with the new window content:
	 bases[bases_logical_start++] = base;

	 // wrap around logical start position
	 if ( bases_logical_start >= W ) bases_logical_start -= W;

	 // if the anchor at the current window position exists,
	 // count in the current window's gc count with the anchor's
	 // multiplicity (if there is no anchor at window_position,
	 // anchor_pos_counts MUST return 0!):
	 GC_counter[gc_cnt] += anchor_pos_counts[window_position];
       } // end of for ( window_position ) loop
       // all windows on the contig are counted now

       delete [] bases; // free memory!
     }
}





/// \brief Accumulates statistics of G and C counts observed across windows of
/// specified length \c W  placed 1) uniformly, one at each
/// position on the specified \c reference sequence and 2) at specified
/// positions with specified multiplicity.
///
/// This is a convenience method: saves onle line of typing and works slightly
/// faster (in most cases, insignificantly) then two separate methods provided
/// for computing GC count statistics for the cases 1) and 2) listed above. See
/// documentation for the other two overloaded implementations for more details.
/// This method simply comp utes both stats in one pass.
/// @param W window size
/// @param reference[in] reference sequence
/// @param GC_counter_ref[out] template argument, an indexable container; upon
/// return it will hold at position \c n the number of windows with
/// GC count equal to n found on the reference. The contract is that
/// 1) CounterOut::operator[unsigned int] is defined and returns (mutable)
/// reference that supports operator ++() and operator +=(int)
/// (in a typical situation it will
/// simply return 'int &'); 2) CounterOut::operator[] call is valid with any
/// argument in the range [0...W] (i.e. the storage is pre-allocated).
/// @param GC_counter_anchor[out] same as \c GC_counter ref but will hold
/// numbers of anchored windows with given gc counts (with multiplicities
/// counted in).
/// @param anchor_pos_counts[in] indexable container that provides multiplicities
///  n(pi) at each anchor position pi. The actual argument (type) must conform
/// to the following contract: 1) it must have ::operator[unsigned int] defined
/// and returning an int (or a type automatically convertible to int),
/// 2) it must be legal to invoke this operator with any argument in the range
/// [ 0, size(reference) ), and 3) ::operator[pi] is expected to return anchor
/// multiplicity n(pi) for all positions pi on the reference contig, where the
/// anchor is present, and it must return 0 for all other positions.
/// @return window GC count statistics on the reference and for the anchored windows;
/// this method will \em increment all the elements of \c GC_counter's at
/// positions [0,W] by the numbers of windows (with multiplicities) having the
/// corresponding GC count. If the output containers have elements
/// at positions W+1 and above legitimately accessible, they will remain
/// accessible and the values at those positions will stay unchanged.
template <class CounterIn, class CounterOut>
  void ComputeWindowGCCounts (CounterOut & GC_counter_ref, ///< GC stats on ref(output)
			      CounterOut & GC_counter_anchor,
			      int W,  ///< window size
			      const basevector & reference, ///< windows are placed on this sequence
			      const CounterIn & anchor_pos_counts) ///< window positions/multiplicities
{
      unsigned int contig_size = reference.size();

     // contig too short to accomodate a single window, nothing to do!
     if ( contig_size < (unsigned int)W ) return;

     // see the comments in the overloaded implementation of the
     // single counter on the reference for the details of the algorithm

       base_t * bases = new base_t[W]; // circular array
       int bases_logical_start = 0;

       // start iterating from the first position on the contig:
       basevector::const_iterator b_iter = reference.Begin();

       int gc = 0 ; // will keep local gc-count within the current window

       unsigned int window_position = 0;

       for ( int pos = 0 ; pos < W ; pos++, ++b_iter ) {
	 base_t base = bases[pos] = *b_iter;
	 if ( base == BASE_C || base == BASE_G ) gc++;
       }

       // don't forget to count in the gc value observed in the first window:
       ++GC_counter_ref[gc];
       GC_counter_anchor[gc]+=anchor_pos_counts[window_position];

       // now walk through all the (remaining) windows on the contig:
       for ( window_position++ ; window_position <= contig_size - W ; window_position++, ++b_iter ) {

         base_t base = *b_iter;
         if ( base == BASE_C || base == BASE_G ) gc++;

	 base_t old_base = bases[bases_logical_start];
	 if ( old_base == BASE_C || old_base == BASE_G ) gc--;

	 bases[bases_logical_start++] = base;
	 if ( bases_logical_start >= W ) bases_logical_start -= W;

	 ++GC_counter_ref[gc]; // count in the occurence of the count on the ref.
	 GC_counter_anchor[gc]+=anchor_pos_counts[window_position];

       } // end of for ( window_position ) loop
       // all windows on the contig are counted now

       delete [] bases; // free memory!
}


/// See ComputeWindowGCCounts(CounterOut &, CounterOut &, int,
/// const basevector &, const CounterIn &) for details. This
/// method is a simple wrapper that calls the above method
/// on each forward and rc contig in the \c reference and
/// accumulates the total counts across all the contigs into the
/// \c GC_counter_ref and \c GC_counter_anchor counters.
///
/// Note: the CounterInCollection is any indexable container,
/// such that \c anchor_pos_counts[2*i+d] is legal and
/// returns a start position counter on the forward (d=0) or
/// reverse complement (d=1) strand of the reference contig
/// \c reference[i] (see the cited above method for the contract
/// each such counter should conform to).

template < class CounterInCollection, class CounterOut >
  void ComputeWindowGCCounts(CounterOut & GC_counter_ref,
			     CounterOut & GC_counter_anchor,
			     int W,
			     const vecbasevector & reference,
			     const CounterInCollection & anchor_pos_counts) {

  for ( size_t i = 0; i < 2*reference.size( ); i++ ) {
      // if i is even, we will look at the forward strand
      // of the reference contig reference[i], otherwise
      // we will count on the reverse strand:
      basevector RC;
      if ( i%2 == 1 ) RC.ReverseComplement(reference[i/2]);
      const basevector & contig = ( i%2 == 0 ? reference[i/2] : RC );
      ComputeWindowGCCounts(GC_counter_ref, GC_counter_anchor, W, contig, anchor_pos_counts[i]);
  }
}


/// This is a shortcut wrapper for
/// ComputeWindowGCCounts(CounterOut &, int, basevector, CounterIn &)
/// (see documentation for that method). The latter method will be called
/// for every forward and reverse complement contig of reference with
/// corresponding element of the <CounterInCollection> as window anchor counter, and
/// the results across all contigs will be accummulated into <GC_counter>.
///
///
/// Note: the CounterInCollection is any indexable container,
/// such that \c anchor_pos_counts[2*i+d] is legal for i=0...reference.size()-1
/// and returns a start position counter on the forward (d=0) or
/// reverse complement (d=1) strand of the reference contig
/// \c reference[i] (see the cited above method for the contract
/// each such counter should conform to).

template <class CounterInCollection, class CounterOut>
  void ComputeWindowGCCounts (CounterOut & GC_counter, ///< GC stats (output)
			      int W,  ///< window size
			      const vecbasevector & reference, ///< windows are placed on this sequence
			      const CounterInCollection & anchor_pos_counts) ///< window positions/multiplicities
{
  for ( int i = 0 ; i < 2*reference.size(); i++ ) {
      // if i is even, we will look at the forward strand
      // of the reference contig reference[i], otherwise
      // we will count on the reverse strand:
    basevector RC;
    if ( i%2 == 1 ) RC.ReverseComplement(reference[i/2]);
    const basevector &contig = ( i%2 == 0 ? reference[i/2] : RC );
    ComputeWindowGCCounts(GC_counter, W, contig, anchor_pos_counts[i]);
  }
}

#endif
//#ifndef __BASEVECTOR_TOOLS_H


