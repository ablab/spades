///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <set>
#include <map>

#include "pairwise_aligners/AlignFromMutmersAndSW.h"
#include "pairwise_aligners/MakeAlignsMethod.h"
#include "pairwise_aligners/Mutmer.h"
#include "PrintAlignment.h"
#include "system/Types.h"
#include "math/Functions.h"
#include "pairwise_aligners/AlignFromMutmers.h"
#include "pairwise_aligners/SmithWatBandedA.h"

// If there's a perfect match of length at least perf but we produced no proper 
// alignment, use banded Smith-Waterman to create one.

void NoProperAligns( vec<align>& aligns, int& aligns_length, 
     const vec<mutmer>& mm, const basevector& rd1, const basevector& rd2,
     int perf, int bandwidth )
{    Bool proper = False;
     for ( int i = 0; i < aligns_length; i++ )
     {    if ( Proper( aligns[i], rd1.size( ), rd2.size( ) ) )
          {    proper = True;
               break;    }    }
     if (proper) return;
     int p1, p2, len, e, max_perfect_match = 0, max_mutmer = 0;
     int best_mutmer = -1;
     for ( unsigned int t = 0; t < mm.size( ); t++ )
     {    mm[t].Unpack( p1, p2, len, e );
          if ( len < perf + e ) continue;
          int last_diff = 0;
          for ( int i = 0; i < len; i++ )
          {    if ( rd1[ p1 + i ] != rd2[ p2 + i ] )
               {    if ( i - last_diff > max_perfect_match )
                    {    max_perfect_match = i - last_diff;
                         best_mutmer = t;    }
                    last_diff = i + 1;    }    }
          max_mutmer = Max( max_mutmer, len );
          if ( len - last_diff > max_perfect_match )
          {    max_perfect_match = len - last_diff;
               best_mutmer = t;    }    }
     if ( max_perfect_match < perf ) return;
     /*
     cout << "see no alignment between sequences of lengths "
          << rd1.size( ) << " and " << rd2.size( )
          << ", but there is a perfect match of length "
          << max_perfect_match << " and a mutmer of length " << max_mutmer
          << "\n";    
     cout << "aligning:\n";
     */
     align a;
     int a_errors;
     mm[best_mutmer].Unpack( p1, p2, len, e );
     SmithWatBandedA( rd1, rd2, p1 - p2, bandwidth, a, a_errors );
     if ( aligns_length < (int) aligns.size( ) ) aligns[aligns_length++] = a;
     // PrintVisualAlignment( True, cout, rd1, rd2, a );    
          }

// --- ORIG METHOD

Bool makealigns_orig_method::MutmersToAlign( const vec<mutmer>& mm, int k,
     const basevector& rd1, const basevector& rd2, vec<align>& aligns, 
     vec<int>& errors, int& aligns_length, int min_mutmer, ostream* log )
{    
  // Old method.

  // Unpack the mutmers.
  
  vector<char> dead(100);
  if ( mm.size( ) > dead.size( ) ) dead.resize( mm.size() );    

  for ( unsigned int t = 0; t < mm.size( ); ++t )
    dead[t] = 0;    

  // for ( unsigned int t = 0; t < mm.size( ); t++ ) // XXX
  //      cout << t << ": pos1 = " << p1[t] << ", pos2 = " << p2[t] // XXX
  //          << ", len = " << l[t] << ", e = " << e[t] << "\n"; // XXX
  
  // Eliminate mutmers which are clearly dominated by another mutmer.
  // We use the following criteria to define when mutmer M dominates mutmer m:
  //
  // (i)   M.pos1 <= m.pos1 and M.Pos1 >= m.pos1
  // (ii)  M.pos2 <= m.pos2 and M.Pos2 >= m.pos2
  // (iii) length(M) >= 10 * length(m)
  // (iv)  the number of errors of M between m.pos1 and m.Pos1 is strictly
  //       less than the number of errors of m.
  
  int min_l = mm[0].Length();
  for ( unsigned int i = 1; i < mm.size( ); ++i )
    min_l = Min( min_l, mm[i].Length() );
  for ( unsigned int i = 0; i < mm.size( ); ++i )
  {
    if ( dead[i] ) continue;
    if ( mm[i].Length() < 10 * min_l ) continue;
    for ( unsigned int j = 0; j < mm.size( ); ++j )
    {
      if ( dead[j] || i == j ) continue;
      
      // Determine if M = mutmer i dominates m = mutmer j.
      
      if ( mm[i].Length() < 10 * mm[j].Length() ) continue;
      if ( mm[i].Pos1() > mm[j].Pos1() || mm[i].End1() < mm[j].End1() ) continue;
      if ( mm[i].Pos2() > mm[j].Pos2() || mm[i].End2() < mm[j].End2() ) continue;
      int errs = 0;
      for ( int h = mm[j].Pos1(); h < mm[j].End1(); ++h )
	if ( rd1[h] != rd2[ h + mm[i].Pos2() - mm[i].Pos1() ] ) ++errs;
      if ( errs < mm[j].Errors() ) dead[j] = 1;  
    }
  }

//      if ( rd1.Length( ) == 220678 )
//      {    cout << "live mutmers:\n";
//           for ( unsigned int ii = 0; ii < mm.size( ); ii++ )
//           {    if ( dead[ii] ) continue;
//                int pos1, pos2, len, e;
//                mm[ii].Unpack(pos1, pos2, len, e);
//                if ( ii == 0 ) *log << "\n" << mm.size( ) << " max-mutmers:\n";
//                *log << "   [" << ii << "] pos1 = " << pos1 << ", pos2 = " << pos2
//                     << ", diff = " << pos2 - pos1 <<
//                     ", length = " << len << ", errors = " << e << "\n" << flush; }  }

  // Identify and eliminate mutmers which appear to consist entirely of
  // tandem repeats.
  
  // The reason for doing this is that tandem repeats can give rise to large
  // numbers of mutmers, many of which have slightly different offsets and
  // are otherwise nearly identical.  In turn, these large numbers of similar
  // mutmers give rise to staggering numbers of alignments, bringing the 
  // program to a halt, unless one caps the number of alignments.  But if one
  // does the latter, good alignments are lost.
  
  // It is clear that one can do a good job of identifying the tandem repeats
  // by looking only at the mutmer data: p1, p2, and l.  However, the approach
  // used here is a primitive and rudimentary implementation, which could
  // easily be improved.
  
  // Here's what we do.  We consider all ordered pairs of mutmers whose offset 
  // differs by at most 10.  For each such pair, suppose that the overlap
  // between the mutmers (as measured by the first read) accounts for >= 70%
  // of the length of one (or both) of the mutmers.  According we tag one
  // (or both) of the mutmers.  Any mutmer which ultimately gets tagged twice
  // is killed.
  
  vector<int> hits(100);
  if ( mm.size( ) > hits.size( ) ) hits.resize( mm.size( ) );

  for ( unsigned int i = 0; i < mm.size( ); ++i )
    hits[i] = 0;
  for ( unsigned int i = 0; i < mm.size( ); ++i )
  {    
    if ( dead[i] ) continue;
    int start1 = mm[i].Pos1(), stop1 = mm[i].End1();
    int offset1 = mm[i].Pos2() - mm[i].Pos1();
    for ( unsigned int j = i + 1; j < mm.size( ); ++j )
    {
      int start2 = mm[j].Pos1(), stop2 = mm[j].End1();
      float overlap = float( Min( stop1, stop2 ) - Max( start1, start2 ) );
      int offset2 = mm[j].Pos2() - mm[j].Pos1();
      if ( Abs( offset1 - offset2 ) > 10 ) continue;
      if ( overlap >= 0.7 * float(mm[i].Length()) ) ++hits[i];
      if ( overlap >= 0.7 * float(mm[j].Length()) ) ++hits[j]; 
    }
  }
  for ( unsigned int i = 0; i < mm.size( ); ++i )
    if ( hits[i] >= 2 ) dead[i] = 1;

  // Try the following new method.  First pick the top 20 mutmers by length.
  // Then use only those mutmers whose offset is within 20 of one of the top
  // 20 mutmers.
  
  if ( mm.size( ) > 20 )
  {
    vec< pair<int, int> > length_offset;
    for ( unsigned int i = 0; i < mm.size( ); ++i )
      if ( !dead[i] )
	length_offset.push_back( pair<int, int>( mm[i].Length(), mm[i].Pos2() - mm[i].Pos1() ) );
    sort( length_offset.rbegin( ), length_offset.rend( ) );
    int max_offset = Min( 20, (int) length_offset.size( ) );
    for ( unsigned int i = 0; i < mm.size( ); ++i )
    {
      if ( dead[i] ) continue;
      int offset = mm[i].Pos2() - mm[i].Pos1();
      int j;
      for ( j = 0; j < max_offset; ++j )
	if ( Abs( offset - length_offset[j].second ) < 20 )
	  break;    
      if ( j == max_offset ) dead[i] = 1;
    }
  }

  // Now create a graph with a vertex for each max-mutmer, and an edge
  // from one max-mutmer to another if one closely follows the other,
  // by which we mean that the gaps on both reads are between
  // -nstretch*k and stretch*k.
  
  // Exception: no edge is created if the two max-mutmers overlap by
  // more than 50%.
  
  // Exception (added 9/4/00): no edge is created if the two
  // max-mutmers overlap and have offsets different by more than 5.
  
  // Changed (11/2/00): no edge is created if the two max-mutmers
  // overlap and have offsets different by more than 5 plus 1% of the
  // shortest of the two max-mutmers.
  
  // Exception (added 9/5/00): no edge is created if the two
  // max-mutmers have offsets different by more than 25.
  
  map< int, vector<int> > edges;		// a sparse matrix
  for ( unsigned int i = 0; i < mm.size( ); ++i )
    for ( unsigned int j = 0; j < mm.size( ); ++j )
    {
      if ( i == j || dead[i] || dead[j] ) continue;
      int offset1 = mm[i].Offset();
      int offset2 = mm[j].Offset();
      int discrep = Abs( offset1 - offset2 );
      if ( discrep > 25 ) continue;
      if ( mm[i].Length() < min_mutmer || mm[j].Length() < min_mutmer ) continue;
      if ( mm[i].End1() > mm[j].Pos1() + nstretch_*k ) continue;
      if ( mm[i].End1() > mm[j].Pos1() + mm[j].Length()/2 ) continue;
      if ( mm[i].End2() > mm[j].Pos2() + nstretch_*k ) continue;
      if ( mm[i].End2() > mm[j].Pos2() + mm[j].Length()/2 ) continue;
      if ( mm[j].Pos1() - mm[i].Pos1() - mm[i].Length() > stretch_*k ) continue;
      if ( mm[j].Pos2() - mm[i].Pos2() - mm[i].Length() > stretch_*k ) continue;
      int max_discrep = 5 + Min( mm[i].Length(), mm[j].Length() ) / 100;
      if ( mm[i].End1() > mm[j].Pos1() && discrep > max_discrep ) continue;
      if ( mm[i].End2() > mm[j].Pos2() && discrep > max_discrep ) continue;
      
      edges[i].push_back(j);
    }

  // Identify edges i --> j which are required, in the sense that if you enter 
  // vertex i, you must then follow the edge to vertex j.  The rule is that 
  // i.length >= 20, j.length >= 20, Abs(i.Pos1 - j.pos1) <= 4, 
  // Abs(i.Pos2 - j.pos2) <= 4, i.error_rate <= 5%, j.error_rate <= 5%, and
  // that there is no j' such that Abs(i.Pos1 - j'.pos1) <= 20, 
  // Abs(i.Pos2 - j.pos2) <= 20.   (Actually 20 is cl.)
  
  for ( map< int, vector<int> >::iterator i = edges.begin(); i != edges.end(); ++i )
  {
    if ( i->second.size() == 1 || mm[i->first].Length() < cl_ || 20 * mm[i->first].Errors() > mm[i->first].Length() ) continue;
    for ( vector<int>::const_iterator j = i->second.begin(); j != i->second.end(); ++j )
    {
      if ( mm[*j].Length() >= cl_ && 20 * mm[*j].Errors() <= mm[*j].Length() )
      {
	if ( Abs( mm[i->first].End1() - mm[*j].Pos1() ) <= 4 &&
	     Abs( mm[i->first].End2() - mm[*j].Pos2() ) <= 4 )
        {
          vector<int>::const_iterator jp;
          for ( jp = i->second.begin(); jp != i->second.end(); ++jp )
	    if ( jp != j
		 && Abs( mm[i->first].End1() - mm[*jp].Pos1() ) <= cl_ 
		 && Abs( mm[i->first].End2() - mm[*jp].Pos2() ) <= cl_ ) break;
	  if ( jp == i->second.end() )
	  { 
	    int x = *j;
	    i->second.clear();
	    i->second.push_back(x);
            break;
	  }
        }
      }
    }
  }

  // Mark as pseudo-dead any vertices in the graph that do have an edge
  // coming into them (i.e., that are not sources)
  
  for ( map< int, vector<int> >::const_iterator i = edges.begin(); i != edges.end(); ++i )
    for ( vector<int>::const_iterator j = i->second.begin(); j != i->second.end(); ++j )
      dead[*j] = 2;

  // Follow all maximal paths in the graph, and create the
  // corresponding alignments.

  vector<int> path(255);
  unsigned int path_max = Min ( 255, (int)mm.size() );
  path[0] = -1;
  unsigned int level = 1;

  int alignment_constructor_calls = 0;
  
  aligns_length = 0;
  while(1)
  {    
    --level;	// Advance to next path.
    map< int, vector<int> >::const_iterator x;
    if ( level == 0 )
    {
      unsigned int i;
      for ( i = path[0] + 1; i < mm.size( ) && dead[i]; ++i ) { }
      if ( i == mm.size( ) ) break;
      path[level++] = i;
    }
    else if ( (x = edges.find( path[level-1] )) != edges.end() )
    {
      vector<int>::const_iterator a;
      for ( a = x->second.begin(); *a != path[level]; ++a ) { }
      if ( ++a == x->second.end() )
        continue;
      path[level++] = *a;
    }
    else
      continue;

    x = edges.find( path[level-1] );
    while ( x != edges.end() && level < path_max ) {
      x = edges.find( path[level++] = x->second[0] );
    }

    // Create an alignment for the current path.
    
    shortvector<mutmer> m(level);
    for ( unsigned int i = 0; i < level; ++i )
      m(i) = mm[ path[i] ];
    if ( m(0).Length() < min_mutmer ) continue;
    /*
      cout << "using"; // XXX
      for ( unsigned int i = 0; i < level; i++ ) // XXX
        cout << " " << path[i]; // XXX
      cout << "\n"; // XXX
      cout << "with offsets"; // XXX
      for ( unsigned int i = 1; i < level; i++ ) // XXX
      {
        int offset1 = p1[path[i]] - p2[path[i]]; // XXX
        int offset2 = p1[path[i-1]] - p2[path[i-1]]; // XXX
        cout << " " << offset1 - offset2;    // XXX
      } // XXX
      cout << "\n"; // XXX
    */

    if ( aligns_length == (int) aligns.size( ) )
    {
      if ( log != 0 )
      {   *log << "MutmersToAlign: over " << aligns.size( ) 
	       << " alignments were found while aligning two sequences "
	       << "of lengths " << rd1.size( ) << " and "
	       << rd2.size( ) << ".\n";
           *log << "(called with min_mutmer = " << min_mutmer << ")\n";
           *log << "Here are the stats on the max-mutmers:\n\n";
           int mutmer_count = 0;
           for ( unsigned int ii = 0; ii < mm.size( ); ++ii )
           {    if ( dead[ii] == 1 ) continue;
	        if ( mm[ii].Length() < min_mutmer ) continue;
	        ++mutmer_count;    }
           *log << "\n" << mutmer_count << " max-mutmers:\n";
           for ( unsigned int ii = 0; ii < mm.size( ); ++ii )
          {    if ( dead[ii] == 1 || mm[ii].Length() < min_mutmer ) continue;
	        *log << "   [" << ii << "] " << mm[ii] << "\n";     }    }
      if (try_if_no_proper_) NoProperAligns( aligns, aligns_length, mm, rd1, rd2,
          try_if_no_proper_perf_, try_if_no_proper_bandwidth_ );
      return False;    
    }

    ++alignment_constructor_calls;
    if ( alignment_constructor_calls > max_alignment_constructor_calls_ )
    {
       if (try_if_no_proper_) NoProperAligns( aligns, aligns_length, mm, rd1, rd2,
          try_if_no_proper_perf_, try_if_no_proper_bandwidth_ );
       return False;    }

    align a;
    int errors_found;
    a.CreateFromMutmers( k, m, rd1, rd2, max_errs_, max_badness_ * 1.5, 
			 local_max_errs_, end_stretch_, 
			 local_max_errors_done_, errors_found );

    if ( a.Nblocks( ) > 0 )
    {    
      ForceAssertLe( a.Pos1( ), (int) rd1.size( ) );
      ForceAssertLe( a.Pos2( ), (int) rd2.size( ) );
      errors[aligns_length] = errors_found;
      aligns[aligns_length++] = a;   
    }
  }
  if (try_if_no_proper_) NoProperAligns( aligns, aligns_length, mm, rd1, rd2,
          try_if_no_proper_perf_, try_if_no_proper_bandwidth_ );
  return True;    
}

// --- ALT METHOD

Bool makealigns_alt_method::MutmersToAlign( const vec<mutmer>& mm, int k,
					    const basevector& rd1, const basevector& rd2, 
					    vec<align>& aligns, vec<int>& errors, 
					    int& aligns_length, 
					    int min_mutmer, ostream* log )
{    
  // Unpack the mutmers.
  
  vector<char> dead(100);
  dead.resize(mm.size( ));    

  for ( unsigned int t = 0; t < mm.size( ); ++t )
    dead[t] = False;    

  // New method.  
  //
  // Find the ten longest max-mutmers.  Then do a banded Smith-Waterman 
  // (bandwidth given as an argument), using the offsets given by the ten 
  // longest max-mutmers.  Obviously this could be refined!!
  //
  // Refinement added: the mutmer must have length at least 10% of the
  // implied overlap between the two basevectors, or have length at least 100.

  aligns_length = 0;
  const int mutmers_to_use = 10;
  ForceAssertLe( mutmers_to_use, (int) aligns.size( ) );
  vec<int> offsets;
  if ( (int) mm.size( ) <= mutmers_to_use )
  {
    for ( unsigned int i = 0; i < mm.size( ); ++i )
      if ( mm[i].Length() >= min_mutmer ) 
      {
	int offset = mm[i].Offset();
	int mutmer_length = mm[i].Length();
	int length1 = rd1.size( ), length2 = rd2.size( );
	int implied_overlap = 
	  IntervalOverlap(0, length1, offset, offset + length2);
	if ( mutmer_length < 100 
	     && mutmer_length < implied_overlap/10 ) continue;
	offsets.push_back(offset); 
      }
  }
  else
  {
    vec< vec<int> > length_offset;
    for ( unsigned int i = 0; i < mm.size( ); ++i )
    {
      vec<int> lo(2);
      lo[0] = mm[i].Length();
      lo[1] = mm[i].Offset();
      if ( mm[i].Length() >= min_mutmer ) length_offset.push_back(lo);    
    }
    sort( length_offset.rbegin( ), length_offset.rend( ) );
    for ( int i = 0; 
	  i < Min( mutmers_to_use, (int) length_offset.size( ) ); ++i )
    {
      int offset = length_offset[i][1];
      int mutmer_length = length_offset[i][0];
      int length1 = rd1.size( ), length2 = rd2.size( );
      int implied_overlap 
	= IntervalOverlap( 0, length1, offset, offset + length2 );
      if ( mutmer_length < 100 && mutmer_length < implied_overlap/10 )
	continue;
      offsets.push_back( length_offset[i][1] ); 
    }
  }
  for ( unsigned int i = 0; i < offsets.size( ); ++i )
  {
    align a;
    int a_errors;
    if ( rd1.size( ) <= 10 * rd2.size( ) )
      SmithWatBandedA( rd1, rd2, offsets[i], bandwidth_, a, a_errors );
    else
    {
      SmithWatBandedA( rd2, rd1, -offsets[i], bandwidth_, a, a_errors );
      a.Flip( );    
    }
    a_errors = ActualErrors(rd1, rd2, a);
    if ( a_errors > max_errs_ ) continue;
    errors[aligns_length] = a_errors;
    aligns[aligns_length++] = a;    
  }
  return True;    
}


// --- PERFECT METHOD

Bool makealigns_perfect_method::MutmersToAlign( const vec<mutmer>& mm, int k,
                                                const basevector& rd1, const basevector& rd2, 
                                                vec<align>& aligns, vec<int>& errors, 
                                                int& aligns_length, 
                                                int min_mutmer, ostream* log )
{
  aligns_length = 0;

  for ( unsigned int mm_idx = 0; mm_idx < mm.size( ); ++mm_idx )
  {
    const mutmer &theMutmer = mm[mm_idx];

    // The mutmer must be perfect.
    if ( theMutmer.Errors() > 0 )
      continue;

    // The mutmer must start at zero on at least one read.
    if ( theMutmer.Pos1() > 0 && theMutmer.Pos2() > 0 )
      continue;

    // The mutmer must end at the end of at least one read.
    if ( theMutmer.End1() < (int) rd1.size() && theMutmer.End2() < (int) rd2.size() )
      continue;

    errors[aligns_length] = 0;
    align &theAlign = aligns[aligns_length++];
    
    theAlign.Setpos1( theMutmer.Pos1() );
    theAlign.Setpos2( theMutmer.Pos2() );
    theAlign.SetNblocks( 1 );
    theAlign.SetGap( 0, 0 );
    theAlign.SetLength( 0, theMutmer.Length() );
  }

  return True;    
}


// --- SWGAP METHOD

// The following class tracks mutmer starts and stops.

class mutmer_terminus {
  
public:

  enum begin_or_end {
    begin,
    end
  };

  mutmer_terminus() {}
  mutmer_terminus( int position, int mutmer_id, begin_or_end which_end )
    : pos_(position), mm_id_(mutmer_id), which_(which_end) { }

  int Position() const { return pos_; }
  int MutmerId() const { return mm_id_; }
  bool IsBegin() const { return which_ == begin; }
  bool IsEnd() const { return ! this->IsBegin(); }

  void SetPosition( int position ) { pos_ = position; }
  void SetMutmerId( int mutmer_id ) { mm_id_ = mutmer_id; }
  void SetIsBegin( begin_or_end which_end ) { which_ = which_end; }

  friend bool operator< ( const mutmer_terminus &lhs, const mutmer_terminus &rhs )
    { return ( ( lhs.pos_ < rhs.pos_ ) ||
	       ( lhs.pos_ == rhs.pos_ &&
		 lhs.mm_id_ < rhs.mm_id_ ) ||
	       ( lhs.pos_ == rhs.pos_ &&
		 lhs.mm_id_ == rhs.mm_id_ &&
		 lhs.which_ == begin && rhs.which_ == end ) ); }

  friend ostream & operator<< ( ostream &out, const mutmer_terminus &mt )
    { return out << mt.mm_id_ 
		 << ( mt.which_ == mutmer_terminus::begin ? "B" : "E" )
		 << "@" << mt.pos_; }

private:
  int pos_;
  int mm_id_;
  begin_or_end which_;
};

void add_or_remove_mutmer(const mutmer_terminus &term, set<int> &overlapping_mutmers) 
  // Adds or removes a mutmer reference from mut
{
  if ( term.IsBegin() )
  {
    pair<set<int>::iterator,bool> insert_result;
    insert_result = overlapping_mutmers.insert( term.MutmerId() );
    ForceAssertEq( insert_result.second, true );
  }
  else
  {
    int erase_result = overlapping_mutmers.erase( term.MutmerId() );
    ForceAssertEq( erase_result, 1 );
  }
}

set<int>::iterator find_maxmut(const vec<mutmer> &path, const set<int> &mut) 
  // Returns i such that path[mut[i]].Length() >= path[mut[j]].Length() for 0 <= j < mut.size()
{
  ForceAssertGt( mut.size(), 0u );

  set<int>::iterator max_iter = mut.begin();
  int maxlen = path[ *max_iter ].Length();

  set<int>::iterator mut_iter = max_iter;
  for ( ++mut_iter; mut_iter != mut.end(); ++mut_iter )
  {
    int mut_len = path[ *mut_iter ].Length();
    if ( mut_len > maxlen) {
      max_iter = mut_iter;
      maxlen = mut_len;
    }
  }
  return max_iter;
}

class progression {
  
public:
  progression()
    : length_( 0 ), offset_( 0 ),
      start_on_seq1_( 0 ), start_on_seq2_( 0 ),
      sum_mutmer_lens_( 0 ), max_mutmer_len_( 0 )
    { }

  progression( const mutmer &start_mm )
    : length_( start_mm.Length() ), offset_( start_mm.Offset() ),
      start_on_seq1_( start_mm.Pos1() ), start_on_seq2_( start_mm.Pos2() ),
      sum_mutmer_lens_( start_mm.Length() ), max_mutmer_len_( start_mm.Length() ),
      path_( 1, start_mm ) 
    { }

  int Length() const { return length_; }
  int Offset() const { return offset_; }
  int BeginOnSeq1() const { return start_on_seq1_; }
  int EndOnSeq1() const { return start_on_seq1_ + length_; }
  int BeginOnSeq2() const { return start_on_seq2_; }
  int EndOnSeq2() const { return start_on_seq2_ + length_; }
  int SumMutmerLengths() const { return sum_mutmer_lens_; }
  int MaxMutmerLength() const { return max_mutmer_len_; }

  void AddMutmer( const mutmer &mut );

  void TrimMutmers( const basevector &rd1, const basevector &rd2, int min_size,
		    ostream *log );

  int CreateAlign( align &a, int k, 
		   const basevector &rd1, const basevector &rd2,
		   int max_errs, int end_stretch, bool affine_penalties ) const;

  friend bool operator< ( const progression &lhs, const progression &rhs )
    {
      if ( lhs.length_ < rhs.length_ ) return true;
      if ( lhs.length_ > rhs.length_ ) return false;
      if ( lhs.offset_ < rhs.offset_ ) return true;
      if ( lhs.offset_ > rhs.offset_ ) return false;
      if ( lhs.start_on_seq1_ < rhs.start_on_seq1_ ) return true;
      if ( lhs.start_on_seq1_ > rhs.start_on_seq1_ ) return false;
      if ( lhs.start_on_seq2_ < rhs.start_on_seq2_ ) return true;
      return false;
    }

  friend ostream& operator<< ( ostream &out, const progression &prog )
    {
      return out << "PROG:\t" 
		 << prog.Length() << "\t" 
		 << prog.Offset() << "\t" 
		 << prog.BeginOnSeq1() << "-" 
		 << prog.EndOnSeq1() << "\t" 
		 << prog.BeginOnSeq2() << "-" 
		 << prog.EndOnSeq2();
    }

private:
  int length_;
  int offset_;
  int start_on_seq1_;
  int start_on_seq2_;
  int sum_mutmer_lens_;
  int max_mutmer_len_;

  vec<mutmer> path_;
};


void progression::AddMutmer( const mutmer &mut )
{
  int new_start1 = Min( this->BeginOnSeq1(), mut.Pos1() );
  int new_start2 = Min( this->BeginOnSeq2(), mut.Pos2() );
  int new_stop1 = Max( this->EndOnSeq1(), mut.End1() );
  int new_stop2 = Max( this->EndOnSeq2(), mut.End2() );
  int new_length = Max( new_stop1 - new_start1, new_stop2 - new_start2 );
  
  path_.push_back( mut );
  
  offset_ += (mut.Offset() - offset_) / (int) path_.size();
  start_on_seq1_ = new_start1;
  start_on_seq2_ = new_start2;
  length_ = new_length;
  sum_mutmer_lens_ += mut.Length();
  max_mutmer_len_ = max( max_mutmer_len_, mut.Length() );
}

void progression::TrimMutmers( const basevector &rd1, const basevector &rd2,
			       int min_size, ostream *log ) 
{
//    *log << *this << "\n";
//    *log << "UNTRIMMED PATH: " << "\n";
//    copy( path_.begin(), path_.end(),
//          ostream_iterator<mutmer>( *log, "\n" ) );
//    *log << "\n";

  // Build vector of endpoints
  vec<mutmer_terminus> ends1;
  ends1.reserve( path_.size() * 2 );
  for ( unsigned int mm_idx = 0; mm_idx < path_.size(); ++mm_idx )
  {
    mutmer &p = path_[mm_idx];
    ends1.push_back( mutmer_terminus( p.Pos1(), mm_idx, mutmer_terminus::begin ) );
    ends1.push_back( mutmer_terminus( p.End1(), mm_idx, mutmer_terminus::end ) );
  }
  sort( ends1.begin(), ends1.end() );

  // Walk through the vector of endpoints and trim overlapping ends.
  set<int> overlapping_mutmers;
  for(unsigned int end_idx = 0; end_idx < ends1.size(); ++end_idx) 
  {
    const mutmer_terminus &end1 = ends1[end_idx];
    
    // Add or remove the mutmer whose endpoint we have reached
    add_or_remove_mutmer(end1, overlapping_mutmers);
    
    // Trim the mutmers on each sequence
    if ( overlapping_mutmers.size() > 1 ) 
    {
      set<int>::iterator max_iter = find_maxmut( path_, overlapping_mutmers );
      mutmer &maxmut = path_[ *max_iter ];

      for ( set<int>::iterator notmax_iter = overlapping_mutmers.begin();
	    notmax_iter != overlapping_mutmers.end(); ++notmax_iter )
      {
	if ( notmax_iter == max_iter )
	  continue;
	
	mutmer &notmaxmut = path_[ *notmax_iter ];
	
	if ( notmaxmut.Length() == 0 )
	  continue;
	
	// If the maxmut starts after the notmaxmut...
	if ( notmaxmut.Pos1() < maxmut.Pos1() ) 
	{
	  int overlap = notmaxmut.End1() - maxmut.Pos1();
	  if ( overlap > 0 )
	    notmaxmut.TrimRight( overlap, rd1, rd2 );
	}
	// If the maxmut starts before the notmaxmut...
	else 
	{
	  int overlap = maxmut.End1() - notmaxmut.Pos1();
	  if ( overlap > 0 )
	    notmaxmut.TrimLeft( overlap, rd1, rd2 );
	}
      }
    }
  }
  ForceAssertEq( overlapping_mutmers.size(), 0u );

//    *log << "MID TRIM: " << *this << "\n";
    
//    *log << "TRIMMED PATH: " << "\n";
//    copy( path_.begin(), path_.end(),
//          ostream_iterator<mutmer>( *log, "\n" ) );
//    *log << "\n";

  // Build vector of endpoints
  vec<mutmer_terminus> ends2;
  ends2.reserve( path_.size() * 2 );
  for ( unsigned int mm_idx = 0; mm_idx < path_.size(); ++mm_idx )
  {
    mutmer &p = path_[mm_idx];
    ends2.push_back( mutmer_terminus( p.Pos2(), mm_idx, mutmer_terminus::begin ) );
    ends2.push_back( mutmer_terminus( p.End2(), mm_idx, mutmer_terminus::end ) );
  }
  sort( ends2.begin(), ends2.end() );
      
  // Walk through the vector of endpoints and trim overlapping ends.
  overlapping_mutmers.clear();
  for(unsigned int end_idx = 0; end_idx < ends2.size(); ++end_idx) 
  {
    const mutmer_terminus &end2 = ends2[end_idx];
    
    // Add or remove the mutmer whose endpoint we have reached
    add_or_remove_mutmer(end2, overlapping_mutmers);
    
    if ( overlapping_mutmers.size() > 1 ) 
    {
      set<int>::iterator max_iter = find_maxmut(path_, overlapping_mutmers);
      mutmer &maxmut = path_[ *max_iter ];
      
      for ( set<int>::iterator notmax_iter = overlapping_mutmers.begin();
	    notmax_iter != overlapping_mutmers.end(); ++notmax_iter )
      {
	if ( notmax_iter == max_iter )
	  continue;
	
	mutmer &notmaxmut = path_[ *notmax_iter ];

	if ( notmaxmut.Length() == 0 )
	  continue;
	
	// If the maxmut starts after the notmaxmut...
	if ( notmaxmut.Pos2() < maxmut.Pos2() ) 
	{
	  int overlap = notmaxmut.End2() - maxmut.Pos2();
	  if ( overlap > 0 )
	    notmaxmut.TrimRight( overlap, rd1, rd2 );
	}
	// If the maxmut starts before the notmaxmut...
	else 
	{
	  int overlap = maxmut.End2() - notmaxmut.Pos2();
	  if ( overlap > 0 )
	    notmaxmut.TrimLeft( overlap, rd1, rd2 );
	}
      }	     
    }
  }
  ForceAssertEq( overlapping_mutmers.size(), 0u );
  
  sort( path_.begin(), path_.end() );
      
  // Check for any remaining mutmers that have inconsistent offsets.
  int last_pos1 = -1, last_pos2 = -1;
  for ( unsigned int mm_idx = 0; mm_idx < path_.size(); ++mm_idx )
  {
    mutmer &this_mm = path_[mm_idx];
    if ( this_mm.Length() > 0 )
    {
      ForceAssert( this_mm.Pos1() > last_pos1 ||
		   this_mm.Pos2() > last_pos2) ;
      last_pos1 = this_mm.Pos1(); 
      last_pos2 = this_mm.Pos2();
    }
  }

  // Erase mutmers with length less than k.
  path_.erase( remove_if( path_.begin(), path_.end(),
			  mutmer_is_short( min_size ) ),
	       path_.end() );

  if ( path_.empty() )
  {
    length_ = 0;
  }
  else
  {
    int min_pos1 = path_.front().Pos1();
    int max_end1 = path_.front().End1();
    int min_pos2 = path_.front().Pos2();
    int max_end2 = path_.front().End2();

    for ( unsigned int mm_idx = 0; mm_idx < path_.size(); ++mm_idx )
    {
      mutmer &this_mm = path_[mm_idx];
      min_pos1 = min<int>( min_pos1, this_mm.Pos1() );
      max_end1 = max<int>( max_end1, this_mm.End1() );
      min_pos2 = min<int>( min_pos2, this_mm.Pos2() );
      max_end2 = max<int>( max_end2, this_mm.End2() );
    }

    start_on_seq1_ = min_pos1;
    start_on_seq2_ = min_pos2;

    length_ = max<int>( max_end1 - min_pos1, max_end2 - min_pos2 );
  }

//    *log << "AFTER TRIM: " << *this << "\n";
  
//    *log << "TRIMMED PATH: " << "\n";
//    copy( path_.begin(), path_.end(),
//          ostream_iterator<mutmer>( *log, "\n" ) );
//    *log << "\n";
}
  

int progression::CreateAlign( align &a, int k, 
			      const basevector& rd1, const basevector& rd2,
			      int max_errs, int end_stretch, bool affine_penalties ) const
{
  shortvector<mutmer> m( path_.size() );
  for( unsigned int m_idx = 0; m_idx < path_.size(); ++m_idx ) 
    m(m_idx).SetFrom( path_[m_idx] );

  int errors_found;
  a.CreateFromMutmersAndSW( k, m, rd1, rd2, max_errs, end_stretch, errors_found, affine_penalties );

  return errors_found;
}


Bool makealigns_sw_gap_method::MutmersToAlign( const vec<mutmer>& const_mm, int k,
					       const basevector& rd1, const basevector& rd2, 
					       vec<align>& aligns, vec<int>& errors, 
					       int& aligns_length, 
					       int min_mutmer, ostream* log )
{    
  // SW-Gap method
  //
  // Create vectors of non-overlapping mutmers with consistent offsets.
  // Then do a Smith-Waterman to fill in the gaps.
  //
  // Currently there is no edge extension.

  vec<mutmer> mm;

  mm = const_mm;
  sort( mm.begin(), mm.end() );

  // Track usable mutmers.
  vector<bool> tagged( mm.size(), false );

  // Remove mutmers with length less than min_mutmer.
  for ( unsigned int mm_idx = 0; mm_idx < mm.size( ); ++mm_idx ) 
    if ( mm[mm_idx].Length() < min_mutmer ) 
      tagged[ mm_idx ] = true;

  // Correct the number of errors (previously limited to 31 for
  // historical reasons).
  for ( unsigned int mm_idx = 0; mm_idx < mm.size( ); ++mm_idx ) 
    mm[mm_idx].CountErrors( rd1, rd2 );

  // Group all mutmers of similar offset into "paths", the total range of the
  // grouped mutmers are called a "prog(ression)". Discard any progression that
  // is shorter than min_prog_length.
  vec<progression> progs;
  for ( unsigned int mm_idx = 0; mm_idx < mm.size(); ++mm_idx ) 
  {
    if ( tagged[ mm_idx ] ) 
      continue;

    const mutmer &this_mm = mm[mm_idx];

    progs.push_back( progression( this_mm ) );

    progression &prog = progs.back();

    tagged[ mm_idx ] = true;

    int last_mm_idx = mm_idx;

    for ( unsigned int other_mm_idx = mm_idx+1; other_mm_idx < mm.size(); ++other_mm_idx ) 
    {
      if ( tagged[other_mm_idx] )
	continue;
      
      const mutmer &other_mm = mm[other_mm_idx]; 
      const mutmer &last_mm  = mm[last_mm_idx]; 

      // If this mutmer is too far away on both sequence 1 and
      // sequence 2, we skip it.
      if ( other_mm.Pos1() - prog.EndOnSeq1() > max_gap_ &&
           other_mm.Pos2() - prog.EndOnSeq2() > max_gap_ )
 	continue;

      // If this mutmer is starts before the last mutmer on sequence 1
      // and sequence 2, we skip it.
      if ( other_mm.Pos1() < last_mm.Pos1() &&
           other_mm.Pos2() < last_mm.Pos2() )
          continue;

      // Add this mutmer to this progression if its offset is close
      // enough.

      int offset_diff = other_mm.Offset() - prog.Offset();

      if ( abs( offset_diff ) <= max_mutmer_offset_diff_ )
      {
	prog.AddMutmer( other_mm );
	tagged[ other_mm_idx ] = true;
        last_mm_idx = other_mm_idx;
      }
    }
    
    // Only use paths with at least one mutmer of length greater
    // than or equal to min_mutmerlength.
    
    if ( prog.MaxMutmerLength() <= min_max_mutmer_length_ )
    {
      if ( verbose_ && log != 0 )
      {
	*log << "Discarding progression due to insufficient maxmutmer length" << "\n";
	*log << prog << "\n";
      }
      progs.pop_back();
      continue;
    }

    // Test prog for inclusion
    if ( prog.Length() < min_prog_length_ ) 
    {
      if ( verbose_ && log != 0 )
      {
	*log << "Discarding progression due to insufficient length: "
	    << prog.Length() << " < " << min_prog_length_ << "\n";
	*log << prog << "\n";
      }
      progs.pop_back();
      continue;
    }
    
    double prog_ratio = ( double( prog.SumMutmerLengths() ) / double( prog.Length() ) );
    if ( prog_ratio < min_prog_ratio_ ) 
    {
      if ( verbose_ && log != 0 )
      {
	*log << "Discarding progression due to insufficient ratio: "
	    << prog_ratio << " < " << min_prog_ratio_ << "\n";
	*log << prog << "\n";
      }
      progs.pop_back();
      continue;
    }

    if ( prog.Length() < ignore_overlap_frac_length_ )
    {
      const int length1 = rd1.size( ), length2 = rd2.size( );
      const int implied_overlap = 
	IntervalOverlap(0, length1, prog.Offset(), prog.Offset() + length2);
      if ( (float) prog.Length() < implied_overlap * min_overlap_frac_ )
      {
	if ( verbose_ && log != 0 ) 
	{
	  *log << "Discarding progression due to insufficient overlap: "
               << prog.Length() << " < " << implied_overlap << " * " << min_overlap_frac_ << "\n";
	  *log << prog << "\n";
	}
	progs.pop_back();
	continue;
      }
    }

    if ( verbose_ && log != 0 )
    {
      *log << "ACCEPTED:" << "\n";
      *log << prog << "\n";
    }
  }

  // For each create two new vectors holding their sorted endpoints. 
  // Trim away, from the shorter mutmers, any overlap between two or more mutmers. 

  const int min_size = k/2;

  vec<progression>::iterator prog_iter;

  for ( prog_iter = progs.begin(); prog_iter != progs.end(); ++prog_iter ) 
    prog_iter->TrimMutmers( rd1, rd2, min_size, log );

  // Feed each path of non-overlapping mutmers into a function that heuristically
  // closes the gaps between them.
  align a;
  int errors_found;
  aligns_length = 0;

  for ( prog_iter = progs.begin(); prog_iter != progs.end(); ++prog_iter ) 
  {
    if ( prog_iter->Length() < min_prog_length_ )
    {
      if ( verbose_ && log != 0 )
      {
	*log << "Discarding progression due to insufficient length: "
	    << prog_iter->Length() << " < " << min_prog_length_ << "\n";
	*log << *prog_iter << "\n";
      }
      continue;
    }

    if ( verbose_ && log != 0 )
    {
      *log << "Making align with progression:" << "\n";
      *log << *prog_iter << "\n";
    }

    errors_found = prog_iter->CreateAlign( a, k, rd1, rd2, max_errs_, 
                                           end_stretch_, affine_penalties_ );

//     PRINT( errors_found );
//     PrintVisualAlignment( False, *log, rd1, rd2, a );

    if ( aligns_length == (int) aligns.size( ) )
    {
      if ( log != 0 )
      {    *log << "MutmersToAlign: over " << aligns.size( ) 
	        << " alignments were found while aligning two sequences "
	        << "of lengths " << rd1.size( ) << " and "
	        << rd2.size( ) << ".\n";
           *log << "(called with min_mutmer = " << min_mutmer << ")\n";
           *log << "Here are the stats on the progressions:\n";
           for ( unsigned int ii = 0; ii < progs.size( ); ++ii )
           {    *log << "\n";
	        *log << "Progression " << ii << ": ";
	        *log << progs[ii] << "\n";    }
           *log << flush;    }
      return False;    
    }

    if ( a.Nblocks( ) == 0 )
    {  
      if( verbose_ && log != 0 ) 
      {
        *log << *prog_iter << "\n";
        *log << "ALIGNMENT DISCARDED, ERRORS:\t" << errors_found << "\n"; 
      }
      continue;
    }

    if ( a.pos1() < 0 || 
	 a.pos2() < 0 ||
	 a.Pos1() > (int) rd1.size() ||
	 a.Pos2() > (int) rd2.size() )
    {
      PRINT2( a.pos1(), a.Pos1() );
      PRINT2( a.pos2(), a.Pos2() );
      cout << "(gaps) lengths: ";
      for ( int b = 0; b < a.Nblocks(); ++b )
	cout << "(" << a.Gaps(b) << ") " << a.Lengths(b) << " ";
      cout << endl;
      
      ForceAssertGe( a.pos1( ), 0 );
      ForceAssertGe( a.pos2( ), 0 );
      ForceAssertLe( a.Pos1( ), (int) rd1.size( ) );
      ForceAssertLe( a.Pos2( ), (int) rd2.size( ) );
    }
    
    errors[aligns_length] = errors_found;
    aligns[aligns_length++] = a;    
  }

  return True;        
}
