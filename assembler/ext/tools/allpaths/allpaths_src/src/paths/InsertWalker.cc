/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This file contains the InsertWalkerLoc and Jaunt class definition and implementation,
// as well as the implementation for the function InsertWalker::WalkUnipaths().

#include "paths/InsertWalker.h"

#include <map>
#include <set> // multiset
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h" // BuildUnipathAdjacencyHyperKmerPath
#include "TaskTimer.h" // TaskTimer




/*******************************************************************************
 *
 *        Jaunt class
 *
 * A jaunt is a walk in progress.  It contains the following data:
 * -- multiset<int> (inherited): A set of unipath IDs that have been visited.
 * -- int _dist: The total distance traveled so far.
 * -- uint _n_repeats: The number of times this Jaunt visits a unipath it's
 *                    already visited.
 *
 * For a further explanation of this class, see the documentation for
 * InsertWalkerLoc below.
 *
 ******************************************************************************/
class Jaunt : public multiset<int>
{
  friend class InsertWalkerLoc;
  
public:
  
  // Constructor: Create an empty jaunt from a separation size.
  Jaunt( const int dist ) : _dist( dist ), _n_repeats(0) {}
  
  
  void addDist( const int x ) { _dist += x; }
  int dist() const { return _dist; }
  
  void tallyRepeat() { _n_repeats++; }
  uint nRepeats() const { return _n_repeats; }
  
  // Combine another Jaunt into this one by taking take the union of the two
  // Jaunts' unipath ID lists.
  // This can only be done if the two Jaunts have the same _dist.
  void mergeJaunt( const Jaunt& j ) {
    ForceAssertEq( _dist, j._dist );
    for ( Jaunt::iterator iter = j.begin();
	  iter != j.end();
	  iter++ ) {
      size_t this_count = count( *iter ), that_count = j.count( *iter );
      
      // This loop criterion may immediately fail
      for ( size_t i = this_count; i < that_count; i++ ) {
	insert( *iter );
	if ( i > 0 ) ++_n_repeats;
      }
    }
  }
  
  // Return true iff:
  // 1. This unipath's distance not greater than j's distance
  // 2. Every unipath visited by j is also visited by this unipath, with equal
  //    or smaller frequency
  // These criteria tend to be satisfied if there is a loop which this unipath
  // traverses once and j traverses more than once.
  bool simplifies( const Jaunt& j ) const {
    if ( _dist > j._dist ) return false;
    
    for ( Jaunt::iterator iter = j.begin();
	  iter != j.end();
	  iter++ ) {
      size_t c = count( *iter );
      if ( c == 0 || c > j.count( *iter ) ) return false;
    }
    return true;
  }
  
  // Comparison operators.
  // Note that these operators are implicitly used for ordering Jaunts in the
  // set<Jaunt> in the InsertWalkerLoc class.
  friend bool operator<( const Jaunt& a, const Jaunt& b ) {
    if ( a._dist  != b._dist  ) return ( a._dist  < b._dist  );
    if ( a.size() != b.size() ) return ( a.size() < b.size() );
    
    Jaunt::iterator ai = a.begin(), bi = b.begin();
    while ( ai != a.end() ) {
      if ( (*ai) != (*bi) ) return ( (*ai) < (*bi) );
      ai++;
      bi++;
    }
    return false; // if control reaches here, the jaunts are equal!
  }
  
  friend bool operator==( const Jaunt& a, const Jaunt& b ) {
    if ( a._dist  != b._dist  ) return false;
    if ( a.size() != b.size() ) return false;
    
    Jaunt::iterator ai = a.begin(), bi = b.begin();
    while ( ai != a.end() ) {
      if ( (*ai) != (*bi) ) return false;
      ai++;
      bi++;
    }
    return true;
  }
  
  
private:
  
  int _dist; // The total distance traveled by this Jaunt so far
  uint _n_repeats; // The number of times this Jaunt visits a unipath it's already visited
  
};







/*******************************************************************************
 *
 *        InsertWalkerLoc (IWL)
 *
 * An InsertWalkerLoc (sometimes abbreviated IWL) is a helper class for the
 * InsertWalker.  Any insert walks that require moving between unipaths (as
 * opposed to inserts that lie entirely on a single unipath) are walked here.
 *
 * THE INSERT WALK
 * In order to walk an insert whose reads fall on different unipaths, we must
 * walk through the unipath graph, from the starting unipath to the target
 * unipath.  This is an insert walk.  An ongoing insert walk will consist of one
 * or more active InsertWalkLocs, each of which represents a particular unipath
 * (_uID) at which the insert walk is proceeding.
 *
 * JAUNTS
 * Each InsertWalker contains one or more Jaunts.  A Jaunt is an incomplete
 * walk: it represents a path that may be taken from the starting unipath to the
 * current unipath.  Each Jaunt consists of an integer amount of distance
 * traveled, along with a multiset of unipath IDs visited (not including the
 * starting and target unipaths.)  A Jaunt may visit the same unipath multiple
 * times; the order of visiting does not matter.  If a Jaunt travels too far
 * to make a reasonable separation (_sep_max), it is deleted.  There are two
 * heuristic parameters (MAX_VISITS_PER_UNIPATH, MAX_REPEAT_VISITS_PER_JAUNT)
 * to prevent the number of Jaunts from exploding due to loops.
 *
 * ADVANCING THE WALK
 * When an InsertWalkerLoc advances (advance()), it progresses from one unipath
 * to that unipath's successors in the adjacency graph.  If there are no
 * successors, the InsertWalkerLoc dies.  If there is more than one successor,
 * the InsertWalkerLoc splits, each child inheriting its complete list of
 * Jaunts.
 *
 * MERGING
 * The InsertWalker maintains a list of InsertWalkerLocs.  It will notice when
 * multiple InsertWalkerLocs converge onto a single unipath (e.g. at the closing
 * of a bubble) and will merge them.  This enables the walk to progress through
 * knotty parts of the unipath graph without the number of InsertWalkLocs
 * blowing up: instead, the number of Jaunts per InsertWalkLoc will grow large.
 * This is kept to a reasonable number by MAX_JAUNTS_PER_LOC.
 *
 * COMPLETING THE WALK
 * When an InsertWalkerLoc reaches a unipath whose successors include the target
 * unipath, it has found a walk.  The Jaunts in this InsertWalkerLoc become
 * completed walks, assuming they have traveled far enough (_sep_min).
 * Note that the InsertWalkerLoc remains live and may advance onto and beyond
 * the target unipath, in case it loops back to the target later.
 *
 ******************************************************************************/
class InsertWalkerLoc {
  
public:
  
  InsertWalkerLoc() {}
  ~InsertWalkerLoc() {}
  
  // Initial constructor.
  // Specifies the parameters for this insert walk.
  // Use this constructor only for the InsertWalkerLoc at the start of a walk.
  InsertWalkerLoc( const vecKmerPath * unipaths,
		   const digraph * unipath_adjs,
		   const int * sep_min, const int * sep_max,
		   const int start_uID, const int start_sep,
		   const int target_uID )
    : _unipaths( unipaths ),
      _adjs( unipath_adjs ),
      _sep_min( sep_min ),
      _sep_max( sep_max ),
      _target_uID( target_uID )
  {
    _uID = start_uID;
    
    // Create a single Jaunt.
    // By convention, the set of unipaths visited by a Jaunt does NOT include
    // the starting or ending unipath.
    Jaunt j( start_sep );
    
    
    // If possible, walk the target unipath backwards until the graph branches.
    // This is purely a time-saving step.
    set<int> visited_uIDs;
    while ( (*_adjs).To( _target_uID ).solo() ) {
      int prev_uID = (*_adjs).To( _target_uID )[0];
      
      // Be sure to avoid getting caught in an infinite loop of backward walking!
      pair< set<int>::iterator, bool > insert_value = visited_uIDs.insert( prev_uID );
      if ( insert_value.second == false ) break; // i.e., if the set already contained it
      if ( _target_uID == prev_uID ) break;
      
      _target_uID = prev_uID;
      j.insert( _target_uID );
      j.addDist( (*_unipaths)[_target_uID].KmerCount() );
    }
    
    // Filter out cases where the target is inaccessible: these inserts can
    // never be walked.
    if ( (*_adjs).To( _target_uID ).empty() ) return;
    if ( (*_adjs).To( _target_uID ).solo() &&
	 (*_adjs).To( _target_uID )[0] == _target_uID ) return;
    
    // Add this Jaunt to the list. 
    _jaunts.insert( make_pair( j.dist(), j ) );
  }
  
  // Splitting constructor.
  // When an InsertWalkerLoc splits into two or more at a unipath junction, each
  // of the new InsertWalks is created using this constructor.
  InsertWalkerLoc( const InsertWalkerLoc & loc, const int new_uID )
    : _unipaths  ( loc._unipaths ),
      _adjs      ( loc._adjs ),
      _sep_min   ( loc._sep_min ),
      _sep_max   ( loc._sep_max ),
      _target_uID( loc._target_uID )
  {
    // Make a copy of the original InsertWalkerLoc, then call advance().
    _uID = loc._uID;
    _jaunts = loc._jaunts;
    
    advance( new_uID );
  }
  
  
  
  // Merging constructor.
  // This is invoked whenever two InsertWalkerLocs advance from different
  // unipaths into the same unipath - e.g., at the closing of a bubble.
  InsertWalkerLoc( const InsertWalkerLoc& loc1, const InsertWalkerLoc & loc2 )
    : _unipaths  ( loc1._unipaths ),
      _adjs      ( loc1._adjs ),
      _sep_min   ( loc1._sep_min ),
      _sep_max   ( loc1._sep_max ),
      _target_uID( loc1._target_uID )
  {
    ForceAssertEq( loc1._uID, loc2._uID );
    _uID = loc1._uID;
    _jaunts = loc1._jaunts;
    
    // Merge the second set of Jaunts into the first set.
    // Find whether any Jaunts in the first set have the same distance as any
    // Jaunts in the second set.  If so, merge them.
    for ( map<int, Jaunt>::const_iterator iter = loc2._jaunts.begin();
	  iter != loc2._jaunts.end();
	  iter++ ) {
      
      map<int, Jaunt>::iterator j = _jaunts.find( (*iter).first );
      if ( j == _jaunts.end() ) // if no match (typical)
	_jaunts.insert( *iter );
      else
	(*j).second.mergeJaunt( (*iter).second );
    }
    
    
    // Look for redundant Jaunts that we can remove.
    // Specifically:  Suppose we have a small local loop in our unipath graph.
    // There will eventually be two Jaunts j1, j2 which traverse this loop 1 and
    // 2 times, respectively.  j2 is longer than j1 and it has no information
    // (no unipaths visited) that j1 lacks.  If j1.dist >= *_sep_min, we have no
    // reason to keep j2.
    map<int, Jaunt>::iterator j1, j2, j2_temp;
    for ( j1 = _jaunts.begin(); j1 != _jaunts.end(); j1++ ) {
      if ( (*j1).first < *_sep_min ) continue;
      
      j2 = j1;
      j2++;
      while ( j2 != _jaunts.end() ) {
	if ( (*j1).second.simplifies( (*j2).second ) ) {
	  j2_temp = j2;
	  j2--;
	  _jaunts.erase( j2_temp );
	}
	j2++;
      }
    }
    
  }
  
  
  /*****************************************************************************
   *
   *                 ACTION FUNCTIONS
   *
   ****************************************************************************/
  
  // Advance this InsertWalkerLoc from its current unipath to a successor.
  // This may decrease the number of jaunts because it will remove jaunts that
  // have traveled too far.
  void advance( int new_uID );
  
  
  /*****************************************************************************
   *
   *                 QUERY FUNCTIONS
   *
   ****************************************************************************/
  
  int loc() const { return _uID; }
  bool hasJaunts() const { return ( !_jaunts.empty() && _jaunts.size() <= _MAX_JAUNTS_PER_LOC ); }
  int minJauntDistance( ) const;
  
  // Get the list of unipaths that succeed this InsertWalkerLoc's unipath.
  // This InsertWalkerLoc may advance to any of these unipaths, although any or
  // all of these advances may leave it with no acceptable jaunts.
  vec<int> getSuccessors() const { return (*_adjs).From( _uID ); }
  
  // Return a list of jaunts that have made it to the target unipath.
  // If the target unipath is not a successor to the current unipath, this
  // returns empty.
  vec<Jaunt> getCompletedWalks() const;
  
  void print( ostream& out = cout ) const;
  
  
  
private:
  
  // Input variables.
  // These are all const, and none of them are locally owned.
  const vecKmerPath * _unipaths;
  const digraph * _adjs; // unipath adjacency graph
  const int * _sep_min, * _sep_max; // allowable range of separatio
  
  // Local variables.
  int _uID; // ID of current unipath
  int _target_uID; // ID of unipath containing read2 - when we get here, the walk is done!
  map<int, Jaunt> _jaunts; // The set of jaunts at this location
  
  // Heuristic parameters.
  static const size_t _MAX_JAUNTS_PER_LOC = 1000;
  static const size_t _MAX_VISITS_PER_UNIPATH = 2;
  static const size_t _MAX_REPEAT_VISITS_PER_JAUNT = 5;
};







// Advance this InsertWalkerLoc from its current unipath to a successor.
// This may decrease the number of jaunts because it will remove jaunts that
// have traveled too far.
void
InsertWalkerLoc::advance( int new_uID ) {
  ForceAssert( (*_adjs).HasEdge( _uID, new_uID ) );
  _uID = new_uID;
  int length = (*_unipaths)[_uID].KmerCount();
  
  // Update all of the jaunts to include this new unipath's ID and length.
  // We must re-create the set of jaunts because a map::iterator is const.
  map<int, Jaunt> new_jaunts;
  for ( map<int, Jaunt>::iterator iter = _jaunts.begin();
	iter != _jaunts.end();
	iter++ ) {
    
    Jaunt j = (*iter).second;
    j.addDist( length );
    // Remove any jaunts that have passed through so many kmers that they are no
    // longer within the tolerated range of the insert.  (We can use 'break'
    // instead of 'continue' because the set of Jaunts is sorted by _dist.)
    if ( j.dist() > *_sep_max ) break;
    
    j.insert( _uID );
    
    // Remove any jaunts that have visited the same unipath too many times.
    // This is necessary to prevent the algorithm from blowing up in the
    // presence of small loops.
    size_t count = j.count( _uID );
    if ( count > _MAX_VISITS_PER_UNIPATH ) continue;
    
    
    // Remove any jaunts that have accrued too many repeat visits across all
    // unipaths.  This prevents the algorithm from blowing up in the presence
    // of repeated small loops.
    if ( count > 1 ) {
      j.tallyRepeat();
      if ( j.nRepeats() >= _MAX_REPEAT_VISITS_PER_JAUNT ) continue;
    }
    
    new_jaunts.insert( make_pair( j._dist, j ) );
  }
  
  _jaunts = new_jaunts;
}




// Return the distance of the Jaunt with the smallest distance.
int
InsertWalkerLoc::minJauntDistance( ) const
{
  int min_dist = INT_MAX;
  for ( map<int, Jaunt>::const_iterator iter = _jaunts.begin();
	iter != _jaunts.end();
	++iter )
    if ( min_dist > (*iter).first )
      min_dist = (*iter).first;
  return min_dist;
}


// Return a list of jaunts that have made it to the target unipath.
vec<Jaunt>
InsertWalkerLoc::getCompletedWalks() const {
  vec<Jaunt> walks;
  
  // If we are not yet approaching the target, there are no complete walks!
  if ( !(*_adjs).HasEdge( _uID, _target_uID ) ) return walks;
  
  for ( map<int, Jaunt>::const_iterator iter = _jaunts.begin();
	iter != _jaunts.end();
	++iter )
    if ( (*iter).first >= *_sep_min )
      walks.push_back( (*iter).second );
  
  return walks;
}



void
InsertWalkerLoc::print( ostream& out ) const
{
  // Print basic info about this walk.
  out << "InsertWalkerLoc at unipath #" << _uID << " (length: " << (*_unipaths)[_uID].KmerCount() << ")" << endl;
  out << "Targeting unipath #" << _target_uID << " and aiming for a separation in the range [" << *_sep_min << "-" << *_sep_max << "]" << endl;
  
  // Print a chart of information about each jaunt in progress.
  out << "Jaunt#\tDist.\tUnipaths visited (not including start/stop)" << endl;
  size_t count = 0;
  for ( map<int, Jaunt>::const_iterator iter = _jaunts.begin();
	iter != _jaunts.end();
	iter++ ) {
    out << count++ << "\t" << (*iter).first;
    for ( Jaunt::const_iterator iter2 = (*iter).second.begin();
	  iter2 != (*iter).second.end();
	  ++iter2 )
      out << "\t" << (*iter2);
    out << endl;
  }
}














/*******************************************************************************
 *
 * InsertWalker constructor
 *
 * Load the data structures (as const pointers) and do some pre-processing.
 *
 ******************************************************************************/
InsertWalker::InsertWalker( const int K, 
			    const vecKmerPath * unipaths,
			    const vec<tagged_rpint> * unipathsdb,
			    const digraph * unipath_adjs )
  : _K( K ),
    _unipaths( unipaths ),
    _unipathsdb( unipathsdb ),
    _unipath_adjs( unipath_adjs )
{
  // Sanity check.
  ForceAssertEq( (int)_unipaths->size(), _unipath_adjs->N() );
  
  // Pre-processing step: Find the component ID of each local unipath.
  // This saves time later.
  _unipath_component.resize( _unipaths->size() );
  vec< vec<int> > comps;
  _unipath_adjs->Components( comps );
  for ( size_t i = 0; i < comps.size(); i++ )
    for ( size_t j = 0; j < comps[i].size(); j++ )
      _unipath_component[ comps[i][j] ] = i;
}



/*******************************************************************************
 *
 * InsertWalker::PlaceReadsOnUnipaths
 *
 * Locate a pair of reads on the unipath graph.
 *
 * Output:
 * u1, u2: The IDs of the unipaths on which read1 and read2 land, respectively.
 * May be the same.
 * sep_offset: The distance between the start of u1 and the start of u2,
 * relative to the (non-inclusive) distance between read1 and read2.  May be
 * positive or negative.  Will be used in walking between u1 and u2.
 *
 *
 ******************************************************************************/

void
InsertWalker::PlaceReadsOnUnipaths( const KmerPath & read1, const KmerPath & read2, int & u1, int & u2, int & sep_offset ) const
{
  u1 = -1, u2 = -1, sep_offset = -1;
  
  // Locate these reads (via their first kmer) on the unipath graph.
  // Each kmer should appear no more than once, because this is a UNIpath graph;
  // and reads should never fail to appear, if the unipaths were made from them.
  vec<longlong> places1, places2;
  Contains( *_unipathsdb, read1.Start(), places1 );
  Contains( *_unipathsdb, read2.Start(), places2 );
  if ( places1.empty() || places2.empty() ) return; // this should never happen
  ForceAssertEq( places1.size(), 1u );
  ForceAssertEq( places2.size(), 1u );
  
  const tagged_rpint &place1 = (*_unipathsdb) [ places1[0] ];
  const tagged_rpint &place2 = (*_unipathsdb) [ places2[0] ];
  u1 = place1.PathId();
  u2 = place2.PathId();
  
  // Determine exactly where each read lands on its unipath.
  int read_start1 = read1.Start() - place1.Start();
  int read_start2 = read2.Start() - place2.Start();
  for ( size_t i = 0; i < place1.PathPos(); i++ )
    read_start1 += (*_unipaths)[u1].Length(i);
  for ( size_t i = 0; i < place2.PathPos(); i++ )
    read_start2 += (*_unipaths)[u2].Length(i);
  
  // Calculate a separation offset for this walk between unipaths u1 and u2.
  // Note that a kmer separation of N is equal to a base separaton of N-K+1.
  sep_offset = -_K + 1;
  sep_offset -= ( read_start1 + read1.KmerCount() ); // subtract read1's offset
  sep_offset += read_start2; // add read2's offset
}


  
/*******************************************************************************
 *
 * InsertWalker::WalkUnipaths
 *
 * This is the main workhorse function in the InsertWalker class.
 *
 * Attempt to bridge the separation between two unipaths u1, u2 by walking along
 * the unipath graph.  The output is a set of unipath IDs that have been walked.
 *
 * The algorithm is as follows:
 * 2. If u1 = u2, and the input sep_offset matches the expected separation,
 *    return the unipath.
 * 3. If the reads are on different unipaths, walk through the unipath graph
 *    from u1 to u2, with the help of an InsertWalkerLoc object.  (See below for
 *    a description of this algorithm.)  Return the set of all unipaths
 *    encountered in all complete walks.
 *
 *
 * INPUT:
 * u1, u2: The IDs of the two unipaths that the reads land on.
 * sep_range: The allowable separation range between the *start* of u1 and the
 * *start* of u2, based on the read placements.
 * min_sep, max_sep: The range of allowable insert separations (in bases)
 * between the reads, not counting the read lengths themselves.
 *
 * Note that the order and orientation of the reads matters: the reads must be
 * pointing in the same direction, with read1 before read2.  If you have a pair
 * of reads (A,B) that are pointing toward each other, you can orient them
 * properly by flipping B and then passing in read1=A, read2=B.
 *
 *                 read1                          read2
 *             ------------->                 ------------->
 *
 *
 *
 * The insert may fail to be walked for any of the following reasons:
 * -- Both reads land on the same unipath (u1=u2) with a distance between them
 *    that is outside of the tolerance.  This read pair is probably chimeric.
 * -- The reads land so far from the edges of their unipaths that, even if u1,u2
 *    were adjacent, the separation would be too big.
 * -- No insert walks can be found by advancing through the unipath graph from
 *    u1 to u2.
 *
 * In any of these events, WalkUnipaths() returns an empty set.
 * Each of these bullets corresponds to a 'return null_set' line below.
 *
 ******************************************************************************/
set<int>
InsertWalker::WalkUnipaths( const int & u1, const int & u2,
			    const ho_interval & sep_range,
			    TaskTimer & timer,
			    const bool verbose ) const
{
  // Refuse to walk insert if there are too many iterations.
  const int MAX_ITERATIONS = 5000;

  set<int> null_set;

  // Log attempted walk.
  if ( verbose )
    cout << "Trying to walk " << u1 << " (" << (*_unipaths)[u1].TotalLength( )
	 << " kmers) -> " << u2 << " (" << (*_unipaths)[u2].TotalLength( )
	 << " kmers). Allowed range: [" << sep_range.Start( )
	 << ", " << sep_range.Stop( ) << ")" << endl;
  
  // CASE: Both reads are on the same unipath.
  // Remember that we expect to see read1 fall BEFORE read2 on the unipath.
  if ( u1 == u2 ) {
    
    // If the allowable separation range includes no separation at all,
    // then we have an immediate solution!
    // If not, then we must ignore this pair; it's probably chimeric.
    if ( ! sep_range.Contains( 0 ) ) {
      if ( verbose )
	cout << " FAIL: Reads are on same unipath but wrong distance" << endl;
      return null_set;
    }
    
    // Report this one-unipath "walk".
    if ( verbose )
      cout << " WALKED: found a one-unipath walk!" << endl;
    set<int> walk;
    walk.insert( u1 );
    return walk;
  }
  
  // CASE: The reads are on different unipaths.  This is the hard case.
  int sep_traveled = (*_unipaths)[u1].KmerCount();
  
  // If u1 and u2 are on different components in the adjacency graph,
  // there can't possibly be a walk between them.
  if ( _unipath_component[u1] != _unipath_component[u2] ) {
    if ( verbose )
      cout << " FAIL: unipaths are on different components" << endl;
    return null_set;
  }
  
  if ( sep_traveled > sep_range.Stop() ) {
    if ( verbose )
      cout << " FAIL: Reads are so far from the ends of their unipaths that they could not possibly form a tolerable insert length (" << sep_traveled << " > " << sep_range.Stop() << ")" << endl;
    return null_set;
  }
  
  // This will be a list of all of the unipaths involved in walks from the end
  // of u1 to the start of u2.  Note that u1 and u2 may themselves be included
  // here, but only if a walk loops through them.
  set<int> unipaths_in_walks;
  
  // Set up an InsertWalkerLoc (IWL) for finding all the possible walks.
  int sep_start = sep_range.Start(), sep_stop = sep_range.Stop();
  InsertWalkerLoc initial_loc( _unipaths, _unipath_adjs, &sep_start, &sep_stop, u1, sep_traveled, u2 );
  vec<InsertWalkerLoc> IWLs( 1, initial_loc );
  bool found_walk = false;
  vec<Bool> IWLs_to_erase(1, False );
  /* Iterative loop for progessively walking along the unipath graph.
   *
   * In each iteration, one InsertWalkerLoc is chosen to advance from its
   * current unipath.  This will cause one of the following events, based on
   * whether the Loc's unipath has any successor unipaths:
   * 1.  No successor unipaths.  The InsertWalkerLoc dies.
   * 2.  One or more successor unipaths.  Perform the following three actions:
   * 2A. The InsertWalkerLoc splits into one or more InsertWalkerLocs, one for
   *     each successor.
   * 2B. Each new InsertWalkLoc is checked for overreach.  If all of the jaunts
   *     (i.e., potential walks) in an InsertWalkLoc have traveled beyond the
   *     range of the insert we are trying to walk, the InsertWalkLoc dies.
   * 2C. The set of InsertWalkLocs is checked for mergers.  If any two
   *     InsertWalkLocs are sitting on the same unipath, they will be merged.
   *
   * We ultimately compile a list of all of the unipaths traversed by all of the
   * completed insert walks.
   *
   ****************************************************************************/
  int n_iterations = 0;
  while( Sum(IWLs_to_erase) < IWLs_to_erase.isize() ) {
    n_iterations++;
    if ( n_iterations > MAX_ITERATIONS ) {
      if ( verbose )
	cout << " FAIL: MAX_ITERATIONS limit reached" << endl;
      return null_set;
    }
    
    // Find the InsertWalkLoc that contains the shortest overall jaunt.  This
    // will be the InsertWalkLoc we attempt to advance.
    int min_jaunt = INT_MAX, min_jaunt_ID = 0;
    for ( size_t i = 0; i < IWLs.size(); i++ ){
      if ( IWLs_to_erase[i] ) continue;
      if ( min_jaunt > IWLs[i].minJauntDistance() ) {
	min_jaunt = IWLs[i].minJauntDistance();
	min_jaunt_ID = i;
      }
    }
    const InsertWalkerLoc IWL_to_advance = IWLs[min_jaunt_ID];
    
    // Check if there are any completed walks coming out of this InsertWalkLoc.
    // This will only happen if the target unipath (u2) is one of the
    // successors to the InsertWalkLoc's current unipath location.
    // Note that if u2 is a direct successor of u1, there will be an empty walk.
    vec<Jaunt> walks = IWL_to_advance.getCompletedWalks();
    if ( walks.nonempty() ) found_walk = true;
    for ( size_t i = 0; i < walks.size(); i++ )
      unipaths_in_walks.insert( walks[i].begin(), walks[i].end() );
    
    // Case 1: No successor unipaths.  The InsertWalkLoc dies.
    vec<int> successors = IWL_to_advance.getSuccessors();
    size_t n_old_locs = IWLs.size(), n_successors = successors.size();
    if ( n_successors == 0 ) {
      IWLs_to_erase[ min_jaunt_ID] = True;
      continue;
    }
    
    // Case 2: One or more successor unipaths.
    
    // 2A. Split the InsertWalkLoc, creating a new InsertWalkLoc for each
    // successor unipath.
    IWLs.resize( n_old_locs + n_successors );
    for ( size_t i = 0; i < n_successors; i++ )
      IWLs[ n_old_locs + i ] = InsertWalkerLoc( IWL_to_advance, successors[i] );
    
    IWLs_to_erase.resize( n_old_locs + n_successors, False );
    IWLs_to_erase[min_jaunt_ID] = True;
    
    // 2B. Check each new InsertWalkLoc for overreach.
    for ( size_t i = n_old_locs; i < n_old_locs + n_successors; i++ )
      if ( !IWLs[i].hasJaunts() )
	IWLs_to_erase[i] = True;
    
    // 2C. Check each new InsertWalkLoc against each old InsertWalkLoc for
    // possible merges.
    // Note that it should never be necessary (or possible) to merge two new
    // IWLs, or a set of more than two IWLs.
    for ( size_t i = 0; i < n_old_locs; i++ ) {
      if ( IWLs_to_erase[i] ) continue;
      for ( size_t j = n_old_locs; j < n_old_locs + n_successors; j++ ) {
	if ( IWLs_to_erase[j] ) continue;
	// Found a merge!  Create a merged InsertWalkLoc, and mark the
	// pre-merge IWLs for deletion.
	if ( IWLs[i].loc() == IWLs[j].loc() ) {
	  IWLs[j] = InsertWalkerLoc( IWLs[i], IWLs[j] );
	  IWLs_to_erase[i] = True;
	  break;
	}
      }
    }
  }
  
  if ( !found_walk ) {
    if ( verbose )
      cout << " FAIL: Didn't find any walks in the unipath graph" << endl;
    return null_set;
  }
  
  // We now have a set of all the unipaths that were traversed in insert walks.
  // Add the starting and target unipaths to the set.
  unipaths_in_walks.insert( u1 );
  unipaths_in_walks.insert( u2 );
  
  if ( verbose ) {
    cout << " WALKED: found a multi-unipath walk!\n"
	 << " Total set of unipaths in walks: " << unipaths_in_walks.size( )
	 << " (unipaths in walks:";
    set<int>::iterator iter;
    for ( iter = unipaths_in_walks.begin();
	  iter != unipaths_in_walks.end();
	  ++iter )
      cout << " " << *iter;
    cout << ")" << endl;
  }
  
  return unipaths_in_walks;
}

