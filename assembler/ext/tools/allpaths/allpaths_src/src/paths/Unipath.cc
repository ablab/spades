/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Equiv.h"
#include "math/Functions.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"


///////////////////////////////////////////////////////////////////////
//
// Stuff for Unipath(...)
//
// NB: explicit instantiation for all TAGGED_RPINTs after the template.
//
///////////////////////////////////////////////////////////////////////


// Step 1: Record what kmer precedes and follows each segment in the pathsdb.
template<class TAGGED_RPINT>
void BuildAdjacentKmers( const vecKmerPath& paths, const vecKmerPath& paths_rc,
			 const vec<TAGGED_RPINT>& pathsdb, 
			 vec<kmer_id_t>& next_read_kmer, vec<kmer_id_t>& prev_read_kmer,
			 const uint n_dots ) {
  
  // If a KmerPathInterval has no predecessor/successor, the value remains at -1.
  next_read_kmer.resize( pathsdb.size(), -1 );
  prev_read_kmer.resize( pathsdb.size(), -1 );
  longlong pid;
  int seg;
  ulonglong i = 0;
  
  // Dotting parameters.
  const size_t dot_n = pathsdb.size();
  size_t dot_i = 0;
  
  for( typename vec<TAGGED_RPINT>::const_iterator db_iter = pathsdb.begin();
       db_iter != pathsdb.end(); db_iter++) {

    // Find the immediate predecessor and successor for this KmerPathInterval.
    pid = db_iter->PathId();
    seg = db_iter->PathPos();
    const KmerPath& path = ( (pid >= 0) ? paths[pid] : paths_rc[-pid-1] );
    if ( seg+1 != path.NSegments() ) next_read_kmer[i] = path[seg+1].Start();
    if ( seg   != 0 )                prev_read_kmer[i] = path[seg-1].Stop ();
    i++;

    if (n_dots) dots_pct(dot_i++, dot_n);
  }
}


// uinterval: An interim structure to keep track of unipath intervals.
// A uinterval is basically a KmerPathInterval tagged with two kmer_ids
// indicating its predecessor and successor (or -1, if at the beginning/end of a
// read.)  There is a 1-to-1 correspondence between uintervals and
// KmerPathIntervals in the dataset.
struct uinterval {
  KmerPathInterval  kpi;
  kmer_id_t          prev_kmer;  // or -1 if beginning of unipath
  kmer_id_t          next_kmer;  // or -1 if end of unipath

  uinterval() {};
  uinterval( const kmer_id_t& start, const kmer_id_t& stop,
	     const kmer_id_t& prev_kmer, const kmer_id_t& next_kmer ) :
    kpi(start,stop), prev_kmer(prev_kmer), next_kmer(next_kmer) { }

  // Also, next_kmer=-2 means this uinterval has been used.
  void SetUsed() { next_kmer = -2; }
  bool WasUsed() { return (next_kmer == -2); }

  friend bool operator<(const uinterval& lhs, const uinterval& rhs)
  { return lhs.kpi.Start() < rhs.kpi.Start(); }

  friend ostream& operator<<(ostream& out, const uinterval u) {
    return out << u.kpi << " betw " << u.prev_kmer << " and " << u.next_kmer;
  }
};


// Step 2: Build unipath intervals
// Returns num_known_path_starts, a poor lower bound on the number of unipaths.
template<class TAGGED_RPINT>
int BuildUnipathIntervals( const vec<TAGGED_RPINT>& pathsdb, 
			   const vec<kmer_id_t>& next_read_kmer, 
			   const vec<kmer_id_t>& prev_read_kmer,
			   vec<uinterval>& uni_intervals,    /* Note idiotic interface!! */
			   const uint n_dots,
			   const String& uni_intervals_file = "",
			   bool verbose = false ) {

  // This is an idiotic interface, but I can't think of a better one:
  // If uni_intervals_file is nonempty, we write the intervals to disk, 
  // instead of pushing them back onto uni_intervals (which we leave untouched)
  Binary3Writer<uinterval> uni_writer;
  bool write_to_file = ! uni_intervals_file.empty();
  if( write_to_file )
    uni_writer.Open( uni_intervals_file );


  ForceAssert( pathsdb.nonempty() );
  kmer_id_t current_kmer = pathsdb.front().Start();
  // In the loop, we build the unipath segment starting with current_kmer
  typename vec<TAGGED_RPINT>::const_iterator db_iter, db_start = pathsdb.begin();
  // all db intervals before db_start end before current_kmer

  int num_known_path_starts = 0;
  // these variables are effectively local to each pass through the loop:
  bool begins_uni, ends_uni;  // true if we know this KPI must start/end unipath
  kmer_id_t ending_kmer;  // kmer ending this KPI
  kmer_id_t next_seg_kmer, prev_seg_kmer;  // kmer before/after this KPI in unipath
  kmer_id_t to_my_left, to_my_right, my_start, my_stop;  // kmers in *db_iter's read.
  kmer_id_t max_stop;  // largest Stop() we saw; for incrementing current_kmer

  // Dotting parameters.
  const size_t dot_n = pathsdb.size();
  size_t dot_i = 0;
  
  while( db_start != pathsdb.end() ) {  // actually the loop always exits via a break

    if (verbose && n_dots > 0) dots_pct(dot_i++, dot_n);

    if(verbose)
      cout << "Starting with kmer " << current_kmer << "...  " << endl;

    // We now construct the KmerPathInterval that appears in a unipath
    // and begins with current_kmer.

    begins_uni = ends_uni = false; // true if we know this KPI begins/ends unipath
    ending_kmer = next_seg_kmer = prev_seg_kmer = max_stop = -1;

    while( db_start != pathsdb.end() && db_start->Stop() < current_kmer )
      db_start++;

    max_stop = db_start->Stop();

    if(verbose)
      cout << "    unipath interval is now [" << current_kmer
	   << "-" << ending_kmer << "]" << endl;

    // step through the TAGGED_RPINTs of all intervals containing current_kmer
    for( db_iter = db_start; 
	 db_iter != pathsdb.end() && 
	   ( db_iter == db_start || db_iter->Start() <= ending_kmer );
	 db_iter++ ) {
      my_start = db_iter->Start();
      my_stop = db_iter->Stop();
      if( my_stop < current_kmer ) continue;
      max_stop = max( max_stop, my_stop );

      if(verbose)
	cout << "  Checking tagged_rpint [" << db_iter->Start()
	     << "-" << db_iter->Stop() << "], between "
	     << prev_read_kmer[ db_iter - pathsdb.begin() ]
	     << " and " 
	     << next_read_kmer[ db_iter - pathsdb.begin() ]
	     << ", my_stop=" << my_stop
	     << endl;

      // CHECK LEFT of current_kmer, if needed
      if( ! begins_uni && my_start <= current_kmer ) {
	// What kmer does this pathsdb entry think comes before current_kmer?
	if( current_kmer == my_start )
	  to_my_left = prev_read_kmer[ db_iter - pathsdb.begin() ];
	else
	  to_my_left = current_kmer - 1;

	// If we didn't know before, record the kmer before this unipath segment.
	// If we *did* know before, and the current read disagrees, this is
	// the beginning of a unipath
	if( prev_seg_kmer == -1 )
	  prev_seg_kmer = to_my_left;  // may still be -1, if at left end of read
	else if( to_my_left != -1 && to_my_left != prev_seg_kmer ) {
	  begins_uni = true;
	  prev_seg_kmer = -1;
	}
      }

      // CHECK LEFT of this interval, to see if the unipath must
      // end because of a join coming in from the left:
      if( my_start > current_kmer &&
	  prev_read_kmer[ db_iter - pathsdb.begin() ] != -1 ) {
	ending_kmer = my_start - 1;
	ends_uni = true;
	next_seg_kmer = -1;
      }

      // CHECK RIGHT.  This is a little trickier since we're both
      // figuring out where the unipath segment ends as well as
      // whether this is the last segment of the unipath or not.
      to_my_right = next_read_kmer[ db_iter - pathsdb.begin() ];  // maybe -1
      if( ending_kmer == -1 ) {  // first time through the loop only
	ending_kmer = my_stop;
	next_seg_kmer = to_my_right;
      }
      else if( ends_uni ) {  // we might need to shorten the interval
	if( my_stop < ending_kmer && to_my_right != -1 )
	  ending_kmer = my_stop;
      }
      else if( my_stop >= ending_kmer ) {
	if( next_seg_kmer == -1 ) {  // lengthen interval
	  ending_kmer = my_stop;
	  next_seg_kmer = to_my_right;
	}
	else if( my_stop > ending_kmer ||
		 (to_my_right != next_seg_kmer && to_my_right != -1) ) {
	  // declare this ends the unipath
	  ends_uni = true;
	  next_seg_kmer = -1;
	}
      }
      else if( my_stop < ending_kmer ) {  // maybe declare this ends the unipath
	if( to_my_right != -1 ) {
	  ending_kmer = my_stop;
	  ends_uni = true;
	  next_seg_kmer = -1;
	}
      }
      
      if(verbose) 
	cout << "    unipath interval is now [" << current_kmer
	     << "-" << ending_kmer 
	     << "], prev_seg_kmer=" << prev_seg_kmer
	     << ", next_seg_kmer=" << next_seg_kmer
	     << (begins_uni ? ", begins unipath" : "")
	     << (ends_uni ? ", ends unipath" : "")
	     << endl;

    }

    // Record the unipath interval I just created.
    uinterval u(current_kmer, ending_kmer, prev_seg_kmer, next_seg_kmer );
    if( write_to_file )
      uni_writer.Write( u );
    else
      uni_intervals.push_back( u );

    if( begins_uni ) num_known_path_starts++;

    if(verbose)
      cout << "  Found unipath interval " << uni_intervals.back().kpi 
	   << ", prev_seg_kmer=" << prev_seg_kmer
	   << ", next_seg_kmer=" << next_seg_kmer << endl;

    // Increment current_kmer for the next interval creation.
    if( max_stop > ending_kmer )
      current_kmer = ending_kmer+1;
    else if( db_iter == pathsdb.end() )
      break;  // actually the loop always exist this way
    else {  // Everything before db_iter is taken care of.
      current_kmer = db_iter->Start();
      db_start = db_iter;
    }
  }



  if( write_to_file )
    uni_writer.Close();

  return num_known_path_starts;

}




// Mark this UI as used, and return the next one, if unused.
// NOTE THAT THIS IS DESTRUCTIVE, since marking as used kills the next_kmer field.
// (But the prev_kmer field is untouched -- walk forwards first, then backwards.)
vec<uinterval>::iterator 
NextUI( vec<uinterval>::iterator iter, vec<uinterval>& uis ) {
  if( iter->next_kmer == -1 ) {
    iter->SetUsed();
    return uis.end();
  }
  vec<uinterval>::iterator ans =
    lower_bound( uis.begin(), uis.end(), uinterval( iter->next_kmer, 
						    iter->next_kmer,0,0 ) );
  iter->SetUsed();

  if( ! ans->WasUsed() && ans->prev_kmer == iter->kpi.Stop() )
    return ans;
  else
    return uis.end();
}


// Mark this UI as used, and return the previous one.
// NOTE THAT THIS IS DESTRUCTIVE, since marking as used kills the next_kmer field.
vec<uinterval>::iterator 
PrevUI( vec<uinterval>::iterator iter, vec<uinterval>& uis ) {
  if( iter->prev_kmer == -1 ) {
    iter->SetUsed();
    return uis.end();
  }
  vec<uinterval>::iterator ans =
    upper_bound( uis.begin(), uis.end(), uinterval( iter->prev_kmer, 
						    iter->prev_kmer,0,0 ) );
  if( ans == uis.end() || ans->kpi.Start() != iter->prev_kmer ) --ans;
  iter->SetUsed();

  if( ! ans->WasUsed() && ans->next_kmer == iter->kpi.Start() )
    return ans;
  else
    return uis.end();
}


// Step 3: string unipath segments together into unipaths.
// If the caller sorted uni_intervals and we used a bitvector to track
// which have been used so far, the uni_intervals argument could be const.
void BuildUnipathsFromIntervals( vec<uinterval>& uni_intervals,
				 vecKmerPath& unipaths,
				 bool verbose = false ) {

  // unipaths should already be reserved as well as possible by the caller.
  // That is, it should be prepared to hold uni_intervals KmerPathIntervals,
  // but we don't know how many inner vectors they will be in yet.
  sort( uni_intervals.begin(), uni_intervals.end() );

  KmerPath this_uni;
  vec<KmerPathInterval> kpis_forward, kpis_backward;
  int intervals_used = 0;
  int intervals_lost = 0;
  
  // Loop over all uintervals (i.e., all KmerPathIntervals) and follow the
  // links between uintervals to form KmerPaths from KmerPathIntervals.
  vec<uinterval>::iterator start_iter, random_iter;
  for( start_iter = uni_intervals.begin(); start_iter != uni_intervals.end(); 
       start_iter++ ) {
    if( start_iter->WasUsed() ) continue;
    
    kpis_forward.clear();
    kpis_backward.clear();

    if(verbose)
      cout << "Unipath " << unipaths.size()
	   << ": starting from unipath interval " << *start_iter << endl;

    // Keep track of the kmer_ids at the beginning and at the end of this
    // unipath.  This allows us to identify self-contained loops.
    kmer_id_t first_kmer_fw = start_iter->next_kmer;
    kmer_id_t  last_kmer_fw = start_iter->kpi.Stop( );
    kmer_id_t first_kmer_bw = start_iter->kpi.Start( );
    kmer_id_t  last_kmer_bw = start_iter->prev_kmer;
    
    // Step forward from the uinterval at start_iter.
    for( random_iter = start_iter;
	 random_iter != uni_intervals.end();
	 random_iter = NextUI( random_iter, uni_intervals ) ) {
      if(verbose)
	cout << "kpis_forward gets " << *random_iter << "\n";
      kpis_forward.push_back( random_iter->kpi );
      
      first_kmer_fw = random_iter->next_kmer;
      last_kmer_fw  = random_iter->kpi.Stop( );
    }

    // Step backward from the uinterval at start_iter.
    for( random_iter = PrevUI( start_iter, uni_intervals ); 
	 random_iter != uni_intervals.end();
	 random_iter = PrevUI( random_iter, uni_intervals ) ) {
      if(verbose)
	cout << "kpis_backward gets " << *random_iter << "\n";
      kpis_backward.push_back( random_iter->kpi );
      
      first_kmer_bw = random_iter->kpi.Start( );
      last_kmer_bw  = random_iter->prev_kmer;
    }

    this_uni.Clear(); 
    for(size_t i=kpis_backward.size(); i>0; i--)
      this_uni.AddSegment( kpis_backward[i-1] );
    for(size_t i=0; i < kpis_forward.size(); i++)
      this_uni.AddSegment( kpis_forward[i] );
    
    // If this unipath is a self-contained loop, we must make sure to rotate the
    // order of the KmerPathIntervals that comprise it, so that it is consistent
    // with its rc.
    bool is_loop = false;
    if ( first_kmer_fw == first_kmer_bw &&  last_kmer_fw == last_kmer_bw ) 
      is_loop = true;
        
    if ( is_loop ) {
      
      // Criterion: Examine all of the KmerPathIntervals and find the one with
      // the lowest-numbered kmer.
      // Note that we consider the rc'd KmerPathIntervals as well (so that this
      // unipath's result will be consistent with its rc's.)
      
      ForceAssertGt( this_uni.NSegments(), 0 );
      kmer_id_t min_kmer_fw_id = this_uni.Start(0);
      int kpi_start = 0;
      for ( int i = 1; i < this_uni.NSegments( ); i++ )
	if ( this_uni.Start( i ) < min_kmer_fw_id ){
	  min_kmer_fw_id = this_uni.Start( i );
	  kpi_start= i;
	}

      KmerPath this_uni_rc = this_uni;
      this_uni_rc.Reverse( );
      ForceAssertGt( this_uni_rc.NSegments(), 0 );
      kmer_id_t min_kmer_rc_id = this_uni_rc.Start(0);
      for ( int i = 1; i < this_uni_rc.NSegments( ); i++ )
	if ( this_uni_rc.Start( i ) < min_kmer_rc_id ) 
	  min_kmer_rc_id = this_uni_rc.Start( i );
	  

      if ( min_kmer_fw_id !=  min_kmer_rc_id ){
	
	if ( min_kmer_rc_id < min_kmer_fw_id ){
	  kmer_id_t last_kmer_id = flip_kmer( min_kmer_rc_id ); 
	  for ( int i = 0; i < this_uni.NSegments( ); i++ )
	    if ( this_uni.Stop( i ) == last_kmer_id )
	      kpi_start = ( i < this_uni.NSegments() -1 ? i + 1 : 0 );
	}
	// Rotate the chosen KmerPathInterval into the front of this unipath.
	KmerPath this_uni_new;
	for ( int i = kpi_start; i < this_uni.NSegments( ); i++ )
	  this_uni_new.AddSegment( this_uni.Segment( i ) );
	for ( int i = 0; i < kpi_start; i++ )
	  this_uni_new.AddSegment( this_uni.Segment( i ) );
      
	this_uni = this_uni_new;
	unipaths.push_back_reserve( this_uni ); 
	
      } else {
	// purely palindromic unipath. Break it into single kmer paths.
	for ( int i = 0; i < this_uni.NSegments(); i++ ){
	  for ( kmer_id_t k = this_uni.Start(i); k <= this_uni.Stop(i); k++ ){
	    KmerPathInterval kpi( k, k );
	    KmerPath this_uni_new;
	    this_uni_new.AddSegment( kpi );
	    unipaths.push_back_reserve( this_uni_new );
	  } 
	}

      }
    
    }else{ // not a loop
    
      unipaths.push_back_reserve( this_uni ); 

      if(verbose)
	cout << "Final path: " << this_uni << endl;
    }

  }

}  




// The old Unipath() couldn't handle a full mammalian-sized genome
// with small kmer size (eg human at k=16) -- it ended up in eternal
// swap, probably due to memory non-locality because of lots of random
// access to pathsdb.  This tries to do better: we pass linearly
// through pathsdb twice, with random access to the (much smaller)
// paths{_rc} the first time.
//
// The idea is to first build vec<kmer_id_t>s parallel to pathsdb which
// tell us what kmer comes right before/after each interval in the db.
// With this data in hand, we can do a second entirely local pass in 
// which we build all the intervals that will appear in the unipaths,
// along with enough extra data to go through them again and string
// them together into the full unipaths.


template<class TAGGED_RPINT>
void Unipath( const vecKmerPath& paths, const vecKmerPath& paths_rc,
	      const vec<TAGGED_RPINT>& pathsdb, vecKmerPath& unipaths, 
	      vec<TAGGED_RPINT>& unipathsdb, Bool log_progress, 
	      const String& unipaths_file, bool verbose ) {

  if ( paths.empty() ) {
    ForceAssert( paths_rc.empty() && pathsdb.empty() );
    unipaths.clear();
  } else {

    if (log_progress) 
      cout << Date( ) << ": starting unipath construction" << endl;




    if (log_progress)
      cout << "\n" << Date( )
	   << ": Step 1 (building next_read_kmer, prev_read_kmer)" << endl;

    // During step 1, we iterate through pathsdb and linearly create two large vectors.
    // We use random read access to paths and paths_rc; after step 1 they are not used.
    // Record, for each pathsdb entry, what kmer follows it, or -1 if end of read.
    vec<kmer_id_t> next_read_kmer, prev_read_kmer;
    BuildAdjacentKmers( paths, paths_rc, pathsdb, 
			next_read_kmer, prev_read_kmer, 
			log_progress ? 100 : 0 );

    // At this point we no longer need paths and paths_rc.
    // If we owned them, we'd destroy them to reclaim their memory.
  
  
    if (log_progress)
      cout << "\n" << Date( )
	   << ": Step 2 (building unipath intervals)" << endl;
    // Find all the kmer path intervals (KPIs) which will appear in the unipaths,
    // and for each, record what kmer follows it (or -1 if it ends the unipath).
    // During step 2, we iterate through pathsdb and linearly create the
    // vector of all the KmerPathIntervals in all unipaths.
    vec<uinterval> uni_intervals;  // a uinterval is a KPI plus some extra data
    int num_known_path_starts =
      BuildUnipathIntervals( pathsdb, next_read_kmer, prev_read_kmer,
			     uni_intervals, 
			     log_progress ? 100 : 0,
			     "", /* No disk temp needed */
			     verbose );

    // reclaim memory
    Destroy( next_read_kmer );
    Destroy( prev_read_kmer );
  
    // At this point we no longer need pathsdb.
    // If we owned it, we'd destroy it to reclaim memory.
  
  
    if (log_progress)
      cout << "\n" << Date( )
	   << ": Step 3 (stringing together unipaths)" << endl;
    // String together the uintervals into unipaths.  This involves
    // sorting the vec<uintervals> and then binary searching it lots,
    // but we no longer use paths{,_rc,db} at all.
    // There may be circular unipaths(!), but we ignore those for now.

    unipaths.clear();
    unipaths.Reserve( uni_intervals.size(), num_known_path_starts );
    // rawsize is right, but there are more unipaths than this.
    // So we'll have to push_back_reserve.
    BuildUnipathsFromIntervals( uni_intervals, unipaths, verbose );



    // Placed here for benchmarking vs old code.
    if (log_progress) cout << Date( ) << ": post-processing" << endl;
    //if (log_progress) PRINT2( unipaths.size(), uni_intervals.size() );

  
    // reclaim memory
    Destroy( uni_intervals );
  
  }  // !paths.empty()


  // If the string unipaths_file is nonempty, we're supposed to save to it.
  // This is vestigial (from the old Unipath(), which used it as a temp file too),
  // but there are callers which use it deliberately.
  if( unipaths_file != "" )
    unipaths.WriteAll(unipaths_file);

  CreateDatabase( unipaths, unipathsdb );

}















// Instantiate the Unipath(...) template for each TAGGED_RPINT type.

#define INSTANTIATE_UNIPATH(T)                                 \
template void Unipath( const vecKmerPath&, const vecKmerPath&, \
                       const vec<T>&, vecKmerPath&, vec<T>&,   \
                       Bool, const String&, bool );            \
template void BuildUnipathAdjacencyGraph( const vecKmerPath& paths, const vecKmerPath& paths_rc, const vec<T>& pathsdb, const vecKmerPath& unipaths, const vec<T>& unipathsdb, digraph& A );

INSTANTIATE_UNIPATH(tagged_rpint)
INSTANTIATE_UNIPATH(big_tagged_rpint)
INSTANTIATE_UNIPATH(new_tagged_rpint)
#undef INSTANTIATE_UNIPATH





///////////////////////////////////////////////////////////////////////
//
// Stuff that isn't part of Unipath(...) below here.
//
///////////////////////////////////////////////////////////////////////


void PrintInColumns( const vec<String>& s )
{    if ( s.empty( ) ) return;
     int colwidth = 0;
     for ( size_t i = 0; i < s.size( ); i++ )
          colwidth = Max( 2 + (int) s[i].size( ), colwidth );
     int ncols = 80 / colwidth;
     for ( size_t i = 0; i < s.size( ); i++ )
     {    if ( i > 0 && ( i % ncols ) == 0 ) cout << "\n";
          cout << s[i];
          for ( int j = 0; j < colwidth - (int) s[i].size( ); j++ )
               cout << " ";    }
     cout << "\n";    }

void DecomposePath( const KmerPath& x, Bool brief )
{    ForceAssertGt( x.KmerCount( ), 0 );
     vecKmerPath paths, paths_rc, unipaths;
     paths.push_back(x);
     vec<tagged_rpint> pathsdb, unipathsdb;
     x.AppendToDatabase( pathsdb, 0 );
     Prepare(pathsdb);
     Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
     longlong seg = 0, kmer = 0;
     while(1)
     {    kmer_id_t next = x.Segment(seg).Start( ) + kmer;
          vec<unipath_interval_id_t> places;
          Contains( unipathsdb, next, places );
          ForceAssertEq( places.size( ), 1u );
          const tagged_rpint& t = unipathsdb[ places[0] ];
          int lead = 0;
          int id = t.PathId( ), pp = t.PathPos( );
          if ( next > t.Start( ) || pp > 0 )
          {    ForceAssert( seg == 0 && kmer == 0 );
               lead = next - t.Start( );    
               for ( int j = 0; j < pp; j++ )
                    lead += unipaths[id].Segment(j).Length( );    }
          if ( seg > 0 || kmer > 0 ) cout << ".";
          cout << BaseAlpha(id);
          int nkmers = unipaths[id].KmerCount( ) - lead;
          if ( lead > 0 ) cout << "[" << nkmers << "]";
          lead = 0;
          for ( int j = 0; j < nkmers; j++ )
          {    if ( seg == x.NSegments( ) )
               {    cout << "[" << j << "]";
                    break;    }
               kmer++;
               if ( kmer == x.Segment(seg).Length( ) )
               {    kmer = 0;
                    ++seg;    }    }    
          if ( seg == x.NSegments( ) ) break;    }

     cout << "\nwhere:\n";
     if ( !brief )
     {    for ( size_t i = 0; i < unipaths.size( ); i++ )
          {    cout << BaseAlpha(i) << ": (" << unipaths[i].KmerCount( ) 
                    << " kmers)\n";
               cout << unipaths[i] << "\n";    }    }
     else
     {    vec<String> s;
          for ( size_t i = 0; i < unipaths.size( ); i++ )
          {    s.push_back( BaseAlpha(i) + "[" 
                    + ToString( unipaths[i].KmerCount( ) ) + "]" );    }
          PrintInColumns(s);    }    }

/**
   Function: BuildUnipathAdjacencyHyperKmerPath

   Create a HyperKmerPath where edges are unipaths, and two edges share a node
   if the two unipaths are adjacent (i.e. their end and start kmers are adjacent
   according to some read).

   *Warning*: This step loses information!  Algorithmically, it is better to work with
   the <unipath adjacency graph>.  The unipath adjacency HyperKmerPath may be easier
   to understand visually.
*/

void BuildUnipathAdjacencyHyperKmerPath(
     int K, const digraph& A, const vecKmerPath& unipaths, HyperKmerPath& h )
{    int nuni = unipaths.size( );
     vec<KmerPath> edges;
     edges.reserve(nuni);
     for ( int i = 0; i < nuni; i++ )
          edges.push_back( unipaths[i] );
     equiv_rel e( 2*nuni );
     for ( int v = 0; v < nuni; v++ )
     {    for ( size_t j = 0; j < A.From(v).size(); j++ )
          {    int w = A.From(v)[j];
               e.Join( 2*v + 1, 2*w );    }    }
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     Sort(reps);
     int N = reps.size( );
     vec< vec<int> > from(N), to(N);
     vec< vec<int> > from_edge_obj(N), to_edge_obj(N);
     vec<int> o, r;
     for ( int i = 0; i < N; i++ )
     {    e.Orbit( reps[i], o );
          for ( size_t j = 0; j < o.size(); j++ )
          {    int x = o[j];
               if ( x % 2 == 0 )
               {    int c = e.ClassId(x+1);
                    e.Orbit( c, r );
                    int v = -1;
                    for ( size_t m = 0; m < r.size(); m++ )
                    {    if ( e.Representative( r[m] ) )
                         {    v = BinPosition( reps, r[m] );
                              break;    }    }
                    from[i].push_back(v);
                    from_edge_obj[i].push_back( x/2 );    }
               else
               {    int c = e.ClassId(x-1);
                    e.Orbit( c, r );
                    int v = -1;
                    for ( size_t m = 0; m < r.size(); m++ )
                    {    if ( e.Representative( r[m] ) )
                         {    v = BinPosition( reps, r[m] );
                              break;    }    }
                    to[i].push_back(v);
                    to_edge_obj[i].push_back( (x-1)/2 );    }    }
          SortSync( from[i], from_edge_obj[i] );
          SortSync( to[i], to_edge_obj[i] );    }
     h.Initialize( K, from, to, edges, from_edge_obj, to_edge_obj );    }

void BuildUnibaseAdjacencyHyperBasevector(
     int K, const digraph& A, const vecbasevector& unibases, HyperBasevector& h )
{    int nuni = unibases.size( );
     vec<basevector> edges;
     edges.reserve(nuni);
     for ( int i = 0; i < nuni; i++ )
          edges.push_back( unibases[i] );
     equiv_rel e( 2*nuni );
     for ( int v = 0; v < nuni; v++ )
     {    for ( size_t j = 0; j < A.From(v).size( ); j++ )
          {    int w = A.From(v)[j];
               e.Join( 2*v + 1, 2*w );    }    }
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     Sort(reps);
     int N = reps.size( );
     vec< vec<int> > from(N), to(N);
     vec< vec<int> > from_edge_obj(N), to_edge_obj(N);
     vec<int> o, r;
     for ( int i = 0; i < N; i++ )
     {    e.Orbit( reps[i], o );
          for ( size_t j = 0; j < o.size(); j++ )
          {    int x = o[j];
               if ( x % 2 == 0 )
               {    int c = e.ClassId(x+1);
                    e.Orbit( c, r );
                    int v = -1;
                    for ( size_t m = 0; m < r.size(); m++ )
                    {    if ( e.Representative( r[m] ) )
                         {    v = BinPosition( reps, r[m] );
                              break;    }    }
                    from[i].push_back(v);
                    from_edge_obj[i].push_back( x/2 );    }
               else
               {    int c = e.ClassId(x-1);
                    e.Orbit( c, r );
                    int v = -1;
                    for ( size_t m = 0; m < r.size(); m++ )
                    {    if ( e.Representative( r[m] ) )
                         {    v = BinPosition( reps, r[m] );
                              break;    }    }
                    to[i].push_back(v);
                    to_edge_obj[i].push_back( (x-1)/2 );    }    }
          SortSync( from[i], from_edge_obj[i] );
          SortSync( to[i], to_edge_obj[i] );    }
     h.Initialize( K, from, to, edges, from_edge_obj, to_edge_obj );    }
