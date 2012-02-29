/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <numeric>

#include "lookup/LookAlign.h"
#include "lookup/FirstLookup.h"
#include "system/Worklist.h" // Worklist (for parallelization)

class FirstLookupWorker {

private:
  static pthread_mutex_t LOCK_FL;

  void processChunk(lookup_table& look, const int chunk) {

    cout << Date() << ": Processing Chunk: " << chunk + 1 << " of " << look.NChunks() << endl;
    
    unsigned int K = look.K( );

    const unsigned int bad_kmer = numeric_limits<unsigned int>::max();
  
    look_align la;
    la.a.SetNblocks(1);
    la.a.SetGap( 0, 0 );

    vec<unsigned int> offsets( 0 );
    look.ReadChunk(chunk);
    
    unsigned int aligns_per_chunk = 0;
    // Go through the query sequences.
    
    for ( size_t id = 0; id < query->size( ); id++) {
      const basevector& s = (*query)[id];
      const unsigned int query_length = s.size();
      const unsigned int min_align_length = (filter.min_size == 0 ? query_length : filter.min_size);

      // Keep track of the best alignment on this query that we have seen.
      // There may be more than one, if their alignment scores are tied.
      vec<look_align> best_aligns( 0 );
      const int max_mismatches = 4;

      // Set batch match to minimum number of matching bases we will accept as an alignment
      int best_match = (filter.min_match ? (query_length - max_mismatches) / 2 : 0);

      // Go through the orientations.  Apply orientation filtering here.
      
      for ( int pass = 0; pass < 2; pass++ ) {
	if ( pass == 0 && filter.orientation == FirstLookupFilter::RC_ONLY ) continue;
	if ( pass == 1 && filter.orientation == FirstLookupFilter::FW_ONLY ) continue;
	
	size_t n = 2*id + pass;

	// If this kmer never appears (freq = 0) or appears too often
	// (freq > FirstLookupFilter::max_kmer_freq), ignore it.
	if ( (*all_index)[n] == bad_kmer ) continue;
	
	basevector src = s;
	if ( pass == 1 ) src.ReverseComplement(s);
	const basevector& S = ( pass == 0 ? s : src );
	vec<unsigned char> SX( S.size( ) );
	for ( int i = 0; i < S.isize( ); i++ )
	  SX[i] = S[i];
	
	// Find the set of alignments for this query.
	unsigned int start = look.StartLocs( (*all_index)[n] );
	unsigned int stop = look.StopLocs( (*all_index)[n] );
	offsets.clear( );
	for ( unsigned int l = start; l < stop; l++ ) {
	  unsigned int offset = look.Locs( l );
	  if ( pass == 1) 
	    if (query_length - K > offset)
	      continue;
	    else
	      offset -= ( query_length - K );

	  unsigned int a_start = offset; // this is where the alignment starts on the reference
	  unsigned int tig, rpos; 
	  look.GetContigPos( a_start, tig, rpos );
	  
	  // Ignore alignments where the query sequence hangs off the
	  // beginning of the target sequence.  (Note that we ALLOW alignments
	  // where the query hangs off the END of the target sequence, though
	  // the min_size filter affects these cases.)
	  // Also note: Right now, we cannot find alignments on contigs which
	  // bridge lookup_table chunks.  So we ignore those contigs.
	  if ( a_start < look.BasesStart( ) ) continue;
	  if ( a_start < look.ContigStart( tig ) ) continue;
	  
	  // Find the alignment length, which is equal to the query length,
	  // except it may be truncated to the end of the target sequence.
	  const unsigned int a_length = Min( query_length, look.ContigStop( tig ) - a_start );
	  
	  // Make sure there are at least enough bases on the target sequence
	  // to satisfy the min_size filter or contain the entire read.
	  if ( a_length < min_align_length ) continue;

	  // Stop if we have too many alignments to extend and check
	  if ( ( filter.max_extend != 0) && ( (offsets.size() + 1) > filter.max_extend) ) break;

	  offsets.push_back( (unsigned int) offset );
	}

	if ( offsets.empty( ) ) continue;

	  // Stop if we have too many alignments to extend and check
	if ( ( filter.max_extend != 0) && (offsets.size()  > filter.max_extend) ) continue;
		
	UniqueSort(offsets);
	
	// Loop over all alignments for this query.  Filter them and score
	// them, and keep the alignment with the highest score.
	for ( int u = 0; u < offsets.isize( ); u++ ) {
	  
	  unsigned int a_start = offsets[u]; // this is where the alignment starts on the reference
	  unsigned int tig, rpos; 
	  look.GetContigPos( a_start, tig, rpos );
	  	  
	  // Find the alignment length, which is equal to the query length,
	  // except it may be truncated to the end of the target sequence.
	  const unsigned int a_length = Min( query_length,
					     look.ContigStop( tig ) - a_start,
					     look.BasesStop( ) - a_start );
	  
	  // Score the alignment.
	  // Step through the query/target, starting at the matching K-mer,
	  // and stopping when you've found <max_mismatches> mismatched bases.
	  // The score of the alignment is equal to the number of correct base
	  // matches found.
	  int matches = K, mismatches = 0, match_length = 0;
	  if ( pass == 0 ) {
	    for ( unsigned int y = K; y < a_length; y++ ) {
	      if ( SX[y] == look.Base( a_start + y ) )
		++matches;
	      else if ( ++mismatches == max_mismatches ) {
		break;
	      }
	    }
	  }
	  else {
	    for ( int y = a_length-K-1; y >= 0; y-- ) {
	      if ( SX[y] == look.Base( a_start + y ) )
		++matches;
	      else if ( ++mismatches == max_mismatches ) {
		break;
	      }
	    }
	  }

	  if ( matches < best_match ) continue;
	  if ( matches > best_match ) best_aligns.clear( );
 
	  int mutations = 0;
	  for ( unsigned int y = 0; y < a_length; y++ ) {
	    if ( SX[y] != look.Base( a_start + y ) )
	      ++mutations;
	  }
	  
	  // Create alignment.
	  
	  la.query_id = id;
	  la.target_id = tig;
	  la.a.Setpos1( 0 ); // the alignment always starts at the beginning of the query!
	  la.a.Setpos2( rpos );
	  la.query_length = query_length;
	  la.target_length = look.ContigSize(tig);
	  la.rc1 = ( pass == 1 );
	  la.a.SetLength( 0, a_length );
	  la.mutations = mutations;
	  
	  // Update.
	  
	  best_aligns.push_back(la);
	  best_match = matches;
	}
      }
      
      // Record the best alignment(s) from this query.
      UniqueSort( best_aligns );


      pthread_mutex_lock( &LOCK_FL );

      // Add aligns to global aligns
      aligns->append( best_aligns );
 
      pthread_mutex_unlock( &LOCK_FL );

     aligns_per_chunk += best_aligns.size();
    }
    cout << Date() << ": Processed Chunk : " << chunk + 1 << " - Found " <<  aligns_per_chunk << " alignments" << endl;
  }


public:

  FirstLookupWorker() {};
  FirstLookupWorker(const FirstLookupWorker &) {};

  void operator() (unsigned int chunk) {

    pthread_mutex_lock( &LOCK_FL );
    
    // Find free table
    int free_table = Position(table_in_use, false);
    
    // There should always be a free table available
    ForceAssertGt(free_table, -1);
    
    // Claim this table
    table_in_use[free_table] = true;
    
    pthread_mutex_unlock( &LOCK_FL );
    
    // Process this chunk of the lookup table
    processChunk(*look_pool[free_table], chunk);
    
    pthread_mutex_lock( &LOCK_FL );
    
    // Release this table for next thread to use
    table_in_use[free_table] = false;
    
    pthread_mutex_unlock( &LOCK_FL );
  }

  static const vecbasevector * query;
  static vec<unsigned int> * all_index;
  static FirstLookupFilter filter;
  static vec<look_align> * aligns;
  static vec<lookup_table*>look_pool;
  static vec<bool> table_in_use;
};


const vecbasevector * FirstLookupWorker::query;
vec<unsigned int> * FirstLookupWorker::all_index;
FirstLookupFilter FirstLookupWorker::filter;
vec<look_align> * FirstLookupWorker::aligns;
vec<lookup_table*> FirstLookupWorker::look_pool;
vec<bool> FirstLookupWorker::table_in_use;
pthread_mutex_t FirstLookupWorker::LOCK_FL = PTHREAD_MUTEX_INITIALIZER;


/**
 * FirstLookup
 */
void FirstLookup( const vecbasevector& query,
		  const String& lookup_file, 
		  vec<look_align>& aligns,
		  const FirstLookupFilter & filter,
		  const int NUM_THREADS)
{
  aligns.clear( );
  
  // Read header information from lookup table.
  
  lookup_table look(lookup_file);
  unsigned int K = look.K( );

  // For each query (or its rc), find the indices of each of its kmers.

  const int nqueries = query.size( );
  aligns.reserve( nqueries );

  // use max unsigned int value as flag for bad kmers
  const unsigned int bad_kmer = numeric_limits<unsigned int>::max();

  vec<unsigned int> all_index( nqueries * 2, bad_kmer );
  for ( int id = 0; id < nqueries; id++ ) {
    const basevector& s = query[id];
    int length = s.isize( ) - (int) K;
    
    // Forward
    unsigned int index = Index( s, 0, K );
    unsigned int freq = look.Freq(index);
    bool bad_freq = ( freq == 0 );
    if ( ( filter.max_kmer_freq != 0 ) && ( freq > filter.max_kmer_freq ) ) bad_freq = true;
    if ( !bad_freq ) all_index[ 2*id ] = index; // if bad_freq, all_index remains at bad_kmer 
    // RC
    basevector src;
    src.ReverseComplement( s );
    index = Index( src, length, K );
    freq = look.Freq(index);
    bad_freq = ( freq == 0 );
    if ( ( filter.max_kmer_freq != 0 ) && ( freq > filter.max_kmer_freq ) ) bad_freq = true;
    if ( !bad_freq ) all_index[ 2*id + 1 ] = index; // if bad_freq, all_index remains at bad_kmer
  }

  // Setup common thread data
  FirstLookupWorker::query = &query;
  FirstLookupWorker::all_index =  &all_index;
  FirstLookupWorker::filter = filter;
  FirstLookupWorker::aligns = &aligns;

  // Create pool of lookup tables (one for each running thread)
  FirstLookupWorker::look_pool.resize(NUM_THREADS);
  FirstLookupWorker::table_in_use.resize(NUM_THREADS, false);
  // Use existing lookup table in the pool as entry [0]
  FirstLookupWorker::look_pool[0] = &look;
  // Create copies of existing table to fill up the pool
  for (int i = 1; i < NUM_THREADS; i++) {
    FirstLookupWorker::look_pool[i] = new lookup_table(look, lookup_file);
  }
      
  cout << Date( ) << ": Begin processing " << look.NChunks() << " Chunks "
       << (NUM_THREADS > 1 ? "(parallelized)" : "") << endl;
  FirstLookupWorker worker;
  Worklist<int,FirstLookupWorker> * worklist = new Worklist<int,FirstLookupWorker>( worker, NUM_THREADS - 1 );
  
  for ( int i = 0; i < (int)look.NChunks(); i++ )
    worklist->add( i );
  
  // The Worklist destructor doesn't return until it has done all its processing
  delete worklist;

  // Destroy lookup table pool (skipping first special entry)
  for (int i = 1; i < NUM_THREADS; i++)
    delete FirstLookupWorker::look_pool[i];


}

