///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTable.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "system/Worklist.h"

void PerfectLookup( const unsigned int K, const vecbasevector& query, 
		    const String& lookup_file, vec<look_align>& aligns, const AlignDir direction,
		    const Bool subsumed_only, const unsigned int target_seq_overlap)
{
  PerfectLookup( K, query, lookup_file, aligns, direction, 0, query.size() - 1,
		 subsumed_only, target_seq_overlap);
}

// TODO: potentially dangerous truncation of indexes by range_ args
void PerfectLookup( const unsigned int K, const vecbasevector& query, 
		    const String& lookup_file, vec<look_align>& aligns,
		    const AlignDir direction, int range_start, int range_end, 
		    const Bool subsumed_only, const unsigned int target_seq_overlap)
{    
     aligns.clear( );
     
     // Check we have something to align
     if (query.empty())
       return;

     // Read header information from lookup table.
     lookup_table look(lookup_file);

     ForceAssertEq( look.K( ), K );
     ForceAssertLt( static_cast<size_t>(range_end), query.size());
     ForceAssertLe( range_start, range_end);

     // For each query (or its rc), find a kmer in it that appears a minimal number 
     // of times in the target.

     const int npasses = ( direction == FW ? 1 : 2 );
     const int nqueries = range_end - range_start + 1;

     vec<int> best_pos( nqueries * 2 );
     vec<unsigned int> best_index( nqueries * 2 );
     for ( int id = 0; id < nqueries; id++ )
     {    int query_id = range_start + id;
          const basevector& s = query[query_id];
          if ( s.size() < K ) {
            for ( int pass = 0; pass < npasses; pass++ ) {
              best_pos[2*id+pass] = -1;
              best_index[2*id+pass] = 0;
            }
            continue;
          }
          for ( int pass = 0; pass < npasses; pass++ )
          {    static basevector src;
	       if ( pass == 1 ) src.ReverseComplement(s);
	       const basevector& S = ( pass == 0 ? s : src );
               int pos = -1;
               unsigned int min_freq = 0, mindex = 0;
	       int length = S.isize() - (int) K;
	       unsigned int index = Index( S, 0, K );
	       unsigned int freq = look.Freq(index);
	       if ( freq != 0 ) 
	       {    mindex = index;
		    pos = 0;
		    if (freq != 1)
		    {    min_freq = freq;
		         for ( int j = 1; j <= length; j++ )
			 {   NextIndex( index, S, j, K );
			     freq = look.Freq(index);
			     if ( freq == 1 )
			     {    mindex = index;
		                  pos = j; 
				  break;    }
			     if ( freq == 0 )
			     {    pos = -1;
                                  break;    }
			     if ( pos < 0 || freq < min_freq )
			     {    min_freq = freq;
			          mindex = index;
				  pos = j;    }    }    }    }
               best_pos[ 2*id + pass ] = pos;
               best_index[ 2*id + pass ] = mindex;    }    }

     // Set up for alignment generation.

     static look_align la;
     la.nhits = la.mutations = la.indels = 0;
     la.a.SetNblocks(1);
     la.a.SetGap( 0, 0 );

     // Go through the lookup table chunks.

     const unsigned int max_query_size = 200 * 1000 * 1000;
     for ( unsigned int i = 0; i < look.NChunks( ); i++ )
     {    look.ReadChunk(i);
          // Go through the query sequences.

          for ( int id = 0; id < nqueries; id++ )
          {    int query_id = range_start + id;
               const basevector& s = query[query_id];
   	       unsigned int query_length = s.size();

               // Go through the orientations.

               for ( int pass = 0; pass < npasses; pass++ )
		 {  int r = best_pos[ 2*id + pass ];  // r is low-freq kmer position in query
                    if ( r < 0 ) continue;
                    static basevector src;
                    if ( pass == 1 ) src.ReverseComplement(s);
                    const basevector& S = ( pass == 0 ? s : src );
                    unsigned int index = best_index[ 2*id + pass ];
                    unsigned int start = look.StartLocs(index);
                    unsigned int stop = look.StopLocs(index);
                    for ( unsigned int l = start; l < stop; l++ )
		    {    unsigned int offset = look.Locs(l) + ( max_query_size - r );

                         unsigned int tig, rpos; // rpos is low-freq kmer position in target
                         look.GetContigPos( look.Locs(l), tig, rpos );
			 
			 unsigned int startx;
			 unsigned int q_start;
			 unsigned int t_start;
			 unsigned int length = query_length;

			 // Determine starting position in target and query
			 if ( offset < look.ContigStart(tig) + max_query_size) {
			   q_start = r - rpos;
			   t_start = 0;
			   length -= q_start;
			   startx = look.ContigStart(tig);
			 } else {
			   q_start = 0;
			   t_start = rpos - r;
			   startx = offset - max_query_size;
			 }

			 if (startx + length > look.ContigStop(tig) ) {
			   length = look.ContigSize(tig) - t_start;
			 }

			 // Do we want subsumed alignments only?
			 if (length != query_length && subsumed_only)
			   continue;

                         // Validate alignment, skipping portion covered by kmer.
                         if ( startx < look.BasesStart( ) ) continue;
                         Bool mismatch = False;
			 unsigned int start_of_kmer = r - q_start;
                         for ( unsigned int y = 0; y < start_of_kmer; y++ )
			 {   if ( S[q_start + y] != look.Base( startx + y ) )
			     {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;  // Alignment failed
                         for ( unsigned int y = start_of_kmer + K; y < length; y++ )
			 {   if ( S[q_start + y] != look.Base( startx + y ) )
			     {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;  // Alignment failed

			 // For partial alignments and overlapping target sequences we try to
			 // filter out all but one of the overlapping alignments.
			 if (length != query_length && target_seq_overlap != 0) {
			   if (length <= target_seq_overlap ) // A longer alignment must also 
			     continue;                        // exist so ignore this one.
			   if (q_start != 0 && rpos - t_start <= target_seq_overlap - K) 
			     continue;  // Another overlapping alignment must exist.
			 }

                         // Record alignment.
                         la.query_id = query_id;
                         la.target_id = tig;
                         la.a.Setpos1(q_start);
                         la.a.Setpos2(t_start);
                         la.query_length = query_length;
                         la.target_length = look.ContigSize(tig);
                         la.rc1 = ( pass == 1 );
                         la.a.SetLength( 0, length );
                         aligns.push_back(la);    }    }    }    }   }

// ---------------------------------------------------------------------------------

void ParallelPerfectLookup( const int max_threads, const unsigned int K, 
     const vecbasevector& query, const String& lookup_file, vec<look_align>& aligns,
     const AlignDir direction, const Bool subsumed_only, 
     const unsigned int target_seq_overlap )
{    ParallelPerfectLookup( max_threads, K, query, lookup_file, aligns, direction, 
          0, query.size() - 1, subsumed_only, target_seq_overlap );    }

// Communication for parallel threads:

vec<int> best_pos;
vec<unsigned int> best_index;
const vecbasevector* queryp;
int Range_start;
int Kx;
int npasses;
vec<int> start_id, stop_id;
unsigned int max_query_size = 200 * 1000 * 1000;
Bool Subsumed_only;
unsigned int Target_seq_overlap;
lookup_table* lookp;
vec<look_align>* alignsp;

class PerfectLookupProcessor1 {

     public:

     void operator( )( int i )
     {    basevector src;
          for ( int id = start_id[i]; id < stop_id[i]; id++ )
          {    int query_id = Range_start + id;
               const basevector& s = (*queryp)[query_id];
               if ( s.isize() < Kx ) 
               {    for ( int pass = 0; pass < npasses; pass++ ) 
                    {    best_pos[2*id+pass] = -1;
                         best_index[2*id+pass] = 0;    }
                    continue;    }
               for ( int pass = 0; pass < npasses; pass++ )
               {    if ( pass == 1 ) src.ReverseComplement(s);
	            const basevector& S = ( pass == 0 ? s : src );
                    int pos = -1;
                    unsigned int min_freq = 0, mindex = 0;
	            int length = S.isize( ) - Kx;
	            unsigned int index = Index( S, 0, Kx );
	            unsigned int freq = lookp->Freq(index);
	            if ( freq != 0 ) 
	            {    mindex = index;
		         pos = 0;
		         if (freq != 1)
		         {    min_freq = freq;
		              for ( int j = 1; j <= length; j++ )
			      {   NextIndex( index, S, j, Kx );
			          freq = lookp->Freq(index);
			          if ( freq == 1 )
			          {    mindex = index;
		                       pos = j; 
				       break;    }
			          if ( freq == 0 )
			          {    pos = -1;
                                       break;    }
			          if ( pos < 0 || freq < min_freq )
			          {    min_freq = freq;
			               mindex = index;
				       pos = j;    }    }    }    }
                    best_pos[ 2*id + pass ] = pos;
                    best_index[ 2*id + pass ] = mindex;    }    }    }
};

class PerfectLookupProcessor2 {
    static pthread_mutex_t LOCK_FL;

     public:

     void operator( )( int i )
     {    
          // Set up for alignment generation.

          look_align la;
          la.nhits = la.mutations = la.indels = 0;
          la.a.SetNblocks(1);
          la.a.SetGap( 0, 0 );

          // Go through the query sequences.

          basevector src;
          for ( int id = start_id[i]; id < stop_id[i]; id++ )
          {    int query_id = Range_start + id;
               const basevector& s = (*queryp)[query_id];
   	       unsigned int query_length = s.size();

               // Go through the orientations.

               for ( int pass = 0; pass < npasses; pass++ )
 	       {    // r is low-freq kmer position in query
                    int r = best_pos[ 2*id + pass ];  
                    if ( r < 0 ) continue;
                    if ( pass == 1 ) src.ReverseComplement(s);
                    const basevector& S = ( pass == 0 ? s : src );
                    unsigned int index = best_index[ 2*id + pass ];
                    unsigned int start = lookp->StartLocs(index);
                    unsigned int stop = lookp->StopLocs(index);
                    for ( unsigned int l = start; l < stop; l++ )
		    {    unsigned int offset = lookp->Locs(l) 
                              + ( max_query_size - r );
                         // rpos is low-freq kmer position in target
                         unsigned int tig, rpos; 
                         lookp->GetContigPos( lookp->Locs(l), tig, rpos );
			 unsigned int startx;
			 unsigned int q_start;
			 unsigned int t_start;
			 unsigned int length = query_length;

			 // Determine starting position in target and query

			 if ( offset < lookp->ContigStart(tig) + max_query_size) 
                         {    q_start = r - rpos;
			      t_start = 0;
			      length -= q_start;
			      startx = lookp->ContigStart(tig);   } 
                         else 
                         {    q_start = 0;
			      t_start = rpos - r;
			      startx = offset - max_query_size;    }
			 if (startx + length > lookp->ContigStop(tig) )
			      length = lookp->ContigSize(tig) - t_start;

			 // Do we want subsumed alignments only?

			 if ( length != query_length && Subsumed_only ) continue;

                         // Validate alignment, skipping portion covered by kmer.

                         if ( startx < lookp->BasesStart( ) ) continue;
                         Bool mismatch = False;
			 unsigned int start_of_kmer = r - q_start;
                         for ( unsigned int y = 0; y < start_of_kmer; y++ )
			 {    if ( S[q_start + y] != lookp->Base( startx + y ) )
			      {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;  // Alignment failed
                         for ( unsigned int y = start_of_kmer + Kx; y < length; y++ )
			 {   if ( S[q_start + y] != lookp->Base( startx + y ) )
			     {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;  // Alignment failed

			 // For partial alignments and overlapping target sequences 
                         // we try to filter out all but one of the overlapping 
                         // alignments.

			 if (length != query_length && Target_seq_overlap != 0) 
                         {    // A longer alignment must also 
                              // exist so ignore this one.
                              if (length <= Target_seq_overlap ) continue;
                              // Another overlapping alignment must exist.
			      if (q_start != 0 && rpos - t_start <= 
                                   Target_seq_overlap - Kx) 
			      {    continue;    }    }

                         // Record alignment.

                         la.query_id = query_id;
                         la.target_id = tig;
                         la.a.Setpos1(q_start);
                         la.a.Setpos2(t_start);
                         la.query_length = query_length;
                         la.target_length = lookp->ContigSize(tig);
                         la.rc1 = ( pass == 1 );
                         la.a.SetLength( 0, length );
                         pthread_mutex_lock( &LOCK_FL );
                         alignsp->push_back(la);    
                         pthread_mutex_unlock( &LOCK_FL );    }    }    }    }

};

pthread_mutex_t PerfectLookupProcessor2::LOCK_FL = PTHREAD_MUTEX_INITIALIZER;

void ParallelPerfectLookup( const int max_threads, const unsigned int K, 
     const vecbasevector& query, const String& lookup_file, vec<look_align>& aligns,
     const AlignDir direction, int range_start, int range_end, 
     const Bool subsumed_only, const unsigned int target_seq_overlap )
{    
     aligns.clear( );

     // Define some global variables for parallelization.

     queryp = &query;
     Range_start = range_start;
     Subsumed_only = subsumed_only;
     Target_seq_overlap = target_seq_overlap;
     Kx = K;
     npasses = ( direction == FW ? 1 : 2 );
     alignsp = &aligns;
     
     // Check we have something to align

     if (query.empty()) return;

     // Read header information from lookup table.

     lookup_table look(lookup_file);
     lookp = &look;
     ForceAssertEq( lookp->K( ), K );
     ForceAssertLt( static_cast<size_t>(range_end), query.size());
     ForceAssertLe( range_start, range_end);

     // Define read ranges for parallelization.

     const int nqueries = range_end - range_start + 1;
     int batch_size = (nqueries+9)/10, count = 0;
     while(1)
     {    start_id.push_back(count);
          count = Min( count + batch_size, nqueries );
          stop_id.push_back(count);
          if ( count == nqueries ) break;    }

     // For each query (or its rc), find a kmer in it that appears a minimal number 
     // of times in the target.

     best_pos.resize( nqueries * 2 ), best_index.resize( nqueries * 2 );
     {    PerfectLookupProcessor1 p;
          Worklist<int,PerfectLookupProcessor1> worklist( p, max_threads );
          for ( int i = 0; i < start_id.isize( ); i++ )
               worklist.add(i);    }

     // Go through the lookup table chunks.

     for ( unsigned int i = 0; i < lookp->NChunks( ); i++ )
     {    lookp->ReadChunk(i);

          // Go through the query sequences.

          PerfectLookupProcessor2 p;
          Worklist<int,PerfectLookupProcessor2> worklist( p, max_threads );
          for ( int i = 0; i < start_id.isize( ); i++ )
               worklist.add(i);    }

     // Sort result, thus insuring that it is well defined.

     Sort( aligns, order_lookalign_QueryIdFullTargetMutations( ) );    }
