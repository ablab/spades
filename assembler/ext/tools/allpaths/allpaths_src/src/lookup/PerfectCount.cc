///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "lookup/LookupTable.h"
#include "lookup/PerfectCount.h"
#include "random/RandomSampleFromStream.h"

// GetBestKmers: For each query (or its rc), find a kmer in it that appears a 
// minimal number of times in the target.

void GetBestKmers( const vecbasevector query, const AlignDir direction,
     const lookup_table& look, vec<int>& best_pos, vec<unsigned int>& best_index )
{    unsigned int K = look.K( );
     const int npasses = ( direction == FW ? 1 : 2 );
     const int nqueries = query.size( );
     best_pos.resize( nqueries * 2 );
     best_index.resize( nqueries * 2 );
     for ( int id = 0; id < nqueries; id++ )
     {    const basevector& s = query[id];
          if ( s.size( ) < K ) continue;
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
               best_index[ 2*id + pass ] = mindex;    }    }    }

int PerfectCount( const vecbasevector& query, const String& lookup_file, 
     const AlignDir direction )
{    vec<Bool> perfect;
     PerfectMark( query, lookup_file, direction, perfect );
     return Sum(perfect);    }

void PerfectMark( const vecbasevector& query, const String& lookup_file, 
     const AlignDir direction, vec<Bool>& perfect )
{    
     // Set up to track the perfects.

     perfect.resize_and_set( query.size( ), False );
     if ( query.empty( ) ) return;

     // Read header information from lookup table.

     lookup_table look(lookup_file);
     unsigned int K = look.K( );

     // For each query (or its rc), find a kmer in it that appears a minimal number 
     // of times in the target.

     vec<int> best_pos;
     vec<unsigned int> best_index;
     GetBestKmers( query, direction, look, best_pos, best_index );

     // Go through the lookup table chunks.

     const int npasses = ( direction == FW ? 1 : 2 );
     const int nqueries = query.size( );
     const unsigned int max_query_size = 200 * 1000 * 1000;
     for ( unsigned int i = 0; i < look.NChunks( ); i++ )
     {    look.ReadChunk(i);

          // Go through the query sequences.

          for ( int id = 0; id < nqueries; id++ )
          {    if ( perfect[id] ) continue;
               const basevector& s = query[id];
   	       unsigned int query_length = s.size();
               if ( query_length < K ) continue;

               // Go through the orientations.

               for ( int pass = 0; pass < npasses; pass++ )
    	       {    
                    if ( perfect[id] ) break;
                    
                    // r is low-freq kmer position in query

                    int r = best_pos[ 2*id + pass ];  
                    if ( r < 0 ) continue;
                    static basevector src;
                    if ( pass == 1 ) src.ReverseComplement(s);
                    const basevector& S = ( pass == 0 ? s : src );
                    unsigned int index = best_index[ 2*id + pass ];
                    unsigned int start = look.StartLocs(index);
                    unsigned int stop = look.StopLocs(index);
                    for ( unsigned int l = start; l < stop; l++ )
		    {    unsigned int offset = look.Locs(l) + ( max_query_size - r );

                         // rpos is low-freq kmer position in target

                         unsigned int tig, rpos; 
                         look.GetContigPos( look.Locs(l), tig, rpos );
			 
			 unsigned int startx;
			 unsigned int q_start;
			 unsigned int t_start;
			 unsigned int length = query_length;

			 // Determine starting position in target and query

			 if ( offset < look.ContigStart(tig) + max_query_size) 
                         {    q_start = r - rpos;
			      t_start = 0;
			      length -= q_start;
			      startx = look.ContigStart(tig);    } 
                         else 
                         {    q_start = 0;
			      t_start = rpos - r;
			      startx = offset - max_query_size;    }

			 if (startx + length > look.ContigStop(tig) ) {
			   length = look.ContigSize(tig) - t_start;
			 }

			 // Keep subsumed alignments only.

			 if ( length != query_length ) continue;

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

                         // Record perfectness.

                         perfect[id] = True;    
                         break;    }    }    }    }    }

void PerfectPick( const vecbasevector& query, const String& lookup_file,
		  const AlignDir direction, vec<placement_mark>& places,
                  vec< pair<placement_mark,int> >& places_plus,
                  const int mode, vec<Bool> & queryHasPerfect)
{    
     // Set up to track the perfects.

     if ( mode == 1 ) places.clear( );
     else places_plus.clear( );
     if ( query.empty( ) ) return;
     vec< RandomSampleFromStreamN1<placement_mark> > placesi( query.size( ) );
     queryHasPerfect.resize(query.size(), false);

     // Read header information from lookup table.

     lookup_table look(lookup_file);
     unsigned int K = look.K( );

     // For each query (or its rc), find a kmer in it that appears a minimal number 
     // of times in the target.

     vec<int> best_pos;
     vec<unsigned int> best_index;
     GetBestKmers( query, direction, look, best_pos, best_index );

     // Go through the lookup table chunks.

     const int npasses = ( direction == FW ? 1 : 2 );
     const int nqueries = query.size( );
     const unsigned int max_query_size = 200 * 1000 * 1000;
     for ( unsigned int i = 0; i < look.NChunks( ); i++ )
     {    look.ReadChunk(i);

          // Go through the query sequences.

          for ( int id = 0; id < nqueries; id++ )
          {    const basevector& s = query[id];
   	       unsigned int query_length = s.size();
               if ( query_length < K ) continue;

               // Go through the orientations.

               for ( int pass = 0; pass < npasses; pass++ )
    	       {    
                    // r is low-freq kmer position in query

                    int r = best_pos[ 2*id + pass ];  
                    if ( r < 0 ) continue;
                    static basevector src;
                    if ( pass == 1 ) src.ReverseComplement(s);
                    const basevector& S = ( pass == 0 ? s : src );
                    unsigned int index = best_index[ 2*id + pass ];
                    unsigned int start = look.StartLocs(index);
                    unsigned int stop = look.StopLocs(index);
                    for ( unsigned int l = start; l < stop; l++ )
		    {    unsigned int offset = look.Locs(l) + ( max_query_size - r );

                         // rpos is low-freq kmer position in target

                         unsigned int tig, rpos; 
                         look.GetContigPos( look.Locs(l), tig, rpos );
			 
			 unsigned int startx;
			 unsigned int q_start;
			 unsigned int t_start;
			 unsigned int length = query_length;

			 // Determine starting position in target and query

			 if ( offset < look.ContigStart(tig) + max_query_size) 
                         {    q_start = r - rpos;
			      t_start = 0;
			      length -= q_start;
			      startx = look.ContigStart(tig);    } 
                         else 
                         {    q_start = 0;
			      t_start = rpos - r;
			      startx = offset - max_query_size;    }

			 if (startx + length > look.ContigStop(tig) ) {
			   length = look.ContigSize(tig) - t_start;
			 }

			 // Keep subsumed alignments only.

			 if ( length != query_length ) continue;

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

                         // Record perfectness.
			 queryHasPerfect[id] = true;
                         static placement_mark P;
                         P.Set( tig, pass == 0, t_start );
                         placesi[id](P);    }    }    }    }

     // Prepare answer.

     int count = 0;
     for ( int id = 0; id < nqueries; id++ )
          if ( placesi[id].Exists( ) ) ++count;
     if ( mode == 1 ) places.reserve(count);
     else places_plus.reserve(count);
     for ( int id = 0; id < nqueries; id++ )
     {    if ( placesi[id].Exists( ) ) 
          {    if ( mode == 1 ) places.push_back( placesi[id].Sample( ) );
               else 
               {    places_plus.push_back( make_pair( 
                         placesi[id].Sample( ), id ) );    }    }    }    }

void PerfectPick( const vecbasevector& query, const String& lookup_file,
		  const AlignDir direction, vec<placement_mark>& places,
		  vec<Bool> & queryHasPerfect)
{    vec< pair<placement_mark,int> > places_plus;
     PerfectPick( query, lookup_file, direction, places, places_plus, 1,
          queryHasPerfect );    }

void PerfectPick( const vecbasevector& query, const String& lookup_file,
                  const AlignDir direction, 
                  vec< pair<placement_mark,int> >& places_plus,
                  vec<Bool> & queryHasPerfect )
{    vec<placement_mark> places;
     PerfectPick( query, lookup_file, direction, places, places_plus, 2,
          queryHasPerfect );    }

void PerfectCountPlaces( const vecbasevector& query, const String& lookup_file,
     const AlignDir direction, vec<int>& places )
{    
     // Set up to track the places.

     places.resize_and_set( query.size( ), 0 );
     if ( query.empty( ) ) return;

     // Read header information from lookup table.

     lookup_table look(lookup_file);
     unsigned int K = look.K( );

     // For each query (or its rc), find a kmer in it that appears a minimal number 
     // of times in the target.

     vec<int> best_pos;
     vec<unsigned int> best_index;
     GetBestKmers( query, direction, look, best_pos, best_index );

     // Go through the lookup table chunks.

     const int npasses = ( direction == FW ? 1 : 2 );
     const int nqueries = query.size( );
     const unsigned int max_query_size = 200 * 1000 * 1000;
     for ( unsigned int i = 0; i < look.NChunks( ); i++ )
     {    look.ReadChunk(i);

          // Go through the query sequences.

          for ( int id = 0; id < nqueries; id++ )
          {    const basevector& s = query[id];
   	       unsigned int query_length = s.size();
               if ( query_length < K ) continue;

               // Go through the orientations.

               for ( int pass = 0; pass < npasses; pass++ )
    	       {    
                    // r is low-freq kmer position in query

                    int r = best_pos[ 2*id + pass ];  
                    if ( r < 0 ) continue;
                    static basevector src;
                    if ( pass == 1 ) src.ReverseComplement(s);
                    const basevector& S = ( pass == 0 ? s : src );
                    unsigned int index = best_index[ 2*id + pass ];
                    unsigned int start = look.StartLocs(index);
                    unsigned int stop = look.StopLocs(index);
                    for ( unsigned int l = start; l < stop; l++ )
		    {    unsigned int offset = look.Locs(l) + ( max_query_size - r );

                         // rpos is low-freq kmer position in target

                         unsigned int tig, rpos; 
                         look.GetContigPos( look.Locs(l), tig, rpos );
			 
			 unsigned int startx;
			 unsigned int q_start;
			 unsigned int t_start;
			 unsigned int length = query_length;

			 // Determine starting position in target and query

			 if ( offset < look.ContigStart(tig) + max_query_size) 
                         {    q_start = r - rpos;
			      t_start = 0;
			      length -= q_start;
			      startx = look.ContigStart(tig);    } 
                         else 
                         {    q_start = 0;
			      t_start = rpos - r;
			      startx = offset - max_query_size;    }

			 if (startx + length > look.ContigStop(tig) ) {
			   length = look.ContigSize(tig) - t_start;
			 }

			 // Keep subsumed alignments only.

			 if ( length != query_length ) continue;

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

                         // Record placement.

                         ++places[id];    }    }    }    }    }
