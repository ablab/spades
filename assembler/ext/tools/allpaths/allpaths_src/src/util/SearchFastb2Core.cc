///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "util/SearchFastb2Core.h"

#define CCC ToStringAddCommas

template<int K> void SFB( const vecbasevector& f1, const vecbasevector& f2,
     const int64_t start2, const int max_placements, 
     vec< triple<int,int,int> >& aligns, BitVec& too_many, const Bool verbose )
{    
     // Build kmer records from F2.

     int64_t nbases = 0;
     vec<int64_t> starts;
     int64_t start = 0;
     for ( size_t id = 0; id < f2.size( ); id++ )
     {    starts.push_back(start);
          if ( f2[id].isize( ) >= K )
          {    start += 2 * ( f2[id].isize( ) - K + 1 );
               nbases += 2 * ( f2[id].isize( ) - K + 1 );    }    }

     vec< kmer_record<K,2> > recs(nbases);
     #pragma omp parallel for 
     for ( longlong id = 0; id < (longlong)f2.size( ); id++ )
     {    const basevector& R2 = f2[id];
          kmer_record<K,2> rec;
          basevector bi(K), bir(K);
          for ( int j = 0; j <= R2.isize( ) - K; j++ )
          {    bi.SetToSubOf( R2, j, K );
               bir.ReverseComplement(bi); 
               rec.Set( bi, id, j );
               recs[ starts[id] + 2*j ] = rec;
               rec.Set( bir, id, -j-1 );
               recs[ starts[id] + 2*j + 1 ] = rec;    }    }

     // Sort kmer records from F2.

     if ( verbose )
     {    cout << Date( ) << ": sorting " << CCC( recs.size( ) ) 
               << " records, mem = " 
               << CCC( MemUsageBytes( ) ) << endl;    }
     ParallelSort( recs, kmer_record<K,2>::lt_basevector_id_pos);

     // Search kmer records from F2 for entries in F1.
     
     if ( verbose ) cout << Date( ) << ": searching records" << endl;
     vec< triple<int,int,int> > alignsi;
     basevector bi;
     kmer_record<K,2> rec;
     int64_t parallel_chunk_size
               = Max( (size_t) 1, f1.size( ) / ( 10 * omp_get_max_threads( ) ) );
     size_t naligns0 = aligns.size( );
     size_t total_aligns = aligns.size();
     #pragma omp parallel private(bi, rec, alignsi)
     {
          #pragma omp for schedule(dynamic, parallel_chunk_size)
       for ( longlong id1 = 0; id1 < (longlong)f1.size( ); id1++ )
          {    int aligns0 = alignsi.size( );
               const basevector& R1 = f1[id1];
               int n1 = R1.size( );
               if ( n1 < K ) continue;
               if ( n1 == K ) bi = R1;
               else bi.SetToSubOf( R1, 0, K );
               rec.Set( bi, 0, 0 );
               size_t low = lower_bound( recs.begin( ), recs.end( ), rec,
                    kmer_record<K,2>::cmp_basevector ) - recs.begin( );
               size_t high = upper_bound( recs.begin( ), recs.end( ), rec,
                    kmer_record<K,2>::cmp_basevector ) - recs.begin( );
               if ( low < high )
               {    for ( size_t j = low; j < high; j++ )
                    {    int id2 = recs[j].GetId( ), pos = recs[j].GetPos( );
                         Bool rc2 = ( pos < 0 );
                         if (rc2) pos = -pos - 1;
                         if ( n1 == K )
                         {    if ( alignsi.isize( ) - aligns0 == max_placements )
                              {    too_many.set(id1,true);
                                   alignsi.resize(aligns0);
                                   break;    }
                              if ( !rc2 ) alignsi.push( id1, start2 + id2, pos );
                              else alignsi.push( id1, start2 + id2, 
                                        -(pos-(n1-K)) - 1 );    }
                         else
                         {    const basevector& R2 = f2[id2];
                              Bool matches = True;
                              if ( !rc2 )
                              {    if ( R1.isize( ) + pos <= R2.isize( ) )
                                   {    for ( int u = K; u < n1; u++ )
                                        {    if ( R1[u] != R2[pos+u] )
                                             {    matches = False;
                                                  break;    }    }
                                   if (matches)
                                   {    if ( alignsi.isize( ) - aligns0
                                             == max_placements )
                                        {    too_many.set(id1,true);
                                             alignsi.resize(aligns0);
                                             break;    }
                                        alignsi.push( id1, 
                                             start2 + id2, pos );    }    }    }
                              else
                              {    if ( n1 - K <= pos )
                                   {    for ( int u = 0; u < n1 - K; u++ )
                                        {    if ( 3-R1[n1-u-1] != R2[pos-(n1-K)+u] )
                                             {    matches = False;
                                                  break;    }    }
                                        if (matches)
                                        {    if ( alignsi.isize( ) - aligns0
                                                  == max_placements )
                                             {    too_many.set(id1,true);
                                                  alignsi.resize(aligns0);
                                                  break;    }
                                             alignsi.push( id1, start2 + id2, 
                                                  -(pos-(n1-K)) - 1 );    }
                                             }    }    }    }    }    }
          #pragma omp critical
          { total_aligns += alignsi.size(); }
          
          #pragma omp barrier

          #pragma omp single
          {
	    if (verbose) cout << Date( ) << ": With records    mem = " << CCC( MemUsageBytes( ) ) << endl; 
	    Destroy(recs);
	    if (verbose) cout << Date( ) << ": Without records mem = " << CCC( MemUsageBytes( ) ) << endl;
	    if (total_aligns > aligns.capacity()) {
	      if (verbose) cout << Date( ) << ": Capacity =  " << aligns.capacity()  
		   << ", Need = " << total_aligns << endl;
	      aligns.reserve(total_aligns + (total_aligns - aligns.size()) * 1.25); 
	      if (verbose) cout << Date( ) << ": Capacity Now = " << aligns.capacity() << endl;
	    }
	  }

          #pragma omp critical
          { aligns.append(alignsi); }   
     }
     if (verbose) cout << Date( ) << ": After parallel  mem = " << CCC( MemUsageBytes( ) ) << endl; 
     if (verbose) cout << Date( ) << ": Capacity Now = " << aligns.capacity() << endl;
   
  // Announce results.

     if ( verbose )
     {    cout << Date( ) << ": found " << CCC( aligns.size( ) - naligns0 )
               << " aligns, mem = " << CCC( MemUsageBytes( ) ) << endl;    }    }

void SearchFastb2( const String& F1, const String& F2, const int K,
     vec< triple<int64_t,int64_t,int> >* pALIGNS, BitVec* pTooMany,
     const int MAX_PLACEMENTS, const double MEM_FRAC_TO_USE,
     const Bool verbose )
{
     // Check arguments.
  
     if ( K != 12 && K != 20 && K != 26 && K != 40 && K != 60 && K != 64 && K != 80 && K != 88 && K != 96 && K != 128 && K != 144 )
     {
       cout << "SearchFastb2: Not implemented for K=" << K << "." << endl;
       cout << "Abort." << endl;
       exit(1);
     }

     // Load F1 in batches.  Use at most 30% of available memory per batch and at 
     // most 2^32 - 1 elements per batch.

     size_t N1 = MastervecFileObjectCount(F1), X1 = FileSize(F1);
     if ( pTooMany )
         pTooMany->clear().reserve(N1);

     uint64_t bytes_available = MEM_FRAC_TO_USE * GetMaxMemory();
     const double U1 = 0.3;
     uint64_t BATCH_SIZE1 = Min( uint64_t(2147483647), uint64_t( floor( 
          ( U1 * double(bytes_available) / double(X1) ) * double(N1) ) ) );
     for ( size_t start1 = 0; start1 < N1; start1 += BATCH_SIZE1 )
     {    size_t stop1 = Min( start1 + BATCH_SIZE1, N1 );
          vecbasevector f1;
          if ( verbose )
          {    cout << Date( ) << ": loading seqs " << CCC(start1) << "-"
                    << CCC(stop1) << " of " << CCC(N1) 
                    << " from F1, mem = " << CCC( MemUsageBytes( ) ) << endl;    }
          f1.ReadRange( F1, start1, stop1 );
          
          // Search.

          uint64_t n2 = MastervecFileObjectCount(F2), start2 = 0, stop2 = 0;
          int ch = 0;
          uint64_t kmer_record_size = 0;

          if ( K == 12 )  kmer_record_size = sizeof( kmer_record<12,2> );
          if ( K == 20 )  kmer_record_size = sizeof( kmer_record<20,2> );
          if ( K == 26 )  kmer_record_size = sizeof( kmer_record<26,2> );
          if ( K == 40 )  kmer_record_size = sizeof( kmer_record<40,2> );
	  if ( K == 60 )  kmer_record_size = sizeof( kmer_record<60,2> );
          if ( K == 64 )  kmer_record_size = sizeof( kmer_record<64,2> );
          if ( K == 80 )  kmer_record_size = sizeof( kmer_record<80,2> );
          if ( K == 88 )  kmer_record_size = sizeof( kmer_record<88,2> );
          if ( K == 96 )  kmer_record_size = sizeof( kmer_record<96,2> );
	  if ( K == 128 ) kmer_record_size = sizeof( kmer_record<128,2> );
	  if ( K == 144 ) kmer_record_size = sizeof( kmer_record<144,2> );

          vec< triple<int,int,int> > aligns;
          BitVec too_many( f1.size(), false );
          while( start2 < n2 )
          {    
               // Define and load chunk.
     
               if ( verbose )
               {    cout << Date( ) << ": loading chunk " << ++ch << " from F2" 
                         << ", mem = " << CCC( MemUsageBytes( ) ) << endl;    }
               uint64_t nbases = 0;
               vecbasevector f2;
               while( stop2 < n2 )
               {    uint64_t chunk = Min( (uint64_t) 100000, n2 - stop2 );
                    f2.ReadRange( F2, stop2, stop2 + chunk );
                    uint64_t bytes_used = (int64_t) MemUsage( ) * (int64_t) 1000;
                    size_t id;
                    for ( id = 0; id < chunk; id++ )
                    {    if ( f2[id+stop2-start2].isize( ) >= K )
                              nbases += 2 * ( f2[id+stop2-start2].isize( ) - K + 1 );
                         if (bytes_used + (2*kmer_record_size*nbases) 
                              > bytes_available) 
                         {    break;    }
                         if ( f2.size( ) == 2147483647 ) break;    }
                    stop2 += id;
                    if ( id < chunk ) 
                    {    f2.resize( f2.size( ) - ( chunk - id ) );
                         break;    }    }

               // Search chunk.

               if ( stop2 == start2 )
               {    cout << "Not enough memory available." << endl;
                    cout << "Abort." << endl;
                    TracebackThisProcess( );    }
               if ( verbose )
               {    cout << Date( ) << ": using seqs " << CCC(start2) << "-" 
                         << CCC(stop2) << " of " << CCC(n2)
                         << " in F2, mem = " << CCC( MemUsageBytes( ) ) << endl;    }
               if (K == 12) SFB<12>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 20) SFB<20>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 26) SFB<26>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 40) SFB<40>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 60) SFB<60>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 64) SFB<64>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 80) SFB<80>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 88) SFB<88>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               if (K == 96) SFB<96>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
	       if ( K == 128 ) 
                    SFB<128>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
	       if ( K == 144 ) 
                    SFB<144>(f1, f2, start2, MAX_PLACEMENTS, aligns, too_many, verbose);
               start2 = stop2;    }

          // Sort and filter alignments.

          if ( verbose )
          {    cout << Date( ) << ": sorting " << CCC( aligns.size( ) ) 
                    << " aligns, mem = " << CCC( MemUsageBytes( ) ) << endl;    }
          ParallelSort(aligns);
          if ( MAX_PLACEMENTS > 0 )
          {    vec< triple<int,int,int> > aligns2;
               for ( size_t i = 0; i < aligns.size( ); i++ )
               {    size_t j;
                    for ( j = i + 1; j < aligns.size( ); j++ )
                         if ( aligns[j].first != aligns[i].first ) break;
                    if ( j - i <= (uint64_t) MAX_PLACEMENTS 
                         && !too_many[ aligns[i].first ] )
                    {    for ( size_t k = i; k < j; k++ )
                              aligns2.push_back( aligns[k] );    }
                    i = j - 1;    }
               aligns = aligns2;    
               if ( verbose )
               {    cout << Date( ) << ": after filtering, have " 
                         << CCC( aligns.size( ) )
                         << " aligns, mem = " 
                         << CCC( MemUsageBytes( ) ) << endl;    }    }

          if ( pALIGNS )
          {
              if ( !start1 )
              {
                  pALIGNS->clear();
                  size_t sz = aligns.size();
                  if ( stop1 != N1 )
                      sz = 1.1*N1/stop1*sz;
                  pALIGNS->reserve(sz);
              }
              for ( size_t i = 0; i < aligns.size( ); i++ )
                  pALIGNS->push( start1 + aligns[i].first,
                                  aligns[i].second, aligns[i].third );
          }
          if ( pTooMany )
              pTooMany->append(too_many);
     }
}
