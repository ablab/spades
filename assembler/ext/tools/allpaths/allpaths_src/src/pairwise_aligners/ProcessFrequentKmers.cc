// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


// Under development!!!  Just a plaything for now.

// ProcessFrequentKmers: Ultimate Purpose -- save information regarding large 
// families of reads, all of which overlap each other.  This information may be
// used later to reconstruct highly repetitive regions of the genome.

// Experimental data set:
// ReadsToAligns1 DATA=projects/016 RUN=run_dir MAXCLIQ=50 MAXBAD=100 
// PROCESS_FREQUENT_KMERS=True | fold -s
// Or use L3191.

#include <iomanip>

#include "Basevector.h"
#include "math/Functions.h"
#include "kmers/KmerRecord.h"
#include "pairwise_aligners/ProcessFrequentKmers.h"
#include "String.h"
#include "system/Types.h"
#include "Vec.h"

class clique_record {
     public:
     int clique_id;
     int read_id;
     clique_record( int clique_id_arg, int read_id_arg ) 
       : clique_id(clique_id_arg), read_id(read_id_arg)
     { }

     friend Bool operator<( const clique_record& cr1, const clique_record& cr2 )
     {    return cr1.read_id < cr2.read_id;    }
};

template<int I, int K>
void ProcessFrequentKmers( vec<bvec const*> const& EE,
                           vec< kmer_record<K,I> > R,
                           unsigned int S,
                           int max_clique,
                           String run_dir )
{
     cout << "There are " << S << " records" << endl;

     int j, l;
     vec<int> clique_start, clique_length;
     vec<clique_record> clique_records;
     for ( int i = 0; i < (int) S; i++ )
     {    for ( j = i+1; j < (int) S; j++ )
          {    for ( l = (K+3)/4 - 1; l >= 0; l-- )
                         if ( R[j].Bytes( )[l] != R[i].Bytes( )[l] ) break;
                    if ( l >= 0 ) break;    }
          if ( j - i > max_clique )
          {    sort( R.begin( ) + i, R.begin( ) + j, kmer_record<K,I>::id_cmp );
               int duplicate_ids = 0;
               for ( int m = i + 1; m < j; m++ )
                    if ( R[m].GetId( ) == R[m-1].GetId( ) ) ++duplicate_ids;

               if ( float(duplicate_ids) / float(j-i) < 0.2 )
               {    cout << "clique " << clique_start.size( ) << ": " 
                         << j-i << " entries, of which "
                         << duplicate_ids << " are duplicates, ";
                    const unsigned char* r = R[i].Bytes( );
                    cout << "kmer = ";
                    for ( int u = 0; u < (K+3)/4; u++ )
                    {    for ( int w = 0; w < 8; w += 2 )
                              cout << as_base(int(((r[u] >> w) & 3)));    }
                    cout << "\n";
                    cout << "   ids:";
                    for ( int m = i; m < j; m++ )
                         cout << " " << R[m].GetId( );
                    cout << "\n";    

                    // Select three reads from cluster: read with leftmost extent,
                    // read with rightmost extent, and read with max extent in both
                    // directions.

                    int left_best = -1, left_best_id = 0;
                    int right_best = -1, right_best_id = 0;
                    int mid_best = -1, mid_best_id = 0;
                    for ( int m = i; m < j; m++ )
                    {    int pos = R[m].GetPos( );
                         int id = R[m].GetId( );
                         int left, right;
                         if ( pos >= 0 )
                         {    left = pos;
                              right = EE[id]->size( ) - pos - K;    }
                         else
                         {    right = - pos - 1;
                              left = EE[id]->size( ) + pos + 1 - K;    }
                         if ( left > left_best ) 
                         {    left_best = left;
                              left_best_id = id;   }
                         if ( right > right_best )
                         {    right_best = right;
                              right_best_id = id;    }
                         if ( left > mid_best && right > mid_best )
                         {    mid_best = Min( left, right );
                              mid_best_id = id;    }    }
                    cout << "left = " << left_best_id 
                              << " (" << left_best << " bases)"
                         << ", mid = " << mid_best_id
                              << " (" << mid_best << " bases)"
                         << ", right = " << right_best_id 
                              << " (" << right_best << " bases)"
                         << "\n";
                    clique_records.push_back( 
                         clique_record( clique_start.size( ), left_best_id ) );
                    clique_records.push_back( 
                         clique_record( clique_start.size( ), mid_best_id ) );
                    clique_records.push_back( 
                         clique_record( clique_start.size( ), right_best_id ) );
                    clique_start.push_back(i);
                    clique_length.push_back(j-i);    }    }

          i = j - 1;    }

     sort( clique_records.begin( ), clique_records.end( ) );
     vec<int> clique_records_reads( clique_records.size( ) );
     for ( unsigned int i = 0; i < clique_records.size( ); i++ )
          clique_records_reads[i] = clique_records[i].read_id;
     for ( unsigned int i = 0; i < clique_start.size( ); i++ )
     {    int cstart = clique_start[i], clength = clique_length[i];
          vec<int> overlapping_cliques;
          for ( int j = 0; j < clength; j++ )
          {    int id = R[ cstart + j ].GetId( );
               int p = BinPosition( clique_records_reads, id );
               if ( p >= 0 )
               {    while(1)
                    {    --p;
                         if ( p < 0 ) break;
                         if ( clique_records_reads[p] != id ) break;    }
                    ++p;
                    while(1)
                    {    int cid = clique_records[p].clique_id;
                         if ( cid != (int) i ) 
                              overlapping_cliques.push_back( cid );    
                         ++p;
                         if ( p == (int) clique_records_reads.size( )
                              || clique_records_reads[p] != id ) 
                              break;    }    }    }
          cout << i << ": ";
          UniqueSort( overlapping_cliques );
          vec<int> our_ids;
          for ( int j = 0; j < clength; j++ )
          {    int id = R[ cstart + j ].GetId( );
               our_ids.push_back(id);    }
          vec<float> hit_ratio;
          for ( int j = 0; j < (int) overlapping_cliques.size( ); j++ )
          {    int hits = 0;
               int m = overlapping_cliques[j];
               int cs = clique_start[m], cl = clique_length[m];
               for ( int l = 0; l < cl; l++ )
               {    int id2 = R[ cs + l ].GetId( );
                    if ( BinPosition( our_ids, id2 ) >= 0 ) ++hits;    }
               hit_ratio.push_back( 100.0 * float(hits) / float(cl) );    }
          for ( int j = 0; j < (int) overlapping_cliques.size( ); j++ )
               if ( hit_ratio[j] >= 75.0 )
                    cout << " " << overlapping_cliques[j] 
                         << " (" << setprecision(3) << hit_ratio[j] << "%)";
          cout << "\n";    }

     // calling exit(0) prevents buffers from flushing in cxx!!! -jbutler 2001/11/08
     cout << flush;
     cerr << flush;
     exit(0);    }

#define INSTANTIATE_PROCESS_FREQUENT_KMERS(K,I )                                  \
     template void ProcessFrequentKmers<I, K>( const vecbasevector& EE,          \
     vec< kmer_record<K,I> > R, unsigned int S, int max_clique, String run_dir )

//#ifndef __DECCXX_VER

FOR_ALL_K(INSTANTIATE_PROCESS_FREQUENT_KMERS, 1);
FOR_ALL_K(INSTANTIATE_PROCESS_FREQUENT_KMERS, 2);

//#endif
