// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// #define USE_TIMERS

#include <stdlib.h>
#include <sys/types.h>

#include <algorithm>

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "kmers/GetNextKmerPair.h"
#include "kmers/KmerRecord.h"
#include "pairwise_aligners/MakeAlignsStatic.h"
#include "pairwise_aligners/MutmerGraph.h"
#include "PackAlign.h"
#include "pairwise_aligners/ProcessFrequentKmers.h"
#include "ShortVector.h"
#include "kmers/SortKmers.h"

// This is like MakeAligns, but has the following differences:
//
// 1. R, rid, multiplicities, mid, aligns, and mm are static.
//
// 2. Rather than generate an output file of alignments, MakeAlignsStatic
// return its result in the vectors answer and answer_counts.  Note that the
// vector answer may contain extra stuff at its end.  Use answer_counts to know
// how much stuff is actually part of the MakeAlignsStatic output.
//
// 3. Deleted arguments: aligns_file, dump_multiplicities.
//
// 4. Silent (logging removed).
//
// 5. Miscellaneous small changes.
//
// There may be other differences, if MakeAligns.cc has been modified in the
// meantime.

static makealigns_orig_method *default_makealigns_method = new makealigns_orig_method;

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic( 
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, 
     int total_passes, vec<bvec const*>const& EE,
     to_compare which_to_compare, int max_clique, int max_badness,
     int max_errs, int max_alignments,
     int min_mutmer, int local_max_errs, int stretch, int end_stretch,
     int nstretch, int local_max_errs_done,
     Bool avoid_promiscuous_kmers, int cl, int bandwidth,
     vec< vec< vec<int> > > *allowed_offsets, int max_offset_discrep,
     Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B,
     const vec< vec<Bool> >* allowed_alignments )
{
    if ( local_max_errs_done == 0 ) local_max_errs_done = 100000000;

    default_makealigns_method->SetMaxBadness( max_badness );
    default_makealigns_method->SetMaxErrs( max_errs );
    default_makealigns_method->SetLocalMaxErrs( local_max_errs );
    default_makealigns_method->SetLocalMaxErrsDone( local_max_errs_done );
    default_makealigns_method->SetStretch( stretch );
    default_makealigns_method->SetEndStretch( end_stretch );
    default_makealigns_method->SetNStretch( nstretch );
    default_makealigns_method->SetCl( cl );
    default_makealigns_method->SetMaxAlignCtorCalls( 1000000 );

    MakeAlignsStatic<I,k,BLOCKS_PER_NODE>( answer, answer_counts,
                                           Passes, total_passes, EE, which_to_compare,
                                           default_makealigns_method,
                                           max_clique, max_alignments, min_mutmer,
                                           avoid_promiscuous_kmers, allowed_offsets,
                                           max_offset_discrep, process_frequent_kmers,
                                           run_dir, offset_A, offset_B,
                                           allowed_alignments );
}

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic(
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes,
     int total_passes, const vecbasevector& EE,  
     to_compare which_to_compare, int max_clique, int max_badness, 
     int max_errs, int max_alignments, 
     int min_mutmer, int local_max_errs, int stretch, int end_stretch, 
     int nstretch, int local_max_errs_done, 
     Bool avoid_promiscuous_kmers, int cl, int bandwidth,
     vec< vec< vec<int> > > *allowed_offsets, int max_offset_discrep, 
     Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B,
     const vec< vec<Bool> >* allowed_alignments )
{
    if ( local_max_errs_done == 0 ) local_max_errs_done = 100000000;

    default_makealigns_method->SetMaxBadness( max_badness );
    default_makealigns_method->SetMaxErrs( max_errs );
    default_makealigns_method->SetLocalMaxErrs( local_max_errs );
    default_makealigns_method->SetLocalMaxErrsDone( local_max_errs_done );
    default_makealigns_method->SetStretch( stretch );
    default_makealigns_method->SetEndStretch( end_stretch );
    default_makealigns_method->SetNStretch( nstretch );
    default_makealigns_method->SetCl( cl );
    default_makealigns_method->SetMaxAlignCtorCalls( 1000000 );
    
    MakeAlignsStatic<I,k,BLOCKS_PER_NODE>( answer, answer_counts,
                                           Passes, total_passes, EE, which_to_compare,
                                           default_makealigns_method,
                                           max_clique, max_alignments, min_mutmer,
                                           avoid_promiscuous_kmers, allowed_offsets,
                                           max_offset_discrep, process_frequent_kmers,
                                           run_dir, offset_A, offset_B,
                                           allowed_alignments );
}



template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic(
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, 
     int total_passes, const vecbasevector& EE, to_compare which_to_compare, 
     makealigns_method *method, int max_clique, int max_alignments, int min_mutmer,
     Bool avoid_promiscuous_kmers, vec< vec< vec<int> > > *allowed_offsets, 
     int max_offset_discrep, Bool process_frequent_kmers,
     String run_dir, int offset_A, int offset_B,
     const vec< vec<Bool> >* allowed_alignments )
{
    vec<bvec const*> vvv; vvv.reserve(EE.size());
    for ( vecbvec::const_iterator end(EE.end()), itr(EE.begin()); itr != end; ++itr )
        vvv.push_back(&*itr);
    MakeAlignsStatic<I,k,BLOCKS_PER_NODE>( answer, answer_counts,
                                           Passes, total_passes, vvv, which_to_compare,
                                           method,
                                           max_clique, max_alignments, min_mutmer,
                                           avoid_promiscuous_kmers, allowed_offsets,
                                           max_offset_discrep, process_frequent_kmers,
                                           run_dir, offset_A, offset_B,
                                           allowed_alignments );
}

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic(
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes,
     int total_passes, vec<bvec const*>const& EE, to_compare which_to_compare,
     makealigns_method *method, int max_clique, int max_alignments, int min_mutmer,
     Bool avoid_promiscuous_kmers, vec< vec< vec<int> > > *allowed_offsets,
     int max_offset_discrep, Bool process_frequent_kmers,
     String run_dir, int offset_A, int offset_B,
     const vec< vec<Bool> >* allowed_alignments )
{    START_TIMER( m1, 200 );
     int N = EE.size( );
     int N0 = which_to_compare.N0;
     ForceAssert( N0 > 0 || which_to_compare.compare_type == ALL_VS_ALL );

     if ( I == 1 )
     {    for ( int i = 0; i < N; i++ )
               if ( EE[i]->size( ) >= 1024 )
                    FatalErr( "I'm sorry, but MakeAlignsStatic<1> has been passed "
                         << "a read of length " << EE[i]->size( )
                         << ".  For sequences of length >= 1024, you need to use "
                         << "MakeAlignsStatic<2>." );    }

     static mutmer_graph<I, BLOCKS_PER_NODE> M;
     M.clear_and_resize(N);

     static vec<int> rid;
     rid.resize(N);
     for ( int i = 0; i < N; i++ )
          rid[i] = i;
     longlong S_init = 0;
     for ( size_t l = 0; l < EE.size( ); l++ )
          if ( EE[l]->size( ) >= k ) S_init += EE[l]->size( ) - k;
     S_init += S_init/4;
    
     ForceAssert( Passes == 1 || Passes == 10 || Passes == 100 );
     if ( Passes == 10 ) S_init /= 7;
     if ( Passes == 100 ) S_init /= 33;

     if ( S_init > 3000000000u )
          FatalErr( "MakeAlignsStatic: the value of S is " << S_init
               << ".  This is dangerously large, because it represents\nthe "
               << "size of a vector which is to be stored in an unsigned int,\n"
               << "and this vector could grow. "
               << "You may need to rewrite some code to get it to work." );

     unsigned int S = S_init, max_S = 0, total_S = 0;
     STOP_TIMER(m1);
     START_TIMER(m2, 200);

     {    static vec< kmer_record<k,I> > R;
          R.resize(S);
          static vec<int> multiplicities;
          multiplicities.resize(0);
          for ( int pass = 0; pass < total_passes; pass++ )
          {    dummy<1> d1;
               dummy<10> d10;
               dummy<100> d100;
               if ( Passes == 1 ) SortKmers( d1, EE, rid, pass, R, S );
               else if ( Passes == 10 ) SortKmers( d10, EE, rid, pass, R, S );
               else SortKmers( d100, EE, rid, pass, R, S );
               Assert( R.size( ) >= S );
               if (process_frequent_kmers) 
                    ProcessFrequentKmers( EE, R, S, max_clique, run_dir );
               max_S = max(S, max_S);
               total_S += S;
	       // gcc 4.2 warns pos1 may be used uninit'd, but in reality, this
	       // will never happen.  So, we init to 0 to suppress the warn.
               int read_id1(-1), read_id2, pos1=0, pos2;
               Bool dump_multiplicities = False;
               while(1)
               {    
		    if ( which_to_compare.compare_type != FIRST_VS_SECOND )
                         GetNextKmerPair<k,I>( R, S, read_id1, read_id2, pos1, pos2, 
                              max_clique, multiplicities, dump_multiplicities,
                              avoid_promiscuous_kmers );
                    else GetNextKmerPair<k,I>( N0, R, S, read_id1, read_id2, pos1, 
                              pos2, max_clique, multiplicities, dump_multiplicities,
                              avoid_promiscuous_kmers );

                    if ( read_id1 < 0 ) break;

                    if ( N0 > 0 )
                    {    
                         if ( which_to_compare.compare_type == FIRST_VS_ALL
                              && read_id1 >= N0 && read_id2 >= N0 ) continue;

			 if ( allowed_offsets )
                         {    int p1 = pos1;
                              int p2 = pos2;
                              
                              int rc = (p1 < 0) ^ (p2 < 0);
                              if ( p1 < 0 ) p1 = -p1;
                              if ( p2 < 0 ) p2 = -p2;
                              --p1;
                              --p2;

                              if ( read_id1 < N0 )
                              {    if (rc) p1 = EE[read_id1]->size( ) - k - p1;    }
                              else
                              {    swap( p1, p2 );
                                   if (rc) p1 = EE[read_id2]->size( ) - k - p1;    }

                              int offset = p1 - p2;
                              vec<int>& v = (read_id1 < N0)
                                   ? (*allowed_offsets)[read_id1][read_id2 - N0]
                                   : (*allowed_offsets)[read_id2][read_id1 - N0];
                              unsigned int i;
                              for ( i = 0; i < v.size( ); i++ )
                                   if ( Abs(offset - v[i]) <= max_offset_discrep )
                                        break;
                              if ( i == v.size( ) ) continue;    
			 }   
		    }

                    if ( allowed_alignments != 0 &&
                         !(*allowed_alignments)[read_id1][read_id2] ) continue;
                    M.MergeKmer( pos1, pos2, k, read_id1, read_id2, 
                         *EE[read_id1], *EE[read_id2] );
	       }
	  }
     }

     STOP_TIMER(m2);
     START_TIMER(m3, 200);

     int answer_size = 0;
     answer_counts.resize(0);

     if ( total_S == 0 )
       return;
     /*
     {    cout << "MakeAlignsStatic: I've been passed sequences of length";
          for ( int l = 0; l < EE.size( ); l++ )
               cout << " " << EE[l].size( );
          cout << "\nI can't work with this data.\n";
          ForceAssert( 0 == 1 );    }
     */

     static vec< mutmer_read_id<I> > mid;
     mid.resize( Max(M.Counts( )) );
     static vec<align> aligns;
     aligns.resize(max_alignments);
     static vec<int> errors;
     errors.resize(max_alignments);

     STOP_TIMER(m3);

     START_TIMER(m4, 200);
     int aligns_length, j;
     basevector rcrd2;
     for ( int l = 0; l < N; l++ )
     {    int count = M.All( l, mid );
          sort( mid.begin( ), mid.begin( ) + count );
          for ( int i = 0; i < count; i++ )
          {    for ( j = i+1; j < count; j++ )
                    if ( mid[i].ReadIdRc( ) != mid[j].ReadIdRc( ) ) break;
               static vec<mutmer> mm;
               mm.resize(j-i);
               for ( int r = 0; r < j-i; r++ )
               {    int pos1, pos2, len, e;
                    mid[i+r].Unpack( pos1, pos2, len, e );
                    mm[r].SetFrom( pos1, pos2, len, e );     }
               Bool RC = mid[i].Rc( );
               int id1 = l, id2 = mid[i].ReadId( );
               Bool succeed = False;
               int cur_min_mutmer = min_mutmer;
               while( !succeed )
               {    if ( !RC ) {
                         succeed =
			   method->MutmersToAlign( mm, k, *EE[id1], *EE[id2], aligns,
                                   errors, aligns_length, cur_min_mutmer, 0 );
	            }
                    else {
		      rcrd2.Setsize(EE[id2]->size());
		      rcrd2.ReverseComplement(*EE[id2]);
  		      succeed =
                         method->MutmersToAlign( mm, k, *EE[id1], rcrd2, aligns,
                                   errors, aligns_length, cur_min_mutmer, 0 );
		    }
                    cur_min_mutmer += 15;    }
               if ( aligns_length > 0 ) 
               {    answer_counts.push_back(aligns_length);
   
                    // Fix id1 and id2.

                    int save_id1 = id1;
	            if ( save_id1 < N0 ) save_id1 += offset_A;
                    else save_id1 += offset_B;
		      
                    int save_id2 = id2;
	            if ( save_id2 < N0 ) save_id2 += offset_A;
  	            else save_id2 += offset_B;	  
		      
	            // Write alignment.

                    int space_needed = answer_size + aligns_length;
                    if ( space_needed > (int) answer.size( ) )
                         answer.resize( space_needed + space_needed/2 );
		    for ( int q = 0; q < aligns_length; q++ )
                         answer[answer_size++].SetFrom( aligns[q], RC, 
                              save_id1, save_id2 );    }
               i = j - 1;    }    }    
     STOP_TIMER(m4);
     }

#ifdef __GNUC__

#define INSTANTIATE_MAKEALIGNS_STATIC(I, K, BLOCKS_PER_NODE)            \
template void MakeAlignsStatic<I, K, BLOCKS_PER_NODE>( \
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, \
     int total_passes, vec<bvec const*>const& EE, to_compare which_to_compare, \
     makealigns_method *method, int max_clique, \
     int max_alignments, int min_mutmer, \
     Bool avoid_promiscuous_kmers, \
     vec< vec< vec<int> > > *allowed_offsets, \
     int max_offset_discrep, Bool process_frequent_kmers, \
     String run_dir, int offset_A, int offset_B, \
     const vec< vec<Bool> >* allowed_alignments ); \
template void MakeAlignsStatic<I, K, BLOCKS_PER_NODE>( \
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, \
     int total_passes, vec<bvec const*>const& EE, \
     to_compare which_to_compare, int max_clique, int max_badness, \
     int max_errs, int max_alignments, \
     int min_mutmer, int local_max_errs, int stretch, int end_stretch, \
     int nstretch, int local_max_errs_done, \
     Bool avoid_promiscuous_kmers, int cl, int bandwidth,\
     vec< vec< vec<int> > > *allowed_offsets, int max_offset_discrep, \
     Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B, \
     const vec< vec<Bool> >* allowed_alignments ); \
template void MakeAlignsStatic<I, K, BLOCKS_PER_NODE>( \
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, \
     int total_passes, const vecbasevector& EE, to_compare which_to_compare, \
     makealigns_method *method, int max_clique, \
     int max_alignments, int min_mutmer, \
     Bool avoid_promiscuous_kmers, \
     vec< vec< vec<int> > > *allowed_offsets, \
     int max_offset_discrep, Bool process_frequent_kmers, \
     String run_dir, int offset_A, int offset_B, \
     const vec< vec<Bool> >* allowed_alignments ); \
template void MakeAlignsStatic<I, K, BLOCKS_PER_NODE>( \
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, \
     int total_passes, const vecbasevector& EE, \
     to_compare which_to_compare, int max_clique, int max_badness, \
     int max_errs, int max_alignments, \
     int min_mutmer, int local_max_errs, int stretch, int end_stretch, \
     int nstretch, int local_max_errs_done, \
     Bool avoid_promiscuous_kmers, int cl, int bandwidth,\
     vec< vec< vec<int> > > *allowed_offsets, int max_offset_discrep, \
     Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B, \
     const vec< vec<Bool> >* allowed_alignments )

     #define INSTANTIATE_MAKEALIGNS_STATIC_FOR_K(K,dummy) \
        INSTANTIATE_MAKEALIGNS_STATIC( 1, K, 50 ); \
        INSTANTIATE_MAKEALIGNS_STATIC( 2, K, 50 )

     FOR_ALL_K(INSTANTIATE_MAKEALIGNS_STATIC_FOR_K,unused);
#endif

#ifdef __DECCXX_VER

#pragma define_template MakeAlignsStatic< 2, 8, 50>
#pragma define_template MakeAlignsStatic< 1, 12, 50>
#pragma define_template MakeAlignsStatic< 2, 12, 50>
#pragma define_template MakeAlignsStatic< 2, 16, 50>
#pragma define_template MakeAlignsStatic< 2, 24, 50>
#pragma define_template MakeAlignsStatic< 2, 48, 50>

#endif
