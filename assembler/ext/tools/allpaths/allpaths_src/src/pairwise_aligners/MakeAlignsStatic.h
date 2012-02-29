// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#ifndef MAKEALIGNS_STATIC_H
#define MAKEALIGNS_STATIC_H

#include "system/Assert.h"
#include "Basevector.h"
#include "pairwise_aligners/MakeAlignsMethod.h"
#include "pairwise_aligners/MakeAlignsToCompare.h"
#include "PackAlign.h"
#include "String.h"
#include "system/Types.h"
#include "Vec.h"

class align_plus {

     public:

     align a;
     int RC;
     int id1, id2;

     align_plus( ) { }
     align_plus( const align& a_arg,
		 int RC_arg,
		 int id1_arg,
		 int id2_arg ) 
       : a(a_arg),
	 RC(RC_arg),
	 id1(id1_arg),
	 id2(id2_arg)
     { }

     void SetFrom( const align& ax, int RCx, int id1x, int id2x ) 
     {    a = ax;
          RC = RCx;
          id1 = id1x;
          id2 = id2x;    }

};

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic( 
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, 
     int total_passes, vec<bvec const*> const& EE,
     to_compare which_to_compare, int max_clique, int max_badness,
     int max_errs = 31, int max_alignments = 255,
     int min_mutmer = 0, int local_max_errs = 31, int stretch = 2,
     int end_stretch = 2, int nstretch = 1,
     int local_max_errs_done = 0, Bool avoid_promiscuous_kmers = False, int cl = 20,
     int bandwidth = 0, vec< vec< vec<int> > >* allowed_offsets = 0,
     int max_offset_discrep = 0, Bool process_frequent_kmers = False,
     String run_dir = "", int offset_A = 0, int offset_B = 0,
     const vec< vec<Bool> >* allowed_alignments = 0 );

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic(
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes,
     int total_passes, vec<bvec const*> const& EE, to_compare which_to_compare,
     makealigns_method *method, int max_clique,
     int max_alignments = 255, int min_mutmer = 0,
     Bool avoid_promiscuous_kmers = False,
     vec< vec< vec<int> > >* allowed_offsets = 0,
     int max_offset_discrep = 0, Bool process_frequent_kmers = False,
     String run_dir = "", int offset_A = 0, int offset_B = 0,
     const vec< vec<Bool> >* allowed_alignments = 0 );

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic(
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes,
     int total_passes, const vecbasevector& EE,
     to_compare which_to_compare, int max_clique, int max_badness, 
     int max_errs = 31, int max_alignments = 255, 
     int min_mutmer = 0, int local_max_errs = 31, int stretch = 2, 
     int end_stretch = 2, int nstretch = 1, 
     int local_max_errs_done = 0, Bool avoid_promiscuous_kmers = False, int cl = 20, 
     int bandwidth = 0, vec< vec< vec<int> > >* allowed_offsets = 0,
     int max_offset_discrep = 0, Bool process_frequent_kmers = False,
     String run_dir = "", int offset_A = 0, int offset_B = 0,
     const vec< vec<Bool> >* allowed_alignments = 0 );

template<int I, int k, int BLOCKS_PER_NODE> void MakeAlignsStatic(
     vec<align_plus>& answer, vec<int>& answer_counts, int Passes, 
     int total_passes, const vecbasevector& EE, to_compare which_to_compare, 
     makealigns_method *method, int max_clique, 
     int max_alignments = 255, int min_mutmer = 0,
     Bool avoid_promiscuous_kmers = False,
     vec< vec< vec<int> > >* allowed_offsets = 0,
     int max_offset_discrep = 0, Bool process_frequent_kmers = False,
     String run_dir = "", int offset_A = 0, int offset_B = 0,
     const vec< vec<Bool> >* allowed_alignments = 0 );

#endif
