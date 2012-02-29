///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REMODEL_GAP_TOOLS_H
#define REMODEL_GAP_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "kmers/KmerRecord.h"

// GapComp
//
// d = possible gap
// a = unordered list of left co-starts
// b = unordered list of right stops
// X = histogram of lengths
// p = pairs (i,j) with i indexing an element of a, j indexing an element of b
// L1 = length of left contig
// L2 = length of right contig

vec<long double> GapComp( const vec<int> D, const vec<int>& a, const vec<int>& b, 
     const vec<int>& x, const vec<double>& X, const vec< pair<int,int> >& p, 
     const int L1, const int L2, const int VERBOSITY );

template<int K> void GetBounds( const vec< kmer<K> >& kmers, kmer<K>& x,
     kmer<K>& xrc, Bool& fw, int64_t& low, int64_t& high )
{    xrc = x;
     xrc.ReverseComplement( );
     fw = ( x < xrc );
     if ( !fw ) x = xrc;
     low = lower_bound( kmers.begin( ), kmers.end( ), x ) - kmers.begin( );
     high = upper_bound( kmers.begin( ), kmers.end( ), x ) - kmers.begin( );    }

void AlignReads( const int K1, const int K2, const Bool use_tail, 
     const vecbasevector& tigs, const vecbasevector& bases, 
     const PairsManager& pairs, vec<Bool>& placed_fw, vec<Bool>& placed_rc, 
     vec< pair<int,int> >& placement, vec< vec<longlong> >& aligns_index, 
     const Bool TIME_STAMPS );

template<int K> void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus );

template<int K> void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus );

void FindTrueGaps(  
     // inputs:
     const vecbasevector& genome,
     const vecbasevector& tigs,  
     const vec<superb>& scaffolds,
     // control:
     const vec<int>& tigs_to_process, const Bool TIME_STAMPS,
     // outputs:
     vec<int>& true_gap, vec<Bool>& true_gap_computed );

void PrintSideBySide( const String& z1, const String& z2, const int N );

// A gap_id represents either a gap between contigs, or else a 
// "gap within a contig", defined by a range of bases on it.

class gap_id {

     public:

     enum GapType { BETWEEN, WITHIN };

     gap_id( ) { }

     gap_id( const GapType t, const int m1, const int m2 ) 
          : t_(t), m1_or_m_(m1), m2_(m2)
     {    ForceAssert( t == BETWEEN );    }

     gap_id( const GapType t, const int m, const int start, const int stop )
          : t_(t), m1_or_m_(m), start_(start), stop_(stop)
     {    ForceAssert( t == WITHIN );    } 

     GapType Type( ) const { return t_; }

     int M1( ) const
     {    ForceAssert( t_ == BETWEEN );
          return m1_or_m_;    }
     int M2( ) const
     {    ForceAssert( t_ == BETWEEN );
          return m2_;    }

     int M( ) const
     {    ForceAssert( t_ == WITHIN );
          return m1_or_m_;    }
     int Start( ) const
     {    ForceAssert( t_ == WITHIN );
          return start_;    }
     int Stop( ) const
     {    ForceAssert( t_ == WITHIN );
          return stop_;    }

     private:

     GapType t_;
     int m1_or_m_;
     int m2_;
     int start_;
     int stop_;

};

template<int K> void PredictGap( 
     // inputs:
     const String& libtype,
     const vec< triple<kmer<K>,int,int> >& kmers_plus, const vec< kmer<K> >& kmers,
     const int max_overlap, const int nlibs, const int lib_offset,
     const vec<int>& libs_to_use, const vec<Bool>& lfail, const vecbasevector& tigs,
     const vec< vec<int> >& DISTS, const vec< vec<double> >& X, 
     const vecbasevector& bases,
     const PairsManager& pairs, const vec<unsigned short>& read_lengths, 
     const vec< vec<longlong> >& aligns_index, const vec<Bool>& placed_fw, 
     const vec<Bool>& placed_rc, const vec< pair<int,int> >& placement, 
     const gap_id& gx, const int VERBOSITY, 
     // outputs:
     Bool& gap_predicted, int& gap, double& dev );

int FirstSumTo( const vec<long double>& x, long double bound );

void DefineDistributions(
     // inputs:
     const int nlibs, const int lib_offset, const vecbasevector& tigs, 
     const vec<unsigned short>& read_lengths, const PairsManager& pairs, 
     vec<Bool>& placed_fw, vec<Bool>& placed_rc, vec< pair<int,int> >& placement, 
     vec< vec<longlong> >& aligns_index, const Bool TIME_STAMPS,
     const Bool HISTOGRAMS,
     // outputs:
     vec<Bool>& lfail, vec<int>& max_dist, vec< vec<int> >& DISTS,
     vec< vec<double> >& X );

void WriteReportAboutGaps( vec<int>& devs, vec<int>& adevs, vec<double>& offbys, 
     vec<double>& aoffbys, const vec< triple<int,int,int> >& results, 
     const vec<int>& true_gap, const vec<Bool>& true_gap_computed, 
     const int fails, const vec<int>& tigs_to_process );

void GetLibCount( const String& run_dir, const String& head, 
		  int& nlibs, const int VERBOSITY, const Bool TIME_STAMPS );


void GetLibCounts( const String& run_dir, const vec<String>& lib_types_to_use,
		   const int VERBOSITY, const Bool TIME_STAMPS, 
		   const String& frags, const String& jumps, const String& long_jumps,
		   int& nlibs_frag, int& nlibs_jump, int& nlibs_long );

template<int K> void ComputeOverlaps( const basevector& M1,
     const basevector& M2, const int m1, const vec<int>& max_dist,
     const vec< triple<kmer<K>,int,int> >& kmers_plus, const vec< kmer<K> >& kmers,
     const int VERBOSITY, vec<int>& accepted_overlaps, int& max_overlap );

#endif
