///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef MAKEALIGNSMETHOD_H
#define MAKEALIGNSMETHOD_H

#include "Basevector.h"
#include "PackAlign.h"
#include "pairwise_aligners/Mutmer.h"
#include "system/Types.h"
#include "Vec.h"

// Allows the specification of the alignment creation method in MakeAligns.

// The orig method creates paths of mutmers and tries to walk between
// them adding small gaps.

// The alt method takes the largest mutmers and performs banded
// smith-waterman alignments at the offset implied by those mutmers.

// The sw_gap method creates paths of non-intersecting mutmers and
// performs banded smith-waterman alignments between them to fill in
// the gaps.

class makealigns_method 
{
 public:
  virtual ~makealigns_method() {};
  
  virtual Bool MutmersToAlign( const vec<mutmer>& mm, int k,
			      const basevector& rd1, const basevector& rd2, 
			       vec<align>& aligns, vec<int>& errors, int& aligns_length, 
			       int min_mutmer, ostream* log ) = 0;
};


class makealigns_orig_method : public makealigns_method 
{
 public:
  makealigns_orig_method() : 
    makealigns_method(), 
    local_max_errors_done_(100000000),
    cl_(20), max_alignment_constructor_calls_( 1000000 ),
    try_if_no_proper_(False)
    { }

  Bool MutmersToAlign( const vec<mutmer>& mm, int k,
		       const basevector& rd1, const basevector& rd2, 
		       vec<align>& aligns, vec<int>& errors, int& aligns_length, 
		       int min_mutmer, ostream* log );

  void SetMaxBadness( float max_badness )    { max_badness_ = max_badness; }
  void SetMaxErrs( int max_errs )            { max_errs_ = max_errs; }
  void SetLocalMaxErrs( int local_max_errs ) { local_max_errs_ = local_max_errs; }
  void SetLocalMaxErrsDone( int lmed )       { if ( lmed == 0 ) lmed = 100000000;
                                               local_max_errors_done_ = lmed; }
  void SetStretch( int stretch )             { stretch_ = stretch; }
  void SetNStretch( int nstretch )           { nstretch_ = nstretch; }
  void SetEndStretch( int end_stretch )      { end_stretch_ = end_stretch; }
  void SetCl( int cl )                       { cl_ = cl; }
  void SetMaxAlignCtorCalls( int max )       { max_alignment_constructor_calls_ = max; }

  void SetTryIfNoProper( int perf, int bandwidth )
  {    try_if_no_proper_ = True;
       try_if_no_proper_perf_ = perf;
       try_if_no_proper_bandwidth_ = bandwidth;    }
  
 private:
  float max_badness_;
  int max_errs_;
  int local_max_errs_;
  int local_max_errors_done_;
  int stretch_;
  int end_stretch_;
  int nstretch_;
  int cl_;
  int max_alignment_constructor_calls_;
  Bool try_if_no_proper_;
  int try_if_no_proper_perf_, try_if_no_proper_bandwidth_;
};


class makealigns_alt_method : public makealigns_method 
{
 public:
  Bool MutmersToAlign( const vec<mutmer>& mm, int k,
		       const basevector& rd1, const basevector& rd2, 
		       vec<align>& aligns, vec<int>& errors, int& aligns_length, 
		       int min_mutmer, ostream* log );
  
  void SetMaxErrs( int max_errs )     { max_errs_ = max_errs; }
  void SetBandwidth( int bandwidth )  { bandwidth_ = bandwidth; }
  
 private:
  int max_errs_;
  int bandwidth_;
};


// Only produce perfect proper alignments.

class makealigns_perfect_method : public makealigns_method 
{
 public:
  Bool MutmersToAlign( const vec<mutmer>& mm, int k,
		       const basevector& rd1, const basevector& rd2, 
		       vec<align>& aligns, vec<int>& errors, int& aligns_length, 
		       int min_mutmer, ostream* log );
};


class makealigns_sw_gap_method : public makealigns_method
{
 public:
  makealigns_sw_gap_method() :
    makealigns_method(),
    min_max_mutmer_length_( 200 ),
    min_prog_length_( 0 ),
    min_prog_ratio_( 0.8 ),
    min_overlap_frac_( 0.2 ),
    ignore_overlap_frac_length_( INT_MAX ),
    max_mutmer_offset_diff_( 1000 ),
    max_gap_( 1000 ),
    verbose_( false ),
    affine_penalties_( false )
    { }

  Bool MutmersToAlign( const vec<mutmer>& mm, int k,
		       const basevector& rd1, const basevector& rd2, 
		       vec<align>& aligns, vec<int>& errors, int& aligns_length, 
		       int min_mutmer, ostream* log );

  void SetMaxErrs( int max_errs )           { max_errs_ = max_errs; }
  void SetEndStretch( int end_stretch )     { end_stretch_ = end_stretch; }

  // Minimum length of largest mutmer in a usable path.
  void SetMinMaxMutmerLength( int min )     { min_max_mutmer_length_ = min; }
  
  // Only progressions longer than this will be considered.
  void SetMinProgressionLength( int min )   { min_prog_length_ = min; }

  // Only progressions longer than this fraction of the predicted
  // overlap will be considered.
  void SetMinOverlapFraction( double min )  { min_overlap_frac_ = min; }
  
  // If the progression is longer than this amount, do not filter it
  // based on minimum overlap fraction.
  void SetIgnoreOverlapFractionLength( int min )
    { ignore_overlap_frac_length_ = min; }

  // Minimum ratio of mutmers to total sequence in a mutmer vector
  // (Though currently overlaps are double counted. This is just to
  // filter out accidentally consistent mutmers with a huge gap
  // between them.
  void SetMinProgressionRatio( double min ) { min_prog_ratio_ = min; }
  
  // Maximum absolute difference in offsets for inclusion in a
  // progression.
  void SetMaxMutmerOffsetDiff( int max )    { max_mutmer_offset_diff_ = max; }

  // Maximum gap between mutmers in a progression.
  void SetMaxGap( int max )                 { max_gap_ = max; }
  
  // Voluminous but occasionally helpful logging is turned on or off with this switch.
  void SetVerbose( bool verbose )           { verbose_ = verbose; }

  // Use affine gap penalties in SW alignments
  void SetAffinePenalties( bool affine_penalties ) { affine_penalties_ = affine_penalties; }

 private:
  int max_errs_;
  int end_stretch_;
  int min_max_mutmer_length_;
  int min_prog_length_;
  double min_prog_ratio_;
  double min_overlap_frac_;
  int ignore_overlap_frac_length_;
  int max_mutmer_offset_diff_;
  int max_gap_;
  bool verbose_;
  bool affine_penalties_;
};


#endif
