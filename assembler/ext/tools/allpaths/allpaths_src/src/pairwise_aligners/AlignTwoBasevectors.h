// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef ALIGNTWOBASEVECTORS
#define ALIGNTWOBASEVECTORS

#include <fstream>

#include "Basevector.h"
#include "pairwise_aligners/MakeAlignsMethod.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "system/Types.h"

void MergeTwoBaseVectors( const basevector& b1, const basevector& b2, 
     const align& a, basevector& c, const qualvector& q1,
     const qualvector& q2, qualvector& q );

// Beware of the default parameter, mode = 1.  This will discard alignments where 
// b2 is RC relative to b1 or where b2 hangs off the front end of b2.  A safer 
// default would be mode = 0, which discards no alignments at all, but it is not 
// known whether this would cause problems, i.e., one should not change defaults lightly.

// If you pass a method_ptr, the parameters stretch, nstretch,
// max_err8, local_max_errors, alt_method, bandwidth, sw_gap_method,
// and min_progression_ratio are ignored, and the specified method is
// used to build aligns instead.

int AlignTwoBasevectors( const basevector& b1, const basevector& b2, align& a, 
     int min_overlap, int max_overlap, float max_error_rate, ostream* log, int& RC,
     int mode = 1, int K = 24, int stretch = 2, int nstretch = 1,
     const qualvector& q1 = qualvector(0), const qualvector& q2 = qualvector(0),
     float max_score = 0, float max_errors = 0, ostream& badlog = cout,
     Bool avoid_promiscuous_kmers = False, int max_cliq8 = 1000,
     int max_aligns8 = 10000, int max_err8 = 1000, int local_max_errors = 50,
     Bool alt_method = False, int bandwidth = 0, int min_mutmer = 0,
     float ambiguity_threshold = 3.0, float max_polyscore_percent = 100.0,
     Bool sw_gap_method = False, double min_progression_ratio = 0.4,
     makealigns_method *method_ptr = 0, Bool assert_if_problem = False,
     vec< vec< vec<int> > >* allowed_offsets = 0, int max_offset_discrep = 0);

#endif
