// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#ifndef FIND_GAPS_H
#define FIND_GAPS_H

#include "AnAssemblyClass.h"

void FindGaps( assembly& A, const vec<int>& simple_reads_orig_index, int s, 
     vec<int>& gaps_dev, int start_contig, int stop_contig );

// GapStats: Compute the approximate gap_ave and gapdev_ave,
// discarding outliers if necessary.

void GapStats( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave );

// Expect: compute the expected gap and s.d. for it, for two linked mtigs.

void Expect( assembly& a, int m1, int m2, int& gap_ave, int& gapdev_ave,
     Bool strong = False, Bool rc1 = False, Bool rc2 = False );

// Expect2: an improved (i.e. corrected) version of Expect, specialized to the
// case where rc1 = False.

void Expect2( assembly& a, int m1, int m2, int& gap_ave, int& gapdev_ave,
     Bool rc1 = False, Bool rc2 = False, Bool exit_if_no_pairings = True );

// ExpectedGap: compute the expected gap between two linked mtigs

int ExpectedGap( assembly& a, int m1, int m2 );
int ExpectedGap2( assembly& a, int m1, int m2 );

// ExpectedGapSuper: compute the expected gap between two supers, implied by 
// the links between two of their mtigs.

int ExpectedGapSuper( assembly& a, int m1, int m2, int& gap_dev, 
     Bool rc1 = False, Bool rc2 = False, Bool strong = False );

// ExpectedGapSuperFull: compute the expected gap between two supers, implied by 
// the links between all of their mtigs.  If it can't find any links, it sets
// gap_dev_ave to -1.  You have to pass it a list of indices of the relevant
// read_location's.

const int huge = 1000 * 1000 * 1000;

int ExpectedGapSuperFull( assembly& a, const super& sup1, const vec<int>& dev1,
     const super& sup2, const vec<int>& dev2, const vec<int>& read_indices, 
     int& gap_dev_ave, Bool rc1, Bool rc2, Bool verbose, Bool force_short, 
     const vec<int>& simple_reads_orig_index, int min_gap = -huge,
     int max_gap = huge );

#endif
