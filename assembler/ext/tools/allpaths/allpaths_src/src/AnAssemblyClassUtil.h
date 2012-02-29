// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef AN_ASSEMBLY_CLASS_UTIL
#define AN_ASSEMBLY_CLASS_UTIL

#include "AnAssemblyClass.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "Qualvector.h"

extern ostream *logp, *elogp;

// Load all_aligns and fill all_aligns_index.
void LoadAllAligns( const String &aligns_file,
		    int n_reads,
		    vec<brief_align> &all_aligns,
		    vec<int> &all_aligns_index );

// Load vec of supers (ann_supers_file: annotated supercontigs).
void LoadVecSupers( const String &ann_supers_file, vec<super> &supers );

// Reorder contigs in a supercontig by starting point.
void ReorderSuper( assembly& A, int s );

// Reorder contigs in a supercontig by minimizing negative gaps.
void ReorderSuperMaxGap( assembly& A, int s );

// Find how many contigs which are declared to overlap, fail to do so.
void CheckSuperOverlaps( ostream& log, assembly& Ass, int scg );

// Connected: determine if there is a link from mtig m1 to mtig m2 which is
// logical in the sense that the orientations of the two reads on the two mtigs
// are reverse to each other.

Bool Connected( assembly& a, int m1, int m2, Bool strong = False );

// ConnectionCount: count the number of links between mtig m1 and mtig
// m2 which are logically oriented.

int ConnectionCount( assembly& a, int m1, int m2 );

//  IllogicalLinksCount: count the number of links between mtig m1 and mtig
//  m2 which are illogically (i.e. same direction) oriented.

int IllogicalLinksCount( assembly& a, int m1, int m2 );

// ReadsOverlap: determine if a read from one mtig has a good overlap with a
// read from another mtig.

Bool ReadsOverlap( const assembly& a, int m1, int m2, Bool RC );

// OverlapUpdate is a wrapper for Overlap. The main extra feature is its
// no_overlap book-keeping: contig's pairs tested for but failing to overlap
// are saved as a vec of pairs.

int OverlapUpdate( vec< set<int> > &no_overlap,
     assembly& a, int m1, int m2, align& al, 
     int kmer_size = 12, int min_mutmer = 0, Bool RC = False, int mode = 1,
     int min_align = 100, int min_align_see = 30 );

// Overlap: return the apparent overlap between two mtigs (else -1).  The RC
// is set, the second mtig is treated as reversed.  If require_read_evidence = True,
// the code is faster but a bit less sensitive.

int Overlap( const assembly& a, int m1, int m2, align& al, 
     int kmer_size = 12, int min_mutmer = 0, Bool RC = False, int mode = 1,
     int min_align = 100, int min_align_see = 30, 
     Bool require_read_evidence = False );

// Reverse: reverse an mtig m.

void ReverseMtig( assembly& a, int m );

// Reverse: reverse a supercontig s.

void Reverse( assembly& a, int s );

// MergeMtigs: merge two mtigs m1, m2 in an assembly, using a given
// alignment, and update the assembly accordingly.  The mtig m2 is
// deleted.

void MergeMtigs( assembly& a, int m1, int m2, align& al );

// TryMergeMtigs: same as MergeMtigs, but does not modify the assembly.
// Produces as output a vec<read_location> new_tig and a basevector
// new_tig_bases.

void TryMergeMtigs( const assembly& a, int m1, int m2, const align& al,
     vec<read_location>& new_tig, basevector& new_tig_bases );

// RecomputeGaps: recompute gaps between mtigs in a supercontig, in a
// rather primitive way.  However, it is not so clear what the
// "correct" way is.

void RecomputeGaps( assembly& a, int supercontig_id, Bool reorder = True,
     Bool verbose = True );

// RecomputeGaps: recompute gaps between mtigs in all supercontigs.

void RecomputeGaps( assembly& a, Bool reorder = True, Bool verbose = True );

// PrintLinkAccuracy: Print link accuracy between contigs within
// supercontigs having more than one contig. Requires that logp (see
// above extern declaration) be defined and open.

void PrintLinkAccuracy( assembly& a );

// SpreadOverHoles: for a given supercontig, compute the spread over each
// gap.  This is what you get by looking at all the links over the gap, and
// on each side, looking at all the read ends.  Return the minimum of the span 
// in bases of the read ends on the left side of the gap, and the span in bases of 
// the read ends on the right side of the gap.  When we compute these spans, we
// don't count gaps between contigs.

void SpreadOverHoles( assembly& A, int s, vec<int>& spread );

// A supercontig is well connected iff for any pair of contigs (a, b) in it,
// it is possible to walk (by using links) from a to b.

bool IsSuperConnected( assembly &A, int super_id );

// A contig is flimsy iff it does not have at least two paired reads with
// opposite orientation on contig.

bool IsContigFlimsy( assembly &A, int contig_id );

// A supercontig is degenerate if:
//  1) its gapped length is < sc_info_min_length (defined in Assemble.cc); and
//  2) all its contigs are flimsy.

bool IsSuperDegenerate( assembly &A, int super_id, int sc_info_min_length );

void MergeMtigsHead( const basevector& b1, const basevector& b2,
     const qualvector& q1, const qualvector& q2, const align& al,
     basevector& c, qualvector& q, vec<int>& to1, vec<int>& to2 );

#endif
