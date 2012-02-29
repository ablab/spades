// Copyright (c) 2003 Whitehead Institute for Biomedical Research

// General utilities for read locations.

#ifndef READ_LOCATION_UTIL
#define READ_LOCATION_UTIL

#include "ReadLocation.h"
#include "ReadPairing.h"
#include "VecAlignmentPlus.h"

// Count number of reads aligning the given loc.
//  loc_id: the index in locs 
//  n_aligns_found: filled with the number of aligns found
//  alignsIndexRdr: reads aligns.index to check if align exists
//  min_align: minimum LocsOverlap bases to say that two reads align
//  max_skip: used as a bound on the number of LocOverlap checks performed
int AlignsCount( const vec<read_location> &locs,
		 int loc_id,
                 int *n_aligns_found,
                 AlignsIndexReader const& alignsIndexRdr,
		 int min_align = 48,
		 int max_skip = 3 );

// How many reads in each contig.
void ReadsInContig( const vec<read_location> &locs,
		    vec<int> &reads_in_contig );

// How many contigs in the given locs file (does not need to be sorted).
int LocsFileContigsCount( const String &locs_file );

// Overlap and Offset based only on locs.

// Return the overlap between locs at id1 and id2 (based only on locs).
int LocsOverlap( const vec<read_location> &locs, int id1, int id2 );

// Like LocsOverlap, but can return a <0 overlap (a gap). Notice however
//  that if id1, id2 belong to different contigs, then LocsOverlap returns
//  0, meanwhile SignedLocsOverlap will ForceAssert.
int SignedLocsOverlap( const vec<read_location> &locs, int id1, int id2 );

// It will assert if locs at id1 and id2 belong to different contigs.
int LocsOffset( const vec<read_location> &locs, int id1, int id2 );

// Merge/Undo merge original read locations.

//  merge original read locations where appropriate
void MergeLocations( const vec<read_location> &locs, 
		     vec<read_location> &mapping_locs,
		     vec<read_location> &merged_locs );

//  go from merged read locations to original locations ( ie unmerge reads )
void UnmergeLocations( const vec<read_location> &locs, 
		       vec<read_location> &mapping_locs,
		       vec<read_location> &unmerged_locs );

// Find insert stretch ( ( observed_sep - given_sep ) / given_stdev ). It
// returns false if link is illogical, or if locs belong to different
// contigs. Notice: pair must be the read_pairing between reads of locA
// and locB, failing this will result in a ForceAssert.
bool Stretch( const read_location &locA,
	      const read_location &locB,
	      const read_pairing &pair,
	      float &stretch );

// Find all logic pairs in the given contig, and fill pair2stretch with
// pairs ( pair_id, pair_stretch ). Remark: locs must be sorted, but for
// efficiency reasons this is not checked.
//
// locs: sorted vector of locs
// flocs: vector of first_locs
// to_pair: sends read_id to the id of the pair with read_id in it
// pairs: read_pairings
// pair2stretch: output, each entry is a pair ( pair_id, pair_stretch )
void AllStretches( int contig_id,
		   const vec<read_location> &locs,
		   const vec<int> &flocs,
		   const vec<int> &to_pair,
		   const vec<read_pairing> &pairs,
		   vec< pair<int,float> > &pair2stretch,
		   bool append = false );

// As above, with the difference that it finds all logic pairs (ie all
// pairs contained in any one single contig), and that it ForceAsserts
// if locs is not sorted.
void AllStretches( const vec<read_location> &locs,
		   const vec<int> &flocs,
		   const vec<int> &to_pair,
		   const vec<read_pairing> &pairs,
		   vec< pair<int,float> > &pair2stretch );

#endif
