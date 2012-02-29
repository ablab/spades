// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef LOADER_H
#define LOADER_H

/*
 * Collection of functions to load commonly used data.
 */
#include "AnAssemblyClass.h"
#include "AnnotatedContig.h"
#include "Basevector.h"
#include "ContigPair.h"
#include "ReadLocation.h"
#include "TrimsFile.h"
#include "Vec.h"
#include "VecInserts.h"


// Load supercontigs (as a vec of super's, defined in AnAssemblyClass).
void LoadNewSupers( const String &asupers_file,
		    vec<super> &new_supers );

// Load supercontigs (as a vec of super's, defined in AnAssemblyClass).
void LoadNewSupers( vec<annotated_supercontig> &asupers,
		    vec<super> &new_supers );

// Load contig lengths.
void LoadContigLengths( const String &contigs_file,
			vec<int> &lengths );

// Load contig lengths.
void LoadContigLengths( const vecbasevector &contigs,
			vec<int> &lengths );

// Load first location with given contig id (it asserts if locs is not sorted).
void LoadFirstLocs( const vec<read_location> &locs,
		    vec<int> &first_locs );

// Load contigs coverages 
void LoadContigCoverages( const vec<read_location> &locs,
			  const vec<int> &first_locs,
			  vec<float> &coverages );

// Load read locations (it is_transposon != 0, then flip transposons).
void LoadLocs( const String &locs_file,
	       vec<read_location> &locs,
	       const vec<Bool> *is_transposon = 0,
	       bool check_sorted = true );

// Load contig pairs for the given supers.
void LoadContigPairs( const vec_inserts &super_ids,
		      const assembly &the_assembly,
		      vec< contig_pair > &cg_pairs );

#endif
