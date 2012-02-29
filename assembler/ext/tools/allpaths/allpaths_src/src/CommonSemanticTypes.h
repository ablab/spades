///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Header: CommonSemanticTypes.h

   Defines some common <semantic types> used throughout the code.
   Putting these definitions into a separate header helps avoid creating some
   spurious dependencies among headers.

   Many of the semantic types seem to be inherently unsigned but should be
   declared as signed integral types, to allow for representing the
   "null" value as -1.
   
   @file
*/

#ifndef __INCLUDE_CommonSemanticTypes_h
#define __INCLUDE_CommonSemanticTypes_h

#include "system/Types.h"
#include "String.h"
#include "SemanticTypes.h"

// TODO: potentially dangerous truncation of indexes
// most of these semantic types defined as ints are now wrong.
// Should be unsigned int or size_t.

// Semantic type: nbases_t
// Represents the length of a DNA region in bases.  Must be signed.
SemanticTypeStd( int, nbases_t );

// Semantic type: nbases_dbl_t
// Represents the length of a DNA region in bases, allowing for
// non-integer lengths.  Useful for representing averages of
// lengths.
SemanticTypeStd( double, nbases_dbl_t );

// Semantic type: nkmers_t
// Represents the length of a DNA region in kmers.  Of course, this
// presupposes there is one fixed kmer size in use throughout.
// Can also represent the count of kmers or kmer occurrences in
// other situations.
SemanticTypeStd( int, nkmers_t );


// Semantic type: basevec_id_t 
// Identifier of a basevector in a set of basevectors.
SemanticTypeStd( int, basevec_id_t );

// Semantic type: basevec_pos_t
// Position in a base vector.
SemanticType( nbases_t, basevec_pos_t );

// Semantic type: read_id_t
// Logical type for a read id -- the index of a read in a <reads_t> array of reads.
// A read can also be denoted by a reads_t::const_iterator.
SemanticType( basevec_id_t, read_id_t );

// Const: NULL_READ_ID
// A value guaranteed to be different from every valid read id.
const read_id_t NULL_READ_ID = -1;


// Semantic type: nreads_t
// Logical type for a count of reads.
typedef read_id_t nreads_t;

// Semantic type: genome_part_id_t
// Logical type for a genome part id -- the index of a <genome part> in a <genome_t> array of genome parts.
// A genome part can also be denoted by a genome_t::const_iterator.
SemanticTypeStd( basevec_id_t, genome_part_id_t );

// Semantic type: read_pos_t
// Logical type for representing position in a read.  Must be signed -- this makes
// coding more convenient.
SemanticType( basevec_pos_t, read_pos_t );

// Semantic type: kmer_pos_t
// Logical type for representing position in a kmer.
typedef int kmer_pos_t;

// Semantic type: genome_part_pos_t
// Logical type for representing position in a genome part.  Must be signed -- this makes
// coding more convenient.
SemanticType( basevec_pos_t, genome_part_pos_t );

// Semantic type: orient_t
// Read orientation, relative to a given strand: forward or reverse complement.
// Use the named constants <ORIENT_FW> and <ORIENT_RC>, rather than the physical
// values False/True, to prevent any confusion.
SemanticTypeStd( Bool, orient_t );

// Constants: orient_t values
//
//    ORIENT_FW - the read or unipath occurs on the <genome part> at the specified
//        position as-is
//    ORIENT_RC - the reverse complement of the read or unipath occurs on the
//        genome part at the specified position.
const orient_t ORIENT_FW = False;
const orient_t ORIENT_RC = True;

// Semantic type: copy_num_t
// Represents the copy number of some DNA string in the genome.
// Can be the copy number of a <kmer> or a <unipath>, for example.
SemanticTypeStd( int, copy_num_t );

// Semantic type: prob_t
// A probability value.
SemanticTypeStd( double, prob_t );

// Semantic type: filename_t
// The name of a file.
SemanticType( String, filename_t );

// Semantic type: filenamepart_t
// Part of a filename.
SemanticType( String, filenamepart_t );

// Semantic type: dirname_t
// The name of a directory.
SemanticType( String, dirname_t );

#endif
// #ifndef __INCLUDE_CommonSemanticTypes_h



