/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: PseudoUnipatherAnnex.h

   Code for PseudoUnipather, moved into a separate file so it can be reused by
   other programs.

   Generally, this code has to do with the problem of scaffolding a set of contigs,
   and then filling in the gaps between some of the contigs in a scaffold
   to increase contig N50.

   @file
*/

#ifndef __INCLUDE_paths_PseudoUnipatherAnnex_h
#define __INCLUDE_paths_PseudoUnipatherAnnex_h

#include "Vec.h"
#include "Intvector.h"
#include "lookup/LookAlign.h"
#include "paths/KmerPath.h"
#include "paths/PseudoUnipathGapCloser.h"
#include "CommonSemanticTypes.h"

// Logical type: line_unipaths_t
// A <line> is represented as a list of <unipath ids> of unipaths
// comprising the line.  Often comes together with line_closers_t,
// which represents the sequence in the gap between each adjacent
// pair of unipaths in the line.
typedef vec< unipath_id_t > line_unipaths_t;

// Type: line_gaps_t
// For one <line>, represents the <gaps closers> between adjacent unipaths in that line.
typedef vec< gapcloser > line_closers_t;

// FuncDecl: GetLines
// 
// Get <lines> in the graph implied by <closers>.  Aborts if it finds a cycle.
//
// Input parameters:
//
//     nuni - total number of unipaths
//     closers - the <gap closers> we have found
//
// Output parameters:
//
//     lines - the lines
//     gaps - the gap closers for each line in 'lines'
void GetLines( const int nuni, const vec<gapcloser>& closers, 
	       vec< line_unipaths_t >& lines, vec< line_closers_t >& gaps );


// FuncDecl: ComputeLineSize
// Compute the size of a <line>
nbases_t ComputeLineSize( const vecKmerPath& unipaths, const line_unipaths_t& line,
			  const line_closers_t& closers );


/**
   Local class: Link

   Represents a linking of two <copy-number-one unipaths> by <paired reads>.  A link
   establishes the distance between two unipaths and their relative orientation; it can
   tell us that the two unipaths (or the lines of which they're currently part) lie
   next to each other on the genome, and we may be able to join them.

   We initially create a number of links between each pair of unipaths, one link for each
   paired read with one read of the pair landing on one unipath and the other read on the
   other unipath; but then we <merge> the links, so that between any given pair of
   unipaths there is at most one link.  These merged links are represented as instances of
   <ink>.
*/
class Link {

     public:

  // Fields: Fields of Link
  //    tig1, tig2 - <unipath ids> of the two unipaths linked by this paired read
  //    sep, dev - separation predicted for these two unipaths by aligning them to this
  //       paired read, and the standard deviation of the separation (since we know the
  //       distance between the two reads in a paired read only approximately)
  //    start1, stop1 - where on the first unipath does the first read of the pair align?
  //    start2, stop2 - where on the second unipath does the second read of the pair align?
  
     basevec_id_t tig1;
     basevec_pos_t start1, stop1;
     nbases_dbl_t sep, dev;
     basevec_id_t tig2;
     basevec_pos_t start2, stop2;

     Link( ) { }

     Link( const basevec_id_t tig1, const basevec_pos_t start1, const basevec_pos_t stop1, const nbases_dbl_t sep,
          const nbases_dbl_t dev, const basevec_id_t tig2, const basevec_pos_t start2, const basevec_pos_t stop2 )
          : tig1(tig1), start1(start1), stop1(stop1), sep(sep), dev(dev),
          tig2(tig2), start2(start2), stop2(stop2) { }

  // Function: operator-less-than
  // Orders the links by the ids of the unipaths they're linking.
  // When an array of <Links> is sorted, multiple links between a given
  // unipath pair will be adjacent in the array, so we can
  // conveniently find them all and merge them into an <ink>.
     friend Bool operator<( const Link& l1, const Link& l2 )
     {    if ( l1.tig1 < l2.tig1 ) return True;
          if ( l1.tig1 > l2.tig1 ) return False;
          if ( l1.tig2 < l2.tig2 ) return True;
          if ( l1.tig2 > l2.tig2 ) return False;
          if ( l1.dev < l2.dev ) return True;
          return False;    }

};  // class Link

// Semantic type: Link_id_t
// The id of a <Link> in a vector of Links.
SemanticTypeStd( int, Link_id_t );
typedef VecIntVec vecLinkIDVec;

void AddLinkBetweenBasevecs( nbases_t n1, nbases_t n2,
			     const look_align& la1, const look_align& la2,
			     nbases_dbl_t sep_mean, nbases_dbl_t sep_dev,
			     vec< Link >& Links );

/**
   FuncDecl: FindLinksBetweenBasevecs

   Given a set of basevecs, and a set of read pairs aligned to the basevecs,
   for each contig pair finds the set of links between these two contigs.
   A Link is a pair of reads where one read lands on one basevec and
   the other read on the other basevec.

   Input params:

      aligns_fw, aligns_bw - all alignments of each fw (bw) read of
         each read pair, to the unipaths.  To find the alignments
         of the fw or bw read of a given read pair, use
         aligns_ind_fw / aligns_ind_rc (defined below).

      aligns_ind_fw, aligns_ind_rc - for the fw and bw read of each
         read pair, indices of aligns of that read to the basevecs
         in aligns_fw / aligns_rc.
       
*/
void FindLinksBetweenBasevecs( // Info about the basevecs:
			       const vec<nbases_t>& basevecSizes,

			       // Info about the pairs:
			       const vec< nbases_dbl_t >& sep_means,
			       const vec< nbases_dbl_t >& sep_devs,

			       // Info about alignment of pairs
			       // to the basevecs:
			       const vec<look_align>& aligns_fw,
			       const vec<look_align>& aligns_rc,
			       const vec< vec<align_id_t> >& aligns_ind_fw,
			       const vec< vec<align_id_t> >& aligns_ind_rc,
			       
			       // Output
			       vec< Link >& Links,

			       // For each basevec, the ids in 'links'
			       // of the links between it and some other
			       // basevec.
			       vecLinkIDVec& Links_index1,
			       vecLinkIDVec& Links_index2
			      
			       );


// FuncDecl: CombineStats
//
// Given a vector of means and stddevs, combine them into one value.
// Must <SortSync>(dev,sep) before calling.
void CombineStats( const vec<nbases_dbl_t>& sep, const vec<nbases_dbl_t>& dev,
		   nbases_dbl_t& Sep, nbases_dbl_t& Dev );
     
/**
   FuncDecl: GetLinkStats

   Given a set of <Links> between a pair of contigs, determine the separation
   between these contigs implied by these links.
*/
void GetLinkStats( const vec<Link>& Links, nbases_dbl_t& Sep, nbases_dbl_t& Dev );



#endif
// #ifndef  __INCLUDE_paths_PseudoUnipatherAnnex_h
