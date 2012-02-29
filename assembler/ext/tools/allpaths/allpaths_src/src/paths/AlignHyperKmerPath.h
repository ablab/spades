/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: AlignHyperKmerPath.h

   Stuff relating to the alignment of a HyperKmerPath to the reference.

   Ideally, each connected <component> of a HyperKmerPath would <capture>
   exactly one contiguous chunk of the reference in a <trusted path> through the
   component.  The trusted path would start at a source, end at a sink,
   follow the graph structure of the component, and each edge of the path
   would have a perfect alignment to the corresponding part of the reference.

   In other words, ideally, each connected component would rougly correspond
   to a "contig" in traditional (linear) assembly: it would represent one
   contiguous piece of the reference.

   Code in this module measures how often this is the case in practice,
   identifies components where this is not the case, and classifies
   violations into common classes so we can understand their source and
   if possible fix them.

   Some common missassembly scenarios that can be seen once we find all
   the trusted paths captured by a component:

       - a component captures two chunks of reference from two different,
       far-away places.  that is a misjoin.

*/

#ifndef ALIGN_HYPER_KMER_PATH_H
#define ALIGN_HYPER_KMER_PATH_H

#include "CoreTools.h"
#include "BinaryIO.h"
#include "lookup/LookAlign.h"
#include "graph/DotUtils.h"
#include "paths/HyperKmerPath.h"
#include "paths/SeqOnHyper.h"

#include <map>

class TrustedPath; // forward declaration

// FuncDecl: AlignHyperKmerPath
//
// AlignHyperKmerPath takes a HyperKmerPath h, whose KmerPath edges are assumed
// not to have any gaps, and aligns each of them to a reference genome (using files 
// GENOME.fastb and GENOME.lookup), then removes some edges alignments which appear 
// to be wrong.  The result is returned as a vec<look_align>, in which the query 
// ids are the edge ids in h.  An index is also generated.  A work directory must 
// be provided for intermediate calculations.
//
// We use two programs to do this, QueryLookupTable and FindAllPerfectAlignments.  
// The first finds alignments, whether perfect or not, but will miss multiple 
// alignments of the same sequence to overlapping places on the reference.  The 
// second is guaranteed to find all perfect alignments.
//
// Output parameters:
//
//    aligns - each element gives the alignment of one edge of the HyperKmerPath h
//             to the reference.
//    aligns_index - for each edge, all alignments of this edge to the reference
//       (the indices in 'aligns' of the alignments of this edge).

void AlignHyperKmerPath( const HyperKmerPath& h, const KmerBaseBroker* kbb,
     const String& GENOME, const String& tmp_dir, vec<look_align>& aligns,
     vec< vec<int> >& aligns_index, bool filter = true );

// FuncDecl: ReorderToFollowReference
//
// Given edge alignments, try to reorder the <components>
// of a HyperKmerPath so that they follow the reference.  Flip them (both components 
// and aligns) if necessary.
void ReorderToFollowReference( HyperKmerPath& h, vec<look_align>& aligns,
     const vec< vec<int> >& aligns_index );

// FuncDecl: PrintAlignedHyperKmerPath
//
// Given edge alignments, print a HyperKmerPath, by
// components, with the edge alignments displayed.
void PrintAlignedHyperKmerPath( ostream& out, const HyperKmerPath& h, 
     const KmerBaseBroker* kbb, const vecbasevector& genome, 
     const vec<look_align>& aligns, const vec< vec<int> >& aligns_index,
     Bool print_component_id_line = True, 
     const vec<TrustedPath>* trusted_pathsp = 0, const Bool brief = False,
     const Bool diploid = False );

// FuncDecl: AlignAndPrintHyperKmerPath
//
// AlignAndPrintHyperKmerPath: one-stop shopping for aligning and printing of a
// HyperKmerPath.
void AlignAndPrintHyperKmerPath( ostream& out, const HyperKmerPath& h, 
     const KmerBaseBroker* kbb, const String& GENOME, const String& tmp_dir, 
     Bool print_component_id_line = True, bool filter = true );

/**
   Class: TrustedPath
   
   A TrustedPath is a path through a <contig>, that corresponds to a
   segment of the reference.  Alternately, we can say that a TrustedPath
   is a segment of the reference that can be "threaded" through
   a path on a contig.

   In general, a contig of the final assembly represents a set of possible
   basevectors (in a factored way).  Trusted paths measure which of these
   basevectors are actually genomic.   Ideally, they're all genomic and all
   come from the same genomic neighborhood.  Computing trusted paths helps
   find deviations from that ideal reality.

   A TrustedPath is a represented as a vector of alignments of
   HyperKmerPath edges to the reference.  The TrustedPath may start in the
   middle of an edge and/or end in the middle of an edge.  The TrustedPath
   follows the graph edges according to the graph structure of the contig.

   Note that the TrustedPath through the contig represents a particular
   basevector; that basevector may occur on either the forward version
   of a genome part (the version represented in the genome.fast file)
   or on the reverse-complement version of the genome part.   In the former
   case, all look_aligns of edges in the path have rc1=False, in the latter
   they all have rc1=True.
   
   Other names for this are "unwound path" or "captured path".
*/
class TrustedPath {
 public:
  TrustedPath() {}

  TrustedPath( int contig,
               int contigTotalEdgeLength,
               const vec<int>& vertexIds,
               const vec<look_align>& aligns );

  int GetContig() const { return m_contig; }

  int GetNumAligns() const { return m_aligns.size(); }

  Float GetFractionEdgeLengthAligned() const { return m_fractionEdgeLengthAligned; }

  int GetVertexIdBefore( int i ) const { return m_vertexIds[i]; }
  int GetVertexIdAfter( int i ) const { return m_vertexIds[i+1]; }
  int GetEdgeId( int i ) const { return m_aligns[i].query_id; }

  int GetFirstVertexId() const { return m_vertexIds.front(); }
  int GetLastVertexId() const { return m_vertexIds.back(); }
  const vec<int>& GetVertexIds() const { return m_vertexIds; }

  genome_part_id_t GetFinishedId() const { return m_aligns.front().target_id; }

  // Where on the forward version of the finished sequence does the trusted path align?
  genome_part_pos_t Begin() const { return m_aligns.front().a.pos2(); }
  genome_part_pos_t End() const { return m_aligns.back().a.Pos2(); }
  int Length() const { return this->End() - this->Begin(); }

  // Returns 0 if the trusted path aligns fw on the reference, 1 if rc.
  int GetRc() const { return m_aligns.front().rc1; }
  bool IsFw() const { return GetRc() == 0; }
  bool IsRc() const { return GetRc() == 1; }

  DEFINE_BINARY_IO_5( TrustedPath, m_contig, m_vertexIds, m_aligns,
		      m_fractionEdgeLengthAligned, m_uniquePerfectEdgeIds );
  
  const look_align& 
  GetAlign( int i ) const { return m_aligns[i]; }

  const vec<look_align>& 
  GetAllAligns() const { return m_aligns; }

  bool DominatesAlignsOf( const TrustedPath& other ) const;

  bool DominatesEdgesOf( const TrustedPath& other ) const;

  bool DominatesVerticesOf( const TrustedPath& other ) const;

  void PrintSummary( ostream& out ) const;
  
  const vec<int>&
    GetUniquePerfectEdgeIds() const;

  void TestValid( ) const;
  void TestValid( const HyperKmerPath& ) const;

 private:
  // Private field: m_contig
  // The id of the connected component within the HyperKmerPath,
  // through which component this TrustedPath passes.
  int m_contig;

  // Private field: m_vertexIds
  // Vertices of this TrustedPath, in order along the path.
  vec<int> m_vertexIds;

  // Private field: m_aligns
  // Alignments of the HyperKmerPath edges in this TrustedPath
  // to the reference.  Note that each edge may have other alignments
  // to the reference, but only one alignment for each edge is used
  // here.
  // The alignments are in order along the path.  In each alignment,
  // the query_id field is the edge id of the aligned HyperKmerPath
  // edge, and the target_id is the <genome part id> of the genome part
  // to which this path aligns.  
  vec<look_align> m_aligns;
  Float m_fractionEdgeLengthAligned;


  mutable vec<int> m_uniquePerfectEdgeIds;
};  // class TrustedPath

ostream& operator<< ( ostream& out, const TrustedPath& path );

DEFINE_SWAP(TrustedPath);

inline bool operator< ( const TrustedPath& lhs, const TrustedPath& rhs ) {
  if ( lhs.GetFinishedId() < rhs.GetFinishedId() ) return true;
  if ( lhs.GetFinishedId() > rhs.GetFinishedId() ) return false;
  if ( lhs.Begin() < rhs.Begin() ) return true;
  if ( lhs.Begin() > rhs.Begin() ) return false;
  if ( lhs.End() < rhs.End() ) return true;
  else return false;
}

// Semantic type: tpid_t
// Identifies an TrustedPath in a vec<TrustedPath>
SemanticTypeStd( int, tpid_t );

// Semantic type: tpleg_t
// Index of a trusted path <leg> within a trusted path.
// Each time a trusted path passes over an edge of the contig,
// that is one leg of the trusted path.
SemanticTypeStd( int, tpleg_t );

// FuncDecl: TrustedPathsToIndexedAligns
// Convert a vec of TrustedPaths into an indexed vec of look_aligns.
void TrustedPathsToIndexedAligns( const vec<TrustedPath>& paths,
                                  const int numEdges,
                                  vec<look_align>& aligns,
                                  vec< vec<int> >& aligns_index,
                                  const bool removeImproper = true );


// FuncDecl: FilterByReference
//
// Find edge alignments and filter them by finding
// the alignments most consistent with the graph structure and the
// finished sequence.  These alignments are then packaged into <trusted paths>.
void FilterByReference( const HyperKmerPath& theGraph, 
                        const int K,
                        const vec<look_align>& aligns,
                        const vec< vec<int> >& aligns_index,
                        vec<TrustedPath>& trustedPaths );

// Remove paths that cover finished sequence that is covered by some
// other, longer path from that contig.

void FilterPathsByAlignDominance( vec<TrustedPath>& trustedPaths );

// Remove paths that are edge-dominated by some other path from that contig. 

void FilterPathsByEdgeDominance( vec<TrustedPath>& trustedPaths, int numEdges );

// Remove paths that are vertex-dominated by some other path from that contig.

void FilterPathsByVertexDominance( vec<TrustedPath>& trustedPaths );

// Remove paths that involve less than some percentage of the contig's
// total edge length.

void FilterPathsByEdgeCoverage( vec<TrustedPath>& trustedPaths, 
                                const Float minFractionCovered = 0.01 );

// Remove paths that involve less than some percentage of the contig's
// total edge length and are shorter than some cutoff.

void FilterPathsByLength( vec<TrustedPath>& trustedPaths, 
                          const int minLength = 100,
                          const int minPercent = 10 );

// ReportMisassemblies.  Look for putative misassemblies:
//
// 1. Report edges that have no end-to-end alignment (or that have no alignment 
// at all).
//
// 2. For each component C, consider the parts of the genome that are covered by it.
// Divide these parts into their connected components, and consider those components
// that include a uniquely anchored edge of length >= MIN_LEN.  If the separation
// between two of these components is >= MIN_SEP, report C.  A given component is
// reported at most once.

void ReportMisassemblies( ostream& out, const HyperKmerPath& h,
     const vec<look_align>& aligns, const vec< vec<int> >& aligns_index,
     const int MIN_LEN = 3000, const int MIN_SEP = 2000, 
     const int MIN_LEN_NO_REPORT = 1000 );



  
// Term: component
//
// One connected component of a HyperKmerPath.

#endif
