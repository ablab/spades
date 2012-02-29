/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_MUXWALKGRAPH_H
#define PATHS_MUXWALKGRAPH_H

#include "paths/OrientedKmerPathId.h"
#include "paths/Mux.h"
#include "paths/MuxToPath.h"
#include "paths/HyperKmerPath.h"
#include "paths/SubsumptionList.h"

#include <map>
#include <set>

#if __GNUC__ > 2
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#else
#include <hash_map>
#endif

// Forward declaration of friend classes which convert MWGs to other things
class MuxWalkGraph_to_vecKmerPath;
class MuxWalkGraph_to_HyperKmerPath;
class MuxWalkGraph_to_HyperKmerPath2;


/// A MuxWalkGraph encodes a graph-like structure describing all ways
/// to walk read-to-read across an insert using minimal extensions.
/// Walks which separate and then reconverge (eg to go around a SNP)
/// are represented efficiently.  (But at the moment, if you go around
/// an indel, everything before the indel will be stored twice, with
/// two different distance_to_end values.)
///
/// For an overview of Mux-based insert walking, see paths/doc/mux.pdf

// A Node in the graph is essentially a pair<okpi,dist_to_end>, which
// represents the datum that a particular oriented kmer path appears
// in a closure in which there are dist_to_end kmers after that read.
// There is a map taking each Node to a vector of all of its successors,
// = all the possible next nodes on its various closures, suffix-array-like.
// The object also knows which nodes are the beginnings of closures, though
// this is redundant as the caller could search for the OKPI of the openers.

// The encoded walks are only guaranteed legal extension by extension:
// if it tells you to walk A->B->C, that does not guarantee that A and
// C align, only that they align on the part of each that overlaps
// with B.  So converting this back to the exact list of legal insert
// walks will still take work.  (See MuxToPath.)



class MuxWalkGraph {
  
public:

  /// Each Node (basically pair<okpi,dist_to_end>) in the graph has a Label.
  //
  // Implementation details:
  //   We give it a label because if we just called it pair<okpi,int>,
  //   we'd need to do a binary search each time we wanted to find it.
  //   We label with ints, and maintain a vec of iterators that the
  //   int indexes into.  This costs memory -- both to hold the vec,
  //   and because each node also needs to remember its own label --
  //   but we don't expect huge MuxWalkGraphs.  We could use an iterator
  //   as the label, but the data structure would be harder to define
  //   and harder to save to disk.

  typedef int Label;
  static const Label LabelDone;

private:

  struct Node {
    OrientedKmerPathId okpi;
    int dist_to_end;
    Bool is_good;   
    // We might want to keep "bad" reads in the graph but not use them later.

    Node(OrientedKmerPathId OKPI, int dist, Bool good=True) 
      : okpi(OKPI), dist_to_end(dist), is_good(good) { }

    friend bool operator<(const Node& lhs, const Node& rhs ) {
      return( lhs.okpi < rhs.okpi || 
	      (lhs.okpi == rhs.okpi && 
	       (lhs.dist_to_end < rhs.dist_to_end ||
		(lhs.dist_to_end == rhs.dist_to_end &&
		 (lhs.is_good < rhs.is_good)))) );
    }

    friend bool operator==(const Node& lhs, const Node& rhs ) {
      return( lhs.okpi == rhs.okpi && 
	      lhs.dist_to_end == rhs.dist_to_end &&
	      lhs.is_good == rhs.is_good );
    }

    size_t hash_val() const {
      return( okpi.GetIdRc() ^ 
	      ( size_t(2*dist_to_end+is_good) << 4*sizeof(size_t) ) );
    }
  };

  struct node_hash {
    size_t operator() ( const Node& node ) const { return node.hash_val(); }
  };


  struct Indices {
    Label self;
    set<Label> successors;

    Indices() : self(-1) { }
    bool Uninitialized() { return( self == -1 ); }
  };

  //  typedef map<Node,Indices> node_map_t;
  typedef hash_map<Node,Indices,node_hash> node_map_t;
  typedef node_map_t::iterator map_iter;
  typedef node_map_t::const_iterator map_const_iter;
  
  typedef hash_map< OrientedKmerPathId,vec<Label> > okpi_map_t;

public:

  longlong size() const { return m_label_lookup.size(); }  // constant time

  /// To get KmerPaths out of a MuxWalkGraph.
  bool AllPaths( const vecKmerPath* pathsFw, 
		 const vecKmerPath* pathsRc,
		 vecKmerPath& ans,
		 unsigned int max_ans = 0 ) const;

  friend class MuxWalkGraph_to_vecKmerPath;
  // As implemented, that class does not need to be a friend;
  // all MWG data is retrieved using public accessors.  But
  // that class makes use of the face that a MuxWalkGraph::Label
  // is LabelDone or an int from 0 to size(), and is therefore
  // suitable for indexing into a vector.  We make it a friend in
  // recognition of its dependence on an internal implementation detail.


public:
 

  /// Make a HyperKmerPath represenation of the graph structure in the
  /// MuxWalkGraph.  Conceptually there are two separate changes:
  /// 1. Translation from mux space to path space -- here we
  ///    throw out closures that can't actually arise from paths.
  /// 2. Taking the dual of the graph, sort of.  In the MuxWalkGraph,
  ///    each vertex holds a read, while in a HyperKmerPath, sequence
  ///    lives on the edges.
  /// We do this by building a HKP with a vertex corresponding to the 
  /// left endpoint of each read, and an edge from there to each
  /// read that can appear to the right of it in a path-space closure.
  /// (So directed edges in the HKP go the opposite direction from those
  /// in the MWG; this takes care of the walking-direction convention.)

  void MakeHyperKmerPath( const vecKmerPath* pathsFw, 
			  const vecKmerPath* pathsRc,
			  const SubsumptionList* subList,
			  HyperKmerPath& ans ) const;
  friend class MuxWalkGraph_to_HyperKmerPath;
  friend class MuxWalkGraph_to_HyperKmerPath2;
  // As implemented, that class does not need to be a friend;
  // all MWG data is retrieved using public accessors.  But
  // that class makes use of the face that a MuxWalkGraph::Label
  // is LabelDone or an int from 0 to size(), and is therefore
  // suitable for indexing into a vector.  We make it a friend in
  // recognition of its dependence on an internal implementation detail.


public:
  // Monte Carlo estimation of the number of cloures (in Mux space)
  int EstimateNumClosures() const;
private:
  int EstimateNumClosuresTrial() const;


public:
  // Pretty much everything else is for graph building.

//  Do not clear() a hash_map!
//  The result is a map with zero objects but possibly many buckets,
//  which trashes future performance.  MuxWalkGraphs are single-use.

//   void clear() { 
//     m_successor_map.clear(); 
//     m_label_lookup.clear();
//     m_openers.clear();
//     m_okpi_lookup.clear();
//   }

  OrientedKmerPathId OKPI( const Label& label ) const {
    return( m_label_lookup[label]->first.okpi );
  }
  
  int DistToEnd( const Label& label ) const {
    return( label == LabelDone ? 0 : m_label_lookup[label]->first.dist_to_end );
  }

  Bool IsGood( const Label& label ) const {
    return( m_label_lookup[label]->first.is_good );
  }

  const set<Label>& Successors( const Label& label ) const {
    return( m_label_lookup[label]->second.successors );
  }
  
  // Find all nodes with this OKPI
  void FindNodes( const OrientedKmerPathId& okpi, vec<Label>& ans ) const;

  // If node is already known, learn this successor.
  // If node is unknown, learn it, with one successor.
  // Return the label of the (maybe new) node acted upon.
  Label AddSuccessor( const Node& node, const Label& successor_label );
  Label AddSuccessor( const OrientedKmerPathId& okpi,
		      const int& dist_to_end,
		      const Label& successor_label,
		      Bool is_good = True ) {
    return AddSuccessor( Node(okpi,dist_to_end,is_good), successor_label );
  }

  // Add in an entire closure.  
  // An offset of 0 corresponds to the end of the original untrimmed read.
  // Positive offsets arise when a fill has been trimmed.
  // Negative offsets can arise when a closing read is subsumed by
  //   another read which extends past the end of the insert, but which
  //   we need to record here anyway because the closer itself doesn't appear
  //   in the path.
  void AddClosure( const vec<Mux>& closure, int offset, 
		   const vec<Bool>* is_good = NULL ) {
    MergeClosure( closure, 
		  AddClosingRead( closure.back().GetPathId(), offset,
				  (is_good == NULL || is_good->back()) ),
		  is_good );
  }

  // This path doesn't reach an end read, it merges into an existing closure.
  // closure.back() must have the same OKPI as the node indicated by the Label.
  void MergeClosure( const vec<Mux>& closure, Label successor_label, 
		     const vec<Bool>* is_good = NULL );

  // Add a closer, a read with successor LabelDone.
  // dist_to_end 0 corresponds to the end of the untrimmed read;
  // if the closer has been trimmed, we must be told its relative offset.
  // This means the graph can correctly store paths ending at different
  // fillings (with potentially different trims) of a single read.
  Label AddClosingRead( const OrientedKmerPathId& okpi, 
			int offset=0, Bool is_good = True ) {
    return AddSuccessor( okpi, offset, LabelDone, is_good );
  }

  // Mark some read as an opener -- the beginning of some closure.
  // Probably this just means it's a filling of the opening read of
  // the insert we're searching at an appropriate dist_to_end.
  void MarkOpener( const Label& opener ) {
    m_openers.insert( opener );
  }

  const set<Label>& GetOpeners() const { return m_openers; }


  // For debugging purposes, a little info about this MuxWalkGraph:
  void Summary(ostream& out) {
    out << "MuxWalkGraph: " << size() << " nodes in "
	<< m_successor_map.bucket_count() << " buckets, "
	<< m_okpi_lookup.size() << " OKPIs in "
	<< m_okpi_lookup.bucket_count() << " buckets, " 
	<< m_openers.size() << " openers"
	<< endl;
  }


  // DOT-language output of the graph structure.
  void Dot(ostream& out) const;

private:

  node_map_t m_successor_map;
  vec<map_iter> m_label_lookup;
  set<Label> m_openers;

  okpi_map_t m_okpi_lookup;

};





/// A friend class of MuxWalkGraph which implements the
/// conversion to the vecKmerPath of all closures.

class MuxWalkGraph_to_vecKmerPath {
public:

  MuxWalkGraph_to_vecKmerPath( const MuxWalkGraph& mwg,
			       const vecKmerPath* pathsFw,
			       const vecKmerPath* pathsRc ) 
    : mwg(mwg), pathsFw(pathsFw), pathsRc(pathsRc) { }

  const MuxWalkGraph& mwg;
  const vecKmerPath* pathsFw;
  const vecKmerPath* pathsRc;

  typedef MuxWalkGraph::Label Label;
  static const Label LabelDone;

  bool AllPaths( vecKmerPath& ans,
		 unsigned int max_ans = 0 ) const;

  // Things used by AllPaths:
  static const int DTE_UNKNOWN, DTE_IMPOSSIBLE;

  void CalcMaxDTE( const Label& label ) const;
  mutable vec<int> max_dte;

  void CalcInDegrees( ) const;
  mutable vec<int> in_degree;
  
  mutable vec< set<ulonglong> > exited_in_state;

  void AllPathsHelper( KmerPath path_so_far,
		       int dte_so_far,
		       Label label,
		       bool label_was_added,
		       const MuxToPath& pathifier,
		       vecKmerPath& ans,
		       unsigned int max_ans ) const;

};





class MuxWalkGraph_to_HyperKmerPath {
public:

  MuxWalkGraph_to_HyperKmerPath( const MuxWalkGraph& mwg,
				 const vecKmerPath* pathsFw,
				 const vecKmerPath* pathsRc ) 
    : mwg(mwg), pathsFw(pathsFw), pathsRc(pathsRc) { }

  const MuxWalkGraph& mwg;
  const vecKmerPath* pathsFw;
  const vecKmerPath* pathsRc;

  typedef MuxWalkGraph::Label Label;
  static const Label LabelDone;

  void MakeHyperKmerPath( HyperKmerPath& ans );

  enum VX_STATUS { VX_UNKNOWN, VX_BAD, VX_DONE };
  vec<VX_STATUS> vertex_status;

  enum COMPAT_STATUS { CP_UNKNOWN, CP_INCOMP, CP_GOOD, CP_DOM };
  typedef map<Label,COMPAT_STATUS> compat_map;
  typedef compat_map::iterator compat_iter;

  void CalcPredecessors( );
  vec< vec<Label> > predecessors;
  vec<Label> closers;


  void AddOpeners( HyperKmerPath& ans ) const;

  void AddEdgesFrom( Label label, HyperKmerPath& ans );

  compat_iter CheckCompatibility( Label from_label, 
				  Label to_label, 
				  const KmerPath& the_path,
				  int max_dte,
				  compat_map& compat,
				  const HyperKmerPath& ans ) const;

  void AddEdgeFromTo( Label from_vx, 
		      Label to_vx, 
		      const KmerPath& the_path, 
		      HyperKmerPath& ans ) const;

  void RemoveUnreachable( HyperKmerPath& ans ) const;

  void AssociateClosers( HyperKmerPath& ans ) const;

};


class MuxWalkGraph_to_HyperKmerPath2 {
public:

  MuxWalkGraph_to_HyperKmerPath2( const MuxWalkGraph& mwg,
				  const vecKmerPath* pathsFw,
				  const vecKmerPath* pathsRc,
				  const SubsumptionList* subList ) 
    : mwg(mwg), pathsFw(pathsFw), pathsRc(pathsRc), subList(subList) 
  { Setup(); }
  
  const MuxWalkGraph& mwg;
  const vecKmerPath* pathsFw;
  const vecKmerPath* pathsRc;
  const SubsumptionList* subList;

  typedef MuxWalkGraph::Label Label;
  static const Label LabelDone;

  // helpful utilities:
  OrientedKmerPathId OKPI(Label label) const { return okpids[label]; }
  const KmerPath& Read(Label label) const 
  { return *(OKPI(label).GetPathPtr(*pathsFw,*pathsRc)); }
  int ReadLength(Label label) const { return Read(label).MinLength(); }
  int LeftEndDTE(Label label) const { return left_end_dte[label]; }
  int RightEndDTE(Label label) const { return right_end_dte[label]; }
  int NumLabels() const { return okpids.size(); }
  bool isCloser(Label label) const { return is_closer[label]; }

  vec<int> left_end_dte;
  vec<int> right_end_dte;
  vec<OrientedKmerPathId> okpids;
  vec<Bool> is_closer;

  void Setup() {
    left_end_dte.resize( mwg.size() );
    right_end_dte.resize( mwg.size() );
    okpids.resize( mwg.size() );
    for( Label v=0; v<mwg.size(); v++ ) {
      okpids[v] = mwg.OKPI(v);
      left_end_dte[v] = mwg.DistToEnd(v);
      right_end_dte[v] = left_end_dte[v] + ReadLength(v) - 1;
    }
  }


  // Find all closers *and* reads which subsume them.
  // We also need to worry about reads which subsume closers but
  // do not appear in the MWG at all!  These may provide justification
  // for a closure which was found in mux-space but which appears
  // illegal in path-space (according to reads in the MWG):
  //
  //            opener:                 ----------
  //        short read:               ----
  //  disagreeing read:          -----------xxxx
  //            closer:     -------
  //    subsuming read: ----------------------  (not in MWG)
  void FindClosers( );
  vec<Label> closers;
  typedef map<pair<OrientedKmerPathId,int>,Label> sub_map_t;
  sub_map_t reads_which_subsume_closers;
  map< Label,vec<Label> > subsumed_closers;

  void FindSubsumedOpeners( );
  vec< vec<Label> > openers_subsumed_by;
  // If label i occurs in openers_subsumed_by[j], then i is an opener
  // and read j subsumes it.  Each vec is sorted by right DTE of i.


  // We need to do things in topologically sorted order -- ie do all
  // the MWG-successors of v before v.  (MWGs are acyclic, since we
  // put such an order on coterminal reads.)  Let's figure out such an
  // order before we begin.
  // Perhaps MuxWalkGraph itself should provide this function?
  vec<Label> topo_order;
  void ReverseTopoSort();



  // For each Label, figure out which other Labels' reads are 
  // MWG-descendants which overlap and match perfectly.
  vec< vec<Label> > matching_descendants;
  void FindMatchingDescendants();


  void MakeHyperKmerPath( HyperKmerPath& ans );
  int vx_left_end_of_insert;

  // label_to_vx[i] contains j if and only if:
  // * The read with Label i appears in the HKP, ending with
  //   an edge into vertex j; AND
  // * there is a globally consistent path of reads from a 
  //   closer up to and including read i which doesn't extend
  //   past the right end of read i.  As a corollary, any read
  //   j such that left(i) <= left(j) <= right(i) <= right(j)
  //   also has a globally consistent path to a closer.
  vec< set<int> > label_to_vxs;
  // The no_guarantee version removes the second point above.
  // This is only used for subsumed openers, where we want to
  // record that the read exists even if we can't promise
  // any transitive global consistency to the right.
  vec< set<int> > label_to_vxs_no_guarantee;


  void ProcessLabel( Label label, HyperKmerPath& ans );

  // Add a path to ans.
  void AttachRead( Label label, KmerPathLoc loc, int attach_vx,
		   const int right_end_vx, HyperKmerPath& ans,
		   set<int>* vxs_seen=NULL );

  void AddSubsumedOpenerVertices( int from_vx, int to_vx,
				  int edge_obj_id,
				  const Label label,
				  HyperKmerPath& ans );

  int AssociateOpeners( HyperKmerPath& ans );

  void RemoveUnreachable( int right_vx, HyperKmerPath& ans ) const;
  
};


#endif
