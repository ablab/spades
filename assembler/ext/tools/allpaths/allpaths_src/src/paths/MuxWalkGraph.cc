/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/MuxWalkGraph.h"
#include "paths/HyperKmerPathCleaner.h"
#include "random/Random.h"
#include "graph/DigraphTemplate.h"

const MuxWalkGraph::Label MuxWalkGraph::LabelDone = -2;
// How can the conversion utilities access this easily?
const MuxWalkGraph_to_vecKmerPath::Label MuxWalkGraph_to_vecKmerPath::LabelDone = -2;
const MuxWalkGraph_to_HyperKmerPath::Label MuxWalkGraph_to_HyperKmerPath::LabelDone = -2;
const MuxWalkGraph_to_HyperKmerPath2::Label MuxWalkGraph_to_HyperKmerPath2::LabelDone = -2;

void MuxWalkGraph::FindNodes( const OrientedKmerPathId& okpi, 
			      vec<MuxWalkGraph::Label>& ans ) const {
  ans.clear();

  okpi_map_t::const_iterator it = m_okpi_lookup.find(okpi);

  if( it != m_okpi_lookup.end() )
    ans = it->second;

  return;
}



// If node is already known, learn this successor.
// If node is unknown, learn it, with one successor.
MuxWalkGraph::Label MuxWalkGraph::AddSuccessor( const Node& node,
				  const MuxWalkGraph::Label& successor_label ) 
{
  map_iter mi = (m_successor_map.insert( make_pair(node,Indices()) )).first;
  
  // mi now points to the thing in the map, which may be new or old.
  // If it's new, we need to initialize it.
  if( mi->second.Uninitialized() ) {
    Label label(m_label_lookup.size());
    m_label_lookup.push_back( mi );
    mi->second.self = label;
    m_okpi_lookup[node.okpi].push_back( label );
  }
  
  // Add the requested successor:
  mi->second.successors.insert( successor_label );
  
  // Return the label of the node just acted upon.
  return mi->second.self;
}



void MuxWalkGraph::MergeClosure( const vec<Mux>& closure,
				 MuxWalkGraph::Label successor_label,
				 const vec<Bool>* is_good ) {
  if( closure.empty() ) return;

  ForceAssert( successor_label != LabelDone );

  // Deal with closure.back():
  // Final read in the vec<Mux> is to be associated with successor_label,
  // so it had better have the same OKPI!
  ForceAssert( OKPI(successor_label) == closure.back().GetPathId() );
  // ... and the same goodness, if we're tracking that.
  ForceAssert( is_good == NULL || IsGood(successor_label)==is_good->back() );

  // Now walk through the vec<Mux> backwards and add all the reads.
  // Each time through the loop we need to examine two of the muxes:
  // one to know which read we're adding, and its successor (forwards) to 
  // know how much further away from the end of the closure that read is.
  for( uint i = closure.size()-1; i>0; i-- ) 
    successor_label = AddSuccessor( closure[i-1].GetPathId(),
				    DistToEnd(successor_label) + closure[i].GetNumKmers(),
				    successor_label,
				    ( is_good == NULL || (*is_good)[i-1] ) );
  MarkOpener( successor_label );
}



bool MuxWalkGraph::AllPaths( const vecKmerPath* pathsFw, 
			     const vecKmerPath* pathsRc,
			     vecKmerPath& ans,
			     unsigned int max_ans ) const {

  MuxWalkGraph_to_vecKmerPath mwg2vkp(*this, pathsFw, pathsRc);

  return mwg2vkp.AllPaths( ans, max_ans );

}



// Things for extracting KmerPaths from the graph:

const int MuxWalkGraph_to_vecKmerPath::DTE_UNKNOWN = -1;
const int MuxWalkGraph_to_vecKmerPath::DTE_IMPOSSIBLE = -2;
// Note: DTE_IMPOSSIBLE is negative, so that it's < all possible values

bool MuxWalkGraph_to_vecKmerPath::AllPaths( vecKmerPath& ans,
					    unsigned int max_ans ) const {
  ans.clear();

  // Calculate the max DTE from each node.
  max_dte.clear();
  max_dte.resize( mwg.size(), DTE_UNKNOWN );
  for( Label i=0; i != Label(max_dte.size()); i++ )
    if( max_dte[i] == DTE_UNKNOWN )
      CalcMaxDTE( i );

  // Calculate the in-degree of each node.
  CalcInDegrees();

  // Scratch space for AllPathsHelper to record which nodes have been
  // exited with which state, to avoid doing redundant computation.  
  // Only used for nodes with out-degree>1, but in some cases that could 
  // still bog down a lot, so this may turn out to be a problem.
  exited_in_state.clear();
  exited_in_state.resize( mwg.size() );



  MuxToPath pathifier( pathsFw, pathsRc );

  // Just start up the AllPathsHelper on each (good) opener.
  for( set<Label>::const_iterator opener = mwg.GetOpeners().begin();
       opener != mwg.GetOpeners().end(); opener++ ) {
    if( mwg.IsGood(*opener) )
      AllPathsHelper( *(mwg.OKPI(*opener).GetPathPtr(*pathsFw,*pathsRc)),
		      mwg.DistToEnd(*opener),
		      *opener,
		      true,
		      pathifier,
		      ans,
		      max_ans );

    if( max_ans && ans.size() >= max_ans )
      break;
  }

  // We might generate the same answer twice, from paths
  // differing on bad reads.  So we unique-sort the answers.
  // Down side: we may surprise the user by returning fewer
  // than max_ans answers even if the search was cut off!
  // We indicate an aborted search by returning false.

  bool aborted = ( max_ans && ans.size() >= max_ans );

  cout << " " << ans.size() << " " << "(aborted=" << aborted << ") " << flush;
  
  // We're actually gonna skip this if the construction is aborted, as
  // it's super expensive.
  if ( !aborted )
    ans.UniqueSort();

  return( !aborted );

}


// For each Label, figure out the maximum DTE_so_far that can possibly
// be extended by this node or any of its children.  This info is used
// by AllPathsHelper, to know when to prune a hopeless search.
void MuxWalkGraph_to_vecKmerPath::CalcMaxDTE( const Label& label ) const {
  int& ans = max_dte[label];

  // First, how far back do I reach:
  if( mwg.IsGood(label) ) {
    ans = mwg.DistToEnd(label) + 
      mwg.OKPI(label).GetPathPtr(*pathsFw,*pathsRc)->MinLength() - 1;
                               // MinLength is what is recorded by the mux.
    // I think this number should never be negative -- that would indicate
    // a read that didn't even start until after the closer ends, so even
    // subsumption shouldn't be an issue.  But I'm too chicken to assert.
  }
  else {
    ans = DTE_IMPOSSIBLE;
  }

  // Now, can any of my successors do better?
  const set<Label>& successors = mwg.Successors(label);
  for( set<Label>::const_iterator succ = successors.begin();
       succ != successors.end(); succ++ ) {
    if( *succ == LabelDone ) continue;
    if( max_dte[*succ] == DTE_UNKNOWN )
      CalcMaxDTE(*succ);
    if( ans < max_dte[*succ] )  // DTE_IMPOSSIBLE is less than everything
      ans = max_dte[*succ];
  }
}


// For each Label, calculate the in-degree.
// This ignores which nodes can be reached from the openers, and
// just calculates the in-degree according to the entire graph.
void MuxWalkGraph_to_vecKmerPath::CalcInDegrees( ) const {
  in_degree.clear();
  in_degree.resize( mwg.size(), 0 );
  for( Label i=0; i != Label(mwg.size()); i++ ) {
    const set<int>& succ = mwg.Successors(i);
    for( set<Label>::const_iterator j = succ.begin(); j != succ.end(); j++ )
      if( *j != LabelDone )
	++in_degree[*j];
  }
}



// path_so_far is a path probably ending with the read indicated by label,
//    though it may not reach quite that far, if the read of label didn't
//    align properly -- in which case, label_was_added is false.
// dte_so_far is the dist_to_end of the last read which did get added.
//
// recursively call on all successors.
void MuxWalkGraph_to_vecKmerPath::AllPathsHelper( KmerPath path_so_far,
				   int dte_so_far,
				   Label label,
				   bool label_was_added,
				   const MuxToPath& pathifier,
				   vecKmerPath& ans,
				   unsigned int max_ans ) const {

  if( max_dte[label] == DTE_IMPOSSIBLE || dte_so_far > max_dte[label] )
    // If every end read is marked good, then no node should have
    // DTE_IMPOSSIBLE as its max_dte.  But who knows?
    return;

  if( in_degree[label]>1 
      && ! exited_in_state[label].insert(path_so_far.GetHash(dte_so_far)).second )
    // We've explored from here on before; don't re-do the work
    return;

  KmerPath scratch_path;

  const set<Label>* successors = &(mwg.Successors(label));

  // Iterate until you hit a branch or the end:
  while( // successors.size() == 1 
	 ! successors->empty() &&
	 ++(successors->begin()) == successors->end() &&
	 *(successors->begin()) != LabelDone ) {
    label = *(successors->begin());
    successors = &(mwg.Successors(label));
    if( mwg.IsGood(label) &&
	pathifier.ExtendByKmersIfMatch( mwg.OKPI(label),
					dte_so_far - mwg.DistToEnd(label),
					path_so_far,
					scratch_path ) ) {
      dte_so_far = mwg.DistToEnd(label);
      path_so_far = scratch_path;
      label_was_added = true;
    }
    else {
      // otherwise dte_so_far and path_so_far stay the same
      label_was_added = false;
    }
    if( max_dte[label] == DTE_IMPOSSIBLE || dte_so_far > max_dte[label] )
      return;

    if( in_degree[label]>1 
	&& ! exited_in_state[label].insert(path_so_far.GetHash(dte_so_far)).second )
      // We've explored from here on before; don't re-do the work
      return;
  }

  // Step through successors of current label:
  for( set<Label>::const_iterator succ = successors->begin();
       succ != successors->end(); succ++ ) {

    if( *succ == LabelDone ) { // base case
      if( ! label_was_added ) continue;
      
      // trim off overhang, if any:
      KmerPath path_to_add;
      if( mwg.DistToEnd(label)==0 ) {
	path_to_add = path_so_far;
      }
      else {
	KmerPathLoc insert_left_end = path_so_far.Begin();
	insert_left_end.IncrementMinGap( - mwg.DistToEnd(label) );
	path_so_far.CopyTail( insert_left_end, path_to_add );
      }

      ans.push_back_reserve( path_to_add, 0, 2.0 );
      // We might generate the same answer twice, from paths
      // differing on bad reads.  So we unique-sort at the end.
    }
    else { // recurse on this successor
      KmerPath new_path;
      if( mwg.IsGood(*succ) &&
	  pathifier.ExtendByKmersIfMatch( mwg.OKPI(*succ),
					  dte_so_far - mwg.DistToEnd(*succ),
					  path_so_far,
					  new_path ) ) {
	AllPathsHelper( new_path, mwg.DistToEnd(*succ), *succ, 
			true, pathifier, ans, max_ans );
      }
      else {
	AllPathsHelper( path_so_far, dte_so_far, *succ,
			false, pathifier, ans, max_ans );
      }
    }

    if( max_ans && ans.size() >= max_ans )
      break;
  }
}



/// Monte Carlo estimation of the number of cloures (in Mux space):
/// the size of a tree is the expected value of the product of the
/// degrees of the vertices, taken over all paths out from the root.
/// Naive "Algorithm E" from Knuth's "Estimating the efficiency of 
/// backtrack programs", Math. Comp. 29 (1975), 121--136.
int MuxWalkGraph::EstimateNumClosures() const {
  int num_trials = 100, total=0;
  for(int trial=0; trial < num_trials; trial++)
    total += EstimateNumClosuresTrial();
  return( int( .5 + (float(total)/float(num_trials)) ) );
}

int MuxWalkGraph::EstimateNumClosuresTrial() const {

  int product=1;
  const set<Label>* choices_p = &GetOpeners();
  while( 1 ) {
    int n = choices_p->size();  // takes linear time
    product *= n;
    set<Label>::iterator i = choices_p->begin();
    if (n>1) advance(i, randomx() % n);  // takes linear time again
    Label next = *i;
    if( next == LabelDone )
      return product;
    choices_p = &Successors(next);
  }
}


/// DOT-language output of the graph structure.
void MuxWalkGraph::Dot(ostream& out) const {
  // Begin graph
  out << "digraph MuxWalkGraph {\n  rankdir=LR;\n";
  out << "  edge [dir=back]; // Mux-walking is right-to-left\n";
  
  for( Label a=0; a != Label(m_label_lookup.size()); a++ ) {
    // Useful info to appear on node for a: OKPI, DistToEnd, IsGood.
    out << "  " << a << " [ "
	<< "label = \"node " << a << "\\n" << OKPI(a) << "\\ndte=" << DistToEnd(a) << '"'
	<< ( IsGood(a) ? "" : ", style=filled,color=lightgrey" )
	<< " ];\n";
    // Successor arrows:
    const set<Label>& succ = Successors(a);
    for( set<Label>::const_iterator bi = succ.begin(); bi != succ.end(); bi++ ) {
      if( *bi == LabelDone )
	out << "  " << a << " [shape=doubleoctagon];\n";
      else
	out << "  " << *bi << " -> " << a << ";\n";  // arrow points backwards
// 	out << "  " << a << " -> " << *bi << ";\n";
    }
  }
  
  // Indicate openers:
  for( set<Label>::iterator opi = m_openers.begin(); opi != m_openers.end(); opi++ )
    out << "  " << *opi << " [shape=doublecircle];\n";
  
  // End graph
  out << "}" << endl;
}




// Things for building a HyperKmerPath representation:

void MuxWalkGraph::MakeHyperKmerPath( const vecKmerPath* pathsFw, 
				      const vecKmerPath* pathsRc,
				      const SubsumptionList* subList,
				      HyperKmerPath& ans ) const {



  MuxWalkGraph_to_HyperKmerPath2 mwg2hkp( *this, pathsFw, pathsRc, subList );
  // To use the old HyperKmerPath generation tool, change the above to:
  // MuxWalkGraph_to_HyperKmerPath mwg2hkp( *this, pathsFw, pathsRc );




  // The old HyperKmerPath conversion could give wrong answers in
  // some situations.  It is still present in the source code.

  // It would be nice to have a (slow) mode where we test the old vs new.
// #define VERIFY_HKP

#ifdef VERIFY_HKP
  MuxWalkGraph_to_HyperKmerPath mwg2hkp_old( *this, pathsFw, pathsRc );
  HyperKmerPath ans_old = ans;
  mwg2hkp_old.MakeHyperKmerPath( ans_old );
#endif


  mwg2hkp.MakeHyperKmerPath( ans );


#ifdef VERIFY_HKP
  // To test isomorphism, we need to have a corresponding pair of vertices.
  // The HyperKmerPaths we generate usually have one source and one sink,
  // so either of those is a good place to start.  But it is possible, 
  // in a repetitive region, for there to be no sources and no sinks!
  // The test below will ignore such cases (as if they were empty HKPs).

  vec<int> sources, sources_old;
  ans.Sources(sources); 
  ans_old.Sources(sources_old);

  if( sources.size() != sources_old.size() ||
      ( ! sources.empty() && 
	! ComponentsAreIsomorphic( ans, sources[0], ans_old, sources_old[0] ) ) ) {
    cout << "=====================================" << endl;
    OrientedKmerPathId op_id = OKPI(*(GetOpeners().begin()));
    String op_name = ToString( op_id.GetId() ) + ( op_id.IsRc() ? "rc" : "fw" ) + ".dot";
    cout << "HyperKmerPaths differ!  Saving in *." << op_name << endl;
    Ofstream(dot_old, String("hkp_old.") + op_name);
    ans_old.PrintSummaryDOT0w(dot_old);
    Ofstream(dot_new, String("hkp_new.") + op_name);
    ans.PrintSummaryDOT0w(dot_new);
    Ofstream(dot_mwg, String("mwg.") + op_name);
    this->Dot( dot_mwg );
    cout << "\nOld hkp:" << endl;
    ans_old.PrintSummary(cout);
    cout << "\nNew hkp:" << endl;
    ans.PrintSummary(cout);
    cout << "=====================================" << endl;
  }
#endif

}


void MuxWalkGraph_to_HyperKmerPath::MakeHyperKmerPath( HyperKmerPath& ans ) {

//   {
//     cout << "Saving MWG in mwg.dot" << endl;
//     Ofstream(dot, "mwg.dot");
//     mwg.Dot(dot);
//   }

  CalcPredecessors();

  vertex_status.clear();
  vertex_status.resize( mwg.size(), VX_UNKNOWN );
  // Change to VX_BAD or VX_DONE as each label is processed.


  // We are going to build a HyperKmerPath which reflects
  // the structure of the MuxWalkGraph.  There is a vertex of
  // the HKP corresponding to each Node of the MWG, which
  // represents the left endpoint of the associated read
  // (plus some extras for the opener right endpts).
  // To minimze bookkeeping, the int Label in the MWG will
  // also be the vertex id number in the HKP (the extra
  // vertices will go at the end).

  ans.Clear();
  ans.AddVertices( mwg.size() );

  AddOpeners( ans );

  // Lay reads into HKP, stepping through backwards (probably openers first)
  for( Label label = mwg.size()-1; label>=0; label-- )
    AddEdgesFrom( label, ans );

  // Immediately, while we still know which HKP vertices come from closers,
  // remove anything that cannot be reached by going right from a closer.
  RemoveUnreachable( ans );

  // Now associate closers.  There may be multiple closer nodes in the MWG,
  // arising from either multiple fills of the closing read, or from multiple
  // reads subsuming a closer, or both.  It's not entirely clear to me how to
  // correctly handle both possibilities simultaneously.  At the moment, we
  //  (1) trim all closer-subsuming reads to dte=0 (in AddEdgeFromTo), and
  //  (2) associate the left ends of all closers with dte <= 0 (here)
  // This is somewhat ham-fisted, and if we are seriously using multiple 
  // fills of reads, we'll need to reexamine this -- probably changing the 
  // Mux searcher to store some subsumption data in the MuxWalkGraph.
  AssociateClosers( ans );

//   {
//     cout << "Saving dirty HKP in dirty.dot" << endl;
//     Ofstream(dot, "dirty.dot");
//     ans.PrintSummaryDOT0w(dot);
//   }

  // Muxes which extend by zero kmers will give us empty edges.
  ans.ContractEmptyEdges();

  // Clean up the resulting graph
  HyperKmerPathCleaner().CleanUpGraph( ans );

//   {
//     cout << "Saving clean HKP in clean.dot" << endl;
//     Ofstream(dot, "clean.dot");
//     ans.PrintSummaryDOT0w(dot);
//   }

  return;

}


void MuxWalkGraph_to_HyperKmerPath::CalcPredecessors( ) {
  predecessors.clear();
  predecessors.resize( mwg.size() );
  closers.clear();

  for( Label i=0; i < mwg.size(); i++ ) {
    const set<Label>& succ = mwg.Successors(i);
    for( set<Label>::const_iterator s = succ.begin(); s!=succ.end(); s++ ) {
      if( *s == LabelDone )
	closers.push_back(i);
      else
	predecessors[*s].push_back(i);
    }
  }
}



// Add one edge for each opener.  The left vertex of the
// edge is the preexisting vertex whose number is the 
// opener's Label.  The right vertex corresponds to the OKPI
// of the read.  In this way, when we zip the graph up at the end,
// all closures starting with the same read will be associated,
// even if they are different lengths (so different MWG nodes).
void MuxWalkGraph_to_HyperKmerPath::AddOpeners( HyperKmerPath& ans ) const {

  map<OrientedKmerPathId,int> okpi_to_vx;
  pair< map<OrientedKmerPathId,int>::iterator, bool > insert_result;

  for( set<Label>::const_iterator op = mwg.GetOpeners().begin();
       op != mwg.GetOpeners().end(); op++ ) {

    insert_result = okpi_to_vx.insert( make_pair(mwg.OKPI(*op),-1) );

    int& okpi_vx = insert_result.first->second;
    if( insert_result.second ) {  
      // this OKPI was previously unknown
      // Create a vertex for it.
      okpi_vx = ans.N();  // replacing the -1 just inserted
      ans.AddVertices(1);
    }

    // Make an edge for this opener.
    ans.AddEdge( *op, okpi_vx, *(mwg.OKPI(*op).GetPathPtr(*pathsFw,*pathsRc)) );
  }
}


// Modify the HKP to add in the read from the mux label.
void MuxWalkGraph_to_HyperKmerPath::AddEdgesFrom( Label label,
						  HyperKmerPath& ans ) {
  
  if( vertex_status[label] != VX_UNKNOWN ) return;

  // Make sure all my predecessors' statuses are known:
  const vec<Label>& preds = predecessors[label];
  for( vec<Label>::const_iterator p = preds.begin(); p != preds.end(); p++ ) {
    if( vertex_status[*p] == VX_UNKNOWN )
      AddEdgesFrom( *p, ans );
  }

  // Never use bad reads:
  if( ! mwg.IsGood(label) ) {
    vertex_status[label] = VX_BAD;
    return;
  }

  // Figure out which nodes we'll need edges to:
  const KmerPath& the_path = *(mwg.OKPI(label).GetPathPtr(*pathsFw,*pathsRc));
  int max_dte = mwg.DistToEnd(label) + the_path.MinLength() - 1;

  compat_map compat;
  for( vec<Label>::const_iterator p = preds.begin(); p != preds.end(); p++ )
    CheckCompatibility( label, *p, the_path, max_dte, compat, ans );

  // Add the edges:
  for( compat_iter it = compat.begin(); it != compat.end(); it++ )
    if( it->second == CP_GOOD )
      AddEdgeFromTo( label, it->first, the_path, ans );

  vertex_status[label] = VX_DONE;
}

MuxWalkGraph_to_HyperKmerPath::compat_iter 
MuxWalkGraph_to_HyperKmerPath::CheckCompatibility( Label from_label, 
				      Label to_label, 
				      const KmerPath& the_path,
				      int max_dte,
				      compat_map& compat,
				      const HyperKmerPath& ans ) const {

  if( to_label >= mwg.size() || mwg.DistToEnd(to_label) > max_dte )
    // to_label >= mwg.size()      : for right endpt vx of opener
    // mwg.DistToEnd(to_label) > max_dte : caller should know better
    return compat.insert( make_pair(to_label,CP_DOM) ).first;

  pair<compat_iter,bool> ins = compat.insert( make_pair(to_label,CP_INCOMP) );
  if( ins.second == false ) // compatibility is already known!
    return ins.first;

  // We need to check compatibility with all previous vertices
  // (back to max_dte) in the MuxWalkGraph.  Yes, all of them -- 
  // it's possible that there will be an edge from from_label to 
  // something very far away in Mux-space, regardless of the
  // statuses of intervening nodes.  Sadly, this is time-consuming.

  const vec<Label>& preds = predecessors[to_label];
  for( vec<Label>::const_iterator p = preds.begin(); p != preds.end(); p++ )
    if( mwg.DistToEnd(*p) <= max_dte )
      CheckCompatibility( from_label, *p, the_path, max_dte, compat, ans );

  // Step through edges in the HKP from to_label, looking for
  // compatibility with the appropriate part of the_path.
  // * If a compatible edge reaches the end of the_path,
  //   then this vertex (to_label) is good.
  // * If a compat edge reaches another vx decleared good, then
  //    this vx is good, AND the far edge is dominated.

  // where each edge out goes to:
  const vec<int>& edge_dest = ans.From(to_label);

  for( int e = 0; e < edge_dest.isize(); e++ ) {  // for each edge:

    // This is the KmerPath on the e'th HKP edge out of to_label:
    const KmerPath& edge_path = ans.EdgeObjectByIndexFrom(to_label,e);
    KmerPathLoc edge_loc = edge_path.Begin();

    // This is the corresponding position on the_path:
    KmerPathLoc path_loc = the_path.Begin();
    path_loc.IncrementMinGap( mwg.DistToEnd(to_label)-mwg.DistToEnd(from_label) );
    
    // If they're incompatible, we don't care about this edge:
    if( ! ScanRightPerfectMatchGaps( edge_loc, path_loc ) )
      continue;

    // If we reached the end of the_path, this vx is good.
    // If we reached the end of edge_path and the far vertex
    // is good, this vx is good and that one is now dominated.
    if( path_loc.atEnd() ) {
      ins.first->second = CP_GOOD;
    }
    else {
      // Now we need to know whether the far end is compatible in turn:
      // This was already computed, *unless* edge_dest[e] is a vertex
      // corresponding to the right endpoint of an opener.
      compat_iter far = CheckCompatibility( from_label, edge_dest[e],
						the_path, max_dte, compat, ans );
      if( far->second == CP_GOOD || far->second == CP_DOM ) {
	far->second = CP_DOM;
	ins.first->second = CP_GOOD;
      }
    }
  }
  // If we exit this loop and haven't found an edge that makes this vx good,
  // then its compatibility remains CP_BAD, per the original insert.
  return ins.first;
}


void MuxWalkGraph_to_HyperKmerPath::AddEdgeFromTo( Label from_vx, 
				      Label to_vx, 
				      const KmerPath& the_path, 
				      HyperKmerPath& ans ) const {
  KmerPath new_edge_path;

  int dte_from = mwg.DistToEnd(from_vx), dte_to = mwg.DistToEnd(to_vx);

  KmerPathLoc copy_start = the_path.Begin();
  if( dte_from < 0 )  // this read subsumed the closer; trim it.
    copy_start.IncrementMinGap( -dte_from );

  KmerPathLoc copy_end = the_path.Begin();
  copy_end.IncrementMinGap( max(0,dte_to) - dte_from );

  the_path.CopySubpathNoLastKmer( copy_start, copy_end, new_edge_path );

  ans.AddEdge( from_vx, to_vx, new_edge_path );
}

 
// We may have inserted edges that can't be reached by
// walking right from a closer.  (These can arise from
// Mux-space closures that are illegal in path-space.)
// We need to check this condition as soon as we finish 
// building the graph, while we still know which hkp
// vertices correspond to closing reads. 
void MuxWalkGraph_to_HyperKmerPath::RemoveUnreachable( HyperKmerPath& ans ) const {
  vec<Bool> reachable( ans.N(), false );
  vec<int> to_check( closers );
  while( ! to_check.empty() ) {
    int vx = to_check.back();
    to_check.pop_back();
    if( ! reachable[vx] )
      to_check.insert( to_check.end(), 
		       ans.From(vx).begin(),
		       ans.From(vx).end() );
    reachable[vx] = true;
  }
  // Delete all edges out of unreachable vertices:
  for( int vx = 0; vx < mwg.size(); vx++ )
    if( ! reachable[vx] )
      ans.DeleteEdgesAtVertex(vx);
  // These edgeless vertices can stick around for now;
  // they are removed by CleanUpGraph  
}


// If there are multiple closers with dte <= 0, associate their vertices.
void MuxWalkGraph_to_HyperKmerPath::AssociateClosers( HyperKmerPath& ans ) const {

  vec<Label>::const_iterator base, other;
  for( base = closers.begin(); base != closers.end(); base++ )
    if( mwg.DistToEnd(*base) <= 0 ) break;

  if( base == closers.end() ) return;

  for( other = base+1; other != closers.end(); other++ )
    if( mwg.DistToEnd(*other) <= 0 )
      ans.TransferEdges( *other, *base );
}


void MuxWalkGraph_to_HyperKmerPath2::FindClosers() {

  closers.clear();
  reads_which_subsume_closers.clear();
  subsumed_closers.clear();

  is_closer.clear();
  is_closer.resize( mwg.size(), false );

  for( Label label=0; label < mwg.size(); label++ )
    if( mwg.Successors(label).count(LabelDone) ) {
      closers.push_back(label);
      is_closer[label] = true;
      reads_which_subsume_closers[ make_pair(OKPI(label),
					     LeftEndDTE(label)) ] = label;
    }

  vec<Label> more_closers;

  // Now look for reads which subsume those:
  for( vec<Label>::iterator cl = closers.begin();
       cl != closers.end(); cl++ ) {
    vec<SubsumptionRecord> this_closer_subs;
    subList->GetFullRecordsFor( OKPI(*cl), this_closer_subs );
    for( vec<SubsumptionRecord>::iterator sub_rec = this_closer_subs.begin();
	 sub_rec != this_closer_subs.end(); sub_rec++ ) {
      pair<OrientedKmerPathId,int> key( sub_rec->GetSuperPathId(),
					LeftEndDTE(*cl) - sub_rec->GetLeftOverhang() );
      pair<sub_map_t::iterator,bool> insert_result =
	reads_which_subsume_closers.insert( make_pair(key,NumLabels()) );
      if( insert_result.second ) {
	// this subsuming read was previously unknown
	// Build a fake record for it, so that everything in this class 
	// thinks the Label new_label refers to this subsuming read.
	Label new_label = NumLabels();
	more_closers.push_back(new_label);
	okpids.push_back( key.first );
	left_end_dte.push_back( key.second );
	right_end_dte.push_back( key.second + ReadLength(new_label) - 1 );
      }
      // In either case, record the subsumption
      subsumed_closers[*cl].push_back( insert_result.first->second );
    }
  }

  // Add the more_closers to closers.
  closers.insert( closers.end(), more_closers.begin(), more_closers.end() );
  // All the new labels are closers.
  is_closer.resize( NumLabels(), true );

}



void MuxWalkGraph_to_HyperKmerPath2::FindSubsumedOpeners() {

  map< pair<OrientedKmerPathId,int>, vec<Label> > subsume_openers;

  for( set<Label>::const_iterator opener = mwg.GetOpeners().begin();
       opener != mwg.GetOpeners().end(); opener++ ) {
    if( ! mwg.IsGood(*opener) ) continue;
    int op_dte = LeftEndDTE(*opener);

    vec<SubsumptionRecord> this_opener_subs;
    subList->GetFullRecordsFor( OKPI(*opener), this_opener_subs );
    for( vec<SubsumptionRecord>::iterator sub_rec = this_opener_subs.begin();
	 sub_rec != this_opener_subs.end(); sub_rec++ )
      subsume_openers[make_pair(sub_rec->GetSuperPathId(),
				op_dte - sub_rec->GetLeftOverhang())].push_back(*opener);
  }

  openers_subsumed_by.clear();
  openers_subsumed_by.resize( NumLabels() );
  map< pair<OrientedKmerPathId,int>, vec<Label> >::iterator it;

  for( Label label=0; label < NumLabels(); label++ )
    if( (it = subsume_openers.find(make_pair( OKPI(label),LeftEndDTE(label) )))
	!= subsume_openers.end() ) {
      sort( it->second.begin(), it->second.end(),
	    indirect_compare< int,less<int> >(right_end_dte) );
      openers_subsumed_by[ label ] = it->second;
    }
  // So if label i occurs in openers_subsumed_by[j], then i is an opener
  // and read j subsumes it.
}



// Builds the vec<Label> topo_order, in which each Label
// appears after all its successors.  Also builds closers.
void MuxWalkGraph_to_HyperKmerPath2::ReverseTopoSort() {

//   closers are now found in FindClosers instead.
//   closers.clear();
//   bool is_closer;

  topo_order.clear();
  topo_order.reserve( NumLabels() );
  vec<Bool> done(NumLabels(), false);

  vec<Label> to_do;

  for(Label i=0; i<mwg.size(); i++) {
    if( done[i] ) continue;

    to_do.push_back(i);
    while( ! to_do.empty() ) {

      Label j = to_do.back();
      const set<Label>& succ = mwg.Successors(j);
      for( set<Label>::const_iterator s = succ.begin(); s != succ.end(); s++ ) {
	if( *s == LabelDone ) {
	  // Need to deal with reads that subsume this closer
	  const vec<Label>& subsuming_reads = subsumed_closers[j];
	  for( vec<Label>::const_iterator subsumer = subsuming_reads.begin();
	       subsumer != subsuming_reads.end(); subsumer++ ) {
	    if( done[*subsumer] )
	      continue;
	    else if( *subsumer >= mwg.size() ) {
	      // read not in the MWG, so has no Successors -- add now
	      topo_order.push_back(*subsumer);
	      done[*subsumer] = true;
	    }
	    else {
	      to_do.push_back(*subsumer);
	    }
	  }
	}
	else if( ! done[*s] ) {
	  to_do.push_back(*s);
	}
      }
      if( j == to_do.back() ) { // succ's all done (since MWG is loop free)
	topo_order.push_back( j );
	done[j] = true;
	to_do.pop_back();
      }
    }
  }
  // Each v appears *after* all of its successors, ie vertices appear l-to-r.
}


// For each Label v, find all the Labels w such that it's possible
// to get from v to w by following the MWG, *and* the reads for v and w
// overlap and match.
// Note that there may be reads which subsume a closer that don't appear
// in the MWG at all, so they aren't found here.
void MuxWalkGraph_to_HyperKmerPath2::FindMatchingDescendants() {
  matching_descendants.clear();
  matching_descendants.resize( NumLabels() );
  
  vec< vec<Label> > overlapping_descendants;
  overlapping_descendants.resize( NumLabels() );

  // We cache the KmerPathLocs of the descendants where their last
  // ancestor(?) began, since each ancestor begins further to the
  // right than the last.  It can be expensive to traverse the
  // descendants from their starts each time if they are composed of
  // many intervals.

  vec<KmerPathLoc> cachedLocs;
  vec<int> cachedDTE;
  for ( Label i = 0; i < NumLabels(); ++i ) {
    cachedLocs.push_back( Read(i).Begin() );
    cachedDTE.push_back( LeftEndDTE(i) );
  }

  for( vec<Label>::const_iterator v = topo_order.begin(); 
       v != topo_order.end(); v++ ) {
    if( *v >= mwg.size() ) // Read not in the MWG (but subsumes a closer)
      continue;            // These have no descendants, by definition.
    int dte = LeftEndDTE(*v);
    set<Label> succ = mwg.Successors(*v); // all of these are already done
    if( succ.find(LabelDone) != succ.end() ) {
      // this is a closer -- to its successors, add reads which subsume it!
      const vec<Label>& subsumers = subsumed_closers[*v];
      succ.insert( subsumers.begin(), subsumers.end() );
    }

    for( set<Label>::const_iterator w = succ.begin(); w != succ.end(); w++ ) {
      if( *w == LabelDone ) continue;
      ForceAssertLe( dte, RightEndDTE(*w) ); // reads overlap their muxes
      overlapping_descendants[*v].push_back(*w);
      for( vec<Label>::const_iterator x = overlapping_descendants[*w].begin();
	   x != overlapping_descendants[*w].end(); x++ ) {
	if( dte <= RightEndDTE(*x) )
	  overlapping_descendants[*v].push_back(*x);
      }
    }
    UniqueSort(overlapping_descendants[*v]);

    // Now check which actually agree with the read:
    const KmerPath& read_v = Read(*v);
    for( vec<Label>::const_iterator w = overlapping_descendants[*v].begin();
	 w != overlapping_descendants[*v].end(); w++ ) {
      if ( cachedDTE[*w] != dte ) {
        //PRINT3( okpids[*v], okpids[*w], dte - cachedDTE[*w] );
        cachedLocs[*w].IncrementMinGap( dte - cachedDTE[*w] );
        cachedDTE[*w] = dte;
      }
      KmerPathLoc v_loc = read_v.Begin(), w_loc( cachedLocs[*w] );
      if( v_loc.GetKmer() == w_loc.GetKmer() && 
	  // IsPerfectMatchToRight asserts loc intervals overlap
	  IsPerfectMatchToRight( v_loc, w_loc ) )
	matching_descendants[*v].push_back(*w);
    }
  }
}

// tells which edge into this vx starts with this kmer (-1 if none)
int WhichEdgeIn( int vx, const HyperKmerPath& hkp, longlong kmer ) {
  for( int e=0; e < hkp.To(vx).isize(); e++ )
    if( hkp.EdgeObjectByIndexTo(vx,e).End().GetKmer() == kmer )
      return e;
  return -1;
}
// tells which edge out of this vx starts with this kmer (-1 if none)
int WhichEdgeOut( int vx, const HyperKmerPath& hkp, longlong kmer ) {
  for( int e=0; e < hkp.From(vx).isize(); e++ )
    if( hkp.EdgeObjectByIndexFrom(vx,e).Start(0) == kmer )
      return e;
  return -1;
}

// split an edge of an HKP just after the kmer point at by loc.
// The last edge added -- the one with index EdgeObjectCount()-1
// -- is the right half of the original edge.
void SplitEdge( int from_vx, int e, int split_vx,
		const KmerPathLoc& loc, HyperKmerPath& ans ) {
  int to_vx = ans.From(from_vx)[e];
  const KmerPath& old_edge = ans.EdgeObjectByIndexFrom(from_vx,e);
  KmerPath left_half, right_half;
  old_edge.CopyHead( loc, left_half );
  old_edge.CopyTailNoFirstKmer( loc, right_half );
  ans.DeleteEdgeFrom(from_vx,e);
  ans.AddEdge(from_vx, split_vx, left_half);
  ans.AddEdge(split_vx, to_vx, right_half);
}

// Add a path (generally a read) to the HyperKmerPath, by attaching
// the loc on the path to a given vertex, where it is already known
// that it matches to the left.  Scan right, and create a new edge
// when the read diverges from what is already in the HKP.
// Returns true if the read gets added.
//
// The flag add_if_subsumed controls what happens if the read is
// already entirely subsumed by the graph: if true, add a vertex for
// the read's right end; if false, do nothing.
//
// If a set<int>* is given, add to the set all the vxs the path passes 
// through, and abort if we get to a vertex already seen.  This is to
// keep us from doing redundant work.
void MuxWalkGraph_to_HyperKmerPath2::AttachRead( Label label,
						 KmerPathLoc path_loc,
						 int attach_vx,
						 const int right_end_vx,
						 HyperKmerPath& ans,
						 set<int>* vxs_seen ) {
  // scan right along perfect matches.
  int e;
  KmerPathLoc hkp_edge_loc;

  while(1) { 

    if( vxs_seen && (vxs_seen->insert(attach_vx).second==false) )
      return;  // we've done this attach before

    e = WhichEdgeOut( attach_vx, ans, path_loc.GetKmer() );
    if( e == -1 ) break;

    hkp_edge_loc = ans.EdgeObjectByIndexFrom(attach_vx,e).Begin();
    ScanRightPerfectMatchGaps( path_loc, hkp_edge_loc );
    
    if( hkp_edge_loc.atEnd() && ! path_loc.atEnd() ) {
      attach_vx = ans.From(attach_vx)[e];
      path_loc.IncrementMinGap(1);
      continue;  // successfully followed an edge.  Keep going.
    }
    else break;
  } // that loop has the most awkward flow control I've written in a while.

  // Now go through the various cases for where the scan can end:
  // 1. The path does not match any outgoing edge, so branches off the HKP.
  //    New edge that goes into the global right_end_vx
  if( e == -1 ) {
    label_to_vxs[label].insert(right_end_vx);
    KmerPath new_edge;
    path_loc.GetPath().CopyTail( path_loc, new_edge );
    ans.AddEdge( attach_vx, right_end_vx, new_edge );
    AddSubsumedOpenerVertices( attach_vx, right_end_vx,
			       ans.EdgeObjectCount()-1,
			       label, ans );
  }
  else if( path_loc.atEnd() ) {
    // 2. The path and edge end at the same time: no HKP change
    //    Perhaps we don't even need to record this in label_to_vxs?
    if( hkp_edge_loc.atEnd() )
      label_to_vxs[label].insert( ans.From(attach_vx)[e] );
    // 3. The path ends strictly inside the edge.  Split the edge.
    else {
      int split_vx = ans.N();
      ans.AddVertices(1);
      SplitEdge( attach_vx, e, split_vx, hkp_edge_loc, ans );
      label_to_vxs[label].insert( split_vx );
    }
  }
  // 4. The path splits off from the edge somewhere in the middle.
  //    Split the edge; branch off a new edge to the global right_end_vx
  else {
    ForceAssert( ! hkp_edge_loc.atEnd() );
    label_to_vxs[label].insert( right_end_vx );
    int branch_vx = ans.N();
    ans.AddVertices(1);
    SplitEdge( attach_vx, e, branch_vx, hkp_edge_loc, ans );
    KmerPath new_edge;
    path_loc.GetPath().CopyTailNoFirstKmer( path_loc, new_edge );
    ans.AddEdge( branch_vx, right_end_vx, new_edge );
    AddSubsumedOpenerVertices( branch_vx, right_end_vx,
			       ans.EdgeObjectCount()-1,
			       label, ans );
  }
}


// An edge (number edge_object_id, from from_vx to to_vx) was just added
// to ans, associated with the given label from the MWG.  Does the newly-added
// edge contain any points which should be right endpoints of subsumed openers?
// If so, we need to split the edge, make new vertices, and record them.
void MuxWalkGraph_to_HyperKmerPath2::AddSubsumedOpenerVertices( int from_vx,
								int to_vx,
								int edge_obj_id,
								const Label label,
								HyperKmerPath& ans ) {
  for( vec<Label>::iterator opener_label = openers_subsumed_by[label].begin();
       opener_label != openers_subsumed_by[label].end(); opener_label++ ) {
    KmerPathLoc loc = ans.EdgeObject(edge_obj_id).End();
    bool found_subsumed = 
      loc.IncrementMinGap( - RightEndDTE(label) + RightEndDTE(*opener_label) );
    // Note that this is a negative number!  We're moving left from the end.
    // If that doesn't fall off the left end, then there is a subsumed opener.
    if( ! found_subsumed ) continue;
    int opener_right_vx = ans.N();
    ans.AddVertices(1);
    label_to_vxs_no_guarantee[ *opener_label ].insert( opener_right_vx );
    SplitEdge( from_vx, ans.EdgeObjectIndexToFromIndex(from_vx,edge_obj_id),
             opener_right_vx, loc, ans );
    // focus attention on the new, smaller edge.  This works because
    // openers_subsumed_by was sorted by increasing right dte of opener.
    edge_obj_id = ans.EdgeObjectCount()-1;  // the right edge just added by SplitEdge
    from_vx = opener_right_vx;
  }
}


void MuxWalkGraph_to_HyperKmerPath2::ProcessLabel( Label label,
						   HyperKmerPath& ans ) {

  int right_end_vx = ans.N();
  ans.AddVertices(1);

  if( isCloser(label) ) {
    KmerPathLoc closer_path_loc = Read(label).Begin();
    ForceAssertLe( LeftEndDTE(label), 0 );
    // If we start using data where reads have been trimmed so that
    // the closers don't reach the end of the insert, we'll have
    // to revisit how to build the left end of the graph.
    closer_path_loc.IncrementMinGap( -LeftEndDTE(label) );
    AttachRead( label, closer_path_loc,
		vx_left_end_of_insert, right_end_vx, ans );
    return;
  }

  set<int> vxs_seen;
  const KmerPath& the_read = Read(label);
  int my_left_dte = LeftEndDTE(label), my_right_dte = RightEndDTE(label);

  const vec<Label>& desc = matching_descendants[label];

  for( vec<Label>::const_iterator d = matching_descendants[label].begin();
       d != matching_descendants[label].end(); d++ ) {
    // If we are coterminal on the right with our descendant,
    // then we'll have no effect on the graph.  For completeness,
    // we'll record this in label_to_vxs, but even that only
    // matters if we are a closing read.
    if( RightEndDTE(*d) == my_right_dte ) {
      label_to_vxs[label].insert( label_to_vxs[*d].begin(),
				  label_to_vxs[*d].end() );
      continue;
    }
    for( set<int>::const_iterator attach_vx = label_to_vxs[*d].begin();
	 attach_vx != label_to_vxs[*d].end(); attach_vx++ ) {
      if( vxs_seen.find(*attach_vx) != vxs_seen.end() )
	continue;
      KmerPathLoc read_loc = the_read.Begin();
      bool no_fall_off_end = 
	read_loc.IncrementMinGap( RightEndDTE(*d) - my_left_dte + 1 );
      // +1 since read_loc points at the kmer *after* the last matching one
      // We fall off the end if the right end of *d is right of the
      // right end of this read.  In that case, we don't want to attach.
      if( no_fall_off_end )
	AttachRead( label, read_loc, *attach_vx, right_end_vx, ans, &vxs_seen );
    }
  }
}



int MuxWalkGraph_to_HyperKmerPath2::AssociateOpeners( HyperKmerPath& ans ) {

  // Create a new vertex to represent the right end of the insert:
  int right_end_of_insert_vx = ans.N();
  ans.AddVertices(1);

  // Associate all right vertices of openers with the new vx:
  for( set<Label>::const_iterator op = mwg.GetOpeners().begin();
       op != mwg.GetOpeners().end(); op++ ) {
    for( set<int>::const_iterator vx = label_to_vxs[*op].begin();
	 vx != label_to_vxs[*op].end(); vx++ )
      ans.TransferEdges( *vx, right_end_of_insert_vx );
    for( set<int>::const_iterator vx = label_to_vxs_no_guarantee[*op].begin();
	 vx != label_to_vxs_no_guarantee[*op].end(); vx++ )
      ans.TransferEdges( *vx, right_end_of_insert_vx );
  }

  // Due to subsumed openers (I think?), there might have been two
  // reads joined by an empty edge which just got associated, making a
  // loop at right_end_of_insert_vx of length zero.  Remove such loops.

  for( int i = ans.From(right_end_of_insert_vx).isize()-1; i >= 0; i-- )
    if( ans.From(right_end_of_insert_vx)[i] == right_end_of_insert_vx &&
	ans.EdgeObjectByIndexFrom(right_end_of_insert_vx,i).IsEmpty() )
      ans.DeleteEdgeFrom(right_end_of_insert_vx,i);

  return right_end_of_insert_vx;
}



// Remove all parts of the graph which can't be reached my moving
// left from the given right_vx.  (We run this just after associating
// all opener right endpoints to a single vertex.)
// Among the things this will delete are edges to the *right* of
// the right_vx, which can arise from reads longer than the opener.
void MuxWalkGraph_to_HyperKmerPath2::RemoveUnreachable( int right_vx, 
							HyperKmerPath& ans ) const {
  vec<Bool> reachable( ans.N(), false );
  vec<int> to_check( 1, right_vx );
  while( ! to_check.empty() ) {
    int vx = to_check.back();
    to_check.pop_back();
    if( ! reachable[vx] )
      to_check.insert( to_check.end(), 
		       ans.To(vx).begin(),
		       ans.To(vx).end() );
    reachable[vx] = true;
  }
  // Delete all edges out of unreachable vertices:
  for( int vx = 0; vx < ans.N(); vx++ )
    if( ! reachable[vx] )
      ans.DeleteEdgesAtVertex(vx);
  // These edgeless vertices can stick around for now;
  // they are removed by CleanUpGraph  
}



void MuxWalkGraph_to_HyperKmerPath2::MakeHyperKmerPath( HyperKmerPath& ans ) {

//   {
//     static int count = 1;
//     String dotfile = "mwg" + ToString(count++) + ".dot";
//     cout << "Saving MWG in " << dotfile << endl;
//     Ofstream(dot, dotfile);
//     mwg.Dot(dot);
//     cout << "Done." << endl;
//   }

  FindClosers();
  ReverseTopoSort();
  FindMatchingDescendants();
  FindSubsumedOpeners();

  label_to_vxs.clear();
  label_to_vxs.resize( NumLabels() );
  label_to_vxs_no_guarantee = label_to_vxs;

  ans.Clear();
  vx_left_end_of_insert = ans.N();  // =0
  ans.AddVertices( 1 );

  for( vec<Label>::const_iterator v = topo_order.begin();
       v != topo_order.end(); v++ )
    ProcessLabel( *v, ans );

  int right_vx = AssociateOpeners(ans);

  RemoveUnreachable( right_vx, ans );

//   {
//     temp_file tempname( "tmp/dirty.XXXXXX" );
//     String dotname = tempname + ".dot";
//     cout << "Saving dirty HKP in " << dotname << endl;
//     Ofstream(dot, dotname);
//     ans.PrintSummaryDOT0w(dot);
//   }

  HyperKmerPathCleaner().CleanUpGraph( ans );

//   {
//     temp_file tempname( "tmp/clean.XXXXXX" );
//     String dotname = tempname + ".dot";
//     cout << "Saving clean HKP in " << dotname << endl;
//     Ofstream(dot, dotname);
//     ans.PrintSummaryDOT0w(dot);
//   }

}  

