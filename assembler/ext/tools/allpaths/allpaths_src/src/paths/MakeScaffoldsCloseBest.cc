 /////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Equiv.h"
#include "Fastavector.h"
#include "Charvector.h"
#include "FastIfstream.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/Alignlet.h"
#include "paths/HyperFastavector.h"
#include "paths/MakeScaffoldsCloseBest.h"
#include "paths/PseudoUnipatherLib.h"
#include "paths/SInsertion.h"
#include "paths/SuckScaffolds.h"
#include "paths/BreakBadScaffolds.h"
#include "paths/Sepdev.h"
#include "paths/reporting/CLinkBundle.h"

// An slink is a link between scaffolds.  The first scaffold is not represented
// in the object.

class slink {

     public:

     slink( ) { }
     slink( const int s2, const Bool fw2, const int sep, const int dev )
          : s2(s2), fw2(fw2), sep(sep), dev(dev) { }

     int s2;
     Bool fw2;
     int sep, dev;

     friend Bool operator<( const slink& l1, const slink& l2 ){    
       if ( l1.s2 < l2.s2 ) return True;
       if ( l1.s2 > l2.s2 ) return False;
       return l1.fw2 < l2.fw2;    
     }

};

void
dump_linkgraph(String const& filename,
	       digraphE<sepdev> const& graph,
	       int const iter) 
{
  // Label each edge with super, orientation, and 
  vec< vec<String> > labels;
  for (int v1 = 0; v1 < graph.N(); ++v1) {
    int s1 = v1 / 2;
    String s1dir = (v1 & 1) ? "rc" : "fw";
    vec<String> v1labels;
    vec<int> from = graph.From(v1);
    for (int i = 0; i < from.isize(); ++i) {
      ostrstream label;
      int v2 = from[i];
      int s2 = v2 / 2;
      String s2dir = (v2 & 1) ? "rc" : "fw";
      sepdev s = graph.EdgeObjectsBetween(v1,v2)[0];
      label << s1 << s1dir << ">" << s2 << s2dir << "(" << s.Sep( ) << ")" << ends;
      //label << s << ends;
      v1labels.push_back(label.str());
    }
    labels.push_back(v1labels);
  }
    
  ostrstream dotfile;
  dotfile << filename;
  dotfile << "." << iter;
  dotfile << ".dot" << ends;
  Ofstream(dot, dotfile.str());
  graph.DOT(dot, labels);
  dot.close();
}

// Clean graph (remove inconsistencies, remove all links except nearest).

template<class E> 
void ShaveGraph( digraphE<E> &G ){ 
  
  vec<int> to_delete;
  // remove incoming
  for ( int v = 0; v < G.N(); v++ ){
    if ( G.To(v).size() == 0 && G.From(v).size() == 1 ){ 
      size_t w = G.From(v)[0];
      if ( G.To(w).size() > 1 ){
	for ( size_t s = 0; s < G.To(w).size(); s++ ){
	  if ( G.To(s).size() > 0 ){ 
	    to_delete.append( G.EdgesBetween( v, w ) );
	    break;
	  }
	}
      }
    }
  }
  int nBadIncoming = to_delete.size();
  // remove outgoing 
  for ( int w = 0; w < G.N(); w++ ){
    if ( G.From(w).size() == 0 && G.To(w).size() == 1 ){
      size_t v = G.To(w)[0];
      if ( G.From(v).size() > 1 ){
	for ( size_t s = 0; s < G.From(v).size(); s++ ){
	  if ( G.From(s).size() > 0 ){
	    to_delete.append( G.EdgesBetween( v, w ) );
	    break;
	  }
	}
      }
    }
  }
  
  // cout << "Shaving function identified " << to_delete.size()
  //      << " edges to delete (" << nBadIncoming
  //      << " incoming and " << to_delete.size() -nBadIncoming
  //      << " outgoing)" << endl;

  G.DeleteEdges( to_delete );
  G.RemoveDeadEdgeObjects();
}

template<class E> 
void CleanupGraph( const int VERBOSITY,
		   const int dev_mult,
		   const vec<superb> &scaffolds,
		   digraphE<E> &G )
{ 
  const int nscaff = scaffolds.isize( );
  
  // Cleanup strand bidiriectionals
  for ( int s1 = 0; s1 < nscaff*2 -1; s1++ ){
    for ( int s2 = s1 + 1; s2 < nscaff*2; s2++ ){
      if ( G.EdgesBetween(s1,s2).size() > 0 ){
	      
	if ( G.EdgesBetween(s2,s1).size() > 0 ){
	  G.DeleteEdges(  G.EdgesBetween(s1,s2) );
	  G.DeleteEdges(  G.EdgesBetween(s2,s1) );
	  if (VERBOSITY >= 2){
	    PRINT( G.EdgesBetween(s1,s2).size() );
	    PRINT( G.EdgesBetween(s2,s1).size() );
	    // cout << "DELETING BIDIRECTIONAL: edges between s1 = " 
	    // 	 << s1 << " and s2 = " << s2 << " were deleted" << endl;
	  }
	}
      }	    
    }
  }
  
  // Cleanup outgoing
  int outConflictCt = 0, outConflictLoc = 0;;
  int delEdgeCt = 0, wantDelEdgeCt = 0;
  vec<int> edgesToDel, allOutEdges;
  if ( VERBOSITY >= 2 ) 
    cout << Date() << " Cleaning up graph (out)" << endl;    
  for ( int s = 0; s < nscaff*2; s++ ){
    if ( G.From(s).size() > 1 ){
      vec< pair<int,sepdev> > links( G.From(s).size() );
      for ( size_t i = 0; i < G.From(s).size(); i++ ){
	int s2 = G.From(s).at(i);	
	ForceAssert( G.EdgeObjectsBetween(s,s2).size() ==1 );
	sepdev sdloc( G.EdgeObjectsBetween(s,s2).at(0).Sep(), G.EdgeObjectsBetween(s,s2).at(0).Dev() );
	links.at(i) = pair<int,sepdev>(s2, sdloc); 
      }
      
      // Find closest
      ForceAssert( links.size() >= 1 );
      size_t closest = 0;
      int sep2 = links.at(0).second.Sep();
      for ( size_t i = 1; i < links.size(); i++ ){
	if ( links.at(i).second.Sep() < sep2 ){
	  sep2    = links.at(i).second.Sep();
	  closest = i;
	}
      }
      
      Bool OK = True;
      int s2 = links.at(closest).first;
      for ( size_t i = 0; i < links.size(); i++ ){
	if ( i == closest ) continue;
	int s2x = links.at(i).first;
	int sep2x = links.at(i).second.Sep();
	int dev2x = links.at(i).second.Dev();
	
	const superb &S2 = scaffolds.at( s2/2 );
	const superb &S2x = scaffolds.at( s2x/2 );
	int s2_len = S2.SubSuperLength( 0, S2.Ntigs( ) - 1 );
	int s2_dev = S2.SubSuperLengthDev( 0, S2.Ntigs( ) - 1 );
	int s2x_len = S2x.SubSuperLength( 0, S2x.Ntigs( ) - 1 );
	int s2x_dev = S2x.SubSuperLengthDev( 0, S2x.Ntigs( ) - 1 );
	if ( sep2x + s2x_len - dev_mult * (s2x_dev + dev2x)
	     <= sep2 + s2_len){    
	  OK = False;
	  if ( VERBOSITY >= 2 ){ 
	    cout << "Outgoing conflict on " << s << " between "; PRINT2(s2,s2x);
	  }
	  outConflictCt++;
	  outConflictLoc++;
	  break; // one conflict is enough to stop
	}
      }
      if ( ! OK ){
	// remove connections
	if ( VERBOSITY >= 2 )
	  cout << " removing all connections from scaffold vertex " << s <<  endl;
	for ( size_t i = 0; i < links.size(); i++ ){
	  int s2 = links.at(i).first;
	  vec<int> linkEdgesFrom = G.EdgesBetween(s, s2);
	  ForceAssert( linkEdgesFrom.size() <= 1u );
	  //G.DeleteEdges( linkEdgesFrom );
	  if ( linkEdgesFrom.size() == 1 )
	    edgesToDel.push_back( linkEdgesFrom.at(0) );
	  delEdgeCt += linkEdgesFrom.size();
	  wantDelEdgeCt++;
	}
      }else{
	// remove all except the nearest
	if ( VERBOSITY >= 2 )
	  cout << " removing distant connections from scaffold vertex " << s <<  endl;
	for ( size_t i = 0; i < links.size(); i++ ){
	  if ( i == closest )
	    continue;
	  int s2 = links.at(i).first;
	  vec<int> linkEdgesFrom = G.EdgesBetween(s, s2);
	  ForceAssert( linkEdgesFrom.size() <= 1u );
	  if ( linkEdgesFrom.size() == 1 )
	    edgesToDel.push_back( linkEdgesFrom.at(0) );
	  //G.DeleteEdges( linkEdgesFrom );
	  delEdgeCt += linkEdgesFrom.size();
	  wantDelEdgeCt++;
	}
	
      } 
    }
  }

  G.DeleteEdges( edgesToDel );
  
  // Symmetrize: check for conflicts, chose closest link
  if ( VERBOSITY >= 2 ) 
    cout << Date() << " Cleaning up graph" << endl;    
  int inConflictCt = 0;
  for ( int s = 0; s < nscaff*2; s++ ){
    if ( G.To(s).size() > 1 ){
      vec< pair<int,sepdev> > links( G.To(s).size() );
      for ( size_t i = 0; i < G.To(s).size(); i++ ){
	int s0 = G.To(s).at(i);
	ForceAssert( G.EdgeObjectsBetween(s0,s).size() ==1 );
	sepdev sdloc( G.EdgeObjectsBetween(s0,s).at(0).Sep(), G.EdgeObjectsBetween(s0,s).at(0).Dev() );
	links[i] = pair<int,sepdev>(s0, sdloc); 
      }
	    
      // Find closest
      ForceAssert( links.size() >= 1 );
      size_t closest = 0;
      int sep0 = links.at(0).second.Sep( );
      for ( size_t i = 1; i < links.size(); i++ ){
	if ( links.at(i).second.Sep( ) < sep0 ){
	  sep0    = links.at(i).second.Sep( );
	  closest = i;
	}
      }
	    
      Bool OK = True;
      int s0 = links.at(closest).first;
      for ( size_t i = 0; i < links.size(); i++ ){
	if ( i == closest ) continue;
	int s0x = links[i].first;
	int sep0x = links[i].second.Sep( );
	int dev0x = links[i].second.Dev( );
	      
	const superb &S0 = scaffolds.at( s0/2 );
	const superb &S0x = scaffolds.at( s0x/2 );
	int s0_len = S0.SubSuperLength( 0, S0.Ntigs( ) - 1 );
	int s0_dev = S0.SubSuperLengthDev( 0, S0.Ntigs( ) - 1 );
	int s0x_len = S0x.SubSuperLength( 0, S0x.Ntigs( ) - 1 );
	int s0x_dev = S0x.SubSuperLengthDev( 0, S0x.Ntigs( ) - 1 );
	if ( sep0x + s0x_len - dev_mult * (s0x_dev + dev0x)
	     <= sep0 + s0_len){    
	  OK = False;
	  if ( VERBOSITY >= 2 ){ 
	    cout << "Incoming conflict on " << s << " between "; PRINT2(s0,s0x);
	  }
	  inConflictCt++;
	  break; // one conflict is enough to stop
	}
      }
      if ( ! OK ){
	// remove connections
	if ( VERBOSITY >= 2 )
	  cout << " removing all connections to scaffold vertex " << s <<  endl;
	for ( size_t i = 0; i < links.size(); i++ ){
	  int s0 = links[i].first;
	  vec<int> linkEdgesTo = G.EdgesBetween(s0, s);
	  ForceAssert( linkEdgesTo.size() == 1 );
	  G.DeleteEdges( linkEdgesTo );
	}
      }else{
	// remove all except the nearest
	if ( VERBOSITY >= 2 )
	  cout << " removing distant connections to scaffold vertex " << s <<  endl;
	for ( size_t i = 0; i < links.size(); i++ ){
	  if ( i == closest ) continue;
	  int s0 = links[i].first;
	  vec<int> linkEdgesTo = G.EdgesBetween(s0, s);
	  ForceAssertEq( 1u, linkEdgesTo.size() );
	  G.DeleteEdges( linkEdgesTo );
	}
	      
      } 
    }
  }
  G.RemoveDeadEdgeObjects();
  if ( VERBOSITY >= 2 )
    PRINT( inConflictCt );
  
}

// Generate a scaffold graph G. Notice scaffolds, aligns0, and
// aligns0_index will be modified.

void ScaffoldGraph( const PairsManager &pairs,
		    const vec<fastavector> &scaffold_contigs,
		    const int dev_mult,
		    const int min_pair_sep,		    
		    const int max_pair_sep,
		    const int min_links,
		    const int max_links,
		    const int min_scaffold_len,
		    const int MIN_LINKS_TO_PRINT,
		    const int VERBOSITY,
		    const size_t ntigs,
		    vec<superb> &scaffolds,
		    vec<alignlet> &aligns0,
		    vec<int> &aligns0_index,
		    digraphE<sepdev> &G,
		    int MAX_OVERLAP )
{
  // Define heuristic constants. These could be made into arguments.
  const int mean_jump_construct_size = 300; // mean size of sheared fragments

  // Build links between scaffolds.  There is one pass for each
  // orientation.  If s1 is a scaffold id, then fw_fw[s1] consists of
  // triples (s2, sep, dev) where s2 is another scaffold id such that
  // there are "enough" read pairs having one read forward on s1 and
  // the other forward on s2.  The construction is symmetric:
  // fw_fw[s2] includes (s1, sep, dev).
  //
  // Note that the link is actually from s1 to s2rc (or from s2
  // to s1rc).
  //
  // The vectors fw_rc, rc_fw, and rc_rc are similarly defined.
  longlong outConflictCt = 0, inConflictCt = 0;
  longlong propAlgndPairsCt = 0, pairsCt = 0;
  int nscaff = scaffolds.size( );
  
  vec< vec< triple<int,int,int> > >
    fw_fw(nscaff), fw_rc(nscaff), rc_fw(nscaff), rc_rc(nscaff);
  for ( int pass = 1; pass <= 2; pass++ ){
    // Reverse stuff.
    for ( int i = 0; i < nscaff; i++ )
      scaffolds[i].Reverse( );
    for ( size_t i = 0; i < aligns0.size(); i++ )
      aligns0[i].Reverse( );
    // Build indices that allow one find the scaffold and position on the
    // scaffold for a given contig.
	  
    vec<int> tig_to_scaffold( ntigs, -1 );
    vec<int> tig_to_scaffold_pos( ntigs, -1 );
    for ( int i = 0; i < nscaff; i++ ){    
      for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ ){    
	tig_to_scaffold[ scaffolds[i].Tig(j) ] = i;
	tig_to_scaffold_pos[ scaffolds[i].Tig(j) ] = j;    
      }    
    }
	  
    // Find forward links from scaffolds.
	  
    vec< vec<slink> > slinks(nscaff);
    for ( ulonglong pair_ID = 0; pair_ID < pairs.nPairs(); pair_ID++ ){
      longlong id1 = pairs.ID1( pair_ID );
      longlong id2 = pairs.ID2( pair_ID );
      // We only consider read pairs for which both ends are uniquely
      // placed.
      if ( aligns0_index[id1] < 0 || aligns0_index[id2] < 0 )
	continue;
      if ( pairs.sep(pair_ID) < min_pair_sep ||
	   pairs.sep(pair_ID) > max_pair_sep)
	continue;
      pairsCt++;
	    
      // Figure out what the link is.  Ignore intrascaffold links, as
      // well as those for which id1 is not forward on the contig.
	    
      const alignlet& la1 = aligns0[ aligns0_index[id1] ];
      if ( ! la1.Fw1() ) continue;
      const alignlet& la2 = aligns0[ aligns0_index[id2] ];
      int t1 = la1.TargetId( ), t2 = la2.TargetId( );
      int s1 = tig_to_scaffold[t1], s2 = tig_to_scaffold[t2];
      if ( s1 < 0 || s2 < 0 ) {
	continue;
	      
      }
      if ( s1 == s2 ) {
	// To speed up the next iteration, indexes were reset to -3 at
	//  this point. Keeping these aligns, however, allows us to
	//  iteratively run regapping code (sante - 2010.11.1).

	// aligns0_index[id1] = -3;  
	// aligns0_index[id2] = -3;
	continue;
      }
      
      int p1 = tig_to_scaffold_pos[t1], p2 = tig_to_scaffold_pos[t2];
      const superb &S1 = scaffolds[s1], &S2 = scaffolds[s2]; 
	    
      propAlgndPairsCt++;

      // Now do some math to decide the gap and its deviation, and
      // whether the gap is within bounds.  There are lots of
      // heuristics here.
      int dist_to_end1 = scaffold_contigs[t1].size() - la1.Pos2( )
	+ S1.SubSuperLength( p1+1, S1.Ntigs( ) - 1 );
      int dist_to_end1_dev
	= S1.SubSuperLengthDev( p1+1, S1.Ntigs( ) - 1 );

      int dist_to_oend1 = la1.Pos2( ) + S1.SubSuperLength( 0, p1 -1 );

      int dist_to_end2, dist_to_end2_dev;
      int dist_to_oend2;
      if ( la2.Fw1( ) ){    
	dist_to_end2 = scaffold_contigs[t2].size() - la2.Pos2( )
	  + S2.SubSuperLength( p2+1, S2.Ntigs( ) - 1 );
	dist_to_end2_dev
	  = S2.SubSuperLengthDev( p2+1, S2.Ntigs( ) - 1 );

	dist_to_oend2 = la2.Pos2( ) + S2.SubSuperLength( 0, p2 -1 );
      }
      else{    
	dist_to_end2 = la2.pos2( ) + S2.SubSuperLength( 0, p2-1 );
	dist_to_end2_dev = S2.SubSuperLengthDev( 0, p2-1 );
	      
	dist_to_oend2 = scaffold_contigs[t2].size() - la2.pos2( ) 
	  + S2.SubSuperLength( p2+1, S2.Ntigs( ) - 1 );
      }
	    
      // do not consider near other edge alignments, possibly avoid some 'innies' 
      // from sheared libraries
      if ( dist_to_oend1 + dist_to_oend2 < mean_jump_construct_size ) continue;
	    
      // disregard alignments that are too far away from edges, they
      // are likely to be just random noise. (They could be used to
      // indicate problems with scaffolds)

      if ( dist_to_end1 - dev_mult * dist_to_end1_dev > 
	   pairs.sep(pair_ID) + dev_mult * pairs.sd(pair_ID) ){
	continue;
      } 
      if ( dist_to_end2 - dev_mult * dist_to_end2_dev  > 
	   pairs.sep(pair_ID) + dev_mult * pairs.sd(pair_ID) ){
	continue;
      } 
	    
      int sep = pairs.sep( pair_ID ) - dist_to_end1 - dist_to_end2;
      int dev = dist_to_end1_dev + dist_to_end2_dev + pairs.sd( pair_ID );
      int min_gap = sep - dev_mult * dev;
      int max_gap = sep + dev_mult * dev;
	    
      if ( sep >= -MAX_OVERLAP )
	slinks[s1].push( s2, la2.Fw1( ), sep, dev );    
      else if ( VERBOSITY >= 2 ) 
	cout << "link "<< s1 << "->" << s2<< " " << "sep "<< sep << " removed." << endl;
    }
	  
    // Condense the links.  We look at all the links between two
    // scaffolds and decide if there are enough.
	  
    vec< vec<slink> > slinks2(nscaff), slinks2_to(nscaff);
    for ( int s1 = 0; s1 < nscaff; s1++ ){    
      Sort( slinks[s1] );

      for ( size_t i = 0; i < slinks[s1].size(); i++ ){    
	int s2 = slinks[s1][i].s2;
	Bool fw2 = True;
	size_t j;
	vec<double> seps, devs, seps_rc, devs_rc;
	for ( j = i; j < slinks[s1].size(); j++ ){ 
	  if ( slinks[s1][j].s2 != s2 ) break;
	  if (slinks[s1][j].fw2) {
	    seps.push_back( slinks[s1][j].sep );
	    devs.push_back( slinks[s1][j].dev );
	  } else {
	    seps_rc.push_back( slinks[s1][j].sep );
	    devs_rc.push_back( slinks[s1][j].dev );
	  }
	}
	if (seps.size() < seps_rc.size()) {
	  fw2 = False;
	  seps = seps_rc;
	  devs = devs_rc;
	}
	
	// Sort pairs (sep, dev), and cluster these into consistent
	// sets, where (s1,d1) is consistent with (s2,d2) iff the two
	// intervals [s1-d1, s1+d1) and [s2-d2, s2+d2) overlap.

	{
	  SortSync( seps, devs );
	  vec< vec<int> > clusters( 1 );
	  clusters[0].push_back( 0 );
	  for (int ii=1; ii<seps.isize( ); ii++) {
	    ho_interval prev( seps[ii-1] - devs[ii-1], seps[ii-1] + devs[ii-1] );
	    ho_interval curr( seps[ii] - devs[ii], seps[ii] + devs[ii] );
	    if ( Overlap( prev, curr ) > 0 ) {
	      clusters[ clusters.size( )-1 ].push_back( ii );
	      continue;
	    }
	    vec<int> new_clusters( 1, ii );
	    clusters.push_back( new_clusters );
	  }
	  int winner_id = 0;
	  for (int ii=1; ii<clusters.isize( ); ii++)
	    if ( clusters[ii].size( ) > clusters[winner_id].size( ) )
	      winner_id = ii;
	  const vec<int> &thec = clusters[winner_id];
	  vec<double> select_seps;
	  vec<double> select_devs;
	  select_seps.reserve( thec.size( ) );
	  select_devs.reserve( thec.size( ) );
	  for (int ii=0; ii<thec.isize( ); ii++) {
	    select_seps.push_back( seps[ thec[ii] ] );
	    select_devs.push_back( devs[ thec[ii] ] );
	  }
	  swap( select_seps, seps );
	  swap( select_devs, devs );
	}

	int count = (int)seps.size();
	if ( count >= MIN_LINKS_TO_PRINT ) {
	  PRINT3( s1, s2, count );
	  for (int ii=0; ii<seps.isize( ); ii++)
	    cout << " SANTEMP sep: " << seps[ii] << " +/- " << devs[ii] << "\n";
	}
	
	// Decide if there are enough links.
	
	if ( count < min_links || count > max_links ){
	  i = j - 1;
	  continue;    
	}
	
	// Skip links between short scffolds.

	if ( scaffolds[s1].TrueLength( ) < min_scaffold_len ) continue;
	if ( scaffolds[s2].TrueLength( ) < min_scaffold_len ) continue;
	
	// Convert to a single link.  The mechanism for defining
	// the separation and deviation is heuristic.
	      
	SortSync( devs, seps );
	double Sep, Dev;
	CombineStats( seps, devs, Sep, Dev );
	      
	// Save the link.

	slinks2[s1].push(s2, fw2, int(round(Sep)), int(round(Dev)) );
	// Re-using this class to record the pointers to each scaffold
	slinks2_to[s2].push(s1, fw2, int(round(Sep)), int(round(Dev)) );
	i = j - 1;    
      }
    }
    
    // Identify the cases where of the scaffolds linked forward to by a
    // given scaffold, one is clearly closest.
	  
    for ( int s1 = 0; s1 < nscaff; s1++ ){    
      if ( slinks2[s1].empty( ) ) continue;
      for ( size_t i = 0; i < slinks2[s1].size(); i++ ){
	int s2 = slinks2[s1][i].s2;
	int sep2 = slinks2[s1][i].sep;
	int dev2 = slinks2[s1][i].dev;
	Bool fw2 = slinks2[s1][i].fw2;
	if ( pass == 2 ){    
	  if (fw2) fw_fw[s1].push( s2, sep2, dev2 );
	  else fw_rc[s1].push( s2, sep2, dev2 );    
	} else{    
	  if (fw2) rc_rc[s1].push( s2, sep2, dev2 );
	  else rc_fw[s1].push( s2, sep2, dev2 );    
	} 
      }
    }     
    
  }
	
  if (VERBOSITY >= 2 )
    PRINT2( pairsCt, propAlgndPairsCt );
  
  if ( VERBOSITY >= 2 ){    
    for ( int i = 0; i < nscaff; i++ ){    
      for ( size_t j = 0; j < fw_fw[i].size(); j++ )
	cout << "fw_fw: " << i << " --> " << fw_fw[i][j] 
	     << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
	     << " " << scaffolds[fw_fw[i][j].first].SubSuperLength(0, scaffolds[fw_fw[i][j].first].Ntigs() -1 ) << "\n";
      for ( size_t j = 0; j < fw_rc[i].size(); j++ )
	cout << "fw_rc: " << i << " --> " << fw_rc[i][j] 
	     << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
	     << " " << scaffolds[fw_rc[i][j].first].SubSuperLength(0, scaffolds[fw_rc[i][j].first].Ntigs() -1 ) << "\n";
      for ( size_t j = 0; j < rc_fw[i].size(); j++ )
	cout << "rc_fw: " << i << " --> " << rc_fw[i][j] 
	     << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
	     << " " << scaffolds[rc_fw[i][j].first].SubSuperLength(0, scaffolds[rc_fw[i][j].first].Ntigs() -1 ) << "\n";
      for ( size_t j = 0; j < rc_rc[i].size(); j++ )
	cout << "rc_rc: " << i << " --> " << rc_rc[i][j] 
	     << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
	     << " " << scaffolds[rc_rc[i][j].first].SubSuperLength(0, scaffolds[rc_rc[i][j].first].Ntigs() -1 ) << "\n";   
    }  
  }
	
  if ( VERBOSITY >= 2 ) 
    PRINT( outConflictCt );
  // Make the scaffold links into a graph.  The vertices are 0fw, 0rc, 1fw,
  // 1rc, etc.  This is painful but purely mechanical.
  //if (VERBOSITY >=2 ) cout << Date() << " buiding graph" << endl;
  vec< vec<int> > from(2*nscaff), to(2*nscaff);
  vec< vec<int> > from_edge_obj(2*nscaff), to_edge_obj(2*nscaff);
  vec<sepdev> edges;
  for ( int i = 0; i < nscaff; i++ ){    
    for ( size_t j = 0; j < fw_fw[i].size(); j++ ){    
      int nedges = edges.size( );
      from[2*i].push_back( 2*fw_fw[i][j].first+1 );
      to[ 2*fw_fw[i][j].first+1 ].push_back(2*i);
      from_edge_obj[2*i].push_back( edges.size( ) );
      to_edge_obj[ 2*fw_fw[i][j].first+1 ].push_back(nedges);
      edges.push( fw_fw[i][j].second, fw_fw[i][j].third );    
    }
    for ( size_t j = 0; j < fw_rc[i].size(); j++ ){    
      int nedges = edges.size( );
      from[2*i].push_back( 2*fw_rc[i][j].first );
      to[ 2*fw_rc[i][j].first ].push_back(2*i);
      from_edge_obj[2*i].push_back( edges.size( ) );
      to_edge_obj[ 2*fw_rc[i][j].first ].push_back(nedges);
      edges.push( fw_rc[i][j].second, fw_rc[i][j].third );    
    }
    for ( size_t j = 0; j < rc_fw[i].size(); j++ ){    
      int nedges = edges.size( );
      from[2*i+1].push_back( 2*rc_fw[i][j].first+1 );
      to[ 2*rc_fw[i][j].first+1 ].push_back(2*i+1);
      from_edge_obj[2*i+1].push_back( edges.size( ) );
      to_edge_obj[ 2*rc_fw[i][j].first+1 ].push_back(nedges);
      edges.push( rc_fw[i][j].second, rc_fw[i][j].third );    
    }
    for ( size_t j = 0; j < rc_rc[i].size(); j++ ){    
      int nedges = edges.size( );
      from[2*i+1].push_back( 2*rc_rc[i][j].first );
      to[ 2*rc_rc[i][j].first ].push_back(2*i+1);
      from_edge_obj[2*i+1].push_back( edges.size( ) );
      to_edge_obj[ 2*rc_rc[i][j].first ].push_back(nedges);
      edges.push( rc_rc[i][j].second, rc_rc[i][j].third );    
    }    
  }

  for ( size_t s = 0; s < from.size(); s++ )
    SortSync( from[s], from_edge_obj[s] );
  
  G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
	
  if ( VERBOSITY >= 2 ){ 
    cout << "number of edges in the graph "; PRINT(edges.size());
  }

}

void MakeScaffoldsCloseBest( vec<superb> &scaffolds,
			     vec<fastavector> &scaffold_contigs,
			     vec<alignlet> &aligns0,
			     vec<int> &aligns0_index,
			     const PairsManager &pairs,
			     const String MIN_PAIR_SEPS,
			     const String MAX_PAIR_SEPS,
			     const String MIN_LINKS,
			     const String MAX_LINKS,
			     const int MAX_OVERLAP,
			     const int MIN_SCAFFOLD_LEN,
			     const Bool SUCK_SCAFFOLDS,
			     const Bool SHAVE_GRAPH,
			     const String SCAFFOLD_GRAPH_OUT,
			     const String SCAFFOLD_GRAPH_DOT,
			     const int VERBOSITY,
			     const int MIN_LINKS_TO_PRINT )
{
  // Define heuristic constants. These could be made into arguments.
  const int dev_mult = 3;
  const int pair_seps_bin_size = 100000; // large value (effectively no binning)
  
  // Read pair separation limits:  With long insert read libraries,
  // we may wish to consider them in sequence to find nearby joins
  // before matches which might jump over other scaffolds. Default
  // uses everthing in one iteration, as before. --bruce 30 Sep 09
  
  vec<int> min_links_set;
  ParseIntSet(MIN_LINKS, min_links_set, False);

  vec<int> min_pair_seps;
  vec<int> max_pair_seps;
  if ( ! MIN_PAIR_SEPS.empty() ){   
    // use separation limits from arguments
    ForceAssert( ! MAX_PAIR_SEPS.empty() );
    ParseIntSet(MIN_PAIR_SEPS, min_pair_seps, False);
    ParseIntSet(MAX_PAIR_SEPS, max_pair_seps, False);
  }else{
    // derive separation limits from library information
    vec<int> libID2Sep( pairs.nLibraries() );
    for (unsigned libID = 0; libID < pairs.nLibraries(); libID++ ){
      libID2Sep[libID] = pairs.getLibrarySep( libID );
    }
    vec<int> sortedLibSeps = libID2Sep; 
    UniqueSort( sortedLibSeps );
    int minimum_sep = sortedLibSeps[0] < 0 ? sortedLibSeps[0] : 0;
    if ( VERBOSITY >= 2 ){
      cout << "sortedLibSeps: "; sortedLibSeps.Print( cout ); cout << endl;
    }
    min_pair_seps.push_back( minimum_sep );
    max_pair_seps.push_back( sortedLibSeps[0] );
    int bin_beg = sortedLibSeps[0];
    for ( unsigned i = 1; i < sortedLibSeps.size(); i++ ){      
      if ( sortedLibSeps[i] > bin_beg + pair_seps_bin_size ){
	min_pair_seps.push_back( minimum_sep );
	max_pair_seps.push_back( sortedLibSeps[i] );
	bin_beg = sortedLibSeps[i];
      }else{
	min_pair_seps.back() = minimum_sep;
	max_pair_seps.back() = sortedLibSeps[i];
      }
    }
    ForceAssertEq( max_pair_seps.back(), sortedLibSeps.back() );
  }
  if ( VERBOSITY >= 1 ){
    cout << "min_pair_seps: "; min_pair_seps.Print(cout); cout << endl;
    cout << "max_pair_seps: "; max_pair_seps.Print(cout); cout << endl;
  }
  ForceAssertEq( min_pair_seps.size(), max_pair_seps.size() );
  
  const int infinite = std::numeric_limits<int>::max( );
  vec<int> max_links_set;
  if ( MAX_LINKS == "" ) max_links_set.resize( min_links_set.size( ), infinite );
  else ParseIntSet( MAX_LINKS, max_links_set, False );
  ForceAssertEq( min_links_set.size( ), max_links_set.size( ) );

  // Track orientations of contigs.
  
  size_t ntigs = scaffold_contigs.size();
  vec<Bool> rctig( ntigs, False );
  
  // Now begin the scaffold joining operation.  This is repeated until there
  // no more joins are made.
  
  if ( VERBOSITY >= 1 ){    
    cout << Date( ) << ": Joining scaffolds... " << endl;
    cout << "MAX_OVERLAP = " << MAX_OVERLAP << endl;
    cout << "\ninitial scaffolds:\n\n";
    for ( size_t i = 0; i < scaffolds.size(); i++ )
      scaffolds[i].Print( cout, "scaffold_" + ToString(i) );    
  }

  int iterationCt = 0, allJoinsCt = 0; // help track global progress

  int MAX_BREAK_ITER = 1;
  for ( int break_iter = 1; break_iter <= MAX_BREAK_ITER; break_iter++ ){    
    for ( unsigned mls = 0; mls < min_links_set.size(); mls++ ){ 
      size_t min_links_size = min_links_set.size();
      int min_links         = min_links_set[ mls ];
      int max_links         = max_links_set[ mls ];
      
      for (size_t mps = 0; mps < max_pair_seps.size(); ++mps){
	int min_pair_sep = min_pair_seps[mps];
	int max_pair_sep = max_pair_seps[mps];
	int joins = 0; int localJoinsCt =  0; // help track local progress
	
	if ( VERBOSITY >= 1 )
	  cout << Date() << ": Using min_links=" << min_links
	       << " separations=[" << min_pair_sep
	       << "," << max_pair_sep
	       << "]" << endl;
	do {
	  
	  int nscaff = scaffolds.size( );
	  iterationCt++;
	  if ( VERBOSITY >= 1 ) 
	    PRINT( iterationCt );
	  joins = 0;
	  
	  // Generate scaffold link graph.
	  
	  digraphE<sepdev> G;
	  digraphE<CLinkBundle> Gbundles;
	  ScaffoldGraph( pairs, scaffold_contigs, dev_mult, min_pair_sep,
			 max_pair_sep, min_links, max_links, MIN_SCAFFOLD_LEN,
			 MIN_LINKS_TO_PRINT, VERBOSITY, ntigs, scaffolds,
			 aligns0, aligns0_index, G, MAX_OVERLAP );
			 
	  if ( SCAFFOLD_GRAPH_DOT != "" ) 
	    dump_linkgraph(SCAFFOLD_GRAPH_DOT, G, iterationCt);
	  
	  if ( SHAVE_GRAPH )
	    ShaveGraph( G );
	  CleanupGraph( VERBOSITY, dev_mult, scaffolds, G );	  

	  if ( SCAFFOLD_GRAPH_OUT != "" ) {
	    String out_file = SCAFFOLD_GRAPH_OUT + "." + ToString(iterationCt);
	    int fd = OpenForWrite( out_file );
	    BinaryWrite( fd, G );
	  }
	  
	  // Now use the links to define lines in the scaffold graph.  There are
	  // heuristics embedded here.
	  
	  vec< vec<int> > lines;
	  vec<Bool> used( nscaff, False );
	  for ( int pass = 1; pass <= 2; pass++ ){    
	    for ( int s0 = 0; s0 < nscaff*2; s0++ ){
	      if ( VERBOSITY >= 2 ){
		if ( G.From(s0).size() > 1 )
		  cout << "s0: G.From(" << s0 << ").size() = " << G.From(s0).size() << endl;
		if ( G.To(s0).size() > 1 )
		  cout << "s0: G.To(" << s0 << ").size() = " << G.To(s0).size() << endl;
	      }
	      int s = s0;
	      if ( used[s/2] || ( pass == 1 && !G.Source(s) ) ) continue;
	      if ( pass == 1 && G.From(s0).size() < 1 ) continue;
	      vec<int> line;
	      line.push_back(s);
	      used[s/2] = True;
	      int extendedCt = 0;
	      while(1){    
		if ( VERBOSITY >= 2 ){
		  if ( G.From(s).size() > 1 )
		    cout << "G.From(" << s << ").size() = " << G.From(s).size() << endl;
		  if ( G.To(s).size() > 1 )
		    cout << "G.To(" << s << ").size() = " << G.To(s).size() << endl;
		}
		if ( G.From(s).size() != 1 ) break;
		s = G.From(s)[0];
		if ( G.To(s).size() > 1 || used[s/2] ) break;
		used[s/2] = True;
		line.push_back(s); 
		extendedCt++;
	      }
	      lines.push_back(line);    
	    }    
	  }
	  
	  // Merge the scaffolds in the lines to form new scaffolds.  This is a
	  // purely formal operation.
	  
	  vec<superb> new_scaffolds;
	  vec<Bool> reversed( ntigs, False );
	  if ( VERBOSITY >= 1 ){    
	    cout << "\nlines now:\n";
	    for ( size_t i = 0; i < lines.size(); i++ ){    
	      for ( size_t j = 0; j < lines[i].size(); j++ ){    
		if ( j > 0 ) cout << " --> ";
		cout << lines[i][j]/2
		     << ( lines[i][j] % 2 == 0 ? "fw" : "rc" ) << " (" << scaffolds[ lines[i][j]/2 ].SubSuperLength(0, scaffolds[ lines[i][j]/2 ].Ntigs() -1 ) << ")";    
	      }
	      cout << "\n";    
	    }    
	  }
	  for ( size_t i = 0; i < lines.size(); i++ ){    
	    joins += lines[i].size( ) - 1;
	    superb x;
	    int ntigs = 0;
	    for ( size_t j = 0; j < lines[i].size(); j++ )
	      ntigs += scaffolds[ lines[i][j]/2 ].Ntigs( );
	    x.SetNtigs(ntigs);
	    int count = 0;
	    for ( size_t j = 0; j < lines[i].size(); j++ ){ 
	      int y = lines[i][j];
	      int s = y/2;
	      superb S = scaffolds[s];
	      
	      // In the rc case, reverse S, and also track the change.
	      
	      if ( y % 2 == 1 ){    
		S.Reverse( );
		for ( int k = 0; k < S.Ntigs( ); k++ ){    
		  rctig[ S.Tig(k) ] = !rctig[ S.Tig(k) ];
		  reversed[ S.Tig(k) ] = True;    
		}    
	      }
	      
	      // Add S to the scaffold.
	      
	      for ( int k = 0; k < S.Ntigs( ); k++ ){    
		x.SetTig( count, S.Tig(k) );
		x.SetLen( count, S.Len(k) );
		if ( k < S.Ntigs( ) - 1 ){    
		  x.SetGap( count, S.Gap(k) );
		  x.SetDev( count, S.Dev(k) );
		  count++;    
		}    
	      }
	      
	      // Add the between scaffold gap.
	      
	      if ( j+1 < lines[i].size() ){    
		x.SetGap( count, G.EdgeObjectByIndexFrom( y, 0 ).Sep( ) );
		x.SetDev( count, G.EdgeObjectByIndexFrom( y, 0 ).Dev( ) );    
	      }
	      count++;    
	    }
	    new_scaffolds.push_back(x);    
	  }
	  
	  // Reverse the alignments lying over contigs that have just
	  // been reversed.
	  
	  for ( size_t i = 0; i < aligns0.size(); i++ )
	    if ( reversed[ aligns0[i].TargetId( ) ] ) aligns0[i].Reverse( );
	  
	  // Check to see if done.
	  
	  scaffolds = new_scaffolds;
	  if ( VERBOSITY >= 1 )
	    cout << Date( ) << ": made " << joins << " joins" << endl;
	  allJoinsCt += joins;
	  localJoinsCt += joins;
	} while (joins > 0);
	
	// Absorb small scaffolds into larger ones.
	
	if ( SUCK_SCAFFOLDS ) {
	  ofstream devnull ( "/dev/null" );
	  ostream *plog = ( VERBOSITY > 0 ? &cout : &devnull );
	  int n_sucked = SuckScaffolds( pairs, aligns0_index, aligns0,
					scaffolds, rctig, plog );
	  if ( VERBOSITY > 0 ) 
	    cout << Date( ) << ": " << n_sucked << " scaffolds sucked" << endl;
	}
	
      }
      
    }

  }
  if ( VERBOSITY >= 1 )
    cout << "total number of joins made = " << allJoinsCt << endl;
  
  // Flip rc contigs (aligns do not need to be flipped, since this was
  // done in the main loop).
  
  if ( VERBOSITY > 0 ) 
    cout << Date() << ": flipping rc contigs in the fastavector" << endl;
  for ( size_t i = 0; i < ntigs; i++ )
    if ( rctig[i] )
      scaffold_contigs[i].ReverseComplement();
  
  // Report some stats.
  
  if ( VERBOSITY > 0 )
    cout << "\nEdges in HyperFastavector: " << ntigs << endl;
  vec<int> contig_sizes, scaffold_sizes;
  for ( size_t i = 0; i < scaffolds.size(); i++ ){    
    for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ )
      contig_sizes.push_back(  scaffolds[i].Len(j) );
    scaffold_sizes.push_back( scaffolds[i].ReducedLength( ) );    
  }
  Sort(contig_sizes), Sort(scaffold_sizes);
  contig_sizes.EraseValue( 0 );
  scaffold_sizes.EraseValue( 0 );
  
  if ( VERBOSITY > 0 ) {
    cout << "N50 contig size = " << N50(contig_sizes) 
	 << " (" << contig_sizes.size() << " contigs)" << endl;
    cout << "N50 scaffold size (without gaps) = " << N50(scaffold_sizes) 
	 << " (" << scaffold_sizes.size() << " scaffolds)" << endl;
    cout << "total bases in contigs = " << BigSum(scaffold_sizes) << endl;
  }
  
}


