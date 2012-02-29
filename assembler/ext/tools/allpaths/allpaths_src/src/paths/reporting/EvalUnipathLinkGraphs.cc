/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/* EvalUnipathLinkGraphs
 *
 * Performs reference-based evaluation on a unipath link graph produced by
 * BuildUnipathLinkGraphs.  Currently only analyzes the unipath graph G -
 * not Gplus or FG.
 *
 * REQUIRED INPUT FILES:
 * /<PRE>/<DATA>/<RUN>/<READS>.unipaths.k<K>
 * /<PRE>/<DATA>/<RUN>/<READS>.unipaths.k<K>.locs
 * /<PRE>/<DATA>/<RUN>/<READS>.unipathgraph.k<K>
 *
 * REPORT:
 * 1. Link correctness, measured as number of standard devations off from truth
 * 2. Number and lengths of isolated unipaths (i.e., with no link)
 *
 * Josh Burton
 * December 2008
 *
 ******************************************************************************/

#include "MainTools.h"
#include "SeqInterval.h" // seq_interval
#include "math/NStatsTools.h" // PrintBasicNStats
#include "paths/KmerPath.h" // vecKmerPath
#include "paths/PdfEntry.h"
#include "paths/Unipath.h" // BuildUnipathAdjacencyGraph
#include "paths/Sepdev.h" // digraphE<fsepdev>
#include "paths/simulation/Placement.h" // placement




/* ReportBasicStats: print basic information about the unipath link graph.
 *
 ******************************************************************************/
void
ReportBasicStats( ostream & log,
		  const size_t n_edges,
		  const vecKmerPath & unipaths,
		  const vec< vec<int> > & component_vertices,
		  const vec<int> & unipath_lengths,
		  const String & paths_file,
		  const int K,
		  const bool find_connectivity )
{
  size_t n_components = component_vertices.size();
  
  
  // Calculate the total size (in kmers) of each component in the link graph
  vec<int> component_kmers( n_components, 0 );
  for ( size_t i = 0; i < n_components; i++ )
    for ( size_t j = 0; j < component_vertices[i].size(); j++ )
      component_kmers[i] += unipath_lengths[ component_vertices[i][j] ];
  
  // Order the component IDs in accordance with the component sizes
  vec<int> component_IDs;
  WhatPermutation( component_kmers, component_IDs, less<int>(), false );
  Sort( component_kmers );
  
  // Report on basic stats of LG
  log << "\n\n\n" << Date() << ": BASIC GRAPH STATS\n\n";
  log << "Graph LG has " << n_components << " components, with a total of " << SizeSum( component_vertices ) << " vertices (unipaths), " << n_edges << " edges (unipath links), and " << Sum( unipath_lengths ) << " kmers." << endl;
  log << "N50 component size (in unipaths) is " << N50_size( component_vertices ) << endl;
  log << "N50 component size (in kmers) is " << N50( component_kmers ) << endl;
  log << endl;
  
  
  
  
  if ( !find_connectivity ) return;
  
  // Find the connectivity of the largest component in the link graph.
  // We do this by gathering the set of unipaths that are in this link graph,
  // and making a unipath *adjacency* graph out of them.
  int c = component_IDs.back( ); // largest component ID
  
  cout << Date() << ": Loading stuff for component analysis..." << endl;
  vecKmerPath c_unipaths( component_vertices[c].size() );
  for ( size_t i = 0; i < component_vertices[c].size(); i++ )
    c_unipaths[i] = unipaths[ component_vertices[c][i] ];
  vec<tagged_rpint> pathsdb, c_unipathsdb;
  CreateDatabase( c_unipaths, c_unipathsdb );
  vecKmerPath paths   ( paths_file +    ".k" + ToString(K) );
  vecKmerPath paths_rc( paths_file + "_rc.k" + ToString(K) );
  CreateDatabase( paths, paths_rc, pathsdb );
  digraph AG;
  HyperKmerPath hkp;
  BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, c_unipaths, c_unipathsdb, AG );
  BuildUnipathAdjacencyHyperKmerPath( K, AG, c_unipaths, hkp );
  
  // Count the number of kmers in each component of the adjacency graph
  vec<vec<int> > hkp_edges;
  hkp.ComponentEdges( hkp_edges );
  vec<size_t> c_kmers( hkp_edges.size(), 1 );
  for ( size_t i = 0 ;  i < hkp_edges.size(); i++ )
    for ( size_t j = 0; j < hkp_edges[i].size(); j++ )
      c_kmers[i] += hkp.EdgeObject( hkp_edges[i][j] ).KmerCount( );
  Sort( c_kmers );
  
  
  // Report
  log << "Largest component (in kmers) is #" << c << " with " << component_vertices[c].size() << " unipaths and " << component_kmers.back( ) << " kmers." << endl;
  log << "Let's make a unipath *adjacency* graph out of this component." << endl;
  log << "This graph has " << c_kmers.size() << " components.  The N50 size in kmers is " << N50( c_kmers ) << ".  The largest has " << c_kmers.back() << " kmers." << endl;
}






/* ReportMissingLinks: Find the coverage of the reference genome by CN-1
 * unipaths, and investigate how many of these unipaths are linked together as
 * they ought to be.
 *
 ******************************************************************************/
void
ReportMissingLinks( ostream& log,
		    const int K,
		    const vecKmerPath & unipaths,
		    const digraphE<fsepdev> & LG,
		    const VecPlacementVec & unipath_POGs )
{
  const size_t n_unipaths = LG.N();
  
  // Find the set of all CN-ploidy unipaths.  (This works fine with diploid
  // genomes, assuming a haploid reference was used for alignment.)
  
  vec<int> CN1_unipaths_orig;
  vec<seq_interval> CN1_unilocs_orig;
  
  for ( size_t i = 0; i < n_unipaths; i++ ) {
    if ( unipath_POGs[i].size() != 1 ) continue;
    CN1_unipaths_orig.push_back( i );
    CN1_unilocs_orig .push_back( unipath_POGs[i][0].SeqInterval() );
  }
  
  // Sort these unipaths by location on the reference.  This brings together
  // the fw and rc versions of each unipath.
  SortSync( CN1_unilocs_orig, CN1_unipaths_orig );
  const size_t n_CN1_unipaths_orig = CN1_unipaths_orig.size();
  
  // Condense the fw and rc versions of each unipath.  (It would be simpler to
  // use UnipathInvolution for this, but then we'd have to load in unipathsdb.)
  // Also drop unipaths shorter than MIN_KMERS.
  vec< pair<int,int> > CN1_unipaths;
  vec<seq_interval> CN1_unilocs;
  const int MIN_KMERS = 100;
  size_t n_long_unipaths = 0;
  
  for ( size_t i = 0; i < CN1_unilocs_orig.size(); i++ ) {
    if ( CN1_unilocs_orig[i].Length() < MIN_KMERS + K - 1 ) continue;
    n_long_unipaths++;
    
    // Is this a second orientation of a unipath we've already seen?
    if ( i > 0 && CN1_unilocs_orig[i] == CN1_unilocs_orig[i-1] )
      CN1_unipaths.back().second = CN1_unipaths_orig[i];
    else {
      CN1_unipaths.push_back( make_pair( CN1_unipaths_orig[i], -1 ) );
      CN1_unilocs .push_back( CN1_unilocs_orig[i] );
    }
  }
  Destroy( CN1_unipaths_orig );
  Destroy( CN1_unilocs_orig );
  const size_t n_CN1_unipaths = CN1_unipaths.size();
  
  
  
  // Look for "linkable upairs" of unipaths.
  // A linkable upair of unipaths (u1,u2) must satisfy the following criteria:
  // -- u1, u2 are both at least MIN_KMERS.  Unipaths that fail this test have
  //    already been filtered out.
  // -- The distance between u1 and u2 is such that a if a jumping insert has
  //    size LINK_SIZE, the reads could fall on u1 and u2.
  const int LINK_SIZE = 2000;
  vec< pair<int,int> > linkable_upairs;
    
  for ( size_t u1 = 0; u1 < n_CN1_unipaths; u1++ ) {
    for ( size_t u2 = u1+1; u2 < n_CN1_unipaths; u2++ ) {
      
      // Find all values of j that fit the distance criterion.
      if ( CN1_unilocs[u1].SeqId() != CN1_unilocs[u2].SeqId() ) break;
      if ( CN1_unilocs[u1].End() + LINK_SIZE < CN1_unilocs[u2].Begin() ) break;
      if ( CN1_unilocs[u1].Begin() + LINK_SIZE > CN1_unilocs[u2].End() ) continue;
      
      // u1,u2 are indices into CN1_unipaths - not the unipath IDs themselves!
      linkable_upairs.push_back( make_pair( u1, u2 ) );
    }
  }
  
  // For each of the linkable upairs, see whether they are actually linked.
  size_t n_linked = 0;
  
  for ( size_t i = 0; i < linkable_upairs.size(); i++ ) {
    size_t u1 = linkable_upairs[i].first;
    size_t u2 = linkable_upairs[i].second;
    
    int u1a = CN1_unipaths[u1].first;
    int u1b = CN1_unipaths[u1].second;
    int u2a = CN1_unipaths[u2].first;
    int u2b = CN1_unipaths[u2].second;
    
    bool link = false;
    if ( LG.HasEdge( u1a, u2a ) || LG.HasEdge( u2a, u1a ) ) link = true;
    if ( u1b != -1 )
      if ( LG.HasEdge( u1b, u2a ) || LG.HasEdge( u2a, u1b ) ) link = true;
    if ( u2b != -1 )
      if ( LG.HasEdge( u1a, u2b ) || LG.HasEdge( u2b, u1a ) ) link = true;
    if ( u1b != -1 && u2b != -1 )
      if ( LG.HasEdge( u1b, u2b ) || LG.HasEdge( u2b, u1b ) ) link = true;
    
    if ( link ) n_linked++;
  }
  
  log << "\n\n\n" << Date() << ": MISSING LINKS\n\n";
  log << "Out of " << n_unipaths << " unipaths..." << endl;
  log << "\t..." << n_CN1_unipaths_orig << " have exactly one placement on the haploid reference" << endl;
  log << "\t..." << n_long_unipaths << " have more than MIN_KMERS=" << MIN_KMERS << " kmers" << endl;
  log << "\tFrom this set, rc unipaths were removed, leading to " << n_CN1_unipaths << " distinct linkable unipaths." << endl;
  log << "These " << n_CN1_unipaths << " linkable unipaths lead to " << linkable_upairs.size() << " linkable upairs, of which " << n_linked << " (" << double( 100.0 * n_linked ) / linkable_upairs.size() << "%) are actually linked in the graph." << endl;
  
}






/* ReportLinkCorrectness: evaluate all the links in the unipath link graph by
 * looking at the truth placements of their unipaths.
 *
 ******************************************************************************/
void
ReportLinkCorrectness( ostream& log,
		       const vecKmerPath & unipaths,
		       const digraphE<fsepdev> & LG,
		       const VecPlacementVec & unipath_POGs,
		       const int N_STDEVS )
{
  size_t n_edges = LG.EdgeObjectCount();
  
  vec<int> to_left, to_right;
  LG.ToLeft( to_left );
  LG.ToRight( to_right );
  
  
  // stdev_margins[N_STDEVS] is the number of links correct to within one stdev.
  // stdev_margins[N_STDEVS-i] is the number of links condensed by i stdevs.
  // stdev_margins[N_STDEVS+i] is the number of links stretched by i stdevs.
  vec<int> stdev_margins( 2*N_STDEVS, 0 );
  int stdev_overflow = 0;
  int no_truth = 0;
  
  // Loop over each link and gather info on the correctness of each one
  for ( size_t i = 0; i < n_edges; i++ ) {
    int v = to_left[i], w = to_right[i];
    
    PlacementVec v_places = unipath_POGs[v];
    PlacementVec w_places = unipath_POGs[w];
    
    // These unipaths must BOTH have truth placements, or we can't do anything.
    if ( v_places.size( ) == 0 || w_places.size( ) == 0 ) {
      no_truth++;
      continue;
    }
    
    // Get the true separation between unipaths v and w
    // If v and/or w are multiply placed, this finds the smallest separation
    int true_sep = numeric_limits<int>::max( );
    for ( PlacementVec::size_type vp = 0; vp < v_places.size( ); vp++ )
      for ( PlacementVec::size_type wp = 0; wp < w_places.size( ); wp++ )
	true_sep = Min( true_sep, Distance( v_places[vp].Interval( ), w_places[wp].Interval( ) ) );
    
    // Get the predicted separation (with stdev) between unipaths v and w
    int predicted_sep = LG.EdgeObject( i ).Sep( );
    int predicted_dev = LG.EdgeObject( i ).Dev( );
    // If the stdev of a link size is negative, something is VERY wrong.
    ForceAssert( predicted_dev >= 0 );
    
    // Measure the incorrectness of this link, in standard deviations 
    // (Note the implicit rounding down in the integer division)
    {
      int margin = predicted_dev == 0 ? INT_MAX :
	( predicted_sep - true_sep ) / predicted_dev;
      
      if ( Abs(margin) < N_STDEVS ) {
	int margin_bin = N_STDEVS + margin;
	stdev_margins[margin_bin]++;
      }
      else stdev_overflow++;
    }
          
  }
  
  double no_truth_pct = ( 100.0 * (double)no_truth       ) / n_edges;
  double overflow_pct = ( 100.0 * (double)stdev_overflow ) / n_edges;
  
  // Report on the correctness of links
  log << "\n\n\n" << Date() << ": LINK CORRECTNESS\n\n";
  log << "EvalUnipathLinkGraphs has compared the links in LG to the true\nseparations between the unipaths they link." << endl;
  log << "The errors are spread as follows:" << endl << endl;
  
  for ( int i = N_STDEVS-1; i > 0; i-- )
    log << "Links condensed by less than " << (i+1) << " standard deviations: " << stdev_margins[N_STDEVS-i] << endl;
  
  log << "Links correct to within +/- one standard deviation: " << stdev_margins[N_STDEVS] << endl;
  
  for ( int i = 1; i < N_STDEVS; i++ )
    log << "Links stretched by less than " << (i+1) << " standard deviations: " << stdev_margins[N_STDEVS+i] << endl;
  
  log << "Links incorrect by N_STDEVS=" << N_STDEVS << " standard deviations or more: " << stdev_overflow << " (" << overflow_pct << "%)" << endl;
  log << "Links between unipaths with no truth placement: " << no_truth << " (" << no_truth_pct << "%)" << endl;
  
}
  




/* ReportIsolatedUnipaths: analyze and print info about isolated unipaths,
 * e.g, unipaths with no links in the unipath link graph.
 *
 ******************************************************************************/
void
ReportIsolatedUnipaths( ostream & log,
			const vec< vec<int> > & component_edges,
			const vec< vec<int> > & component_vertices,
			const vec<int> & unipath_lengths,
			const VecPlacementVec & unipath_POGs )
{
  size_t n_unipaths = unipath_lengths.size();
  size_t n_components = component_vertices.size();
  
  
  // "Isos" = "isolated unipaths"
  vec<int> isos;
  isos.reserve( n_components );
  
  for ( size_t i = 0; i < n_components; i++ ) {
    
    // A component with no edges is an isolated unipath
    if ( component_edges   [i].size() == 0 &&
	 component_vertices[i].size() == 1 )
      isos.push_back( component_vertices[i][0] );
    
    
  }
  
  size_t n_isos = isos.size();
  
  
  // Calculate the sum of isolated unipaths' lengths
  vec<int> iso_lengths;
  iso_lengths.reserve( n_isos );
  for ( size_t i = 0; i < n_isos; i++ )
    iso_lengths.push_back( unipath_lengths[ isos[i] ] );
  
  // Find the copy number of isolated unipaths.
  static const uint MAX_CN = 10; // max copy number we track (for output only)
  vec<int> iso_CNs( MAX_CN + 1, 0 );
  for ( size_t i = 0; i < n_isos; i++ ) {
    uint CN = Min( unipath_POGs[ isos[i] ].size(), MAX_CN );
    iso_CNs[CN]++;
  }
  
  
  // Report!
  log << "\n\n\n" << Date() << ": ISOLATED UNIPATHS\n\n";
  log << "An 'isolated unipath' is defined as a unipath with no links." << endl;
  log << "Isolated unipaths are rejected from consideration as seeds." << endl;
  log << "\nThere are " << n_isos << " isolated unipaths ("
      << double( 100.0 * n_isos ) / n_unipaths
      << "% of total), each constituting its own graph compoment." << endl;
  log << "The total length of isolated unipaths is " << Sum( iso_lengths ) << " kmers." << endl << endl;
  log << "Length statistics for isolated unipaths (in kmers):" << endl;
  PrintBasicNStats( "Length", iso_lengths, log );
  log << endl;
  log << "True copy numbers of isolated unipaths (i.e., number of truth placements, with 0 indicating nongenomic):" << endl;
  for ( size_t i = 0; i < iso_CNs.size(); i++ )
    log << "CN-" << i << ( i == MAX_CN ? "+" : "" ) << ":\t" << iso_CNs[i] << endl;
  log << endl;
  
}








int main( int argc, char *argv[] )
{
  RunTime( );
  
  // Command-line arguments
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( GRAPH, "seeds" );
  CommandArgument_Int_OrDefault_Doc( N_STDEVS, 9, "Number of standard deviations to count separately" );
  CommandArgument_Int( K );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  
  EndCommandArguments;
  
  // Filenames
  String data_dir = PRE + "/" + DATA;
  String  run_dir = PRE + "/" + DATA + "/" + RUN;
  String  sub_dir = PRE + "/" + DATA + "/" + RUN + "/ASSEMBLIES/" + SUBDIR;
  String kK = ".k" + ToString( K );
  
  // Open an output file.
  String logfile = sub_dir + "/EvalUnipathLinkGraphs.out";
  Ofstream( log, logfile );
  command.PrintTheCommandPretty( log );
  
  
     
  cout << Date( ) << ": Beginning EvalUnipathLinkGraphs..." << endl;
  cout << Date( ) << ": Sending output to " << logfile << endl;
  
  
  /*****************************************************************************
   *         LOAD INPUT FILES
   ****************************************************************************/
  
  // Find ploidy.
  const int ploidy = FirstLineOfFile( data_dir + "/ploidy" ).Int( );
  
  // Unipaths
  vecKmerPath unipaths( run_dir + "/" + READS + ".unipaths" + kK );
  
  // Find unipath lengths
  size_t n_unipaths = unipaths.size( );
  vec<int> lengths( n_unipaths, 0 );
  for ( size_t i = 0; i < n_unipaths; i++ )
    lengths[i] = unipaths[i].KmerCount( );
  
  // Unipath link graph (currently, just the 'seeds' graph is analyzed)
  digraphE<fsepdev> LG;
  BinaryRead( sub_dir + "/unipath_link_graph." + GRAPH + kK, LG );
  
  /*
  // Unipath copy numbers (predicted)
  VecPdfEntryVec cp;
  cp.ReadAll( (run_dir + "/" + READS + ".unipaths.predicted_count" + kK).c_str() );
  
  // Mark unipaths as 'high_CN' if there is a CHANCE they have CN > ploidy
  vec<bool> high_CN( n_unipaths, false );
  for ( size_t i = 0; i < n_unipaths; i++ )
    for ( PdfEntryVec::size_type j = 0; j < cp[i].size( ); j++ )
      if ( cp[i][j].first > ploidy && cp[i][j].second > 0.1 )
	high_CN[i] = true;
  ForceAssertEq( n_unipaths, cp.size() );
  */
  
  
  // Unipath placements on genome (truth)
  VecPlacementVec unipath_POGs;
  unipath_POGs.ReadAll( run_dir + "/" + READS + ".unipaths" + kK + ".locs" );
  
  
  // Sanity checks!  If any of these asserts fails, the files are inconsistent.
  ForceAssertEq( n_unipaths, (size_t) LG.N() );
  ForceAssertEq( n_unipaths, unipath_POGs.size() );
  
  
  // Gather info on LG
  vec<vec<int> > component_vertices;
  vec<vec<int> > component_edges;
  LG.Components( component_vertices );
  LG.ComponentEdges( component_edges );
  size_t n_edges = LG.EdgeObjectCount();
  
  
  /*****************************************************************************
   *         INVESTIGATE MISSING LINKS
   ****************************************************************************/
  ReportMissingLinks( log, K, unipaths, LG, unipath_POGs );
  
  
  
  /*****************************************************************************
   *         BASIC GRAPH STATS
   ****************************************************************************/
  bool find_link_connectivity = true;
  String paths_file = run_dir + "/" + READS + ".paths";
  ReportBasicStats( log, n_edges, unipaths, component_vertices, lengths, paths_file, K, find_link_connectivity );
  
  
  
  /*****************************************************************************
   *         FIND LINK CORRECTNESS
   ****************************************************************************/
  ReportLinkCorrectness( log, unipaths, LG, unipath_POGs, N_STDEVS );
  
  
  
  /*****************************************************************************
   *         INVESTIGATE ISOLATED UNIPATHS
   ****************************************************************************/
  ReportIsolatedUnipaths( log, component_edges, component_vertices, lengths, unipath_POGs );
  
  
  
  cout << Date( ) << ": Done!" << endl;
  return 0;
}
