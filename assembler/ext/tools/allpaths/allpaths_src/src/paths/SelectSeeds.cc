///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* SelectSeeds
 *
 * Select seed unipaths.  A seed is a unipath that will eventually be used to
 * create a localized assembly, by finding nearby unipaths (nearby on the
 * unipath link graph) and the reads that align to those unipaths.  The reads
 * and unipaths that get pulled in are called the "cloud" that is created by
 * the seed.
 *
 * Criteria for seed selection include:
 * -- Copy number: A unipath with predicted copy number > 1 shouldn't be a seed.
 * -- Length: A unipath with length < MIN_KMERS_IN_SEED will not be chosen.
 * -- Links in graph: A unipath will not be chosen if it is isolated in the
 *    unipath link grapoh.
 * -- Spacing: If seed unipaths appear to be less than SEED_DIST apart on the
 *    unipath graph, they are pruned out.
 * -- Orientation: If a reference genome is supplied, and if USE_TRUTH=True,
 *    unipaths that align in RC to the reference are screened out.
 *
 * Seed unipaths are stored in a vec<int> and are written to the file seeds.ids
 * in SUBDIR.  This file is loaded by LocalizeReadsLG.
 *
 *
 * Josh Burton
 * December 2008
 *
 ******************************************************************************/



#include "MainTools.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "graph/Digraph.h" // digraph
#include "paths/PdfEntry.h" // pdf_entry
#include "paths/KmerPath.h" // vecKmerPath
#include "paths/Unipath.h" // UnipathInvolution
#include "paths/UnipathNhoodLG.h" // sepdev, fsepdev
#include "paths/FindUnipathSeedsLG.h" // FindUnipathSeeds, EvalUnipathSeeds
#include "paths/simulation/Placement.h" // placement











int main( int argc, char *argv[] )
{
  RunTime( );
  
  
  // Command-line arguments
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_Int( K );
  CommandArgument_String_OrDefault( GRAPH, "unipath_link_graph" );
  CommandArgument_Bool_OrDefault_Doc( USE_THEO_GRAPHS, False, "Use the file unipath_link_graph.theo, created by TheoreticalUnipathLinkGraph" );
  CommandArgument_Bool_OrDefault( EVAL, True );
  
  // Heuristic parameters
  CommandArgument_Int_OrDefault( SEED_DIST, 4000 );
  CommandArgument_Int_OrDefault( MIN_KMERS_IN_SEED, 100 );
  CommandArgument_Int_OrDefault( MIN_KMERS_IN_NEIGHBORS, 1000 );
  CommandArgument_Bool_OrDefault_Doc(CHEAT_RC, False, 
       "cheat by removing seeds that are rc to the reference");
  CommandArgument_Bool_OrDefault(WRITE_RESULTS, True);  
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;
  
  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS); 

  cout << Date( ) << ": Beginning SelectSeeds..." << endl;
   
  
  /*****************************************************************************
   *
   *        LOAD INPUT FILES
   *
   ****************************************************************************/
  
  cout << Date( ) << ": Loading input files" << endl;
  
  // Filenames
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String file_head = run_dir + "/" + READS;
  String kK = ".k" + ToString( K );
  String data_dir = PRE + "/" + DATA;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  // Read paths
  vecKmerPath paths   ( file_head + ".paths" + kK );
  vecKmerPath paths_rc( file_head + ".paths_rc" + kK );
  vec<tagged_rpint> pathsdb;
  BinaryRead3(file_head + ".pathsdb" + kK, pathsdb);
  
  // Calculate read lengths
  vec<int> read_lengths = KmerPathSeqLength( paths, K );
  
  // Unipaths [db]
  vecKmerPath unipaths( file_head + ".unipaths" + kK );
  vec<tagged_rpint> unipathsdb;
  BinaryRead3(file_head + ".unipathsdb" + kK, unipathsdb);
  
  // Unipath link graph
  String graph_file_infix = USE_THEO_GRAPHS ? ".theo" : "";
  digraphE<fsepdev> LG;
  BinaryRead( sub_dir + "/" + GRAPH + graph_file_infix + ".cloud" + kK, LG );
  
  // Unipath copy numbers (predicted)
  VecPdfEntryVec CNs( (file_head + ".unipaths.predicted_count" + kK).c_str() );
  
  // Ploidy
  const int ploidy = FirstLineOfFile( data_dir + "/ploidy" ).Int( );
  
  // Variables containing genome information - may be empty!
  // These are filled only if needed by EVAL or CHEAT_RC.
  vecbasevector genome;
  VecPlacementVec unipath_POGs; // unipath Placements On Genome
  
  if ( EVAL )
    genome.ReadAll( data_dir + "/genome.fastb" );
  
  if ( EVAL || CHEAT_RC )
    unipath_POGs.ReadAll( file_head + ".unipaths" + kK + ".locs" );

  // Sanity checks!
  // If any of these asserts fails, the input files are inconsistent.
  size_t n_unipaths = unipaths.size();
  ForceAssertEq( paths.size(), paths_rc.size() );
  ForceAssertEq( paths.size(), read_lengths.size() );
  ForceAssertEq( n_unipaths, (size_t) LG.N() );
  ForceAssertEq( n_unipaths, CNs.size() );
  if ( EVAL || CHEAT_RC ) ForceAssertEq( n_unipaths, unipath_POGs.size() );
  
  cout << "\t\tDataset contains " << paths.size() << " reads and " << n_unipaths << " unipaths." << endl;
  
  
  
  /*****************************************************************************
   *
   *        DATASET PRE-PROCESSING: Build some useful auxiliary data structures
   *
   ****************************************************************************/
  
  cout << Date( ) << ": Pre-processing dataset" << endl;
  
  // Find unipath lengths
  vec<int> unipath_lengths( n_unipaths );
  for ( size_t i = 0; i < n_unipaths; i++ )
    unipath_lengths[i] = unipaths[i].KmerCount( );
  
  
  // Read the CNs and find the predicted copy number (CN) of each unipath
  vec<int> predicted_CNs( n_unipaths, -1 );
  
  for ( size_t i = 0; i < n_unipaths; i++ )
    GetMostLikelyValue( predicted_CNs[i], CNs[i] );
  
  
  // Find unipath involution
  vec<int> to_rc;
  UnipathInvolution( unipaths, unipathsdb, to_rc );

  // Find unipaths that are branches in the graph.
  vec<Bool> branches( n_unipaths );
  {
       digraph AG;
       BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, 
            AG );
       for ( size_t i = 0; i < n_unipaths; i++ )
            branches[i] = ( AG.From( i ).size( ) > 1 || AG.To( i ).size( ) > 1 );   }

  // Unilocs (i.e., read locations on unipaths) and index.

  vec<ReadLocationLG> unilocs;
  BinaryRead2( file_head + ".unilocs." + ToString( K ), unilocs );
  vec<longlong> unilocs_index( n_unipaths + 1, -1 );
  unilocs_index[0] = 0;
  for ( int i = 0; i < unilocs.isize( ); i++ )
       unilocs_index[ unilocs[i].Contig( ) + 1 ] = i + 1;
  for ( size_t i = 1; i <= n_unipaths; i++ )
       if ( unilocs_index[i] < 0 ) unilocs_index[i] = unilocs_index[i-1];

  // Read pairing information
  PairsManager pairs( file_head + ".pairs" );
  pairs.makeCache( );
  
  /*****************************************************************************
   *
   *        SELECT SEED UNIPATHS
   *
   ****************************************************************************/
  
  
  // Initialize the logfile and write some header information to it.
  ofstream seed_log( (sub_dir + "/SelectSeeds.log").c_str() );
  cout << "Logging to " << sub_dir << "/SelectSeeds.log" << endl;
  seed_log << "# SeedStatus.log" << endl;
  seed_log << "#" << endl;
  seed_log << "# This file was automatically generated by the following command:" << endl;
  command.PrintTheCommandPretty( seed_log, "# " );
  seed_log << "# Total number of unipaths: " << n_unipaths << endl;
  seed_log << "#" << endl;
  seed_log << "# Fields:" << endl;
  seed_log << "# Unipath_ID      Seed_status      Reason" << endl;
  
  
  
  // Select seed unipaths
  cout << Date( ) << ": Selecting seed unipaths..." << endl;
  
  vec<int> seeds;
  vec<SeedStatus> seed_status;
  FindUnipathSeeds( K, ploidy, MIN_KMERS_IN_SEED, MIN_KMERS_IN_NEIGHBORS,
		    unipath_lengths, to_rc,
		    predicted_CNs, LG, CHEAT_RC, unipath_POGs, SEED_DIST,
		    branches, &paths,
		    &paths_rc, &pathsdb, &read_lengths, &unilocs,
                    &unilocs_index, &pairs, seeds, seed_status, seed_log,
                    NUM_THREADS );

  cout << Date( ) << ": Done!  " << seeds.isize( ) << " seeds selected."
       << endl << endl;
  
  if ( EVAL )
    EvalUnipathSeeds( seeds, unipath_lengths, genome, unipath_POGs,
		      SEED_DIST, predicted_CNs, LG, seed_status );
  
  // Write output

  if (WRITE_RESULTS) {
    cout << "Writing seeds to file " << sub_dir << "/seeds.ids" << endl;
    WRITE( sub_dir + "/seeds.ids", seeds );
    // Force LocalizeReadsLG to process this newly created list of seeds. 
    Remove( sub_dir + "/seeds.status" );
  }
  
  cout << Date( ) << ": Done!" << endl;
}
