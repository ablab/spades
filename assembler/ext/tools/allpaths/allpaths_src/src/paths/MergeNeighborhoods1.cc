///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



/* MergeNeighborhoods1
 *
 *
 * 1. Load in the HyperBasevectors produced by LocalizeReads3G.
 * 2. Combine the HyperBasevectors with the unibases and path them all together.
 * 3. Convert the HyperBasevectors to HyperKmerPaths in the same kmer space.
 * 4. Save the HyperKmerPaths.  Also save the reads/paths for the future use of
 * KmerBaseBrokers.
 *
 * All of these steps are purely mechanical.  There are no assembly algorithms
 * here, only data-processing and conversion algorithms.
 *
 *
 * INPUTS
 * -- Local neighborhood HyperBasevectors (SUBDIR/seed/.../?.hbv)
 * -- Unibases (RUN/<READS>.unibases)
 * OUTPUTS
 * -- Local neighborhood HyperKmerPaths (SUBDIR/nhood.hypers)
 * -- Global reads/paths (SUBDIR/reads.{fastb,paths,paths_rc,pathsdb})
 *
 *
 ******************************************************************************/


#include "MainTools.h"
#include "ParseSet.h"
#include "feudal/BinaryStream.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/ReadHBVs.h"
#include "paths/ReadsToPathsCoreX.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_Int( K );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault( SEEDS, "" );
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  cout << Date( ) << ": Beginning MergeNeighborhoods..." << endl;

  /*****************************************************************************
   *
   *        LOAD INPUT FILES
   *
   ****************************************************************************/

  cout << Date( ) << ": Loading input files" << endl;

  // Filenames
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String kK = ".k" + ToString( K );
  String unibases_file = run_dir + "/" + READS + ".unibases" + kK;
  String data_dir = PRE + "/" + DATA;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  // Only load selected seeds.
  int n_seeds = FirstLineOfFile( sub_dir + "/seeds.ids" ).Int( );
  vec<bool> selected( n_seeds, ( SEEDS == "" ? True : False ) );
  size_t n_selected = n_seeds;
  if ( SEEDS != "" ) {
    vec<int> sel;
    ParseIntSet( SEEDS, sel );
    cout << Date( ) << ": selected "
         << sel.size( ) << " seeds (out of a set of "
         << n_seeds << ")"
         << endl;
    n_selected = sel.size();
    for (size_t ii=0; ii<sel.size( ) ; ii++)
      selected[ sel[ii] ] = True;
  }

  // HyperBasevectors representing local assemblies
  cout << Date( ) << ": Loading " << n_selected
       << " HyperBasevectors from files in " << sub_dir << "..." << endl;

  vec<HyperBasevector> local_HBVs;
  readHBVs(sub_dir,selected,&local_HBVs);
  // Note:  unpopulated seeds have hbvs with K == 0.

  // Make a vecbasevector out of all of the basevectors in all of the local
  // HyperBasevector objects.

  cout << Date( ) << ": Creating a vecbasevector for the entire assembly" << endl;
  vecbasevector bases;
  for ( size_t i = 0; i < local_HBVs.size(); i++ )
      if ( local_HBVs[i].K() )
          for ( int j = 0; j < local_HBVs[i].EdgeObjectCount( ); j++ )
              bases.push_back( local_HBVs[i].EdgeObject(j) );


  // Add to this basevector all of the global unibases.  Now it contains the
  // entire dataset of basevectors known to be a part of the assembly.
  // (The "true" flag means append to bases.)
  uint64_t n_bases = bases.size();
  bases.ReadAll( unibases_file, true );

  // Re-path.  This is a runtime bottleneck.

  cout << Date( ) << ": Pathing this vecbasevector" << endl;
  vecKmerPath spaths, spaths_rc;
  vec<tagged_rpint> spathsdb;
  ReadsToPathsCoreY(bases, K, spaths, spaths_rc, spathsdb,
                    run_dir + "/MergeNeighborhoods", NUM_THREADS );

  // Write the bases and paths to file.  EvalHyper needs this data, and it's
  // smarter to output it here than to force EvalHyper to re-path everything.

  cout << Date( ) << ": Writing bases/paths to sub_dir" << endl;
  bases.WriteAll( sub_dir + "/reads.fastb" );
  spaths.WriteAll( sub_dir + "/reads.paths" + kK );
  spaths_rc.WriteAll( sub_dir + "/reads.paths_rc" + kK );
  BinaryWrite2( sub_dir + "/reads.pathsdb" + kK, spathsdb );

  // Write Unipaths
  bases.WriteRange( sub_dir + "/reads.unibases" + kK, n_bases, bases.size() );
  spaths.WriteRange( sub_dir + "/reads.unipaths" + kK, n_bases, bases.size() );

  // Clear out some data structures we no longer need.
  Destroy( bases );
  Destroy( spaths_rc );
  Destroy( spathsdb );

  // Convert the local HyperBasevectors into HyperKmerPaths that use the
  // newly generated kmer-numbering system.

  cout << Date( ) << ": Converting local HyperBasevectors to HyperKmerPaths" << endl;
  int count = 0;
  vec<HyperKmerPath> local_HKPs;
  for ( size_t i = 0; i < local_HBVs.size(); i++ )
  {
      if ( local_HBVs[i].K() )
      {
          vec<KmerPath> these_paths(0);
          for ( int j = 0; j < local_HBVs[i].EdgeObjectCount(); j++ )
              these_paths.push_back(spaths[count++]);
          HyperKmerPath h(K, local_HBVs[i], these_paths);
          local_HKPs.push_back(h);
      }
  }
  BinaryOverwrite( sub_dir + "/nhood.hypers", local_HKPs );
  cout << Date( ) << ": done!" << endl;
}
