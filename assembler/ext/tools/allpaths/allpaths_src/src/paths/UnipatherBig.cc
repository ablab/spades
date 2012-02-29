/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Unipather.  Consider the reads, always including the reverse complements in the
// set.  We assume that the reads have no gaps.  Divide all kmers in the read set 
// into kmer paths ("unipaths"), such that:
// - each kmer is in exactly one unipath;
// - if two kmers are adjacent in a unipath then they are adjacent in some read; 
// - unipaths contain no branch points (meaning that if two kmers xy are adjacent
//   in a unipath, then we never have xy' or x'y in a read for some other kmers 
//   x' or y');
// - otherwise, unipaths are maximal.
//
// This construction appears closely related to the construction of unitigs 
// described in Myers et. al. 2000.
//
// Generates: reads.unipaths.k*, reads.unipathsdb.k*, where * is K.
//
// The logic of how the unipather works could probably be improved substantially,
// so as to speed it up.  In particular, two things to look at are 
// (a) redundancy: the code is designed to rebuild a given unipath over and over,
// discarding duplicates on the fly;
// (b) multiple searches of pathsdb, some of which are unnecessary.
//
// Complication.  At least in artificial data sets, there can be circles without
// branch points.  These are arbitrarily broken.  Examples of this type have not
// been tested with VALIDATE=True.

#include "MainTools.h"
#include "graph/Digraph.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

#include "paths/Unipath.h"  // Unipath.h is templatized over tagged_rpint type

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS, "reads");
  CommandArgument_Bool_OrDefault_Doc(BUILD_UNIBASES, True, "Build and write unibases");
  CommandArgument_Bool_OrDefault_Doc(BUILD_UNIGRAPH, True, "Build and write the unipath adjacency graph");
  EndCommandArguments;
  
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String head = run_dir + "/" + READS;
  String KS = ToString(K);
  
  // Read in data.
  
  cout << Date( ) << ": loading data" << endl;
  vecKmerPath paths( head + ".paths.k" + KS );
  vecKmerPath paths_rc( head + ".paths_rc.k" + KS );
  BREAD2( head + ".pathsdb_big.k" + KS, vec<big_tagged_rpint>, pathsdb );
  
  // Build unipaths.
  
  cout << Date( ) << ": building unipaths (with big_tagged_rpint)" << endl;
  vecKmerPath unipaths;
  vec<big_tagged_rpint> unipathsdb;
  String unipaths_file = head + ".unipaths.k" + KS;
  Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb, True, unipaths_file );
  BinaryWrite3( head + ".unipathsdb_big.k" + KS, unipathsdb );
  size_t n_unipaths = unipaths.size();
  cout << "Found " << n_unipaths << " unipaths.\n";
  
  
  // Build unipath adjacency graph.
  if (BUILD_UNIGRAPH) {
    cout << Date() << ": building unipath adjacency graph" << endl;
    digraph unigraph;
    BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, unigraph );
    BinaryWrite( run_dir + "/" + READS + ".unipath_adjgraph.k" + KS, unigraph );
  }
  Destroy( unipathsdb );
  
  
  
  // Build unibases.
  // This requires loading the reads' bases into a KmerBaseBroker, which uses a
  // lot of memory.
  if (BUILD_UNIBASES) {
    cout << Date() << ": building KmerBaseBroker (needed for unibases)" << endl;
    KmerBaseBrokerBig kbb( run_dir, K, paths, paths_rc, pathsdb, READS );
    
    cout << Date() << ": building unibases" << endl;
    vecbasevector unibases;
    unibases.reserve(n_unipaths);
    for ( size_t i = 0; i < n_unipaths; i++ )
      unibases.push_back( kbb.Seq( unipaths[i] ) );
    
    unibases.WriteAll( head + ".unibases.k" + KS );
  }
  
  
  cout << Date() << ": Done with UnipatherBig!" << endl;
  return 0;
}
