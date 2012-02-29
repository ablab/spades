/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: Unipather

   
   Consider the reads, always including the reverse complements in the
   set.  We assume that the reads have no gaps.  Divide all kmers in the read set 
   into kmer paths ("unipaths"), such that:
   - each kmer is in exactly one unipath;
   - if two kmers are adjacent in a unipath then they are adjacent in some read; 
   - unipaths contain no branch points (meaning that if two kmers xy are adjacent
     in a unipath, then we never have xy' or x'y in a read for some other kmers 
     x' or y');
   - otherwise, unipaths are maximal.
  
   This construction appears closely related to the construction of unitigs 
   described in Myers et. al. 2000.
  
   Generates: reads.unipaths.k*, reads.unipathsdb.k*, where * is K.
  
   The logic of how the unipather works could probably be improved substantially,
   so as to speed it up.  In particular, two things to look at are 
   (a) redundancy: the code is designed to rebuild a given unipath over and over,
   discarding duplicates on the fly;
   (b) multiple searches of pathsdb, some of which are unnecessary.
  
   Complication.  At least in artificial data sets, there can be circles without
   branch points.  These are arbitrarily broken.  Examples of this type have not
   been tested with VALIDATE=True.

   @file
*/

#include "Alignment.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "MainTools.h"
#include "paths/ImproperMerge.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "paths/simulation/Placement.h"
#include "paths/UnibaseUtils.h"


static inline 
String Tag(String S = "UP") { return Date() + " (" + S + "): "; } 

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS, "reads");
  CommandArgument_String_OrDefault(PATHS, "paths");
  CommandArgument_String_OrDefault_Doc(UNIPATHS, "unipaths",
	 "Output unipaths name - READS.UNIPATHS.kK");
  CommandArgument_String_OrDefault_Doc(UNIBASES, "unibases",
	 "Output unibases name - READS.UNIBASES.kK");
  CommandArgument_String_OrDefault_Doc(UNIGRAPH, "unipath_adjgraph",
	 "Output unipath adjacency graph name - READS.UNIGRAPH.kK");
  CommandArgument_String_OrDefault_Doc(UNIBASE_GRAPH, "ovrlp_adjgraph",
	 "Output adjacency graph name - READS.UNIBASES.UNIBASE_GRAPH.kK");
  CommandArgument_Bool_OrDefault_Doc(BUILD_UNIBASES, True,
	 "Build and write unibases");
  CommandArgument_Bool_OrDefault_Doc(BUILD_UNIGRAPH, True,
	 "Build and write the unipath adjacency graph");
  CommandArgument_Bool_OrDefault_Doc(BUILD_UNIBASE_OVRLP_GRAPH, False,
	 "Build and write unibase adjacency graph based purely on overlapping unibase data");
  EndCommandArguments;
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  
  String KS = ToString(K);
  
  // Read in data.
  
  cout << Tag() << "loading data" << endl;
  vecKmerPath paths( run_dir + "/" + READS + "." + PATHS + ".k" + ToString(K) );
  vecKmerPath paths_rc( run_dir + "/" + READS + "." + PATHS + "_rc.k" + ToString(K) );    
  BREAD2( run_dir + "/" + READS + "." + PATHS + "db.k" + KS, vec<tagged_rpint>, pathsdb );
  
  // Build unipaths.
  
  cout << Tag() << "building unipaths" << endl;
  vecKmerPath unipaths;
  {
    vec<tagged_rpint> unipathsdb;
    const String unipaths_file = run_dir + "/" + READS + "." + UNIPATHS + ".k" + KS;
    Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb, True, unipaths_file );
    
    cout << Tag() << "writing unipaths db" << endl;
    BinaryWrite3( run_dir + "/" + READS + "." + UNIPATHS + "db.k" + KS, unipathsdb );    

  
    // Output some statistics on unipath sizes.
    {
      const size_t n_unipaths = unipaths.size();
      size_t n_KPIs = 0;
      vec<size_t> ulen(n_unipaths);
      for ( size_t i = 0; i < n_unipaths; i++ ) {
        n_KPIs += unipaths[i].NSegments();
        ulen[i] = unipaths[i].KmerCount();
      }
      Sort(ulen);
    
      cout << Tag() << setw(14) << ToStringAddCommas(n_unipaths) << "  unipaths found." << endl;
      cout << Tag() << setw(14) << ToStringAddCommas(n_KPIs) << "  kmer path intervals." << endl;
      cout << Tag() << setw(14) << ToStringAddCommas(Sum(ulen)) << "  total unipath length." << endl;
      cout << Tag() << setw(14) << ToStringAddCommas(N50(ulen)) << "  N50 unipath length." << endl;
    }
  
    // Build unipath adjacency graph.
    if (BUILD_UNIGRAPH) {
      cout << Tag() << "building unipath adjacency graph" << endl;
      digraph unigraph;
      BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, unigraph );

      cout << Tag() << "writing unipath adjacency graph" << endl;
      BinaryWrite( run_dir + "/" + READS + "." + UNIGRAPH + ".k" + KS, unigraph );
    }
  }
  
  
  // Build unibases.
  // This requires loading the reads' bases into a KmerBaseBroker, which uses a
  // lot of memory.
  if ( BUILD_UNIBASES || BUILD_UNIBASE_OVRLP_GRAPH ) {
    const size_t n_unipaths = unipaths.size();

    cout << Tag() << "building KmerBaseBroker (needed for unibases)" << endl;
    KmerBaseBroker kbb( run_dir, K, paths, paths_rc, pathsdb, READS );
    
    cout << Tag() << "building unibases" << endl;
    vecbasevector unibases;
    unibases.reserve(n_unipaths);
    for ( size_t i = 0; i < n_unipaths; i++ )
      unibases.push_back( kbb.Seq( unipaths[i] ) );
    
    cout << Tag() << "writing unibases" << endl;
    unibases.WriteAll( run_dir + "/" + READS + "." + UNIBASES + ".k" + KS );
    
    if ( BUILD_UNIBASE_OVRLP_GRAPH ) {
      cout << Tag() << "building unibase adjacency graph" << endl;
      digraph unibasegraph;
      BuildUnibaseAdjacencyGraph( unibases, unibasegraph, K );

      cout << Tag() << "writing unibase adjacency graph" << endl;
      BinaryWrite( run_dir + "/" + READS + "." + UNIBASES + "." + UNIBASE_GRAPH + ".k" + KS, unibasegraph );
    }
  }
  
  
  cout << Tag() << "Done with Unipather!" << endl;
  return 0;
}
