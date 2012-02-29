/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2011) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// LittleHelpsBigPart2.  Use unipaths from little K1 to enlarge unipaths from big K_LG.

#include "MainTools.h"
#include "Basevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "paths/UnibaseUtils.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"


static inline 
String Tag(String S = "LHB2") { return Date() + " (" + S + "): "; } 


int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(IN_HEAD_LG);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_UnsignedInt(K_LG);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
  EndCommandArguments;
  
  RunTime();  

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Memory control 

  SetMaxMemory(MAX_MEMORY_GB<<30);
 
  // Load unibases.

  cout << Tag() << "Loading unibases to compute size" << endl;

  size_t nub_lg,  nkmer_lg;
  {
    BaseVecVec unibases_lg(IN_HEAD_LG + ".unibases.k" + ToString(K_LG));
    nub_lg = unibases_lg.size();
    nkmer_lg = unibases_lg.SizeSum() - unibases_lg.size() * (K_LG - 1);
  }

  // ****  Experimental pathing code with hard coded K ****

  vecbasevector K_1_unibases;

  {
    unsigned const K = 96;       // Output unibase K
    unsigned const K_1 = K + 1;  // Internal K for pathing

    ForceAssert(K == K_LG);

    const String fastb_in = OUT_HEAD + ".part1.fastb";

    // Build internal unipath graph using ReadPather at K+1

    cout << Tag() << " Building K=" << K_1 << " unipaths " << endl;

    String graphInfoFilename = UnipathGraph<K_1>::getInfoFilename(fastb_in);
    String countsFile = UnipathGraph<K_1>::getCountsFilename(graphInfoFilename);
    size_t nKmers = nkmer_lg ; // original kmer count at K_LG
    KmerDict<K_1>* pDict;
    pDict = new KmerDict<K_1>(nKmers);
    pDict->process(fastb_in, NUM_THREADS);
    UnipathGraph<K_1>::create(graphInfoFilename,*pDict);
  
    // Generate K+1 unibases
    
    cout << Tag() << " Building K=" << K_1 << " unibases " << endl;
    
    UnipathGraph<K_1> graph(graphInfoFilename);
    UnipathEdgeVec const& edges = graph.getAllEdges();
    
    K_1_unibases.reserve(edges.size());
   
    typedef HugeBVec::const_iterator Itr;
    for ( size_t idx = 0; idx < edges.size(); ++idx ) {
      basevector bv;
      UnipathEdge const& edge = edges[idx];
      Itr bItr = graph.getBases(edge.getInitialKmerID());
      Itr bEnd = graph.getBases(edge.getFinalKmerID()) + K_1;
      bv.append(bItr,bEnd);
      K_1_unibases.push_back(bv);
    }
  }

  // Now build the new unipaths.

  cout << Tag() << "Pathing everything" << endl;
  vecKmerPath       paths_new;
  vecKmerPath       paths_new_rc;
  vec<tagged_rpint> pathsdb_new;
  ReadsToPathsCoreY(K_1_unibases, K_LG, paths_new, paths_new_rc, pathsdb_new,
		    OUT_HEAD + ".LittleHelpsBig.pseudo_reads", NUM_THREADS);

  cout << Tag() << "Unipathing everything" << endl;
  vecKmerPath       unipaths_new;
  vec<tagged_rpint> unipathsdb_new;
  Unipath(paths_new, paths_new_rc, pathsdb_new, unipaths_new, unipathsdb_new);

  cout << Tag() << "Brokering kmers and bases" << endl;
  KmerBaseBroker    kbb_new(K_LG, paths_new, paths_new_rc, pathsdb_new, K_1_unibases);

  cout << Tag() << "Building adjacencies" << endl;
  digraph           adj_graph_new;
  BuildUnipathAdjacencyGraph(paths_new, paths_new_rc, pathsdb_new, 
                             unipaths_new, unipathsdb_new, adj_graph_new);

  // Write output files.
  cout << Tag() << "Writing output files" << endl;
  const String KS_LG = ToString(K_LG);
  unipaths_new.WriteAll(OUT_HEAD + ".unipaths.k" + KS_LG);
  BinaryWrite3(OUT_HEAD + ".unipathsdb.k" + KS_LG, unipathsdb_new);
  BinaryWrite(OUT_HEAD + ".unipath_adjgraph.k" + KS_LG, adj_graph_new);
  IncrementalWriter<basevector> writer(OUT_HEAD + ".unibases.k" + KS_LG, unipaths_new.size());
  for ( size_t i = 0; i < unipaths_new.size( ); i++ )
    writer.add(kbb_new.Seq( unipaths_new[i] ));
  writer.close();
  
  cout << Tag() << "Done!" << endl;
  return 0;
}
