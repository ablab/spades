/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// LittleHelpsBig.  Use unipaths from little K1 to enlarge unipaths from big K_LG.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/GetNexts.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "paths/UnibaseUtils.h"


static inline 
String Tag(String S = "LHB") { return Date() + " (" + S + "): "; } 



// Build unipaths and find terminal unibases.
// Terminal unibases, or "leaves", are unibases that have only a single
// adjacency on the unipath adjacency graph.  Leaves that are anchored at
// their start are called "sinks", while leaves anchored at their end are
// called "sources".

void build_sinks_sources(const size_t K_LG,
			 const String & IN_HEAD_LG,
			 vec<longlong> * iub_sinks_lg_p,
			 vec<longlong> * iub_sources_lg_p)    
{ 
  cout << Tag() << "Loading unipaths" << endl;

  vecKmerPath unipaths(IN_HEAD_LG + ".unipaths.k" + ToString(K_LG));
  digraph Adj;
  BinaryRead(IN_HEAD_LG + ".unipath_adjgraph.k" + ToString(K_LG), Adj);
  HyperKmerPath hkp;
  BuildUnipathAdjacencyHyperKmerPath(K_LG, Adj, unipaths, hkp);
    
  // Find terminal unibases.
  cout << Tag() << "Finding sorces and sinks" << endl;
  for (int v = 0; v < hkp.N(); v++) {
    if (hkp.To(v).empty() && 
	hkp.From(v).size() == 1)
      iub_sources_lg_p->push_back(hkp.EdgeObjectIndexByIndexFrom(v, 0));

    if (hkp.From(v).empty() && 
	hkp.To(v).size() == 1)
      iub_sinks_lg_p->push_back(hkp.EdgeObjectIndexByIndexTo(v, 0));
  }
}



// Pre-processing step:
// For each sink, and each source, determine whether or not an alignment is
// possible using that unibase.  We do this by examining the unibases_sm.

void evaluate_sources_sinks(const BaseVecVec & unibases_sm,
                            const BaseVecVec & unibases_lg,
                            const size_t K_SM,
                            const String & OUT_HEAD,
                            const size_t NUM_THREADS,
                            const vec<longlong> & iub_sinks_lg,
                            const vec<longlong> & iub_sources_lg,
                            vec<bool> * good_sink_p,
                            vec<bool> * good_source_p)
{
  const size_t nub_sm = unibases_sm.size();
  const size_t n_sinks_sources = iub_sinks_lg.size();

  good_sink_p->clear();
  good_sink_p->resize(n_sinks_sources, true);
  good_source_p->clear();
  good_source_p->resize(n_sinks_sources, true);

  // Build paths, using both sets of unibases and K=K_SM.
  cout << Tag() << " Pathing" << endl;
  vecKmerPath paths_all_sm;
  {
    BaseVecVec unibases_all = unibases_sm;
    unibases_all.Append(unibases_lg);
    ReadsToPathsCoreY(unibases_all, K_SM, paths_all_sm, 
		      OUT_HEAD + ".LittleHelpsBig.unibases", NUM_THREADS);
  }


  // Form paths database for for unibases_sm.
  cout << Tag() << " Paths db" << endl;
  vec<tagged_rpint> pathsdb_sm;
  {
    vecKmerPath paths_sm(nub_sm);
    for (size_t iub = 0; iub < nub_sm; iub++)
      paths_sm[iub] = paths_all_sm[iub];
    CreateDatabase(paths_sm, pathsdb_sm);
  }

  // Write output files.
  // cout << Tag() << "Writing paths" << endl;
  // const String KS_SM = ToString(K_SM);
  // paths_all_sm.WriteAll(OUT_HEAD + ".all.paths.k" + KS_SM);
  // BinaryWrite3(OUT_HEAD + ".sm.pathsdb.k" + KS_SM, pathsdb_sm);
    
  cout << Tag() << " Good sinks" << endl;
  size_t n_good_sinks = n_sinks_sources;
  for (size_t iis = 0; iis < n_sinks_sources; iis++) {
    // ---- The last kmer in this unipath must appear EXACTLY once in the pathsdb_sm.
    const size_t iub_sink_all = nub_sm + iub_sinks_lg[iis];
    const kmer_id_t kmer = paths_all_sm[iub_sink_all].Stop();
    vec<longlong> rpints;
    Contains(pathsdb_sm, kmer, rpints);
    if (rpints.size() != 1) {
      (*good_sink_p)[iis] = false;
      n_good_sinks--;
    }
  }
  cout << Tag() << " n_good_sinks=   " << n_good_sinks << endl;
  

  cout << Tag() << " Good sources" << endl;
  size_t n_good_sources = n_sinks_sources;
  for (size_t iis = 0; iis < n_sinks_sources; iis++) {
    // ---- The first kmer in this unipath must appear EXACTLY once in the pathsdb_sm.
    const size_t iub_source_all = nub_sm + iub_sources_lg[iis];
    const kmer_id_t kmer = paths_all_sm[iub_source_all].Start();
    vec<longlong> rpints;
    Contains(pathsdb_sm, kmer, rpints);
    if (rpints.size() != 1) {
      (*good_source_p)[iis] = false;
      n_good_sources--;
    }
  }
    cout << Tag() << " n_good_sources= " << n_good_sources << endl;
  
}



// Create a BaseVecVec containing the last K_LG-1 bases of each of the unibases
// in sinks, and the first K_LG-1 bases of each of the unibases in sources.

void build_leaves(const size_t K_SM,
		  const size_t K_LG,
		  const BaseVecVec & unibases_lg,
		  const vec<longlong> & iub_sinks_lg,
		  const vec<longlong> & iub_sources_lg,
		  BaseVecVec * bvv_leaf_p)
{
  const size_t n_sinks   = iub_sinks_lg.size();
  const size_t n_sources = iub_sources_lg.size();
  const size_t n_total_leaves = n_sources + n_sinks;

  const size_t nb_leaf = K_LG - 1;

  bvv_leaf_p->reserve(n_total_leaves);

  for (size_t iis = 0; iis < n_sinks; iis++) {
    const BaseVec & bv_sink = unibases_lg[iub_sinks_lg[iis]];
    bvv_leaf_p->push_back(BaseVec(bv_sink, bv_sink.size() - nb_leaf, nb_leaf));
  }
  for (size_t iis = 0; iis < n_sources; iis++) {
    const BaseVec & bv_source = unibases_lg[iub_sources_lg[iis]];
    bvv_leaf_p->push_back(BaseVec(bv_source, 0, nb_leaf));
  }
}






// Now, gather a set of eligible pairs of sink and source unibases.
// A pair of sink and source unibases (call them S1 and S2) is considered
// "eligible" if the last K_SM-mer of S1 appears in S2 *and* the first K_SM-mer
// of S2 appears in S1.  If this criterion is not met, these unibases cannot
// possibly overlap by K_SM or more bases.

void build_eligibles(const vecKmerPath & paths_leaf,
		     const vec<tagged_rpint> & pathsdb_leaf,
                     const vec<bool> & good_sink,
                     const vec<bool> & good_source,
		     vec<pair<size_t, size_t> > * eligibles_p,
		     vec<unsigned> * counts_eligibles_p)
{
  const size_t n_sinks   = good_sink.size();
  const size_t n_sources = good_source.size();
  
  // Loop over all sink unipaths, and examine the last kmer in each.
  cout << Tag() << "Loop over sinks" << endl;
  for (size_t ii_sink = 0; ii_sink < n_sinks; dots_pct(ii_sink++, n_sinks)) {
    if (good_sink[ii_sink]) {

      // Find all other appearances of this kmer in pathsdb_leaf.
      vec<longlong> places;
      kmer_id_t kmer = paths_leaf[ii_sink].Stop();
      Contains(pathsdb_leaf, kmer, places);
      const size_t n_places = places.size();
      
      for (size_t ipl = 0; ipl < n_places; ipl++) {
        const longlong i_path = pathsdb_leaf[places[ipl]].PathId();
        if (i_path >= 0 &&                 // ignore RC alignments
            i_path >= longlong(n_sinks)) { // ignore alignments to other sink unipaths
          const size_t ii_source = i_path - n_sinks;
          if (good_source[ii_source])
            eligibles_p->push_back(make_pair(ii_sink, ii_source));
        }
      }
    }
  }

    
  // Loop over all source unipaths, and examine the first kmer in each.
  cout << Tag() << "Loop over sources" << endl;
  for (size_t ii_source = 0; ii_source < n_sources; dots_pct(ii_source++, n_sources)) {
    if (good_source[ii_source]) {
 
      // Find all other appearances of this kmer in pathsdb_leaf.
      vec<longlong> places;
      kmer_id_t kmer = paths_leaf[n_sinks + ii_source].Start();
      Contains(pathsdb_leaf, kmer, places);
      const size_t n_places = places.size();
      
      for (size_t ipl = 0; ipl < n_places; ipl++) {
        const longlong i_path = pathsdb_leaf[places[ipl]].PathId();
        if (i_path >= 0 &&                // ignore RC alignments
            i_path < longlong(n_sinks)) { // ignore alignments to other source unipaths
          const size_t ii_sink = i_path;
          if (good_sink[ii_sink]) 
            eligibles_p->push_back(make_pair(ii_sink, ii_source));
        }
      }
    }
  }
  
  // Every eligible sink/source pair should now appear TWICE in eligibles.
  // Sort the list of eligibles, and keep track of how many times each s/s
  // pair appeared - later, we'll drop all s/s pairs that appeared only once.
  cout << Tag() << "Unique sort" << endl;
  UniqueSortAndCount(*eligibles_p, *counts_eligibles_p);

  cout << Tag() << "n_eligibles = " << eligibles_p->size() << endl;
}  








void add_bridges(const size_t K_LG,
                 const size_t K_SM,
                 const BaseVecVec & unibases_lg,
                 const vec<longlong> & iub_sinks_lg,
                 const vec<longlong> & iub_sources_lg,
                 const vec<bool> & good_sink,
                 const vec<bool> & good_source,
                 const vec< pair<size_t, size_t> > & eligibles,
                 const vec<unsigned> & counts_eligibles,
                 BaseVecVec * bvv_new_p)
{  
  const size_t n_eligibles = eligibles.size();
  size_t n_bridges = 0;
      
  map<size_t, size_t> overlaps;
  for (size_t iel = 0; iel < n_eligibles; dots_pct(iel++, n_eligibles)) {
    const size_t ii_sink   = eligibles[iel].first;
    const size_t ii_source = eligibles[iel].second;
        
    if (good_sink  [ii_sink]   &&    // unnecessary
        good_source[ii_source] &&    // unnecessary
        counts_eligibles[iel] != 1) {
          
      const BaseVec & ubv_sink   = unibases_lg[iub_sinks_lg[ii_sink]];
      const BaseVec & ubv_source = unibases_lg[iub_sources_lg[ii_source]];
          

      if (true) {

        // ---- NEW VERSION: only need last K_LG bases to check for overlaps

        const BaseVec bv_sink(ubv_sink, ubv_sink.size() - K_LG, K_LG);
        const BaseVec bv_source(ubv_source, 0, K_LG);
        
        // Check for direct overlap between the paths. 
        
        const size_t overlap = LargestOverlap(bv_sink, bv_source, K_LG-1, K_SM);
        if (overlap > 0) {
          overlaps[overlap]++;
          
          // We've found an overlap! Add bridge between unibases to bvv_new.
          bvv_new_p->push_back(Cat(BaseVec(bv_sink, 0, bv_sink.size() - overlap), 
                                   bv_source));
          n_bridges++;
        }
        
      }
      else {

        // ---- OLD VERSION: join two full unibases as a bridge and add it
        // Check for direct overlap between the paths. 
        
        const size_t overlap = LargestOverlap(ubv_sink, ubv_source, K_LG-1, K_SM);
        if (overlap > 0) {
          overlaps[overlap]++;
	  
          // We've found an overlap! Add bridge between unibases to bvv_new.
          bvv_new_p->push_back(Cat(BaseVec(ubv_sink, 0, ubv_sink.size() - overlap), 
                                   ubv_source));
          n_bridges++;
        }
      }
    }
  }
  if (false)
    for (map<size_t, size_t>::const_iterator it = overlaps.begin();
	 it != overlaps.end(); it++)
      cout << "overlap: " << it->first << " " << it->second << endl;


  cout << Tag() << "n_bridges= " << n_bridges << endl;
}  








void add_adjacencies(const BaseVecVec & unibases_lg,
                     const size_t K_LG,
                     BaseVecVec * bvv_new_p)
{
  const size_t nub = unibases_lg.size();
  vec< vec<int> > iub_nexts(nub);
  GetNexts(K_LG, unibases_lg, iub_nexts);
  for (size_t iub1 = 0; iub1 < nub; iub1++) {
    const BaseVec & ubv1 = unibases_lg[iub1];

    if (true) {
      // ---- NEW VERSION: only K_LG + 1 common bases

      const BaseVec bv1(ubv1, ubv1.size() - K_LG, K_LG);
      for (size_t ii = 0; ii < iub_nexts[iub1].size(); ii++) {
        const size_t iub2 = iub_nexts[iub1][ii];
        const BaseVec & ubv2 = unibases_lg[iub2];
        bvv_new_p->push_back(Cat(bv1, BaseVec(ubv2, K_LG - 1, 1)));
      }
    }
    else {
      // ---- OLD VERSION: full unibase + 1 base

      for (size_t ii = 0; ii < iub_nexts[iub1].size(); ii++) {
        const size_t iub2 = iub_nexts[iub1][ii];
        const BaseVec & ubv2 = unibases_lg[iub2];
        bvv_new_p->push_back(Cat(ubv1, BaseVec(ubv2, K_LG - 1, 1)));
      }
    }
  }
}










int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(IN_HEAD_SM);
  CommandArgument_String(IN_HEAD_LG);
  CommandArgument_Int(K_SM);
  CommandArgument_Int(K_LG);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_Bool_OrDefault(REPATH, True);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;
  
  RunTime();  

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Load unibases.

  cout << Tag() << "Loading unibases" << endl;
  
  BaseVecVec unibases_lg(IN_HEAD_LG + ".unibases.k" + ToString(K_LG));
  const size_t nub_lg = unibases_lg.size();
  const size_t nkmer_lg = unibases_lg.SizeSum() - unibases_lg.size() * (K_LG - 1);
  cout << Tag() << "n_unibases(K=" << K_LG << ")= " << nub_lg << endl;
  cout << Tag() << "n_kmers(K=" << K_LG << ")= " << nkmer_lg << endl;

  // This will store all new base vectors.
  BaseVecVec bvv_new = unibases_lg;
  
  if (true) {  // Scoping out intermediate data structures


    BaseVecVec unibases_sm(IN_HEAD_SM + ".unibases.k" + ToString(K_SM));
    const size_t nub_sm = unibases_sm.size();
    cout << Tag() << "n_unibases(K=" << K_SM << ")= " << nub_sm << endl;

    // Make output directory, if necessary.
    if (OUT_HEAD.Contains("/")) Mkpath(OUT_HEAD.RevBefore("/"));

  
    // Build unipaths and find terminal unibases ('leaves') which are
    // unibases that have only a single adjacency on the unipath adjacency graph.  
    // Leaves at the unibase start are 'sources' and at the end are 'sinks'.

    vec<longlong> iub_sinks_lg;
    vec<longlong> iub_sources_lg;
    build_sinks_sources(K_LG, IN_HEAD_LG, & iub_sinks_lg, & iub_sources_lg);
    const size_t n_sinks   = iub_sinks_lg.size();
    const size_t n_sources = iub_sources_lg.size();
    ForceAssertEq(n_sinks, n_sources);
    cout << Tag() << "n_sinks=   " << n_sinks << endl;
    cout << Tag() << "n_sources= " << n_sources << endl;
    
    

    // Determine which unibases_sm can follow which.  This is not quite right
    // because it assumes that K_SM-1 overlap is enough for adjacency.
  
    // Pre-processing step:
    // For each sink, and each source, determine whether or not an alignment is
    // possible using that unibase.  We do this by examining the unibases_sm.

    cout << Tag() << "Evaluating sinks and sources" << endl;

    vec<bool> good_sink;
    vec<bool> good_source;
    evaluate_sources_sinks(unibases_sm, unibases_lg, 
                           K_SM, OUT_HEAD, NUM_THREADS, 
                           iub_sinks_lg, iub_sources_lg,
                           & good_sink, & good_source);
  

    vec< pair<size_t, size_t> > eligibles;
    vec<unsigned> counts_eligibles;
    if (true) {

      // Create a BaseVecVec containing the last K_LG bases of each of the unibases
      // in sinks, and the first K_LG bases of each of the unibases in sources.
      // Then, path these bases using K_SM, and create a pathsdb.  Now we can look up
      // a K_SM-mer in the pathsdb to quickly determine whether or not an overlap of
      // size > K_SM is possible.
      vecKmerPath paths_leaf;
      vec<tagged_rpint> pathsdb_leaf;
      {
        cout << Tag() << "Building leaves" << endl;
        BaseVecVec bvv_leaf;
        build_leaves(K_SM, K_LG, unibases_lg, iub_sinks_lg, iub_sources_lg, 
                     & bvv_leaf);
        cout << Tag() << "Pathing leaves" << endl;
        vecKmerPath paths_leaf_rc;
        ReadsToPathsCoreY(bvv_leaf, K_SM, paths_leaf, paths_leaf_rc, pathsdb_leaf);
      }
      
      
      // Now, gather a set of eligible pairs of sink and source unibases.
      // A pair of sink and source unibases (call them S1 and S2) is considered
      // "eligible" iff the last K_SM-mer of S1 appears in S2 *and* the first K_SM-mer
      // of S2 appears in S1.  If this criterion is not met, these unibases cannot
      // possibly overlap by K_SM or more bases.
      cout << Tag() << "Building list of eligible pairs" << endl;
      
      build_eligibles(paths_leaf, pathsdb_leaf, good_sink, good_source,
                      & eligibles, & counts_eligibles);
    }
  
    // MAIN ALGORITHM
    // For each eligible sink/source pair, try to bridge.

    cout << Tag() << "Looking for bridges" << endl;
    add_bridges(K_LG, K_SM, unibases_lg, iub_sinks_lg, iub_sources_lg, 
                good_sink, good_source, eligibles, counts_eligibles, 
                & bvv_new);

    
    // Make sure we keep adjacency information

    cout << Tag() << "Adding adjacencies" << endl;
    add_adjacencies(unibases_lg, K_LG, & bvv_new);

  } // end scoping of intermediate data structures

  if (REPATH) {

    // Now build the new unipaths.
    
    cout << Tag() << "Pathing everything" << endl;
    vecKmerPath       paths_new;
    vecKmerPath       paths_new_rc;
    vec<tagged_rpint> pathsdb_new;
    ReadsToPathsCoreY(bvv_new, K_LG, paths_new, paths_new_rc, pathsdb_new,
		      OUT_HEAD + ".LittleHelpsBig.pseudo_reads", NUM_THREADS);
    
    cout << Tag() << "Unipathing everything" << endl;
    vecKmerPath       unipaths_new;
    vec<tagged_rpint> unipathsdb_new;
    Unipath(paths_new, paths_new_rc, pathsdb_new, unipaths_new, unipathsdb_new);
    
    cout << Tag() << "Brokering kmers and bases" << endl;
    KmerBaseBroker    kbb_new(K_LG, paths_new, paths_new_rc, pathsdb_new, bvv_new);
    
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
    
  } else {

    // Write out new basevectors
    bvv_new.WriteAll(OUT_HEAD + ".part1.fastb");
    
  } 
  
  cout << Tag() << "Done!" << endl;
  return 0;
}
