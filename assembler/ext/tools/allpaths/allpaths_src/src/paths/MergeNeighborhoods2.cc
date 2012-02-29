///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MergeNeighborhoods2.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/InternalMerge.h"
#include "paths/InternalMergeImpl.h" 
#include "paths/KmerBaseBroker.h" 
#include "paths/PdfEntry.h" 
#include "paths/Unipath.h"
#include "feudal/BinaryStream.h"

#include "system/WorklistN.h"  // worklistN


static inline 
String Tag(String S = "MN2") { return Date() + " (" + S + "): "; } 

class overlapit 
{
public:

  int c1;
  int c2;
  int overlap;
  int j1;
  int j2;

  overlapit() { }

  overlapit(const int c1, 
            const int c2, 
            const int overlap, 
            const int j1,
            const int j2) 
    : c1(c1), 
      c2(c2), 
      overlap(overlap), 
      j1(j1), 
      j2(j2) 
  {}
          
  friend bool operator<(const overlapit& x, const overlapit& y)
  {
    if (x.c1 < y.c1) return true;
    if (x.c1 > y.c1) return false;
    if (x.c2 < y.c2) return true;
    if (x.c2 > y.c2) return false;
    if (x.overlap < y.overlap) return true;
    if (x.overlap > y.overlap) return false;
    if (x.j1 < y.j1) return true;
    if (x.j1 > y.j1) return false;
    if (x.j2 < y.j2) return true;
    return false;    
  }

  friend bool operator==(const overlapit& x, const overlapit& y)
  {
    return x.c1 == y.c1 && x.c2 == y.c2 && x.overlap == y.overlap
      && x.j1 == y.j1 && x.j2 == y.j2;     
  }

};

TRIVIALLY_SERIALIZABLE(overlapit);







bool Nice(const HyperKmerPath& h)
{
  vec<int> sources, sinks;
  h.Sources(sources), h.Sinks(sinks);
  return sources.solo() && sinks.solo() /* && h.Acyclic() */;    
}











// A LabeledHyperKmerPath is a HyperKmerPath together with an integer label
// assigned to each of its edges.

class LabeledHyperKmerPath : public HyperKmerPath 
{
private:
  vec<int> L_;

public:

  LabeledHyperKmerPath() { }

  LabeledHyperKmerPath(const HyperKmerPath& h, const vec<int>& L)
    : HyperKmerPath(h), L_(L) { }

  int EdgeId(const int i) const { return L_[i]; }

  vec<int> EdgeIds() const { return L_; }

  // Components: takes as input a possibly disconnected object and returns
  // its components.  Returns true if output != input.

  bool Components(vec<LabeledHyperKmerPath>& H)
  {
    H.clear();
    vec< vec<int> > comp;
    HyperKmerPath(*this).Components(comp);
    ForceAssert(comp.nonempty());
    if (comp.solo()) {
      H.resize(1);
      H[0] = *this;
      return false;
    }
    for (size_t i = 0; i < comp.size(); i++) {
      vec<int> L;
      HyperKmerPath h(*this, comp[i]);
      for (size_t v = 0; v < comp[i].size(); v++) {
        for (size_t j = 0; j < From(comp[i][v]).size(); j++) {
          int e = EdgeObjectIndexByIndexFrom(comp[i][v], j);
          L.push_back(EdgeId(e));
        }
      }
      H.push(h, L); 
    }
    return true; 
  }

  // Disintegrate: takes as input a connected object and splits it into
  // individual edges.  Returns true if output != input.

  bool Disintegrate(vec<LabeledHyperKmerPath>& H)
  {
    H.clear();
    if (EdgeObjectCount() == 1) {
      H.resize(1);
      H[0] = *this;
      return false;
    }
    for (int i = 0; i < EdgeObjectCount(); i++) {
      vec<KmerPath> p(1);
      p[0] = EdgeObject(i);
      HyperKmerPath h(K(), p);
      vec<int> L(1);
      L[0] = EdgeId(i);
      H.push(h, L);
    }
    return true;
  }

  // DivideMultipleEnds: takes as input a connected object and returns
  // a collection of connected objects with multiple ends (sources or sinks)
  // divided up if possible.  Returns true if output != input.

  bool DivideMultipleEnds(vec<LabeledHyperKmerPath>& H)
  {
    H.clear();
    vec<int> sources, sinks;
    Sources(sources), Sinks(sinks);
    if (sources.size() <= 1 && sinks.size() <= 1) {
      H.resize(1);
      H[0] = *this;
      return false;
    }
    vec<int> verts_to_splay;
    if (sinks.size() > 1) {
      
      // Find the 'maximal predecessors' of all sinks.

      vec< vec<int> > pred(sinks.size());
      for (size_t l = 0; l < sinks.size(); l++)
        GetPredecessors1(sinks[l], pred[l]);
      vec<int> p;
      Intersection(pred, p);
      for (size_t l = 0; l < p.size(); l++) {
        bool maximal = true;
        int v = p[l];
        for (size_t m = 0; m < From(v).size(); m++) {
          int w = From(v)[m];
          if (w != v && BinMember(p, w)) maximal = false;    
        }
        if (maximal) verts_to_splay.push_back(v);
      }
    }
    if (sources.size() > 1) {
        // Find the 'minimal successors' of all sources.
      
      vec< vec<int> > succ(sources.size());
      for (size_t l = 0; l < sources.size(); l++)
        GetSuccessors1(sources[l], succ[l]);
      vec<int> s;
      Intersection(succ, s);
      for (size_t l = 0; l < s.size(); l++) {
        bool minimal = true;
        int v = s[l];
        for (size_t m = 0; m < To(v).size(); m++) {
          int w = To(v)[m];
          if (w != v && BinMember(s, w)) minimal = false;
        }
        if (minimal) verts_to_splay.push_back(v);
      }
    }
    if (verts_to_splay.nonempty()) {
      UniqueSort(verts_to_splay);
      for (size_t j = 0; j < verts_to_splay.size(); j++)
        SplayVertex(verts_to_splay[j]);
      Components(H);
      return true;
    }
    else {
      H.resize(1);
      H[0] = *this;
      return false;
    }
  }

};









void RunScaffoldAccuracy(const HyperKmerPath &hyper,
                         const KmerBaseBroker &kbb,
                         const String &ass_base,
                         const String &ref_base,
                         String &info)
{
  String ass_fasta_file = ass_base + ".fasta";
  hyper.DumpFasta(ass_fasta_file, kbb);

  String n_comps = ToString(hyper.ConnectedComponents());
  info += n_comps + " comps, ";

  // Run ScaffoldAccuracy with different D (distances).
  for (int did=0; did<3; did++) {
    int dist = pow(10, did) * 1000;
    String str_dist = ToString(dist);
    String log_file = ass_base + "_d" + str_dist + ".log";
    String loc_info = "d_" + str_dist + "_";
    
// MakeDepend: dependency ScaffoldAccuracy
    String theCommand
      = "ScaffoldAccuracy ASSEMBLY=" + ToString(ass_fasta_file)
      + " D=" + str_dist
      + " REFHEAD=" + ref_base
      + " >& " + log_file;
    System(theCommand);
    
    String aline;
    ifstream in(log_file.c_str());
    while (in) {
      getline(in, aline);
      if (! in) break;
      if (aline.Contains("no contig", 0)) {
	loc_info += "[na] ";
	break;
      }
      if (aline.Contains("invalids", 0)) {
	String strInv = aline.Before(",").After("= ");
	String strTot = aline.After("total = ");
	if (strTot == "0") {
	  loc_info += "[na (0/0)] ";
	  break;
	}
	if (strInv == "0") loc_info += "[ok (";
	else loc_info += "[BAD (";
	loc_info += strInv + "/" + strTot + ")] ";
	break;
      }
    }
    in.close();

    info += loc_info;
  }
  
}


















size_t FindComponents(size_t const K,
                      String const & sub_dir,
                      vec< vec<int> > * comp, 
                      vec<HyperKmerPath> * nhcomp)
{    
  // Load neighborhoods.
  cout << Tag() << "finding components in all nhoods" << endl;
  
  cout << Tag() << "loading local HyperKmerPaths" << endl;
  vec<HyperKmerPath> nhoods;
  BinaryRead(sub_dir + "/nhood.hypers", nhoods);
  
  // Find components.
  
  size_t n_edges = 0;

  vec<LabeledHyperKmerPath> ALL;
  for (size_t i = 0; i < nhoods.size(); i++) {
    if (i % 100000 == 0) DPRINT2(i, nhoods.size());

    const HyperKmerPath& nh = nhoods[i];

    if (nh.EdgeObjectCount() == 0) continue;

    vec<int> ids(nh.EdgeObjectCount());
    for (size_t j = 0; j < ids.size(); j++)
      ids[j] = j + n_edges;
    LabeledHyperKmerPath X0(nh, ids);
    vec<LabeledHyperKmerPath> X;
    X0.Components(X);
    for (size_t j = 0; j < X.size(); j++) {

      if (Nice(X[j])) continue;

      vec<LabeledHyperKmerPath> H;
      bool changed = X[j].DivideMultipleEnds(H);
      if (changed) {
        X[j] = H[0];
        for (size_t l = 1; l < H.size(); l++)
          X.push_back(H[l]);
        j--;
      }
      else {
        changed = X[j].Disintegrate(H);
        if (changed) {
          X[j] = H[0];
          for (size_t l = 1; l < H.size(); l++)
            X.push_back(H[l]);
        }
      }
    }

    bool verbose = false;
    if (verbose) {
      cout << "\nNHOOD " << i << " INPUT:\n";
      nh.PrintSummaryDOT0w(cout);
      cout << "\nNHOOD " << i << " OUTPUT:\n";
      vec<HyperKmerPath> V;
      for (size_t j = 0; j < X.size(); j++)
        V.push_back(X[j]);
      HyperKmerPath H(K, V);
      H.PrintSummaryDOT0w(cout);
      cout << "\n";
    }


    ALL.append(X);
    n_edges += nh.EdgeObjectCount();
  }

  for (size_t i = 0; i < ALL.size(); i++) {
    nhcomp->push_back(ALL[i]);
    comp->push_back(ALL[i].EdgeIds());
  }


  int disint = 0;

/*
  vec< vec<int> > one_comp;
  nh.ComponentsE(one_comp);
  vec<int> to_left, to_right;
  nh.ToLeft(to_left), nh.ToRight(to_right);
  for (size_t j = 0; j < one_comp.size(); j++) {
    vec<int> compv;
    for (size_t l = 0; l < one_comp[j].size(); l++) {
      compv.push_back(to_left[ one_comp[j][l] ]);
      compv.push_back(to_right[ one_comp[j][l] ]);
    }
    UniqueSort(compv);
    HyperKmerPath hlet(nh, compv);
       
    // Decide if the component topology is nice.
    
    if (!hlet.Nice()) {    
      // The component topology is not nice.  For now, we just 
      // disintegrate it into edges.
      
      if (disint == 0) PRINT(i);
      ++disint;
      for (size_t l = 0; l < one_comp[j].size(); l++) {
        int c = one_comp[j][l];
        vec<KmerPath> p;
        p.push_back(nh.EdgeObject(c));
        nhcomp.push(K, p);
        vec<int> C;
        C.push_back(c + n_edges);
        comp.push_back(C);
      }
    }
    else {    
      // The component topology is nice.
      
      for (size_t l = 0; l < one_comp[j].size(); l++)
        one_comp[j][l] += n_edges;
      vec<int> C = one_comp[j];
      
      // Push back the component.
      
      comp.push_back(C);
      nhcomp.push_back(hlet);
    }
  }
  
  n_edges += nh.EdgeObjectCount();
  
*/

  cout << Tag() << "n_comp = " << comp->size() << endl;
  cout << Tag() << "disint = " << disint << endl;

  return n_edges;
}








void ConstructInvolution(const size_t n_uni, 
                         const size_t n_edges, 
                         const vecKmerPath & spaths,
                         const vec<tagged_rpint> & spathsdb,
                         vec<int> * to_rc)
{
  ForceAssertEq(to_rc->size(), n_uni);
  cout << Tag() << "constructing involution" << endl;
  cout << Tag() << "n_edges = " << n_edges << endl;

  for (size_t iu = 0; iu < n_uni; iu++) {
    KmerPath U = spaths[n_edges + iu];
    U.Reverse();
    vec<int64_t> locs;
    Contains(spathsdb, U.Segment(0), locs);
    size_t n_locs = locs.size();
    for (size_t il = 0; il < n_locs; il++) {

      const tagged_rpint& t = spathsdb[ locs[il] ];

      if (t.Fw() && (size_t)t.ReadId() >= n_edges) {
        (*to_rc)[iu] = (size_t)t.ReadId() - n_edges;
        break;
      }
    }
  }
  for (size_t iu = 0; iu < n_uni; iu++)
    ForceAssertEq(iu, (size_t)(*to_rc)[ (*to_rc)[iu] ]);
}







// Find all cases where two components overlap.  Scan the paths database.
void FindAllComponentOverlaps(const size_t ploidy,
                              const vec<tagged_rpint> & spathsdb,
                              const vec<int> & to_rc,
                              const vec<int> & to_comp,
                              const vec<int> & predicted_CNs,
			      const String & working_dir,
                              vec<overlapit> * overlaps,
			      const vec<bool>& edge_repeat)
{
  // write overlaps to disk for memory efficiency on large genomes
  String tmp_file = working_dir + "/" + "MN2_overlaps.tmp";
  BinaryIterativeWriter<vec<overlapit> > overlaps_out(tmp_file.c_str());
  
  const size_t n_edges = to_comp.size();
    
    cout << Tag() << "looking for overlaps" << endl;
  
  int repeat_counter = 0;
  for (size_t i = 0; i < spathsdb.size(); i++) {
    const tagged_rpint& t = spathsdb[i];
    size_t id = t.ReadId();
    
    // We're looking first for a unipath edge of copy number "one".
    int u = id - n_edges;
    if (t.Rc() || 
        u < 0 || 
        predicted_CNs[u] > (int)ploidy || 
        to_rc[u] < u)
      continue;
        
    // Now find all the neighborhood edges that it overlaps.
        
    size_t ndb = spathsdb.size();
    ForceAssertGe((int64_t)i, (int64_t)t.Lookback());
    size_t jlo = i - t.Lookback();
    size_t jhi = i + 1;
    while (jhi < ndb && spathsdb[jhi].Start() <= t.Stop()) 
      jhi++;
        
    // Find all overlaps between different components.
       
    for (size_t j1 = jlo; j1 < jhi; j1++) {
      const tagged_rpint& t1 = spathsdb[j1];
      const size_t id1 = t1.ReadId();

      if (id1 < n_edges) {
        int c1 = to_comp[id1];

        for (size_t j2 = j1 + 1; j2 < jhi; j2++) {
          const tagged_rpint& t2 = spathsdb[j2];
          const size_t id2 = t2.ReadId();
          if (id2 < n_edges) {
            int c2 = to_comp[id2];
            if (c1 != c2) {
              int overlap = IntervalOverlap(t1.Start(), t1.Stop(),
                                            t2.Start(), t2.Stop());
              if (overlap > 0) {
		// skip overlaps of repetative region
		if ( edge_repeat[id1] && edge_repeat[id2] ) { repeat_counter++;  continue; }
                if (c1 <= c2)
		  overlaps_out.write(overlapit(c1, c2, overlap, j1, j2));
                else 
		  overlaps_out.write(overlapit(c2, c1, overlap, j2, j1));
              }
            }
          }
        }

      }

    }


  }
    cout << Tag() << "skipped overlaps " << repeat_counter<< endl;

  overlaps_out.close();

  // load entire set of overlaps and erase tmp file
  BinaryReader::readFile<vec<overlapit> >(tmp_file.c_str(), overlaps);
  Remove(tmp_file);
}





// Mine the overlaps.  Note that the way this works now, if sequence x
// appears once in c1 and twice in c2, it counts double (etc.), which
// doesn't really make sense.



void MineComponentOverlaps(const size_t K, 
                           const vecKmerPath & spaths,
                           const vec< vec<int> > & comp,
                           const vec<HyperKmerPath> & nhcomp,
                           const vec<overlapit> & overlaps,
                           const size_t min_component, 
                           const size_t MIN_OVERLAP_TO_CLUSTER, 
                           vec<HyperKmerPath> * hcomps,
                           size_t * n_kmers_max)
{
  cout << Tag() << "mining overlaps" << endl;
  size_t n_comp = comp.size();

  vec< vec<int> > from(n_comp);
  vec< vec<int> > to(n_comp);
  
  size_t n_overlaps = overlaps.size();

  for (size_t i1 = 0; i1 < n_overlaps; i1++) {
    int c1 = overlaps[i1].c1;
    int c2 = overlaps[i1].c2;
    size_t over = 0;

    size_t i2 = i1 + 1;
    while (i2 != n_overlaps && 
           overlaps[i2].c1 == c1 &&
           overlaps[i2].c2 == c2) 
      i2++;

    /*
    for (i2 = i1 + 1; i2 < overlaps.size(); i2++) {
      if (overlaps[i2].c1 != c1) break;
      if (overlaps[i2].c2 != c2) break;
    }
    */

    for (size_t j = i1; j != i2; j++)
      over += overlaps[j].overlap;
    
    if (over >= MIN_OVERLAP_TO_CLUSTER) {
      from[c1].push_back(c2);
      from[c2].push_back(c1);
      to[c1].push_back(c2);
      to[c2].push_back(c1);
      // PRINT3(c1, c2, over);
    }

    i1 = i2 - 1;
  }


  cout << Tag() << "sorting" << endl;
  for (size_t i = 0; i < n_comp; i++) {
    Sort(from[i]);
    Sort(to[i]);
  }

  
  digraph G(from, to);
  cout << Tag() << "finding components" << endl;
  vec< vec<int> > Gcomp;
  G.Components(Gcomp);
  
  cout << Tag() << "finding component sizes" << endl;
  for (size_t i = 0; i < Gcomp.size(); i++) {

    const vec<int>& o = Gcomp[i];
    size_t n_kmers = 0;
    vec<HyperKmerPath> hs;

    for (size_t j = 0; j < o.size(); j++) {
      for (size_t l = 0; l < comp[ o[j] ].size(); l++)
        n_kmers += spaths[ comp[ o[j] ][ l ] ].KmerCount();
      hs.push_back(nhcomp[ o[j] ]);
    }

    if (n_kmers >= min_component) {
      hcomps->push(K, hs);
      *n_kmers_max = Max(*n_kmers_max, n_kmers);
      // PRINT2(i, n_kmers);   
    } 
  }

}






class MergeProcessor
{
private:
  const size_t _n_threads;
  const size_t _min_overlap;
  const size_t _min_proper_overlap;
  vec<HyperKmerPath> * _hcomps_p;
  const bool _accuracy;
  const String & _base_file;
  const String & _base_ref;
  const KmerBaseBroker & _kbb;
  const size_t _dotter;

public:
  MergeProcessor(const size_t n_threads,
                 const size_t min_overlap,
                 const size_t min_proper_overlap,
                 vec<HyperKmerPath> * hcomps_p,
                 const bool accuracy,
                 const String & base_file,
                 const String & base_ref,
                 const KmerBaseBroker & kbb,
                 const size_t dotter)
    : 
    _n_threads(n_threads),
    _min_overlap(min_overlap),
    _min_proper_overlap(min_proper_overlap),
    _hcomps_p(hcomps_p),
    _accuracy(accuracy),
    _base_file(base_file),
    _base_ref(base_ref),
    _kbb(kbb),
    _dotter(dotter)
  {}

  // Parallelize over n_threads instead of over the number of objects.
  // When the number of objects is too large and processing each object is 
  // very fast, you're better off partitioning the data yourself.
  void operator() (const size_t thread_ID)
  {
    vec<tagged_rpint> uniqdb_null;
    size_t n_hcomps = _hcomps_p->size();
    
    size_t n_dotter = _dotter / _n_threads;
    
    size_t IDlo = (thread_ID * n_hcomps) / _n_threads;
    size_t IDhi = ((thread_ID + 1) * n_hcomps) / _n_threads;
    
    for (size_t ID = IDlo; ID != IDhi; ID++) {
      
      if ((!_accuracy) && 
          thread_ID == 0 && 
          ID % n_dotter == 0) 
        Dot(cout, ID / n_dotter);
      
      
      HyperKmerPath * hkpp = &((*_hcomps_p)[ID]);
      
      String base_ass = _base_file + ToString(ID);
      String info;
      if (_accuracy) {
        info = ToString(ID) + "." + ToString(n_hcomps - 1) + "   ";
        RunScaffoldAccuracy(*hkpp, _kbb, base_ass, _base_ref, info);
      }
    
      if (hkpp->ConnectedComponents() > 1) {
        
        if (_accuracy) info += "- merging - ";
      
        uniqdb_null.clear();
        InternalMergeImpl(*hkpp, NULL, 
                          _min_overlap, _min_proper_overlap, 
                          20, false, false, uniqdb_null);
      }
      else {
        if (_accuracy) info += "- zipping only - ";
      }
      
      hkpp->Zipper();
      RemoveHangingEnds(*hkpp, &KmerPath::KmerCount, 250, 5.0);
      hkpp->RemoveSmallComponents(1000);
      hkpp->RemoveUnneededVertices();
      hkpp->RemoveDeadEdgeObjects();
      
      
      if (_accuracy) {
        RunScaffoldAccuracy(*hkpp, _kbb, base_ass, _base_ref, info);
        info += "\n";
        cout << info;
      }
    }
  
  }

};















int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_String_OrDefault(READS, "all_reads");
  CommandArgument_Int(K);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Int_OrDefault(MIN_OVERLAP, 10000);
  CommandArgument_Int_OrDefault(MIN_PROPER_OVERLAP, 2000);
  CommandArgument_Int_OrDefault(MIN_OVERLAP_TO_CLUSTER, 2000);
  CommandArgument_Bool_OrDefault(SKIP_REPEAT, True);
  CommandArgument_Int_OrDefault(MAX_REPEAT, 10000);
  CommandArgument_Bool_OrDefault(WRITE, true);
  CommandArgument_String_OrDefault(OUT_SUFFIX, "");

  // If ACCURACY is true, run ScaffoldAccuracy on each component,
  // before and after merge (heavy impact on performance!) Notice
  // also that some intermediate files will not be removed.
  CommandArgument_Bool_OrDefault(ACCURACY, false);

  EndCommandArguments;

  // Thread control (OMP used by ParallelSort)
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Define heuristic constants.

  const size_t min_component = 1000;
  
  // Define filenames.

  cout << Tag() << "beginning MergeNeighborhoods2" << endl;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String file_head = run_dir + "/" + READS;
  String kK = ".k" + ToString(K);
  String data_dir = PRE + "/" + DATA;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  // Load some stuff.

  cout << Tag() << "load some stuff" << endl;
  int ploidy = FirstLineOfFile(data_dir + "/ploidy").Int();
  size_t n_uni = MastervecFileObjectCount(file_head + ".unibases" + kK);
  vecKmerPath spaths;
  vec<tagged_rpint> spathsdb;
  spaths.ReadAll(sub_dir + "/reads.paths" + kK);
  BinaryRead2(sub_dir + "/reads.pathsdb" + kK, spathsdb);






  // Find all the components in all the neighborhoods.

  vec< vec<int> > comp;
  vec<HyperKmerPath> nhcomp;
  size_t n_edges = FindComponents(K, sub_dir, 
                                  & comp, & nhcomp);
    
  ForceAssertEq(n_edges + n_uni, spaths.size());
  size_t n_comp = comp.size();





  // Set up to compute HyperKmerPath components.  Scope.

  vec<HyperKmerPath> hcomps;
  size_t n_kmers_max = 0;
  {
    vec<int> to_rc(n_uni, -1);
    ConstructInvolution(n_uni, n_edges, spaths, spathsdb, & to_rc);


    // Index the component edges.
    vec<int> to_comp(n_edges);
    for (size_t i = 0; i < comp.size(); i++) {
      for (size_t j = 0; j < comp[i].size(); j++)
        to_comp[ comp[i][j] ] = i;
    }

    // Predict copy numbers.
    vec<int> predicted_CNs(n_uni, -1); 
    { 
      VecPdfEntryVec CNs((file_head + ".unipaths.predicted_count" + kK).c_str());
      for (size_t i = 0; i < n_uni; i++)
        GetMostLikelyValue(predicted_CNs[i], CNs[i]);
    }


    // identify high-copy-number edges
    vec<int> edge_CN(n_edges, 0);
    for(size_t i=0;i< spathsdb.size();i++)
    {
      tagged_rpint& t = spathsdb[i];
      uint id = t.ReadId();
      if (id < n_edges) continue; // only unipaths 
      if (predicted_CNs[ id - n_edges] <= ploidy) continue; // only HCN paths
      // find the edges that contain this unipaths
      size_t i1  = i - t.Lookback();
      for(size_t j=i1; j<i; j++)
      {
	tagged_rpint& t1 = spathsdb[j];
	uint id1 = t1.ReadId();
	if (id1 >= n_edges ) continue; // looking for the edge paths
	int overlap = IntervalOverlap(t1.Start(), t1.Stop(), t.Start(), t.Stop());
	if (overlap > 0) edge_CN[id1] += overlap;
      }
    }
    vec<bool> edge_repeat(n_edges, False);
    if (SKIP_REPEAT)
      for(size_t i=0;i<edge_CN.size();i++)
	if ( edge_CN[i] > MAX_REPEAT) edge_repeat[i] = True;


    // Find all cases where two components overlap.  Scan the paths database.
    vec<overlapit> overlaps;
    FindAllComponentOverlaps(ploidy, spathsdb, to_rc, to_comp, predicted_CNs, 
	                     sub_dir, & overlaps, edge_repeat);


    // Sort overlaps
    cout << Tag() << "sorting overlaps" << endl;
    ParallelUniqueSort(overlaps);



    // Mine the overlaps.  Note that the way this works now, if sequence x
    // appears once in c1 and twice in c2, it counts double (etc.), which
    // doesn't really make sense.
    MineComponentOverlaps(K, spaths, comp, nhcomp, overlaps, 
                          min_component, MIN_OVERLAP_TO_CLUSTER, 
                          & hcomps, & n_kmers_max);    

  }
  cout << Tag() << "n_kmers_max = " << ToStringAddCommas(n_kmers_max) << endl;








  KmerBaseBroker kbb;
  if (ACCURACY) {
    vecKmerPath spaths_rc;
    spaths_rc.ReadAll(sub_dir + "/reads.paths_rc" + kK);
    
    vecbvec sbases;
    sbases.ReadAll(sub_dir + "/reads.fastb");

    kbb.Initialize(K, sbases, spaths, spaths_rc, spathsdb);
  }





  // Merge HyperKmerPaths. 

  Mkdir777(sub_dir + "/MergeNeighborhoods2");
  String base_file = sub_dir + "/MergeNeighborhoods2/hyper_";
  String base_ref = PRE + "/" + DATA + "/genome";     
     
  cout << Tag() << "starting merger";

  size_t dotter = 1000;
  if (ACCURACY) 
    cout << endl;
  else 
    cout << " ("
            << hcomps.size() << " clusters, . = "
            << dotter << " clusters):"
            << endl;

  // Scope. Worklist destructor waits for all threads to complete.  
  if (true) {
    MergeProcessor proc(NUM_THREADS, MIN_OVERLAP, MIN_PROPER_OVERLAP, 
                        &hcomps, ACCURACY, base_file, base_ref, kbb, dotter);

    // Parallelize over n_threads instead of over the number of objects.
    // When the number of objects is too large and processing each object is 
    // very fast, you're better off partitioning the data yourself.
    WorklistN<MergeProcessor> worklist(proc, NUM_THREADS);
  }

  if (!ACCURACY) cout << endl;

  // Merge HyperKmerpaths and write.

  if (WRITE) {
    HyperKmerPath h(K, hcomps);
    cout << Tag() << "write assembly to file <sub_dir>/hyper.prelim" + 
      OUT_SUFFIX << endl;
    BinaryOverwrite(sub_dir + "/hyper.prelim" + OUT_SUFFIX, h);
    Ofstream(dot, sub_dir + "/hyper.prelim" + OUT_SUFFIX + ".dot");
    h.PrintSummaryDOT0w(dot, true, false, true, 0, true);
  }
     
  cout << Tag() << "done" << endl;
  _exit(0);
}
