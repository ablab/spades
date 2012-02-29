///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CloseUnipathGaps.  For each terminal unipath, from consensus sequence
// for the reads that extend it, and assign quality scores to the consensus bases,
// as Min(40, sum[0] - sum[1]), where sum[0] is the sum of quality scores for the
// winning base, and sum[1] is the sum of quality scores for the runner up base.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "lookup/FirstLookup.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/CloseUnipathGapsCore.h"
#include "paths/GetNexts.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"

// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

class extension {

     public:

     extension() { }

     extension(const int u, const basevector& bases, const qualvector& quals,
          const vec<longlong>& ids) : u(u), bases(bases), quals(quals), ids(ids) { }

     int u;
     basevector bases;
     qualvector quals;
     vec<longlong> ids;

};

void BuildExtenders(const vecbasevector& bases,
                     const vecqualvector& quals,
                     const vecbasevector& unibases,
                     const vec< vec<int> >& nexts,
                     const vec<int>& to_rc,
                     const VecULongVec& seeds,
                     const VecULongVec& seeds_rc,
                     const int K,
                     const int UNIBASES_K,
                     vec< triple<int,int,longlong> >& extenders,
                     vec<Bool>& fw, const int VERBOSITY, ostream& log)
{
     extenders.clear(), fw.clear();
     vec< triple<int,int,longlong> > extenders0;
     CloseUnipathGapsCore(bases, quals, unibases, nexts, to_rc, K, UNIBASES_K,
			     seeds, seeds_rc, extenders0, VERBOSITY, log);
     const int max_dist = 100;
     extenders.reserve(extenders0.isize());
     fw       .reserve(extenders0.isize());
     for (int i = 0; i < extenders0.isize(); i++)
     {    int u = extenders0[i].first, pos = extenders0[i].second;
          longlong id = extenders0[i].third;
	  int nu = unibases[u].size();
          int len = bases[id].size();
          int extends_by = len - (nu - pos);
          if (extends_by >= -max_dist)
          {    extenders.push_back(extenders0[i]);
               fw.push_back(True);    }
          extends_by = -pos;
          if (extends_by >= -max_dist)
          {    extenders.push(to_rc[u], nu - pos - len, id);
               fw.push_back(False);    }    }
     SortSync(extenders, fw);    }

void BuildExtensions(const vec< triple<int,int,longlong> >& extenders,
                      const vec<Bool>& fw,
                      const vecbasevector& bases,
                      const vecqualvector& quals,
                      const PairsManager& pairs,
                      const vecbasevector& unibases,
                      const vec<int>& to_rc,
                      vec<extension>& E,
                      ostream& log,
                      const int VERBOSITY)
{
     E.clear();
     E.reserve(extenders.isize());
     for (int i = 0; i < extenders.isize(); i++)
     {    int u = extenders[i].first;

          // Find the extensions of unipath u.

          int j;
          for (j = i + 1; j < extenders.isize(); j++)
               if (extenders[j].first != u) break;

	  if (VERBOSITY >= 1) {
	    log << "\n------------------------------------------------------------"
		<< "------------------------\n\n";
	    log << "looking off " << u << "[" << to_rc[u] << "]\n";
	  }
          vec< triple<int,int,longlong> > extl;
          vec<int> extends_by;
          vec<Bool> fw_by;
          for (int k = i; k < j; k++)
	  {    int pos = extenders[k].second;
	       longlong id = extenders[k].third;
               int len = bases[id].size(), nu = unibases[u].size();
               extends_by.push_back(len - (nu - pos));
               fw_by.push_back(fw[k]);
               extl.push_back(extenders[k]);    }
          ReverseSortSync(extends_by, fw_by, extl);
	  if (VERBOSITY >= 2)
          for (int k = 0; k < extl.isize(); k++)
	  {    int pos = extl[k].second;
	       longlong id = extl[k].third;
               PRINT5_TO(log, pos, id, pairs.getPartnerID(id),
                    extends_by[k], int(fw_by[k]));    }
          int nrows, ncols = Max(0, Max(extends_by));
          for (nrows = 0; nrows < extl.isize(); nrows++)
               if (extends_by[nrows] <= 0) break;
          vec< vec<char> > B(nrows), Q(nrows);
          basevector e(ncols);
          qualvector q(ncols);
          vec<longlong> ids;
          for (int k = i; k < j; k++)
               if (fw[k]) ids.push_back(extenders[k].third);
          UniqueSort(ids);
          if (ncols == 0)
	  {    if (VERBOSITY >= 1) log << "There are no extensions.\n";
               E.push(u, e, q, ids);
               i = j - 1;
               continue;    }
          for (int x = 0; x < nrows; x++)
          {    B[x].resize(ncols, ' ');
               Q[x].resize(ncols, ' ');    }
          for (int k = 0; k < nrows; k++)
	  {    int pos = extl[k].second;
	       longlong id = extl[k].third;
	       ForceAssertEq(extl[k].first, u);
               int len = bases[id].size(), l0 = unibases[u].isize() - pos;
               for (int l = l0; l < len; l++)
               {    if (fw_by[k])
                    {    B[k][l-l0] = bases[id][l];
                         Q[k][l-l0] = quals[id][l];    }
                    else
                    {    B[k][l-l0] = 3 - bases[id][ len - l - 1 ];
                         Q[k][l-l0] = quals[id][ len - l - 1 ];    }
		 if (VERBOSITY >= 2) {
                    if (Q[k][l-l0] >= 0) log << as_base(B[k][l-l0]);
                    else log << " ";
		 }
	       }
               if (VERBOSITY >= 2) log << "\n";
	  }
          vec<char> ext, extq;
          for (int l = 0; l < ncols; l++)
          {    vec<int> calls(4, 0);
               for (int r = 0; r < nrows; r++)
               {    if (B[r][l] == ' ') continue;
                    calls[ B[r][l] ] += Q[r][l];    }
               vec<char> bs;
               bs.push_back('A', 'C', 'G', 'T');
               ReverseSortSync(calls, bs);
	       if (VERBOSITY >= 2) {
		 for (int x = 0; x < bs.isize(); x++)
		   {    if (x > 0) log << " ";
		     log << bs[x] << "[" << calls[x] << "]";    }
		 log << "\n";
	       }
               ext.push_back(bs[0]);
               const int max_qual = 40;
               extq.push_back(Min(max_qual, calls[0] - calls[1]));    }
          for (int m = 0; m < ncols; m++)
          {    e.Set(m, as_char(ext[m]));
               q[m] = extq[m];    }

          // If the reads are highly inconsistent with the consensus, bail.

          vec<int> errsums;
          for (int k = 0; k < nrows; k++)
	  {    int pos = extl[k].second;
	       longlong id = extl[k].third;
	       ForceAssertEq(extl[k].first, u);
               int len = bases[id].size(), l0 = unibases[u].isize() - pos;
               int errsum = 0;
               for (int l = l0; l < len; l++)
               {    if (B[k][l-l0] != e[l-l0])
                         errsum += Min(Q[k][l-l0], extq[l-l0]);    }
               errsums.push_back(errsum);
               // PRINT2(id, errsum); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    }
          Sort(errsums);
          if (errsums[nrows/2] > 100)
          {    if (VERBOSITY >= 1)
               {    cout << "median errsum is " << errsums[nrows/2] << "\n";
                    cout << "rejecting\n";    }
               i = j - 1;
               continue;    }

          E.push(u, e, q, ids);
	  if (VERBOSITY >= 2) {
          for (int m = 0; m < ext.isize(); m++)
               log << ext[m];
          log << "\n";
          for (int m = 0; m < ext.isize(); m++)
          {    if (extq[m] < 10) log << int(extq[m]);
               else log << int(extq[m]) / 10;    }
          log << "\n";
          for (int m = 0; m < ext.isize(); m++)
          {    if (extq[m] >= 10) log << int(extq[m]) % 10;
               else log << " ";    }
          log << "\n";
	  }
          i = j - 1;    }    }


int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(DIR);
  CommandArgument_String_Doc(IN_HEAD,
			     "looks for DIR/IN_HEAD.{fastb,qualb} = reads");
  CommandArgument_String_Doc(UNIBASES, "looks for DIR/UNIBASES");
  CommandArgument_Int(UNIBASES_K);
  CommandArgument_Int_OrDefault_Doc(K, 20, "seed size");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_Doc(WORKDIR, "looks for DIR/WORKDIR");
  CommandArgument_String_OrDefault_Doc(EXT, "",
					"No effect; just here for backward compatibility.");
  CommandArgument_Int_OrDefault(VERBOSITY, 0);
  CommandArgument_String_OrDefault_Doc(REF, "",
				       "DIR/REF is lookup table for reference sequence, to assess joins");
  CommandArgument_String_OrDefault(OUT_HEAD, "");
  CommandArgument_Bool_OrDefault_Doc(WRITE, True,
				     "set this to false if you don't want any output to be written to files");
   CommandArgument_Bool_OrDefault_Doc(BUILD_UNIGRAPH, False,
	 "Build and write the unipath adjacency graph");
  CommandArgument_Bool_OrDefault(USE_LINKS, False);
  EndCommandArguments;



  // Thread control (use OMP)
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Check arguments.

  ForceAssertEq(K, 20);
  if (WRITE && OUT_HEAD == "") {
    cout << "If WRITE is True, OUT_HEAD " << "must also be specified." << endl;
    cout << "Abort." << endl;
    exit(1);
  }

  // Set up logging.

  String outhead = DIR + "/" + OUT_HEAD;
  String logfile = outhead + ".CloseUnipathGaps.log";

  String work = DIR + "/" + WORKDIR;
  Mkdir777(work);

  ofstream& log
    = (WRITE ? *(new ofstream(logfile.c_str())) : (ofstream&) cout);
  if (WRITE) Mkdir777(outhead.RevBefore("/"));
  if (WRITE) cout << "logging to " << logfile << endl;

  cout << Date() << ": Loading input files..." << endl;

  // Get the reads and ancillary info.

  vecbasevector bases(DIR + "/" + IN_HEAD + ".fastb");
  int nreads = bases.size();
  PairsManager pairs(DIR + "/" + IN_HEAD + ".pairs");
  
  // Get unibases and ancillary info.

  vecbasevector unibases(DIR + "/" + UNIBASES);
  int nuni = unibases.size();
  vec< vec<int> > nexts;
  GetNexts(UNIBASES_K, unibases, nexts);
  vec<int> to_rc;
  UnibaseInvolution(unibases, to_rc, UNIBASES_K);

  // Get links between unibases.  

  vec< vec<int> > jlinked;
  if (USE_LINKS)
  {    jlinked.resize(nuni);
       String line;
       fast_ifstream in( DIR + "/" + UNIBASES +  ".predicted_gaps.txt" );
       while(1)
       {    getline( in, line );
            if ( in.fail( ) ) break;
            if ( line.Contains( "#", 0 ) ) continue;
            int u1, u2, sep, dev, nlinks;
            istringstream iline( line.c_str( ) );
            iline >> u1 >> u2 >> sep >> dev >> nlinks;
            jlinked[u1].push_back(u2);    }    }

  // Find seeds for unipath gaps.

  VecULongVec seeds, seeds_rc;
  cout << Date() << ": Finding seeds for unipath gaps" << endl;
  FindUnipathGapSeeds(bases, unibases, K, NUM_THREADS, work, & seeds);
  ReverseComplement(bases);
  FindUnipathGapSeeds(bases, unibases, K, NUM_THREADS, work, & seeds_rc);
  ReverseComplement(bases);

  // Align to unipath gaps and build extensions.

  vec< triple<int,int,longlong> > extenders;
  vec<Bool> fw;
  cout << Date() << ": Loading read quality scores so we can build extenders"
       << endl;
  vecqualvector quals(DIR + "/" + IN_HEAD + ".qualb");
  cout << Date() << ": Building those extenders" << endl;
  BuildExtenders(bases, quals, unibases, nexts, to_rc, seeds, seeds_rc,
		  K, UNIBASES_K, extenders, fw, VERBOSITY, log);
  vec<extension> E;
  BuildExtensions(extenders, fw, bases, quals, pairs, unibases, to_rc, E, log, VERBOSITY);

  // Extend the unibases by Q20 stretches.

  vec<Bool> touched(nuni, False);
  for (int i = 0; i < E.isize(); i++) {
    extension& e = E[i];
    if (touched[e.u]) continue;
    qvec::size_type add;
    for (add = 0; add < e.quals.size(); add++)
      if (e.quals[add] < 20) break;
    unibases[e.u] = Cat(unibases[e.u], basevector(e.bases, 0, add));
    unibases[ to_rc[e.u] ] = unibases[e.u];
    unibases[ to_rc[e.u] ].ReverseComplement();
    touched[e.u] = touched[ to_rc[e.u] ] = True;
  }

  // Align to unipath gaps and build extensions.

  BuildExtenders(bases, quals, unibases, nexts, to_rc, seeds, seeds_rc, K,
		  UNIBASES_K, extenders, fw, VERBOSITY, log);
  BuildExtensions(extenders, fw, bases, quals, pairs, unibases, to_rc, E, log,
		   VERBOSITY);

  cout << Date() << ": Prepare for merging paths" << endl;

  // Pre-processing step:
  // Create a lookup table of <read ID> -> <extensions containing that read>.
  // This allows us to use a loop over read IDs [O(n)] rather than pairs of
  // extensions [O(n^2)].
  // TODO:  potentially dangerous truncation of IDs, unless n_exts is known to be < 2^31
  VecIntVec reads_to_exts(nreads);
     
  for (size_t i = 0; i < E.size(); i++)
    for (size_t j = 0; j < E[i].ids.size(); j++)
      reads_to_exts[ E[i].ids[j] ].push_back(i);

  // Look for pairs of extensions.  Merge where possible.
  cout << Date() << ": Merging paths, where possible" << endl;

  // Loop over all read IDs.  For each read, find the set of extensions that
  // contain this read, and also the set that contains this read's partner;
  // then attempt to merge the extensions pairwise.
  vecbasevector merged_paths;

  for (longlong read_ID = 0; read_ID < nreads; read_ID++) {
    longlong partner_ID = pairs.getPartnerID(read_ID);
    if (partner_ID < read_ID) continue; // this catches cases of partner = -1
    IntVec const& E1 = reads_to_exts[read_ID];
    IntVec const& E2 = reads_to_exts[partner_ID];
    for (IntVec::size_type i1 = 0; i1 < E1.size(); i1++) {
      for (IntVec::size_type i2 = 0; i2 < E2.size(); i2++) {
	const extension &e1 = E[ E1[i1] ], &e2 = E[ E2[i2] ];
	if (e2.u == e1.u || e2.u == to_rc[e1.u]) continue;
	if (VERBOSITY >= 1) {
	  log << "\n----------------------------------------------------------"
	      << "--------------------------\n";
	  log << "\ntrying to join " << e1.u << " to " << e2.u << "\n";
	}

        if ( USE_LINKS && !Member( jlinked[e2.u], to_rc[e1.u] ) )
        {    log << "There are not linked!\n";
             const int min_to_require_link = 1500;
             if ( unibases[e1.u].isize( ) + unibases[e2.u].isize( ) 
                  >= min_to_require_link )
             {    log << "Giving up, lengths = " << unibases[e1.u].size( )
                       << " and " << unibases[e2.u].size( ) << "\n";    
                  continue;    }    }

	if (nexts[e1.u].nonempty() &&
	     nexts[ to_rc[e2.u] ].nonempty()) {
	  if (VERBOSITY >= 1) {
	    log << "but " << e1.u << " is followed by another "
		<< "unipath, so this join would make no sense\n";
	    log << "Not sure how we got here!\n";
	  }
	  continue;
	}
	/*
	  if (nexts[ to_rc[e2.u] ].nonempty())
	  {    if (VERBOSITY >= 1)
	  {    log << "but " << e2.u << " is preceded by another "
	  << "unipath, so this join would make no sense\n";
	  log << "Not sure how we got here!\n";    }
	  continue;    }
	*/
	basevector u1plus = Cat(unibases[e1.u], e1.bases);
	qualvector q1plus(unibases[e1.u].size() + e1.quals.size());
	for (int i = 0; i < unibases[e1.u].isize(); i++)
	  q1plus[i] = 40;
	for (qvec::size_type i = 0; i < e1.quals.size(); i++)
	  q1plus[ i + unibases[e1.u].size() ] = e1.quals[i];
	basevector u2plus = Cat(unibases[e2.u], e2.bases);
	qualvector q2plus(unibases[e2.u].size() + e2.quals.size());
	for (int i = 0; i < unibases[e2.u].isize(); i++)
	  q2plus[i] = 40;
	for (qvec::size_type i = 0; i < e2.quals.size(); i++)
	  q2plus[ i + unibases[e2.u].size() ] = e2.quals[i];
	u2plus.ReverseComplement();
	q2plus.ReverseMe();
	int n1 = u1plus.size(), n2 = u2plus.size();
	int max_overlap = UNIBASES_K + e1.bases.size() + e2.bases.size();
	max_overlap = Min(max_overlap, n1, n2);
	basevector u1plust, u2plust;
	u1plust.SetToSubOf(u1plus, n1 - max_overlap, max_overlap);
	u2plust.SetToSubOf(u2plus, 0, max_overlap);
	if (VERBOSITY >= 2) {
	  u1plust.Print(log, "u1plus_trunc");
	  u2plust.Print(log, "u2plus_trunc");
	}
	int maxscore = -1000000000;
	vec<int> bests;
	for (int pass = 1; pass <= 2; pass++) {
	  for (int o = 0; o <= max_overlap; o++) {
	    int plus = 0, minus = 0;
	    for (int j = 0; j < o; j++) {
	      int s = Min(q1plus[n1-o+j], q2plus[j]);
	      if (u1plus[n1-o+j] == u2plus[j]) plus += s;
	      else minus += s;    }
	    int score = plus - minus;
	    if (pass == 1) {
	      if (score == maxscore) bests.push_back(o);
	      if (score > maxscore) {
		maxscore = score;
		bests.clear();
		bests.push_back(o);
	      }
	    }
	    if (pass == 2 && score >= maxscore - 200 && VERBOSITY >= 2)
	      PRINT4_TO(log, o, plus, minus, score);
	  }
	}
	if (VERBOSITY >= 1) {
	  log << "\nbest score = " << maxscore << ", o = {";
	  for (int j = 0; j < bests.isize(); j++) {
	    if (j > 0) log << ",";
	    log << bests[j];
	  }
	  log << "}\n";
	}
	
	// Join u1_plus to u2_plus along overlap of size o.
	
	int o = bests[0];
	basevector join(n1+n2-o);
	for (int i = 0; i < n1; i++)
	  join.Set(i, u1plus[i]);
	for (int i = 0; i < n2; i++)
	  join.Set(i+n1-o, u2plus[i]);
	for (int i = 0; i < o; i++) {
	  if (q1plus[n1-o+i] >= q2plus[i])
	    join.Set(n1-o+i, u1plus[n1-o+i]);
	  else join.Set(n1-o+i, u2plus[i]);    }
	if (REF != "") {
	  Ofstream(out, work + "/join.fasta");
	  join.Print(out, "join");
	  out.close();
	  log << "\nalignment of best join:\n";
	  SystemSucceed("QueryLookupTable" + ARG(K, 12) + ARG(MM, 12)
			 + ARG(MC, 0.15) + ARG(NH, True) + ARG(VISUAL, True)
			 + ARG(SEQS, work + "/join.fasta")
			 + ARG(QUIET, True) + ARG(L, DIR + "/" + REF)
			 + " > " + work + "/join.log");
	  CpAppend(work + "/join.log", log);
	}
	
	// Align extended unipaths.
	
	if (REF != "") {
	  log << "alignment of extended unipaths:\n";
	  Ofstream(out, work + "/test.fasta");
	  u1plus.Print(out, "u1plus");
	  u2plus.Print(out, "u2plus");
	  out.close();
	  SystemSucceed("QueryLookupTable" + ARG(K, 12) + ARG(MM, 12)
			 + ARG(MC, 0.15) + ARG(NH, True) + ARG(VISUAL, True)
			 + ARG(SEQS, work + "/test.fasta")
			 + ARG(QUIET, True) + ARG(PARSEABLE, True)
			 + ARG(L, DIR + "/" + REF) + " > " + work + "/test.qltout");
	  SystemSucceed("cat " + work + "/test.qltout | grep -v QUERY"
			 + " > " + work + ".qltout.filt");
	  CpAppend(work + ".qltout.filt", log);
	  vec<look_align> aligns;
	  LoadLookAligns(work + "/test.qltout", aligns);
	  if (aligns.size() == 2 && aligns[0].query_id == 0
	       && aligns[1].query_id == 1) {
	    log << "actual overlap = "
		<< Overlap(aligns[0].Extent2(), aligns[1].Extent2())
		<< "\n";
	  }
	}
	
	// Keep join if trusted.

	const int min_score = 50;
	if (maxscore >= min_score) {
	  if (VERBOSITY >= 1) {
	    log << "\nJOINING!\n";
	    join.Print(log, "join");
	  }
	  merged_paths.push_back_reserve(join);
	}
      }
    }
  }

  // Now build the new unipaths.
  // We can skip this step (to save runtime) if not writing any output files.

  if (WRITE) {
    Destroy(bases), Destroy(quals);
    String KS2 = ToString(UNIBASES_K);

    cout << Date() << ": Building new unipaths" << endl;
    vecbasevector all(unibases);
    all.Append(merged_paths);
    for (int id1 = 0; id1 < nuni; id1++) {
      for (int j = 0; j < nexts[id1].isize(); j++) {
	int id2 = nexts[id1][j];
	basevector b = unibases[id1];
	b.resize(b.size() + 1);
	b.Set(b.size() - 1, unibases[id2][UNIBASES_K-1]);
	all.push_back_reserve(b);
      }
    }
    vecKmerPath newpaths, newpathsrc, newunipaths;
    vec<tagged_rpint> newpathsdb, newunipathsdb;
    digraph unigraph;
    ReadsToPathsCoreY(all, UNIBASES_K, newpaths, newpathsrc, newpathsdb,
		      work + "/CUGaps2", NUM_THREADS);
    Unipath(newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb);
    KmerBaseBroker newkbb(UNIBASES_K, newpaths, newpathsrc, newpathsdb, all);
    vecbasevector newunibases;
    for (size_t i = 0; i < newunipaths.size(); i++)
      newunibases.push_back_reserve(newkbb.Seq(newunipaths[i]));
    BuildUnipathAdjacencyGraph( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb, unigraph );
 
    // Write output files.
    cout << Date() << ": Writing output files" << endl;
    newunipaths.WriteAll(outhead + ".unipaths.k" + KS2);
    BinaryWrite( outhead + ".unipath_adjgraph.k" + KS2, unigraph );
    BinaryWrite3(outhead + ".unipathsdb.k" + KS2, newunipathsdb);
    newunibases.WriteAll(outhead + ".unibases.k" + KS2);

    /*
      cout << "\nNEW UNIBASES:\n\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      for (int i = 0; i < newunipaths.size(); i++) // XXXXXXXXXXXXXXXXXXXXXXX
      newunibases[i].Print(cout, "new_unibases_" + ToString(i)); // XXXXX
    */

    nuni = newunibases.size();
    vecbasevector all2(newunibases);
    GetNexts(UNIBASES_K, newunibases, nexts);
    for (int id1 = 0; id1 < nuni; id1++) {
      for (int j = 0; j < nexts[id1].isize(); j++) {
	int id2 = nexts[id1][j];
	basevector b = newunibases[id1];
	b.resize(b.size() + 1);
	b.Set(b.size() - 1, newunibases[id2][UNIBASES_K-1]);
	all2.push_back_reserve(b);
      }
    }
    all2.WriteAll(outhead + ".unibases.plus_junctions.k" + KS2);
    
  }

  cout << Date() << ": Done!" << endl;    
}
