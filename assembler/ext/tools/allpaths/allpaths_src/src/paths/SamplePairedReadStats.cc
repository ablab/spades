///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: SamplePairedReadStats

   Determines separation stats, orientation, etc. of paired reads.
   This is done by aligning a random sample of pairs to a reference

   INPUT FILES:
     reference.fastb
     (paired reads).pairs

   OUTPUT FILES:
     READS_HEAD.outies (stats for incorrectly oriented pairs)
     if OUT_HEAD nonempty then following files are written
     OUT_HEAD.pairs
     if WRTE_STATS writes READS_HEAD.sample_stats containing detailed stats
     
   @file

*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/FirstLookupFinderECJ.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTabBuilder.h"
#include "random/Shuffle.h"
#include "paths/UnibaseUtils.h"
#include "VecUtilities.h"
#include <map>
#include "reporting/PerfStat.h"
#include "paths/GetNexts.h"
#include "paths/Unipath.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"

// auxiliary routines    -------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, const VecQualNibbleVec &quals,
			  FirstLookupFinderECJ& lfinder, FirstLookupFilterECJ& lfilter,
			  const longlong pi, const int LookupK,
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, int TRIM, int64_t& sep , Bool verbose);
double gm( vec<float> data, double& mean, double& sigma, const Bool verbose);
void actual_gm( const vec<float>& data, vec<double>& means, vec<double>& sigmas, vec<double>& weights, const Bool verbose );
//----------------------------------------------------------

// Fraction of seps falling within N sigma of mean. --bruce
double FractionInDist(const vec<float> &data, const double mean, const double sigma, const double N=3)
{
  ForceAssertGt(sigma, 0.0);
  int count = 0;
  for (u_int i = 0; i < data.size(); ++i) {
    double x = (data[i]-mean) / sigma;
    if (abs(x) < N) ++count;
  }
  return (double)count / (double)data.size();
}

int N50hbv( HyperBasevector& hbv ){
  vec<int> lens( hbv.Edges().size() );
  for ( size_t u = 0; u < hbv.Edges().size(); u++ )
    lens[u] = hbv.Edges()[u].isize();
  PRINT( Max(lens) );
  return N50(lens);
}

void convert_to_unibases( vecbasevector& seqs, const size_t K, int n_threads ){
  
  cout << Date() << ": pathing" << endl;
  vecKmerPath paths, paths_rc;
  vec<tagged_rpint> pathsdb;
  ReadsToPathsCoreY( seqs, K, paths, paths_rc, pathsdb, "", n_threads);
  
  cout << Date() << ": unipathing" << endl;
  vecKmerPath unipaths;
  vec<tagged_rpint> unipathsdb;
  Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb, True);
  
  cout << Date() << ": building KmerBaseBroker (needed for unibases)" << endl;
  KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, seqs );
  
  
  cout << Date() << ": building unibases" << endl;
  vecbasevector unibases;
  unibases.reserve(unipaths.size());
  for ( size_t i = 0; i < unipaths.size(); i++ )
    unibases.push_back( kbb.Seq( unipaths[i] ) );
  seqs = unibases;
  
}

void flatten_unibases( vecbasevector& unibases, const size_t K ){
  cout << Date() << ": Building unibaseses adjacency graph" << endl;
  PRINT( unibases.size() );


  digraph AG;
  BuildUnibaseAdjacencyGraph( unibases, AG, K );
  cout << Date() << ": building hyperbasevector" << endl;
  HyperBasevector hbv;
  BuildUnibaseAdjacencyHyperBasevector( K -1, AG, unibases, hbv );
  AG.Clear();
  PRINT2( hbv.N(), hbv.TotalEdgeLength() );
  PRINT2( N50hbv( hbv ), hbv.N() );
 
  for ( size_t u = 0; u < unibases.size(); u++ ){
    if ( unibases[u] != hbv.Edges().at(u) )
      FatalErr(" unibases do not correspond to edges");

  }

  vec<int> toRc;
  UnibaseInvolution( unibases, toRc, K );


  ForceAssertEq( unibases.size(), hbv.Edges().size() );
  vec<Bool> uniUsed( unibases.size(), False );
  vec< vec<int> > comps;
  hbv.Components(comps);
  PRINT( comps.size() );
  
  vec<int> keep;
  for ( size_t ic = 0; ic < comps.size(); ic++ ){
    const vec<int> c = comps[ic];
    Bool anyUsed = False;
    for ( size_t j = 0; j < c.size(); j++ ){
      int v = c[j];
      for ( size_t t = 0; t < hbv.FromEdgeObj(v).size(); t++ ){
	ForceAssert( hbv.FromEdgeObj(v)[t] >= 0);
	ForceAssert( hbv.FromEdgeObj(v)[t] < hbv.Edges().isize() );
	if ( uniUsed.at( hbv.FromEdgeObj(v)[t] ) ){
	  anyUsed = True;
	  break;
	}
      }
    }
    if ( ! anyUsed )
      keep.append(c);
    
    for ( size_t j = 0; j < c.size(); j++ ){
      int v = c[j];
      for ( size_t t = 0; t < hbv.FromEdgeObj(v).size(); t++ ){
	//PRINT3( v, t, hbv.FromEdgeObj(v)[t] );
	uniUsed.at( hbv.FromEdgeObj(v)[t] ) = True;
	uniUsed.at( toRc.at( hbv.FromEdgeObj(v)[t] ) ) = True;
      }      
    }
  }
  HyperBasevector h( hbv, keep );
  hbv = h; 
 


  cout << Date() << ": removing small components\n" << endl;
  hbv.RemoveSmallComponents( 100 );

  hbv.RemoveUnneededVertices();
  hbv.RemoveDeadEdgeObjects();
  hbv.RemoveEdgelessVertices();

  PRINT2( hbv.N(), hbv.TotalEdgeLength() );
  PRINT2( N50hbv( hbv ), hbv.N() );
  
  cout << Date() << ": popping bubbles" << endl;
  int npopped = 0;
  while( 1 ) {
    vec<int> bubbleLocs; bubbleLocs.reserve( hbv.N() );
    for ( int v = 0; v < hbv.N(); v++ ){
      if ( hbv.From(v).size() == 2 && hbv.From(v)[0] == hbv.From(v)[1] ){
	int w = hbv.From(v)[0];
	if ( hbv.To(v).size() > 2 || hbv.From(w).size() > 2 )
	  continue;
	int u0 = hbv.EdgeObjectIndexByIndexFrom( v, 0 );
	int u1 = hbv.EdgeObjectIndexByIndexFrom( v, 1 );
	if ( unibases[u0].size() == unibases[u1].size() )
	  bubbleLocs.push_back( v );
      }
    }
    if ( bubbleLocs.size() == 0 )
      break;
    npopped += bubbleLocs.size();
    hbv.PopBubbles( bubbleLocs );
  }
  PRINT( npopped );
  
  PRINT2( N50hbv( hbv ), hbv.N() );
  
  {
    vec<int> to_delete;
    // remove incoming
    for ( int v = 0; v < hbv.N(); v++ ){
      if ( hbv.To(v).size() == 0 && hbv.From(v).size() == 1 ){ 
	size_t w = hbv.From(v)[0];
	if ( hbv.To(w).size() > 1 ){
	  for ( size_t s = 0; s < hbv.To(w).size(); s++ ){
	    if ( hbv.To(s).size() > 0 ){ 
	      for ( size_t i = 0; i < hbv.EdgesBetween( v, w ).size(); i++ ){
		size_t uid = hbv.EdgesBetween(v,w)[i];
		if ( hbv.Edges()[uid].size() < 200 )
		  to_delete.push_back( uid );
	      }
	      break;
	    }
	  }
	}
      }
    }
    int nBadIncoming = to_delete.size();
    // remove outgoing 
    for ( int w = 0; w < hbv.N(); w++ ){
      if ( hbv.From(w).size() == 0 && hbv.To(w).size() == 1 ){
	size_t v = hbv.To(w)[0];
	if ( hbv.From(v).size() > 1 ){
	  for ( size_t s = 0; s < hbv.From(v).size(); s++ ){
	    if ( hbv.From(s).size() > 0 ){
	      for ( size_t i = 0; i < hbv.EdgesBetween( v, w ).size(); i++ ){
		size_t uid = hbv.EdgesBetween(v,w)[i];
		if ( hbv.Edges()[uid].size() < 200 )
		  to_delete.push_back( uid );
	      }
	      break;
	    }
	  }
	}
      }
    }
    
    cout << "Shaving function identified " << to_delete.size()
	 << " edges to delete (" << nBadIncoming
	 << " incoming and " << to_delete.size() -nBadIncoming
	 << " outgoing)" << endl;
    
    hbv.DeleteEdges( to_delete );
    hbv.RemoveUnneededVertices();
    hbv.RemoveDeadEdgeObjects();
    hbv.RemoveEdgelessVertices();
  }
  PRINT3( hbv.N(), hbv.TotalEdgeLength(), hbv.EdgeObjectCount() );
  PRINT2( N50hbv( hbv ), hbv.N() );
  unibases.resize( 2 * hbv.Edges().size() );
  for ( size_t i = 0; i < hbv.Edges().size(); i++ ){
    unibases[i] = hbv.Edges()[i];
    unibases[ hbv.Edges().size() + i ] = hbv.Edges()[i];
    unibases[ hbv.Edges().size() + i ].ReverseComplement();
  }
  PRINT( unibases.size() );
  unibases.UniqueSort();
  PRINT2( "after sort", unibases.size() );
  cout << "Number of flattened edges = " << unibases.size() << endl;
}

int main( int argc, char *argv[] ){
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault_Doc(REF,"", 
    "looks for RUN_DIR/REF.fastb");
  CommandArgument_String_OrDefault_Doc(UNIBASES,"", 
    "looks for RUN_DIR/UNIBASES as a reference");
  CommandArgument_UnsignedInt_OrDefault_Doc(UNIBASES_K,0, 
    "needed by involution algorithm");
  CommandArgument_String_Doc(READS_HEAD, 
    "looks for RUN_DIR/READS_HEAD.{fastb, pairs/pairs}");
  CommandArgument_Bool_OrDefault_Doc( FLIP, True,
    "reverse-complement reads for alignment; use for jumps, not frags.");
  CommandArgument_String_OrDefault_Doc(OUT_HEAD, "", 
    "Name for output files. No output files written if not specified");
  CommandArgument_Bool_OrDefault_Doc(WRITE_STATS, True, 
    "Write sample stats into READS_HEAD.sample_stats file.");
  CommandArgument_Bool_OrDefault_Doc(WRITE_SEPS, True, 
    "Write sampled separations into READS_HEAD.seps.{in,out} files.");
  CommandArgument_Int_OrDefault_Doc(ORIG_DEV_BRACKET, 5, 
    "Radius (in std dev) of data censoring around original (lab specified) mean paired read separation.");

  CommandArgument_Bool_OrDefault_Doc( NEW_LOOKUP_TABLE, True,
    "build new lookup table even if a file already exists.");
  CommandArgument_Bool_OrDefault_Doc( KEEP_LOOKUP_TABLE, False,
    "erase lookup table after use.");

  CommandArgument_Int_OrDefault_Doc(TRIM, 0, 
    "Trim that many bases from the beginning of a read (if negative then from the end) before aligning.");
  CommandArgument_Int_OrDefault_Doc(TARGET_SAMPLE_SIZE, 2000, 
    "sample size");
  CommandArgument_Double_OrDefault_Doc(TARGET_UNIBASE_COVERAGE, .2,
    "unibases covering this fraction of total are enough");
  CommandArgument_Int_OrDefault_Doc(MIN_SAMPLE_SIZE, 100, 
    "minimum sample size");
  CommandArgument_Int_OrDefault_Doc(MAX_SEPARATION, 100000, 
    "maximum aligned pair separation considered");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_UNIBASE_LENGTH, 0, 
    "only include unibases of this length or longer for alignments");
  CommandArgument_Bool_OrDefault_Doc(FLATTEN, False, 
    "Flatten unibase graph to get longer unibases");
  CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 133333,
    "Seed value for the random generator." );
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_Bool_OrDefault_Doc(WRITE_NOTHING, False, 
    "Don't write any ouput files.");
  EndCommandArguments;
  
  ForceAssertGe( TARGET_SAMPLE_SIZE, MIN_SAMPLE_SIZE );
  ForceAssertGt( TARGET_SAMPLE_SIZE, 0 );

  ForceAssert( !REF.empty() || !UNIBASES.empty() );
  if ( !UNIBASES.empty() )
    ForceAssertNe( UNIBASES_K, 0u );

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  int LookupK = 12;

  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  
  String full_ref_head          = run_dir + "/" + REF;
  String ref_fastb_file         = full_ref_head + ".fastb";
  
  String full_reads_head        = run_dir + "/" + READS_HEAD;
  String reads_fastb            = full_reads_head + ".fastb";
  String quals_file             = full_reads_head + ".qualb";
  String pairs_file             = full_reads_head + ".pairs";
  String pairto_file            = full_reads_head + ".pairto";

  String results_file           = full_reads_head + ".sample_stats";
  String seps_file              = full_reads_head + ".sample_seps";
  String outies_file            = full_reads_head + ".outies";

  vecbasevector ref;
  String out_lookup_file; 
  if ( UNIBASES.empty() ){
    cout << "Using reference file" << endl;
    out_lookup_file = full_ref_head + ".lookuptab";
    if ( ! IsRegularFile( out_lookup_file ) || NEW_LOOKUP_TABLE ) {
      cout << Date( ) << ": creating and saving lookup table " << out_lookup_file << endl;
      LookupTabBuilder lookup( LookupK );
      lookup.addFastb( ref_fastb_file.c_str( ) );
      lookup.write( out_lookup_file.c_str( ) );
    }
    cout << Date() << ": reading reference" << endl;
    ref.ReadAll( ref_fastb_file );
  } else {
    cout << "Using unibases as reference" << endl;
    String uni_fastb_file = run_dir + "/" + UNIBASES;
    String strand_uni_fastb_file = uni_fastb_file + ".oneway.tmp2";
    out_lookup_file = strand_uni_fastb_file + ".lookuptab";
    if ( ! IsRegularFile( uni_fastb_file ) || NEW_LOOKUP_TABLE ) {
      cout << Date() << ": loading unibases" << flush;
      vecbasevector unibases(uni_fastb_file);
      cout << unibases.size() << " loaded" << endl;

      
      if ( FLATTEN) {
	flatten_unibases( unibases, UNIBASES_K );
	cout << unibases.size() << " unibases after flattening" << endl;
	convert_to_unibases( unibases, UNIBASES_K, NUM_THREADS );
	cout << unibases.size() << " after repathing" << endl;
      }

      vec<int> toRc;
      cout << Date() << ": removing rc copies" << endl;
      UnibaseInvolution( unibases, toRc, UNIBASES_K );
      
      cout << Date() << ": building unibase reference file" << endl;
      uint64_t ulensSum = 0; int ulongest = 0;
      vec<int> suids( unibases.size(), vec<int>::IDENTITY );
      vec<int> ulens( unibases.size() );
      vec<int> ustatus( unibases.size(), 1 );
      for ( size_t i = 0; i < unibases.size(); i++ ){
	ulens[i] = unibases[i].size();
	ulensSum += ulens[i];
	if ( ulongest < ulens[i] )
	  ulongest = ulens[i];
	if ( ustatus[i] == 1 )
	  ustatus[ toRc[i] ] = -1; //don't use rc, we will not worry about palindromes
      }
      vec<int> sulens( ulens );
      
      PRINT3( N50(ulens), Max(ulens), Min(ulens) );

      ReverseSortSync( sulens, suids ); 
      
      uint64_t approxLength = (ulensSum - unibases.size() * UNIBASES_K) / 2;
      uint64_t targetLength = round(approxLength * TARGET_UNIBASE_COVERAGE);
      
      cout << Date() << ": total unibase length " << approxLength
	   << ", targetting " << (ToString(100*TARGET_UNIBASE_COVERAGE))
	   << "% (" << targetLength << ")" << endl;
      
      uint64_t targetLensSum = 0;
      vecbasevector strand_unibases;
      for ( size_t i = 0; i < suids.size(); i++ ){
	size_t uid = suids[i];
	if ( ustatus[uid] > 0 && unibases[uid].size() > (max(MIN_UNIBASE_LENGTH, 2 * UNIBASES_K))) {
	  strand_unibases.push_back( unibases[uid] );
	  targetLensSum += ulens[uid];
	  if ( targetLensSum >= targetLength)
	    break;
	}
      }
      cout << Date() << ": selected " << strand_unibases.size() << " unibases of size range "
	   << strand_unibases.back().size() << "-" << strand_unibases[0].size()
	   << " covering " << targetLensSum << " bases" << endl;
      
      strand_unibases.WriteAll( strand_uni_fastb_file );
    }
    if ( ! IsRegularFile( out_lookup_file ) || NEW_LOOKUP_TABLE ){
      cout << Date() << ": building lookup table" << endl;
      LookupTabBuilder lookup( LookupK );
      lookup.addFastb( strand_uni_fastb_file.c_str( ) );
      cout << Date() << ": writing lookup table" << endl;
      lookup.write( out_lookup_file.c_str( ) );
    }
    cout << Date() << ": reading reference" << endl;
    ref.ReadAll( strand_uni_fastb_file );
    // Erase tmp file
    Remove(strand_uni_fastb_file);
  }

  // find longest contig/unipath in the reference
  int maxRefTigLen = 0;
  for ( size_t u = 0; u < ref.size(); u++ )
    maxRefTigLen = ref[u].isize() > maxRefTigLen ? ref[u].isize() : maxRefTigLen;
  
  cout << Date( ) << ": loading lookup table " << out_lookup_file << endl;
  LookupTab lookup_tab( out_lookup_file.c_str( ) );

  size_t nreads = MastervecFileObjectCount( reads_fastb );
   
  cout << Date() << ": loading pairs" << endl;
  PairsManager pairs;
  if ( IsRegularFile( pairs_file ) )
    pairs.Read( pairs_file );
  else if ( IsRegularFile( pairto_file ) )
    pairs.ReadFromPairtoFile( pairto_file, nreads );
  else
    FatalErr(" ERROR: did not find pairs file ");

  uint64_t npairs   = pairs.nPairs();
  size_t nLibraries = pairs.nLibraries();
  cout << Date( ) << ": found " << npairs << " pairs in " << nLibraries << " libraries" << endl;
  vec<String> titles;
  titles.push_back("libNo", "origSep", "origDev");
  titles.push_back("newInSep", "newInDev", "3sigma%");
  titles.push_back("newOutSep", "newOutDev");
  titles.push_back("nLibraryPairs", "nPairsSampled", "nRead1Algnd", "nRead2Algnd");
  titles.push_back("nPairsAlgnd", "nPairsUniq", "nPairsSameContig", "nPairsLogical", "nPairsInnies", "nPairsOuties");
  titles.push_back("nSampleInnies", "nSampleOuties");
  vec< map<String,longlong> > libStats( nLibraries ); // for storing results
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    libStats[lib]["readLen"] = 0;
    for ( size_t ti = 0; ti < titles.size(); ti++ )
      libStats[ lib ][ titles[ti] ] = 0;
  }

  vec<String> libID2libName( nLibraries );
  for (size_t lib = 0; lib < nLibraries; lib++ ){
    libStats[ lib ]["origSep"] = pairs.getLibrarySep( lib );
    libStats[ lib ]["origDev"] = pairs.getLibrarySD( lib );
    libID2libName[ lib ]       = pairs.getLibraryName( lib );
    libStats[ lib ]["libNo"]   = lib;
  }
  cout << Date() << ": computing pairs in each library" << endl;
  vec< vec<size_t> > lib_pairs(nLibraries);
  for (size_t pi = 0; pi < npairs; pi++ ) {
    libStats[pairs.libraryID(pi)]["nLibraryPairs"]++;
    lib_pairs[pairs.libraryID(pi)].push_back(pi);
  }
  
  srand(RANDOM_SEED);
  for (size_t lib = 0; lib < nLibraries; lib++ ){
    random_shuffle(lib_pairs[lib].begin(), lib_pairs[lib].end());
  }

  cout << Date() << ": reading reads from " << reads_fastb << endl;
  vecbasevector reads( reads_fastb );

  cout << Date() << ": reading quals from " << quals_file << endl;
  VecQualNibbleVec quals;
  LoadQualNibbleVec(quals_file, &quals);
  ForceAssertEq(nreads, quals.size());

  // set alignment filtering parameters
  FirstLookupFilterECJ lfilter;
  lfilter.orientation = FirstLookupFilterECJ::ALL;
  lfilter.min_size = 20;
  lfilter.max_kmer_freq = 10000;
  lfilter.max_extend = lfilter.max_kmer_freq;
  lfilter.score_delta = 20;
  lfilter.score_max = 100;
  lfilter.min_match = 20;
  lfilter.mismatch_threshhold = 3;
  lfilter.mismatch_neighborhood = 8;
  lfilter.mismatch_backoff = 3;
  lfilter.max_placements = 2; // we only care if there are multiple plausible placements
  FirstLookupFinderECJ lfinder(lfilter, lookup_tab, ref, UNIBASES_K);

  cout << Date() << " Computing separations for libraries:" << endl;
  vec< vec< float> >   pairSepDataInnies( nLibraries );
  vec< vec< float> >   pairSepDataOuties( nLibraries );
  vec< vec< float> >   readLengthData( nLibraries );
  
  vec< Bool > found_data( nLibraries, False );
  
#pragma omp parallel for
  for (size_t lib = 0; lib < nLibraries; lib++ ) {
    if ( libStats[lib]["origSep"] >  maxRefTigLen )
      continue;
    for ( size_t  mi= 0;  mi < lib_pairs[lib].size() && !found_data[lib]; mi++ ) {
        size_t pi = lib_pairs[lib][mi];
      ForceAssertEq(lib, (size_t)pairs.libraryID(pi));
    
      //if ( mi % (perm.size()/100) == 0 )
      //Dot( cout, mi / (double)perm.size() * 100.0);
      //if ( mi == perm.size() -1 ) cout << endl;


      if ( found_data[ lib ] ) continue;

      if ( pairSepDataInnies[ lib ].isize() >= TARGET_SAMPLE_SIZE ){
	double outiePercLoc = (double)libStats[lib]["nPairsOuties"] / 
	  ( double)(libStats[lib]["nPairsOuties"] + 
		    libStats[lib]["nPairsInnies"] ) * 100.0;
      
	if ( pairSepDataOuties[ lib ].isize() >= MIN_SAMPLE_SIZE ||
	     outiePercLoc < 1.0){
	  if (VERBOSE) {
	    cout << "\nfound enough data for "; PRINT2( lib, libID2libName[lib] );
	    PRINT3(  pairSepDataInnies[ lib ].isize(), 
		     pairSepDataOuties[ lib ].isize(), outiePercLoc );
	  } else cout << " " << lib << flush;
	}
	found_data[ lib ] = True; 
      }

      libStats[ lib ]["nPairsSampled"]++;
    
      //lfilter.min_size = Min( reads[ pairs.ID1(pi) ].size(), reads[ pairs.ID2(pi) ].size() );
      Bool Read1Algnd = False, Read2Algnd = False, PairAlgnd = False, PairUnique = False,
	PairSameContig = False, PairLogical = False, PairInnie = False, PairOutie = False;
      int64_t sep = MAX_SEPARATION + 1;  
      pair_alignment_data( pairs, reads, quals, lfinder, lfilter,
			   pi, LookupK,
			   Read1Algnd, Read2Algnd, PairAlgnd, PairUnique, 
			   PairSameContig, PairLogical, PairInnie, PairOutie,
			   FLIP, TRIM, sep, VERBOSE );

      if ( Read1Algnd )     libStats[ lib ]["nRead1Algnd"]++;
      if ( Read2Algnd )     libStats[ lib ]["nRead2Algnd"]++;
      if ( PairAlgnd )      libStats[ lib ]["nPairsAlgnd"]++;
      if ( PairUnique )     libStats[ lib ]["nPairsUniq"]++;
      if ( PairSameContig ) libStats[ lib ]["nPairsSameContig"]++;
      if ( PairLogical )    libStats[ lib ]["nPairsLogical"]++;
      if ( PairInnie )      libStats[ lib ]["nPairsInnies"]++;
      if ( PairOutie )      libStats[ lib ]["nPairsOuties"]++;


      if ( abs(sep) <= MAX_SEPARATION ){
	if ( PairInnie )
	  pairSepDataInnies[lib].push_back( sep );
	else if ( PairOutie )
	  pairSepDataOuties[lib].push_back( sep );
	
	float rlen = 
	  (float)(reads[ pairs.ID1(pi) ].size() + reads[ pairs.ID2(pi) ].size())/2.0;
	readLengthData[lib].push_back( rlen );
      }

    }
  }

  cout << endl;

  // compute library stats
  cout << Date() << " computing stats from data" << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    if ( VERBOSE )
      PRINT2( lib, libID2libName[lib] );
    if (  libStats[ lib ]["nPairsInnies"] < MIN_SAMPLE_SIZE ){
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": unable to compute stats; sample too small (" << 
	pairSepDataInnies[lib].isize() << " < " << MIN_SAMPLE_SIZE << 
	"), preserving original sep = " << libStats[lib].at("origSep") << ", dev = " <<
	libStats[lib].at("origDev") << endl;
      
      // return original values if computation not possible
      libStats[ lib ]["newInSep"]  = libStats[ lib ]["origSep"]; 
      libStats[ lib ]["newInDev"]  = libStats[ lib ]["origDev"];
      libStats[ lib ]["newOutSep"] = 0; 
      libStats[ lib ]["newOutDev"] = 0;
      continue;
    }
    
    
    libStats[ lib ]["nSampleInnies"] = pairSepDataInnies[lib].size();
    libStats[ lib ]["nSampleOuties"] = pairSepDataOuties[lib].size();
    libStats[ lib ]["readLen"] = Mean( readLengthData[lib] );
    
    double meanInnies        = (double) libStats[lib].at("origSep");
    double devInnies         = (double) libStats[lib].at("origDev");

    double max_weight_innies = gm(pairSepDataInnies[lib], meanInnies, devInnies, VERBOSE );

    double sigma3 = FractionInDist(pairSepDataInnies[lib], meanInnies, devInnies, 3);
    libStats[ lib ]["3sigma%"]   = round(sigma3*100.0);

    // sanity check samples against mu,sigma
    if ( sigma3 < .80)
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": only " << libStats[lib]["3sigma%"] << "% of pairs are within 3 sigma of mean. Possible problem with data." << endl;
    
    libStats[ lib ]["newInSep"]   = round(meanInnies);
    libStats[ lib ]["newInDev"]   = round(devInnies);

    // update pairs library information
    pairs.changeLibrarySepSd( lib, round(meanInnies), round(devInnies) );
    
    if ( pairSepDataOuties[lib].size() >= (unsigned) MIN_SAMPLE_SIZE ) {
      Sort( pairSepDataOuties[lib] );
      double meanOuties        = (double) pairSepDataOuties[lib].front() - 1.0;
      double devOuties         = 1.0;
      double max_weight_outies = gm(pairSepDataOuties[lib], meanOuties, devOuties, VERBOSE );
	
      double sigma3 = FractionInDist(pairSepDataOuties[lib], meanOuties, devOuties, 3);
      // sanity check samples against mu,sigma
      if ( sigma3 < .80)
	cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	  ": only " << round(sigma3*100) << "% of outie pairs are within 3 sigma of mean. Possible problem with data." << endl;

      libStats[ lib ]["newOutSep"]    = round(meanOuties);
      libStats[ lib ]["newOutDev"]    = round(devOuties);
    } else {
      libStats[ lib ]["outSkewness"]   = 0;
      libStats[ lib ]["outKurtosis"]   = 0;
    }

    
    if ( abs( libStats[lib].at("newInSep") - libStats[lib].at("origSep") ) > 
	 2.0 * libStats[lib].at("origDev") )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": computed separation (" << libStats[lib].at("newInSep") << " +/- " << 
	libStats[lib].at("newInDev") << ") significantly different from the original one (" << libStats[lib].at("origSep") << " +/- 2 * " << 
	libStats[lib].at("origDev") << ")" << endl;

  }

  uint64_t refLength = 0;
  vec<int64_t> srefLens( ref.size(), 0 );
  for (size_t i = 0; i < ref.size(); i++){
    srefLens[i] = ref[i].size();
    refLength += srefLens[i];
  }
  ReverseSort( srefLens );
  
  vec<double> refInWeights( nLibraries, 0.0 );
  vec<double> refOutWeights( nLibraries, 0.0 );
  for ( size_t lib = 0; lib < libStats.size(); lib++ ){
    for ( size_t sri = 0; sri < srefLens.size(); sri++ ){
      if ( srefLens[sri] >= libStats[lib]["newInSep"] + 2 * libStats[lib]["readLen"])
	refInWeights[lib] += 
	  srefLens[sri]  - libStats[lib]["newInSep"] - 2 * libStats[lib]["readLen"];
      if ( srefLens[sri] >= libStats[lib]["newOutSep"] + 2 * libStats[lib]["readLen"])
	refOutWeights[lib] += 
	  srefLens[sri] - libStats[lib]["newOutSep"] - 2 * libStats[lib]["readLen"];
    }
    refInWeights[lib] /= (double)refLength;
    refOutWeights[lib] /= (double)refLength;
  }
  if ( VERBOSE ){
    cout << "reference weights for innies:" << endl;
    refInWeights.Print(cout); cout << endl;
    cout << "reference weights for outies:" << endl;
    refOutWeights.Print(cout); cout << endl;
  }

  for ( size_t lib = 0; lib < libStats.size(); lib++ ){
    libStats[ lib ]["percInnies"]   = 0;
    libStats[ lib ]["percOuties"]   = 0;
    if ( refInWeights[lib] > 0.0 && refOutWeights[lib] > 0.0 ){
      double sumBoth                       = 
	(double)libStats[ lib ].at("nPairsInnies")/refInWeights[lib] + 
	(double)libStats[ lib ].at("nPairsOuties")/refOutWeights[lib];
      if ( sumBoth > 0.0 ){
	libStats[ lib ]["percInnies"]   =  
	  round( (double)libStats[ lib ].at("nPairsInnies")/refInWeights[lib]/sumBoth * 100.0 );
	libStats[ lib ]["percOuties"]   =  
	  round( (double)libStats[ lib ].at("nPairsOuties")/refOutWeights[lib]/sumBoth * 100.0 );
      }
    }

    if (libStats[ lib ]["percOuties"] < 1 ){
      // very small sample size, likely there are no real outies present
      libStats[ lib ]["newOutSep"] = libStats[ lib ]["newOutDev"] = 0;
    }

    longlong nPairsInniesWghtd = refInWeights[lib] > 0.0 ? libStats[ lib ].at("nPairsInnies")/refInWeights[lib] : libStats[ lib ].at("nPairsInnies");
    longlong nPairsOutiesWghtd = refOutWeights[lib] > 0.0 ? libStats[ lib ].at("nPairsOuties")/refOutWeights[lib] : libStats[ lib ].at("nPairsOuties");
    if ( nPairsInniesWghtd < nPairsOutiesWghtd )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": found more outies (" << nPairsOutiesWghtd << ") than innies (" 
	   << nPairsInniesWghtd << "); possible problem with reads  orientation." 
	   << endl;
    
  }


  if ( ! WRITE_NOTHING ){
    // Write .outies output file
    ofstream dotout( outies_file.c_str() );
    for ( size_t lib = 0; lib < nLibraries; lib++ )
      dotout << libStats[ lib ].at("newOutSep") << " " << libStats[lib].at("newOutDev") << " " << libStats[ lib ].at("percOuties") << endl;
    
    
    // WRITING NEW PAIRS FILE
    if ( ! OUT_HEAD.empty() ){
      cout << Date() << ": writing pairs" << endl;
      pairs.Write( run_dir + "/" + OUT_HEAD + ".pairs" );
    }
  }
  
  if (!KEEP_LOOKUP_TABLE) {
    cout << "File: " << ToString (IsRegularFile(out_lookup_file)) << endl;
    Remove(out_lookup_file);
    cout << "File: " << ToString (IsRegularFile(out_lookup_file)) << endl;
  }
  
  cout << Date() << ": preparing stats output" << endl;
  // compute the character length of longest library name for formatting purposes
  string longestLibName = "libraryName";
  size_t longestLibNameSize = longestLibName.size();
  for ( size_t lib = 0; lib < nLibraries; lib++ )
    if ( libID2libName[lib].size() > longestLibNameSize )
      longestLibNameSize = libID2libName[lib].size();

  cout << Date() << ": preparing short stats output" << endl;
  /// ----------------- BUILDING SHORT STATS OUTPUT ----------------------------------
  stringstream sout;
  vec<String> short_titles; 
  short_titles.push_back("libNo", "origSep", "origDev");
  short_titles.push_back("newInSep", "newInDev");
  short_titles.push_back("3sigma%");
  short_titles.push_back("newOutSep", "newOutDev");
  short_titles.push_back("percInnies", "percOuties");
  sout << "                              ------ SHORT STATS ------\n";
  String libNameTag = "libraryName";
  sout.width( longestLibNameSize ); sout.fill(' '); 
  sout << libNameTag; sout << " ";
  short_titles.Print( sout ); sout << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    sout.width( libNameTag.size() ); sout.fill(' '); sout << libID2libName[ lib ]; sout << " ";
    for (unsigned it = 0; it < short_titles.size(); it++ ){
      String value = 
      //libStats[lib].at( short_titles.at(it) ) != -1 ? 
      //ToString( libStats[lib].at( short_titles.at(it) ) ) : "NA";
	ToString( libStats[lib].at( short_titles.at(it) ) );
      sout.width( short_titles[it].size() ); sout.fill(' '); 
      sout << value << " ";
    }
    sout << endl;
    
  }
  sout << "--------------------------------------------------------------\n\n\n";
  
  string shold;
  if (!VERBOSE) {
    shold = sout.str();
    cout << shold;
  }

  cout << Date() << ": preparing long stats output" << endl;
  /// ----------------- BUILDING LONG STATS OUTPUT----------------------------------  
  sout << "                         ------ MORE DETAILED STATS ------\n";
  sout.width( longestLibNameSize );  sout.fill(' '); sout << libNameTag; sout << " ";
  titles.Print( sout ); sout << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    sout.width( libNameTag.size() ); sout.fill(' '); sout << libID2libName[ lib ];  sout << " ";
    for (unsigned it = 0; it < titles.size(); it++ ){
      String value = 
	//libStats[lib].at( titles[it] ) != -1 ? ToString( libStats[lib][ titles[it] ] ) : "NA";
	ToString( libStats[lib][ titles[it] ] );
      sout.width( titles[it].size() ); sout.fill(' '); 
      sout << value << " ";
    }
    sout << endl;
  }
  sout << "--------------------------------------------------------------\n\n";

  // PRINTING AND WRITING STATS
  cout << Date() << ": writing stats" << endl;
  shold = sout.str();
  if (VERBOSE)
    cout << shold;
  if ( ! WRITE_NOTHING && WRITE_STATS ){ 
    ofstream rout( results_file.c_str() );
    rout << shold;
  }
  
  if (! WRITE_NOTHING && WRITE_SEPS) {
    cout << Date() << ": writing seps" << endl;
    ofstream in_seps((seps_file+".in").c_str());
    ofstream out_seps((seps_file+".out").c_str());
    for ( size_t lib = 0; lib < nLibraries; lib++ ) {
      for (u_int i = 0; i < pairSepDataInnies[lib].size(); ++i)
	in_seps << libStats[lib]["libNo"] << " " << pairSepDataInnies[lib][i] << endl;
      for (u_int i = 0; i < pairSepDataOuties[lib].size(); ++i)
	out_seps << libStats[lib]["libNo"] << " " << pairSepDataOuties[lib][i] << endl;
    }
    in_seps.close();
    out_seps.close();
  }
      
  vec< vec<String> > perfTable;
  PerfStat::log( ) << PerfStatBlockStart( "Paired Read Separation Stats" );
  perfTable.push_back( 
		      MkVec( String("Lib"),String("OrigSep"),String("NewSep"),String("NewDev"),String("3sigma%"),String("%NonJumps"),String("%ReadsAlgnd") )   );
  
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    int PercReadsAlgnd = libStats[lib]["nPairsSampled"] > 0 ? round( (libStats[lib]["nRead1Algnd"] + libStats[lib]["nRead2Algnd"]) / 2.0 / libStats[lib]["nPairsSampled"] * 100.0 ) : 0;
    String sNewSep = libStats[ lib ]["nPairsInnies"] < MIN_SAMPLE_SIZE ? "NA" : ToString(libStats[lib]["newInSep"]);
    String sNewDev = libStats[ lib ]["nPairsInnies"] < MIN_SAMPLE_SIZE ? "NA" : ToString(libStats[lib]["newInDev"]);
    perfTable.push_back(
			MkVec( ToString(libID2libName[ lib ]), ToString(libStats[lib]["origSep"]), sNewSep, sNewDev, ToString(libStats[lib]["3sigma%"]), ToString(libStats[lib]["percOuties"]), ToString(PercReadsAlgnd) )
			  );
  }
  PrintTabular( cout, perfTable, 2, "rrrrrr");
  PerfStat::log( ) << PerfStatBlockStop( );
  cout << Date() << ": Done!" << endl;

} // main() ends here.







//-------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, const VecQualNibbleVec& quals,
			  FirstLookupFinderECJ& lfinder, FirstLookupFilterECJ& lfilter,
			  const longlong pi, const int LookupK,
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, int TRIM, int64_t& sep, Bool verbose ){


  ulonglong id1 = pairs.ID1(pi);
  bvec read1 = reads[ id1 ];
  QualNibbleVec qual1 = quals[id1];
  if ( FLIP ) {
    read1.ReverseComplement();
    qual1.ReverseMe();
  }
  if ( TRIM > 0 ) read1.SetToSubOf( read1, TRIM, -1 );
  else if ( TRIM < 0 ) read1.SetToSubOf( read1, 0, read1.size() + TRIM );
  

  //lfilter.orientation = FirstLookupFilterECJ::FW_ONLY;
  list<first_look_align> hits1;
  lfinder.getAlignments( read1, qual1, id1, &hits1);
  
  ulonglong id2 = pairs.ID2(pi);
  bvec    read2 = reads[ id2 ];
  QualNibbleVec qual2 = quals[id2];
  if ( FLIP ) {
    read2.ReverseComplement();
    qual2.ReverseMe();
  }
  if ( TRIM > 0 )read2.SetToSubOf( read2, TRIM, -1 );
   else if ( TRIM < 0 ) read2.SetToSubOf( read2, 0, read2.size() + TRIM );
  
  //lfilter.orientation = FirstLookupFilterECJ::RC_ONLY;
  list<first_look_align> hits2;
  lfinder.getAlignments( read2, qual2, id2, &hits2);
 
  
  if ( hits1.size() > 0 ) Read1Algnd = True;
  if ( hits2.size() > 0 ) Read2Algnd = True;

  if ( hits1.size() > 0 && hits2.size() > 0 )
    PairAlgnd = True;
  
  if ( verbose ){
    if ( hits1.size() > 1 ){
      cout << "multiple hits for read1:" << endl;
      for ( list<first_look_align>::iterator it = hits1.begin(); it != hits1.end(); it++ ){
	(*it).PrintReadableBrief( cout, read1.size(), LookupK );
      }
      cout << endl;
    }
    if ( hits2.size() > 1 ){
      cout << "multiple hits for read2:" << endl;
      for ( list<first_look_align>::iterator it = hits2.begin(); it != hits2.end(); it++ ){
	(*it).PrintReadableBrief( cout, read2.size(), LookupK );
      }
      cout << endl;
    }
  }

  if ( hits1.size() != 1 || hits2.size() != 1 )
    return;  // take only unique alignments
  
  PairUnique = True;
 
  first_look_align &fla1 = *hits1.begin();
  first_look_align &fla2 = *hits2.begin();
  
  if ( fla1.target_loc.getContig() != fla2.target_loc.getContig() ) 
    return;  // not the same contig
  
  PairSameContig = True;
  
  // orientations
  Bool fd1 = fla1.is_FW();
  Bool fd2 = fla2.is_FW();
  
  if ( (fd1 && fd2) || (!fd1 && !fd2) ) // require correct (logical) orientation
    return;
  PairLogical = True;
  
  // if we flipped for alinment, reverse sense of forwardness
  if (FLIP) {
    fd1 = !fd1;
    fd2 = !fd2;
  }

  // use convention [s1,e1] refer to forward read, [s2,e2] to reverse
  int64_t s1,e1,s2,e2;
  if (fd1) {
    s1 = fla1.get_start_on_target(read1.isize(),LookupK);
    e1 = s1 + read1.isize();
    e2 = fla2.get_start_on_target(read2.isize(),LookupK);
    s2 = e2 + read2.isize();
  } else {
    e2 = fla1.get_start_on_target(read1.isize(),LookupK);
    s2 = e2 + read1.isize();
    s1 = fla2.get_start_on_target(read2.isize(),LookupK);
    e1 = s1 + read2.isize();
  }

  // positive->innie, negative->outie
  int64_t dist = s2 - s1;

  // fragment length
  int64_t length;

  // This is a bit of a hack which assumes we will be handling jump
  // libraries with FLIP=True, but not regular frag libraries.  Jump
  // library "innies" appear as outies here, because they are flipped
  // earlier in the pipeline.  However, those fragments smaller than the
  // read length cause confusion. --bruce
  if (FLIP) {
    int64_t readlengths = read1.isize() + read2.isize();
    dist -= readlengths;
    length = (dist > 0) ? (dist+readlengths) : (e1 - e2);
  } else {
    length = (dist > 0) ? dist : (e1 - e2);
  }

#ifdef notdef
  // show bases of reads beyond end of fragment
  if (length < read1.isize()) {
    bvec b;
    b.SetToSubOf(read1,length,read1.size()-length);
    b.Print(cout);
    b.SetToSubOf(read2,length,read2.size()-length);
    b.Print(cout);
  }
#endif

  // we define separation as fragment length minus read portion,
  // regardless of orientation (per David). --bruce
  sep = length - (read1.isize() + read2.isize());

  if ( dist >= 0 )
    PairInnie = True;
  else
    PairOutie = True;
}




// ----------------------
double gm( vec<float> data, double& emean, double& esigma, const Bool verbose ){
  Sort( data );
  size_t N = data.size();
  int i1 = round( N * 1.0/100.0 ), i2 = round( N * ( 1.0 - 1.0/100.0 ) );
  data.SetToSubOf( data, i1, i2 - i1 + 1 );

  double sigma = 1.0;
  ForceAssertGt( sigma, 0 );
  size_t K = 3;
  vec<double> means(K,0.0), sigmas(K, sigma), weights(K, 1.0/K );
  
  means[0] = data.front(); means[K -1 ] = data.back();
  
  for ( size_t i = 1; i < K -1; i++ ){
    int j = i * (N - 1.0) / (K -1.0);
    means[i] = data[ j ];
  }

  if (verbose){
    cout << "meansOrig:   "; means.Print(cout); cout << endl;
    cout << "sigmasOrig:  "; sigmas.Print(cout); cout << endl;
    cout << "weightsOrig: "; weights.Print(cout); cout << endl;
  }
  
  Bool RemovedGaussian = True;
  while( K > 1 && RemovedGaussian ){
    RemovedGaussian = False;
    actual_gm( data, means, sigmas, weights, verbose );    
    vec<int> indices( K );
    for ( size_t i = 0; i < K; i++ )
      indices[i] = i;
    vec<double> weights_sorted( weights );
    ReverseSortSync( weights_sorted, indices );
    for ( size_t i = 0; i < indices.size() -1; i++ ){
      int c = indices[i];
      int n = indices[i+1];
      if ( abs( means[c] - means[n] ) <= 3.0 * sigmas[c] || 
	   abs( means[c] - means[n] ) < 3.0 * sigmas[n] ){
	K--;
	RemovedGaussian = True;
	means.Erase(n,n+1); sigmas.Erase(n,n+1); weights.Erase(n,n+1); indices.Erase(i,i+1);
	break;
      }
    }

  }
  
  vec<int> indices( K );
  for ( size_t i = 0; i < K; i++ )
    indices[i] = i;
  vec<double> weights_sorted( weights );
  ReverseSortSync( weights_sorted, indices );
  double mean0 = means[indices[0]], sigma0 = sigmas[indices[0]];
  double maxweight = weights[ indices[0] ];
  if (maxweight > 0) {
    emean  = mean0;
    esigma = sigma0;
  }
  if ( verbose ){
    PRINT4( K, emean, esigma, maxweight );
    cout << endl;
  }

  return maxweight;
}

void actual_gm( const vec<float>& data, vec<double>& means, vec<double>& sigmas, vec<double>& weights, const Bool verbose ){
  size_t N = data.size();
  size_t K = means.size();
  double ratio = 9999;
  size_t nsteps = 0;
  double oldloglike = 1.0, loglike = 1.0;
  while( nsteps < 2  || abs(1.0 - ratio) > 0.001 && nsteps < 100){
    nsteps++;
    oldloglike = loglike; 
    vec< vec<double> > ddat( N, vec<double> (K, 0.0 ) );
    for ( size_t k = 0; k < K; k++ )
      for ( size_t n = 0; n < N; n++ )
	if ( weights[k] > 0.0 )
	  ddat[n][k] = -0.5 * pow( (data[n] - means[k]) / sigmas[k], 2.0 ) +
	    log( weights[k] )  + log( 1.0 / sqrt(2.0 * M_PI ) ) 
	    - log( sigmas[k] );
    
    
    loglike = 0.0; 
    for ( size_t n = 0; n < N; n++ ){
      double max = ddat[n][0];
      for ( size_t k = 0; k < K; k++ )
	if ( weights[k] > 0.0 )
	  if ( ddat[n][k] > max ) max = ddat[n][k];

      double sum = 0.0;
      for ( size_t k = 0; k < K; k++ )
	if ( weights[k] > 0.0 )
	  sum += exp( ddat[n][k] - max );
      
      double tmp = max + log( sum );
      for ( size_t k = 0; k < K; k++ )
	if ( weights[k] > 0.0 )
	  ddat[n][k] = exp( ddat[n][k] - tmp );
      
      loglike += tmp;
    }
    
    for ( size_t k = 0; k < K; k++ ){
      if ( weights[k] > 0.0 ){
	double wgt = 0.0;
	for ( size_t n = 0; n < N; n++ )
	  wgt += ddat[n][k];
	
	weights[k] = wgt / N; 
	
	double sum = 0.0;
	for ( size_t n = 0; n < N; n++ )
	  sum += ddat[n][k] * data[n];
	
	means[k] = sum / wgt;
	double sum2 = 0.0;
	for ( size_t n = 0; n < N; n++ )
	  sum2 += ddat[n][k] * pow( data[n] - means[k], 2.0 );
	// since we are dealing with integer seps, it's not improbable to get a peak
	// with only one value (sigma = 0)...this causes all sorts of trouble, so
	// putting floor of 1.0 on sigma. --bruce
	sigmas[k] = max(sqrt( sum2 / wgt ), 1.0);
      }else{
	means[k] = data.front() - 100000;
	sigmas[k] = 1e-99;
      }
    }
    
    ratio = abs( loglike / oldloglike );
    if (verbose){
      PRINT2(ratio, nsteps);
      cout << "means:   "; means.Print(cout); cout << endl;
      cout << "sigmas:  "; sigmas.Print(cout); cout << endl;
      cout << "weights: "; weights.Print(cout); cout << endl;
    }
  }
  if (verbose) cout << "\n\n";
}
