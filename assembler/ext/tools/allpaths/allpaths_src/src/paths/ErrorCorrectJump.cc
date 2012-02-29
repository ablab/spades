///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "feudal/IncrementalWriter.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "STLExtensions.h"
#include "lookup/FirstLookupFinderECJ.h"
#include "lookup/FirstLookupPerfect.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTab.h"
#include "lookup/LookupTabBuilder.h"
#include "math/Functions.h"
#include "paths/GetNexts.h"
#include "paths/KmerPath.h"
#include "reporting/PerfStat.h"
#include "util/ReadTracker.h"
#include "util/RunCommand.h"

static inline 
String Tag(String S = "ECJ") { return Date() + " (" + S + "): "; } 

void WriteSummaryStats( ostream & out, const vec<PM_LibraryStats> & stats_in,
			const vec<PM_LibraryStats> & stats_out );


void SaveAlignments( const String & fn, const BaseVecVec & queries, 
		     const BaseVecVec & ref, const size_t K_lookup_tab,
		     const vec<first_look_align> & first_aligns,
		     const bool FLIP, const bool PERFECT );

size_t FindAlignments( const size_t K, const String & queries_quals_fn, 
		       const String & ref_tab_fn, const BaseVecVec & queries, 
		       const BaseVecVec & ref,
		       const FirstLookupFilterECJ & lookup_filter,
		       vec<first_look_align> * first_aligns, 
		       size_t * K_lookup_tab, const bool PERFECT,
		       const bool FLIP, const size_t NUM_THREADS);

/**
 * ErrorCorrectJump
 *
 * Do error correction of jumping reads by aligning them onto a set of
 * unibases (plus adjacency graph), and replacing the bases of the
 * reads with those of the underlying unibase(s). If FLIP is true, run
 * error correction by seeding (and extending) from the last K-mer in
 * the read, rather than from the first.
 *
 * Input files:
 *   <REF_HEAD>.unibases.k<K>
 *   <REF_HEAD>.unipaths.k<K>
 *   <REF_HEAD>.unipath_adjgraph.k<K>
 *   <QUERY_HEAD>.fastb
 *   <QUERY_HEAD>.qualb
 *   <QUERY_HEAD>.pairs
 *
 * Output files (if SAVE_ALIGNS_AND_QUIT=False):
 *   <OUT_HEAD>.orig_id
 *   <OUT_HEAD>.fastb
 *   <OUT_HEAD>.paths.k<K>
 *   <OUT_HEAD>.pairs*
 *   <OUT_HEAD>.qltout       (if SAVE_ALIGNS=True)
 *
 * Output fles (if SAVE_ALIGNS_AND_QUIT=True):
 *   <OUT_HEAD>.qltout
 *
 */
int main( int argc, char **argv)
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_UnsignedInt(K);
  CommandArgument_String_Doc(REF_HEAD,
    "Base name for unibases and unipaths" );
  CommandArgument_String_Doc(QUERY_HEAD,
    "Base name for query reads" );
  CommandArgument_String_Doc(OUT_HEAD,
    "Base name for output reads" );
  CommandArgument_Int_OrDefault_Doc(TARGET_READ_LEN, -1,
    "Ensure all corrected reads are at least this length (default is K bases)");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Int_OrDefault(LOOKUP_K, 12);
  CommandArgument_Bool_OrDefault_Doc(FLIP, True,
    "If true, rc reads, error correct, and then rc back" );
  CommandArgument_Bool_OrDefault_Doc(PERFECT, False,
    "Require perfect alignments (useful for very short reads)");
  CommandArgument_Bool_OrDefault_Doc(SAVE_ALIGNS, False,
    "Save alignments (as look_align objects)" );
  CommandArgument_Bool_OrDefault_Doc(SAVE_ALIGNS_AND_QUIT, False,
    "Same as SAVE_ALIGNS, but quit after dumping aligns" );
  CommandArgument_Int_OrDefault(MAX_KMER_FREQ, 10000);
  CommandArgument_Int_OrDefault(MIN_MATCH, 20);
  CommandArgument_Int_OrDefault(MISMATCH_THRESHOLD, 3);
  CommandArgument_Int_OrDefault(MISMATCH_NEIGHBORHOOD, 8);
  CommandArgument_Int_OrDefault(MISMATCH_BACKOFF, 3);
  CommandArgument_Int_OrDefault(SCORE_DELTA, 20);
  CommandArgument_Int_OrDefault(SCORE_MAX, 100);
  CommandArgument_Bool_OrDefault(TRACK_READS, True);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_Bool_OrDefault(WRITE_EXTRA, False);
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_String_OrDefault(DEBUG_READS, "");
  EndCommandArguments;

  // Thread control
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Check arguments.
  if (FLIP && PERFECT) {
    cout << "WARNING! WARNING! WARNING! Since you asked for PERFECT=True, "
	 << "I am resetting FLIP to False.\n" << endl;
    FLIP = False;
  }
  
  // File names.
  String strK = ToString(K);
  String lK = (LOOKUP_K == 12) ? "" : ToString(LOOKUP_K);

  String ref_bases_fn = REF_HEAD + ".unibases.k" + strK;
  String ref_tab_fn = REF_HEAD + ".unibases.k" + strK + ".lookup" + lK;
  String ref_unipaths_fn = REF_HEAD + ".unipaths.k" + strK;

  String queries_bases_fn = QUERY_HEAD + ".fastb";
  String queries_quals_fn = QUERY_HEAD + ".qualb";
  String queries_pairs_fn = QUERY_HEAD + ".pairs";

  String out_source_fn = OUT_HEAD + ".orig_id";
  String out_bases_fn = OUT_HEAD + ".fastb";
  String out_bases_extra_fn = OUT_HEAD + "_extra.fastb";
  String out_pairs_extra_fn = OUT_HEAD + "_extra.pairs";
  String out_kpaths_fn = OUT_HEAD + ".paths.k" + strK;
  String out_qlt_fn = OUT_HEAD + ".qltout";

  // Generate a LookupTab from the unibases and write it to file.
  // (Skip this if the file already exists.)
  if (!IsRegularFile(ref_tab_fn)) {
    cout << Tag() << "Creating LookupTab file" << endl;
    LookupTabBuilder lookup_tab_builder(LOOKUP_K);
    lookup_tab_builder.addFastb(ref_bases_fn.c_str());
    lookup_tab_builder.write(ref_tab_fn.c_str());
  }

  // Load.
  cout << Tag() << "loading reference bases" << endl;
  BaseVecVec ref(ref_bases_fn);

  // Set up filtering options for call to FirstLookupFinedECJ.
  FirstLookupFilterECJ lookup_filter;
  lookup_filter.orientation = FLIP ?
    FirstLookupFilterECJ::RC_ONLY :
    FirstLookupFilterECJ::FW_ONLY;
  lookup_filter.min_size = 20;
  lookup_filter.max_kmer_freq = MAX_KMER_FREQ;
  lookup_filter.max_extend = MAX_KMER_FREQ;
  lookup_filter.score_delta = SCORE_DELTA;
  lookup_filter.score_max = SCORE_MAX;
  lookup_filter.min_match = MIN_MATCH;
  lookup_filter.mismatch_threshhold = MISMATCH_THRESHOLD;
  lookup_filter.mismatch_neighborhood = MISMATCH_NEIGHBORHOOD;
  lookup_filter.mismatch_backoff = MISMATCH_BACKOFF;
  lookup_filter.max_placements = 2;   // if there are multiple placements
  ParseIntSet(DEBUG_READS, lookup_filter.debug_reads, False);

  // Load queries bases/quals for alignment.
  cout << Tag() << "loading query bases" << endl;
  BaseVecVec queries(queries_bases_fn);
  size_t n_queries = queries.size();
  
  // Flip the bases and quals if necessary.
  if (FLIP) {
    cout << Tag() << "rc-ing query bases" << endl;
    for (size_t ii = 0; ii < n_queries; ii++)
      queries[ii].ReverseComplement();
  }    

  // Run FirstLookup.
  vec<first_look_align> first_aligns;
  size_t K_lookup_tab;
  const size_t n_aligns = FindAlignments(K, queries_quals_fn, ref_tab_fn, 
                                   queries, ref, lookup_filter,
                                   & first_aligns,
                                   & K_lookup_tab, 
                                   PERFECT, FLIP, NUM_THREADS);
  ForceAssertEq( (int64_t)LOOKUP_K, (int64_t)K_lookup_tab );

  // Save aligns.
  if ( SAVE_ALIGNS || SAVE_ALIGNS_AND_QUIT ) {
    ForceAssertLe(n_queries, (1ul << 30));   // too many queries!
    SaveAlignments(out_qlt_fn, queries, ref, K_lookup_tab,
                   first_aligns, FLIP, PERFECT);
    if ( SAVE_ALIGNS_AND_QUIT ) {
      cout << Tag( ) << "Done with ErrorCorrectJump" << endl;
      return 0;
    }
  }
  
  // Generate placement map.
  vec< vec<size_t> > multi_align_index(n_queries);
  for (size_t ia = 0; ia < n_aligns; ia++)
     multi_align_index[ first_aligns[ia].query_ID ].push_back( ia );

  size_t n_uniquely = 0, n_non_uniquely = 0, n_unplaced = 0;
  for ( size_t qid = 0; qid < n_queries; qid++ ) {
    if ( multi_align_index[qid].size() == 1 ) n_uniquely++;
    else if ( multi_align_index[qid].size() > 1 ) n_non_uniquely++;
    else n_unplaced++;
  }
  
  if (true) {
    cout.width(10); 
    cout << n_aligns << " alignments"  << endl << endl;
    cout.width(10); 
    cout << n_queries << " reads" << endl;
    cout.width(10); 
    cout << n_uniquely << " reads aligned uniquely" << endl;
    cout.width(10); 
    cout << n_non_uniquely << " reads aligned non_uniquely" << endl;
    cout.width(10); 
    cout << n_unplaced << " reads unaligned" << endl;
  }
  
  // Output.
  IncrementalWriter<BaseVec> b_out(out_bases_fn.c_str());
  IncrementalWriter<BaseVec> *b_extra_out = 0;
  if (WRITE_EXTRA)
    b_extra_out = new IncrementalWriter<BaseVec>(out_bases_extra_fn.c_str());
  PairsManager pairs_extra_out;
  IncrementalWriter<KmerPath> k_out(out_kpaths_fn.c_str());
  PairsManager pairs_out;
  ofstream orig_out(out_source_fn.c_str());
  ReadTracker rt;
  size_t rt_source = 0;
  size_t n_good_pairs = 0, n_good_extra_pairs = 0;
  size_t n_shorts = 0, n_multiples = 0;

  if (TRACK_READS) {
    rt_source = rt.AddSource(QUERY_HEAD);
  }
  
  cout << Tag() << "loading reference unipaths" << endl;
  vecKmerPath unipaths(ref_unipaths_fn);

  // Get next unibases.
  vec< vec<int> > nexts;
  GetNexts(K, ref, nexts);

  // Load PairsManager.
  cout << Tag() << "loading pairs" << endl;
  PairsManager pairs(queries_pairs_fn);
  uint64_t n_pairs = pairs.nPairs();
  ForceAssertEq(n_queries, pairs.nReads());

  // Parse all pairs.
  size_t dotter = 1e6;

  cout << Tag() << "parsing "
       << n_pairs << " pairs (.="
       << dotter << " pairs):"
       << endl;

  // To fix the crash when inconsistent sep are found for variable
  // length reads save the new sep information for the library
  int nLib= pairs.nLibraries();
  vec<int> newSep(nLib, 0);
  vec<bool> hasNewSep(nLib, False);
  uint64_t max_uint64 = numeric_limits<uint64_t>::max();
  for (size_t pair_ID = 0; pair_ID < n_pairs; pair_ID++) {
    DotMod(cout, pair_ID, dotter);
    
    size_t q1_ID = pairs.ID1(pair_ID);
    size_t q2_ID = pairs.ID2(pair_ID);
    
    // only accept pairs with both reads aligned
    if ( multi_align_index[q1_ID].size() >= 1 &&
	 multi_align_index[q2_ID].size() >= 1 ) { 

      // check if alternative alignment regions provide consistent correction
      bool IsConsistent1 = true;     
      const uint64_t q1_len = queries[q1_ID].size();
      basevector b10;
      uint64_t s1 = -1, t1_ID = queries.size();
      for ( int alt1 = multi_align_index[q1_ID].isize() -1; alt1 >= 0; alt1-- ){
	size_t ia1 = multi_align_index[q1_ID][alt1];
	first_look_align & hit1 = first_aligns[ia1];
	ForceAssertEq(hit1.query_ID, q1_ID);
	s1 = hit1.get_start_on_target(q1_len, K_lookup_tab);
	t1_ID = hit1.get_target_ID();
	size_t aliLen = min( q1_len, ref[t1_ID].size() - s1 );
	basevector b1(ref[t1_ID], s1, aliLen);
	if ( alt1 == multi_align_index[q1_ID].isize() -1 )
	  b10 = b1;
	else if ( b1 != b10 ){
	  IsConsistent1 = false;
	  break;
	}
      }
      ForceAssert( s1 < max_uint64 && t1_ID < max_uint64 );
      if ( ! IsConsistent1 ) continue; // correction not consistent, skip pair
      
      bool IsConsistent2 = true;     
      const uint64_t q2_len = queries[q2_ID].size();
      basevector b20;
      uint64_t s2 = -1, t2_ID = -1;
      for ( int alt2 = multi_align_index[q2_ID].isize() -1; alt2 >= 0; alt2-- ){
	size_t ia2 = multi_align_index[q2_ID][alt2];
	first_look_align & hit2 = first_aligns[ia2];
	ForceAssertEq(hit2.query_ID, q2_ID);
	s2 = hit2.get_start_on_target(q2_len, K_lookup_tab);
	t2_ID = hit2.get_target_ID();
	size_t aliLen = min( q2_len, ref[t2_ID].size() - s2 );
	basevector b2(ref[t2_ID], s2, aliLen);
	if ( alt2 == multi_align_index[q2_ID].isize() -1 )
	  b20 = b2;
	else if ( b2 != b20 ){
	  IsConsistent2 = false;
	  break;
	}
      }
      ForceAssert( s2 < max_uint64 && t2_ID < max_uint64 );
      if ( ! IsConsistent2 ) continue;   // correction not consistent, skip pair

      // multiple alignmnents found.
      if (multi_align_index[q1_ID].size() != 1 ||
	  multi_align_index[q2_ID].size() != 1) {
	if (WRITE_EXTRA) {
	  b_extra_out->add( b10 );
	  b_extra_out->add( b20 );
	  pairs_extra_out.addPair(2*n_good_extra_pairs, 2*n_good_extra_pairs + 1,
				  pairs.sep(pair_ID), pairs.sd(pair_ID),
				  pairs.libraryName(pair_ID)+"_extra", True);
	}
	n_good_extra_pairs++;
	n_multiples++;
	continue;
      }

      // Figure out if we can extend read to be at least K bases.
      uint64_t maxlen1 = ref[t1_ID].size() - s1;      
      uint64_t n1 = t1_ID;
      while (maxlen1 < q1_len && nexts[n1].size() == 1) {
	  n1 = nexts[n1][0];
	  maxlen1 += ref[n1].size() - K + 1;
      }

      uint64_t maxlen2 = ref[t2_ID].size() - s2;
      uint64_t n2 = t2_ID;
      while (maxlen2 < q2_len && nexts[n2].size() == 1) {
	n2 = nexts[n2][0];
	maxlen2 += ref[n2].size() - K + 1;
      }
      
      if ( maxlen1 < K || maxlen2 < K ) { // reads cannot be both extended
	if (WRITE_EXTRA) {
	  b_extra_out->add( b10 );
	  b_extra_out->add( b20 );
	  pairs_extra_out.addPair(2*n_good_extra_pairs, 2*n_good_extra_pairs + 1,
				  pairs.sep(pair_ID), pairs.sd(pair_ID),
				  pairs.libraryName(pair_ID)+"_extra", True);
	}
	n_good_extra_pairs++;
	n_shorts++;
	continue;
      }
      else { // both reads can be extended
        orig_out << q1_ID << "\n" << q2_ID << "\n";
	if (TRACK_READS) {
	  rt.AddRead(rt_source, q1_ID);
	  rt.AddRead(rt_source, q2_ID);
	}
	  
        // Define the error-corrected reads.
	uint64_t N1 = (TARGET_READ_LEN > 0 ? TARGET_READ_LEN : MAX(q1_len, K));
	uint64_t N2 = (TARGET_READ_LEN > 0 ? TARGET_READ_LEN : MAX(q2_len, K));

        basevector b1(ref[t1_ID], s1, Min(N1, ref[t1_ID].size() - s1));
        basevector b2(ref[t2_ID], s2, Min(N2, ref[t2_ID].size() - s2));
        vec<longlong> kmers1, kmers2;

        for (size_t j = 0; j < b1.size()-K+1 && kmers1.size() != N1-K+1; j++)
          kmers1.push_back(unipaths[t1_ID].GetKmer(s1 + j));
	
        for (size_t j = 0; j < b2.size()-K+1 && kmers2.size() != N2-K+1; j++)
          kmers2.push_back(unipaths[t2_ID].GetKmer(s2 + j));
	
        n1 = t1_ID;
        n2 = t2_ID;
        while (b1.size() < N1 && nexts[n1].size() == 1) {
          n1 = nexts[n1][0];
          for (size_t j = K-1; j < ref[n1].size() && b1.size() != N1; j++)
            b1.AppendBase(ref[n1][j]);

          const KmerPath& U = unipaths[n1];
          for (size_t j=0; j<(size_t)U.NSegments() && kmers1.size() != N1-K+1; j++)
            for (longlong r=U.Start(j); r<=U.Stop(j) && kmers1.size() != N1-K+1; r++)
              kmers1.push_back(r);
        }
	
        while (b2.size() < N2 && nexts[n2].size() == 1) {
          n2 = nexts[n2][0];
          for (size_t j = K-1; j < ref[n2].size() && b2.size() != N2; j++)
            b2.AppendBase(ref[n2][j]);
	  
          const KmerPath& U = unipaths[n2];
          for (size_t j = 0; j < (size_t)U.NSegments() && kmers2.size() != N2 - K + 1; j++)
            for (longlong r = U.Start(j); r <= U.Stop(j) && kmers2.size() != N2 - K + 1; r++)
              kmers2.push_back(r);
	  
        }

        b_out.add(b1);
        b_out.add(b2);

        KmerPath kp1(kmers1), kp2(kmers2);
        k_out.add(kp1);
        k_out.add(kp2);
	
	// Adjust the separation between the reads to account for the
	// change we have made in the read length.  Hopefully all
	// pairs in a library will get their separations adjusted by
	// the same amount; if not, there may be a crash in addPair,
	// below.
	int sep = pairs.sep(pair_ID);
	sep += ( q1_len - N1 );
	sep += ( q2_len - N2 );

	// To fix the crash when inconsistent sep are found for
	// variable length reads save the new sep information for the
	// library
	int libID = pairs.libraryID(pair_ID);
	if ( hasNewSep[libID] ){ 
	  // if the sep is not consistent with the previous reads
	  if(sep != newSep[libID]) sep = newSep[libID];
	}else{
	  hasNewSep[libID] = True;
	  newSep[libID] = sep;
	}

        if (VERBOSE) PRINT2(pair_ID, pairs_out.nPairs());
	pairs_out.addPair(2*n_good_pairs, 2*n_good_pairs + 1,
                          sep, pairs.sd(pair_ID),
                          pairs.libraryName(pair_ID), True);
	n_good_pairs++;
      }
    }
  }
  

  orig_out.close();
  cout << endl;

  if (TRACK_READS) rt.Dump(OUT_HEAD);


  if (WRITE) {
    // Save.
    cout << Tag() << "merging output fastb and KmerPaths" << endl;
    b_out.close();
    k_out.close();

    cout << Tag() << "saving pairing info" << endl;
    pairs_out.Write(OUT_HEAD + ".pairs");
    
    if (WRITE_EXTRA) {
      b_extra_out->close();
      pairs_extra_out.Write( out_pairs_extra_fn );
    }
  }

  WriteSummaryStats(cout, pairs.getLibraryStats(), pairs_out.getLibraryStats());


  
  PerfStat::log( ) << std::fixed << std::setprecision(2) 
		   << PerfStat( "frac_jumps_ec", "% of jump reads pairs that are error corrected",
				100.0 * n_good_pairs / n_pairs );
  // Done.
  cout << Tag() << "Done with ErrorCorrectJump" << endl;
  return 0;
}



void WriteSummaryStats( ostream & out, const vec<PM_LibraryStats> & stats_in, const vec<PM_LibraryStats> & stats_out )
{

  // Make a chart header line.
  vec< vec<String> > rows;

  vec<String> row1, row2, row3;
  row1.push_back( "lib_name", "lib_ID", "sep", "dev", "orig_reads", "corr_reads", "frac_corr" );
  row2.push_back( "--------", "------", "---", "---", "----------", "----------", "---------" );
  rows.push_back( row1, row2 );
  
  uint64_t in_count = 0;
  uint64_t out_count = 0;

  // For each library, make a line of info.
  for ( size_t i = 0; i < stats_out.size(); i++ ) {
    vec<String> row;
    if (stats_out[i].name == "Unpaired" )
      row.push_back( stats_out[i].name, "-", "-", "-" );
    else
      row.push_back( stats_out[i].name,
		     ToString(i),
		     ToString( stats_out[i].sep ),
		     ToString( stats_out[i].sd ) );
    row.push_back( ToString( stats_in[i].n_reads ) );
    row.push_back( ToString( stats_out[i].n_reads ) );
    row.push_back( ToString( stats_out[i].n_reads / static_cast<double>(stats_in[i].n_reads) , 2 ) );
    rows.push_back(row);
    in_count += stats_in[i].n_reads;
    out_count += stats_out[i].n_reads;
  }

  row3.push_back( "Total", "", "", "", ToString(in_count), ToString(out_count) ,
		  ToString(out_count / static_cast<double>(in_count), 2) );
  rows.push_back( row2, row3 );
  
  // Convert the lines into a chart.
  out << endl;
  PrintTabular( out, rows, 2, "lrrrrrr" );
  out << endl;
}



/****************************/
/*   FUNCTIONS START HERE   */
/****************************/



void SaveAlignments( const String & fn, 
		     const BaseVecVec & queries, 
		     const BaseVecVec & ref, 
		     const size_t K_lookup_tab,
		     const vec<first_look_align> & first_aligns,
		     const bool FLIP, 
		     const bool PERFECT )
{
  cout << Tag() << "saving aligns" << endl;
  
  ofstream os_la(fn.c_str());
  SetWritePrettyLookAligns(os_la);

  size_t n_aligns = first_aligns.size();
  for (size_t ia = 0; ia < n_aligns; ia++) {
    
    // Convert the first_look_aligns into a look_align
    look_align align;
    first_aligns[ia].convert_to_look_align(align, queries, ref, K_lookup_tab);
    if (FLIP) align.rc1 = false;
    
    os_la << align;
  }

  os_la.close();
}

size_t FindAlignments(const size_t K, 
                      const String & queries_quals_fn, 
                      const String & ref_tab_fn, 
                      const BaseVecVec & queries, 
                      const BaseVecVec & ref, 
                      const FirstLookupFilterECJ & lookup_filter,
                      vec<first_look_align> * first_aligns,
                      size_t * K_lookup_tab,
                      const bool PERFECT,
                      const bool FLIP, 
                      const size_t NUM_THREADS)
{
  // Load LookupTab for alignment.
  cout << Tag() << "Loading LookupTab from file " << ref_tab_fn << endl;
  LookupTab lookup_tab(ref_tab_fn.c_str());
  *K_lookup_tab = lookup_tab.getK();
    
  if (PERFECT) {   // Run the alignment using FirstLookupPerfect
    cout << Tag() << "looking for perfect aligns" << endl;
    FirstLookupPerfect lperfect(lookup_tab, ref);
    lperfect.getAllAlignments(queries, first_aligns, NUM_THREADS);

  }
  else {   // Run the alignment using FirstLookupFinderECJ
    size_t n_queries = queries.size();

    cout << Tag() << "loading query quals" << endl;
    VecQualNibbleVec quals;
    LoadQualNibbleVec(queries_quals_fn, &quals);
    ForceAssertEq(n_queries, quals.size());
    if (FLIP) {
      cout << Tag() << "rc-ing query quals" << endl;
      for (size_t ii = 0; ii < n_queries; ii++)
        quals[ii].ReverseMe();
    }
    
    cout << Tag() << "running FirstLookupFinderECJ" << endl;
    FirstLookupFinderECJ lfinder(lookup_filter, lookup_tab, ref, K);
    lfinder.getAllAlignments(queries, quals, first_aligns, NUM_THREADS);
  }
  
  return first_aligns->size();
}

