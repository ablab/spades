///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FindErrors

#include "MainTools.h"
#include "reporting/PerfStat.h"

#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/QualNibbleVec.h"

#include "kmers/KmerSpectrumCore.h"

#include "paths/FindErrorsCore.h"


static inline 
String Tag(String S = "FE") { return Date() + " (" + S + "): "; } 




template<class QVV_t>
void assert_bases_quals_consistency(const BaseVecVec & bvv,
                                    const QVV_t & qvv)
{
  ForceAssertEq(bvv.size(), qvv.size());
  const size_t nv = bvv.size();
  size_t n_diff = 0;
  for (size_t iv = 0; iv < nv; iv++)
    if (bvv[iv].size() != qvv[iv].size())
      n_diff++;
  ForceAssertEq(n_diff, 0u);
}






float gc_content_compute(const BaseVecVec & bvv)
{
  const size_t nbv = bvv.size();
  size_t nb = 0;
  size_t nb_gc = 0;
  for (size_t ibv = 0; ibv < nbv; ibv++) {
    const BaseVec & bv = bvv[ibv];
    nb += bv.size();
    for (size_t ib = 0; ib < bv.size(); ib++)
      if (bv[ib] == 1 || bv[ib] == 2) 
        nb_gc++;
  }
  
  return float(nb_gc) / float(nb);
}





void heuristics_estimate(const u_int ESTIMATED_COVERAGE,
                         const unsigned ploidy, 
                         const BaseVecVec & bases,
                         const KmerSpectrum & kspec,
                         u_int & MIN_READSTACK_DEPTH,
                         u_int & MIN_BASESTACK_DEPTH,
                         int & MIN_QSS_TO_CONFIRM,
                         int & MIN_QSS_TO_SUPPORT,
                         int & MIN_QSS_TO_WIN)
{
  double kf_lo = 0;
  double kf_mode = 0;  // actually, not used currently

  if (ESTIMATED_COVERAGE) {
    kf_lo = round(ESTIMATED_COVERAGE / 7.0);
    kf_mode = round(4.0 * ESTIMATED_COVERAGE / 7.0);
    cout << Tag() << "Setting heuristic parameters based on estimated coverage." << endl;
  } 
  else {
    // ---- Analyze kmer spectrum (needs mean read length)
    const size_t n_reads = bases.size();
    size_t n_bases = 0;
    for (size_t i = 0; i < n_reads; i++) 
      n_bases += bases[i].size();
    
    kspec.analyze(ploidy, n_bases / n_reads);

    if (kspec.genome_size() == 0) {
      cout << Tag() << "Can't determine scaling parameters." << endl;
      return;
    }
    cout << Tag() << "Setting heuristic parameters based on kmer spectrum." << endl;
    kf_lo = kspec.kf_min1();
    kf_mode = kspec.kf_max2();
  }
  cout << Tag() << "(kf_min= " << kf_lo << ", kf_mode= " << kf_mode << ")" << endl;

  // Set min read stack 1/4 of the way off the minimum toward the mode (log scale);
  // reasoning is that we only want to error correct off stacks which we believe to
  // consist of valid kmers, and it's about 50/50 at the minimum.
  //double good_stack = pow((double)kf_lo * kf_lo * kf_lo * kf_mode, .25);
  // Scale base determined by looking at examples at our nominal default of 50x coverage

  double good_stack = kf_lo;
  double scale = good_stack / 7.0;

  MIN_QSS_TO_CONFIRM = round(scale * MIN_QSS_TO_CONFIRM);
  cout << Tag() << "Setting MIN_QSS_TO_CONFIRM to " << MIN_QSS_TO_CONFIRM << endl;

  MIN_QSS_TO_SUPPORT = round(scale * MIN_QSS_TO_SUPPORT);
  cout << Tag() << "Setting MIN_QSS_TO_SUPPORT to " << MIN_QSS_TO_SUPPORT << endl;

  MIN_QSS_TO_WIN = round(scale * MIN_QSS_TO_WIN);
  cout << Tag() << "Setting MIN_QSS_TO_WIN to " << MIN_QSS_TO_WIN << endl;

  MIN_READSTACK_DEPTH = round(good_stack);
  cout << Tag() << "Setting MIN_READSTACK_DEPTH to " << MIN_READSTACK_DEPTH << endl;

  // Don't go too far into the fringes off the base kmer stack; since we're scaling MIN_QSS...
  // parameters linearly, this should probalby be linear with the stack size as well.
  MIN_BASESTACK_DEPTH = round(good_stack / 4.0);
  cout << Tag() << "Setting MIN_BASESTACK_DEPTH to " << MIN_BASESTACK_DEPTH << endl;
}




// ---- temporary class to store locations of quals of 0 
//      when we change to only QualNibbleVec this will be obsolete

class BVLocFE
{
public:
  uint64_t iv : 40;
  uint64_t i  : 24;
  BVLocFE(const size_t _iv, const size_t _i) : iv(_iv), i(_i) {}
};












int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;

  // File locations:

  CommandArgument_String_Doc(HEAD_IN, "Looks for 'HEAD_IN.{fastb,qualb}'.");
  CommandArgument_String_Doc(HEAD_OUT, "Generates 'HEAD_OUT.{fastb,keep,log}'.");

  CommandArgument_String_OrDefault_Doc(PLOIDY_FILE, "", "If file exists, reads ploidy.");
  
  // Computational performance:

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0).");
  CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB, 0);



  // General control:

  CommandArgument_Bool_OrDefault_Doc(PRE_CORRECT, True, 
     "Pre correct reads before error correction.");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_CYCLES, 2, 
     "Number of correction cycles.");
  CommandArgument_Bool_OrDefault_Doc(USE_KMER_SPECTRUM, False, 
     "Use kmer spectrum to adjust heuristics.");
  CommandArgument_UnsignedInt_OrDefault_Doc(ESTIMATED_COVERAGE, 0, 
     "If positive, use this as estimated coverage for adjusting heuristics instead of kmer spectrum.");
  CommandArgument_Bool_OrDefault_Doc(WRITE, True, "write results");


  // Heuristics:

  CommandArgument_Int_OrDefault_Doc(K, 24, 
     "Kmer size for error correction.");

  CommandArgument_Int_OrDefault_Doc(MIN_Q_TO_SUPPORT, efp_default.min_Q_to_support, 
     "Minimum base quality score (Q) to mark it a 'supporting' or 'confirming' base");
  CommandArgument_Int_OrDefault_Doc(MIN_N_BASES_TO_SUPPORT, efp_default.min_n_bases_to_support, 
     "Minimum number supporting bases of a nucleotide to mark it 'supported'");
  CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CORRECT, efp_default.max_QSsum_ratio_to_correct, 
     "Maximum ratio of QSS (Quality score sum) of existing to best base to mark "
     "it for 'correction'");
  CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CORRECT2, efp_default.max_QSsum_ratio_to_correct2, 
     "Maximum ratio of QSS (Quality score sum) of second best to best base to mark "
     "any base for 'correction'");
  CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CONFIRM, efp_default.max_QSsum_ratio_to_confirm, 
     "Maximum ratio of QSS (Quality score sum) of second best to best base to "
     "confirm any base");
  CommandArgument_Int_OrDefault_Doc(MIN_QSS_TO_CONFIRM, efp_default.min_QSsum_to_confirm, 
     "Minimum QSS (Quality score sum) of bases of a nucleotide to mark it 'base_locked'; "
     "overrides MIN_N_BASES_TO_CONFIRM & MIN_Q_TO_SUPPORT");
  CommandArgument_Int_OrDefault_Doc(MIN_QSS_TO_WIN, efp_default.min_QSsum_to_win, 
     "Ignore columns for which the quality score sum for the winning base is not "
     "at least this value.");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_READSTACK_DEPTH, efp_default.min_readstack_depth, 
     "Minimum number of same-kmer containing reads to consider building a ReadStack");
  CommandArgument_UnsignedInt_OrDefault_Doc(MAX_READSTACK_DEPTH, efp_default.max_readstack_depth, 
     "Maximum number of same-kmer containing reads to consider building a ReadStack");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_BASESTACK_DEPTH, efp_default.min_basestack_depth, 
     "Minimum number of bases in a column to make inferences from it; "
     "if zero, set to same value as MIN_READSTACK_DEPTH");
  CommandArgument_Bool_OrDefault_Doc(SKIP_UNDER, efp_default.skip_under, 
     "Don't do anything with bases under kmers");
  CommandArgument_Int_OrDefault_Doc(AUTO_CONFIRM, efp_default.auto_confirm, 
     "if set, automatically confirm bases whose quality score sum is at least this");

  CommandArgument_Bool_OrDefault_Doc(DO_PALINDROMES, efp_default.do_palindromes, 
     "Includes read stacks built from palidromic kmers.");
  CommandArgument_Bool_OrDefault_Doc(DO_BRANCHES, efp_default.do_branches, 
     "Continues error correction after finding a suspected branching.");
  CommandArgument_Int_OrDefault_Doc(MIN_QSS_TO_SUPPORT, efp_default.min_QSsum_to_support, 
     "If positive, minimum quality score sum of supporting bases of a nucleotide to "
     "mark it 'supported'; overrides MIN_N_BASES_TO_SUPPORT & MIN_Q_TO_SUPPORT");
  CommandArgument_Int_OrDefault_Doc(MIN_N_BASES_TO_CONFIRM, efp_default.min_n_bases_to_confirm,
     "Minimum number supporting bases of a nucleotide to mark it 'confirmed'");
  CommandArgument_Int_OrDefault_Doc(QUAL_CEIL_RADIUS, efp_default.qual_ceil_radius, 
     "Replace each quality score by the minimum of the quality scores " 
     "over a window of this radius.");
  CommandArgument_Bool_OrDefault_Doc(QCR_ALWAYS, efp_default.qcr_always, 
     "Apply QUAL_CEIL_RADIUS every cycle. Doesn't really make sense to do "
     "it more than once but may improve results.");



  // Logging control:

  CommandArgument_UnsignedInt_OrDefault_Doc(VERBOSITY, 1, "Verbosity of logging.");
  CommandArgument_LongLongSet_OrDefault_Doc(I_READS_PRINT, "",
     "List of read indexes for stack printing.");

  // These options are broken at this point.

  CommandArgument_Int_OrDefault_Doc(MAX_PRINT, -1, "controls verbose logging: "
     "if set, terminate execution after printing this many columns");
  CommandArgument_Int_OrDefault_Doc(PRINT_STACK_N, efp_default.print_stack_n, 
     "Number of read stacks to print per thread");
  CommandArgument_Int_OrDefault_Doc(PRINT_STACK_N0, efp_default.print_stack_n0, 
     "Starting index of read stacks to print");
  CommandArgument_Int_OrDefault_Doc(PRINT_STACK_MIN_DEPTH, efp_default.print_stack_min_depth, 
     "Minimum size of a read stack to print");
  CommandArgument_Int_OrDefault_Doc(PRINT_STACK_MAX_DEPTH, efp_default.print_stack_max_depth, 
     "Maximum size of a read stack to print");
  CommandArgument_Int_OrDefault_Doc(PRINT_STACK_STYLE, efp_default.print_stack_style, 
     "Read stack print style: 0:BW 1:16col. 2:256col 3:256col{^(-.} "
     "(positive | negative:with | without error recommendations)");

  CommandArgument_String_OrDefault_Doc(TRACE_IDS, "", "filter processing: if "
     "specified (in ParseIntSet format), only process stacks intersecting these "
                                       "read ids - note that this affects output!");


  EndCommandArguments;
  



  // ---- Thread control and memory limits
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  SetMaxMemory(MAX_MEMORY_GB << 30);


  // ---- Some parameter checking

  if (MIN_BASESTACK_DEPTH == 0)
    MIN_BASESTACK_DEPTH = MIN_READSTACK_DEPTH;
  
  if (MIN_READSTACK_DEPTH < MIN_BASESTACK_DEPTH)
    MIN_READSTACK_DEPTH = MIN_BASESTACK_DEPTH;

  const unsigned K_PC = (K & 1u) ? K : K + 1; // used for pre-correction 

  cout << Tag() << " K_FE    = " << setw(12) << K << endl;
  cout << Tag() << " K_PC    = " << setw(12) << K_PC << endl;
  

  if (PRE_CORRECT && K_PC > 29) {
    cout << Tag() << "Currently PRE_CORRECT=True only works for K <= 29." << endl; 
    exit(1);
  }   



  // ---- Getting ploidy
  
  const unsigned ploidy = (PLOIDY_FILE != "" && IsRegularFile(PLOIDY_FILE)) ?
    StringOfFile(PLOIDY_FILE, 1).Int() : 1;
  
  cout << Tag() << " ploidy  = " << setw(12) << ploidy << endl;


  // ---- load reads

  const String fastb_fn = HEAD_IN + ".fastb";
  cout << Tag() << "Reading bases from '" << fastb_fn << "'." << endl;
  BaseVecVec bases(fastb_fn);

  const size_t n_reads = bases.size();
  size_t n_bases = 0;
  for (size_t i = 0; i != n_reads; i++) 
    n_bases += bases[i].size();

  float gc_mu = gc_content_compute(bases);

  cout << Tag() << setw(12) << n_reads << "  reads read." << endl;

  // ---- write some stats
  { 
    PerfStat::log() << std::fixed << std::setprecision(0)
                    << PerfStat("n_frag_reads", 
                                "total number of original fragment reads", 
                                n_reads);
    PerfStat::log() << std::fixed << std::setprecision(1)
                    << PerfStat("mean_frag_read_len", 
                                "mean length of original fragment reads in bases", 
                                float(n_bases) / float(n_reads));
    PerfStat::log() << std::fixed << std::setprecision(1)
                    << PerfStat("gc_content_frags", 
                                "% gc content of fragment reads", 
                                100.0 * gc_mu);
  }


  // ---- quals definition

  const String qualb_fn = HEAD_IN + ".qualb";
  vec<BVLocFE> quals0; // hack to store quals of 0 (TO CHANGE)


  // ---- Kmer spectrum definition

  KmerSpectrum kspec(K_PC);







  // ==== FIRST TASK pre-correct reads
  // 
  //   note that it uses 8 bit quals; should change it at some point.
  //






  if (PRE_CORRECT) {

    // ---- Reading 8 bit quals (TO CHANGE)
    
    cout << Tag() << "Reading quals from '" << qualb_fn << "'." << endl;
    QualVecVec quals8(qualb_fn);
    assert_bases_quals_consistency(bases, quals8);
    
    // ---- Heuristics for pre-correction algorithm
    
    PC_Params pcp; // defaults are OK
    pcp.i_reads_print.insert(I_READS_PRINT.begin(), I_READS_PRINT.end());
    
    // ---- Do the precorrection
    
    cout << Tag() << "Pre-correcting." << endl;
    pre_correct_parallel(pcp, K_PC, & bases, & quals8, & kspec, 
                         VERBOSITY, NUM_THREADS);
    
    // ---- HACK! find the quals that were set to zero during pre-correction
    //      to add them back after reading the quals as nibble 

    for (size_t iv = 0; iv < n_reads; iv++) {
      const QualVec & qv = quals8[iv];
      for (size_t i = 0; i < qv.size(); i++) 
        if (qv[i] == 0) quals0.push_back(BVLocFE(iv, i));
    }
  } 
  else {

    // ---- compute kmer spectrum of input reads

    kmer_spectrum_compute(bases, & kspec, 
                          VERBOSITY, NUM_THREADS);
  }




  // ---- kmer spectrum analysis

  cout << Tag() << "Analysing kmer spectrum." << endl;

  kspec.to_text_file(HEAD_IN);

  const size_t KF_LOW = 10; // under-estimate of kmer spec minimum
  genome_analysis_report(kspec, n_bases / n_reads, ploidy, KF_LOW, VERBOSITY);
  genome_size_perf_stats(kspec, n_bases / n_reads, ploidy, KF_LOW, VERBOSITY);


  // ==== 2nd TASK: Error Correction ====




  if (true) {

    // ---- Loading quals as nibble.  

    cout << Tag() << "Reading quals as nibble from '" << qualb_fn << "'." << endl;
    QualNibbleVecVec quals;
    LoadQualNibbleVec(qualb_fn, &quals);
    assert_bases_quals_consistency(bases, quals);

    // ---- HACK! reset quals that were touched by pre-correct

    {
      const size_t n0 = quals0.size();
      for (size_t i = 0; i < n0; i++) 
        quals[quals0[i].iv].set(quals0[i].i, 0);
    }

    // ---- Estimate heuristics 

    if (USE_KMER_SPECTRUM || ESTIMATED_COVERAGE > 0) {
      heuristics_estimate(ESTIMATED_COVERAGE,
                          ploidy,
                          bases, 
                          kspec,                        
                          MIN_READSTACK_DEPTH,
                          MIN_BASESTACK_DEPTH,
                          MIN_QSS_TO_CONFIRM,
                          MIN_QSS_TO_SUPPORT,
                          MIN_QSS_TO_WIN);
    }
    
    // ---- The error correction parameters
    
    EF_Params efp(MIN_Q_TO_SUPPORT, 
                  MIN_N_BASES_TO_SUPPORT, 
                  MAX_QSS_RATIO_TO_CORRECT, 
                  MAX_QSS_RATIO_TO_CORRECT2, 
                  MAX_QSS_RATIO_TO_CONFIRM,
                  MIN_QSS_TO_CONFIRM,
                  MIN_QSS_TO_WIN,
                  MIN_READSTACK_DEPTH, 
                  MAX_READSTACK_DEPTH,
                  MIN_BASESTACK_DEPTH,
                  SKIP_UNDER,
                  AUTO_CONFIRM,
                  PRINT_STACK_N, 
                  PRINT_STACK_N0, 
                  PRINT_STACK_MIN_DEPTH, 
                  PRINT_STACK_MAX_DEPTH,
                  PRINT_STACK_STYLE,
                  DO_PALINDROMES,  // remove this
                  DO_BRANCHES,
                  MIN_QSS_TO_SUPPORT,
                  MIN_N_BASES_TO_CONFIRM,
                  QUAL_CEIL_RADIUS,
                  QCR_ALWAYS,
                  I_READS_PRINT);
    
    // ==== error correct the reads
    
    
    find_errors_parallel(efp, K, &bases, &quals, NUM_CYCLES, 
                         VERBOSITY, NUM_THREADS); 
    
    // ---- save reads
    
    if (WRITE) {
      const String fastb_fn = HEAD_OUT + ".fastb";
      cout << Tag() << "Writing bases to '" << fastb_fn << "'." << endl;
      bases.WriteAll(fastb_fn);
      
      const String qualb_fn = HEAD_OUT + ".qualb";
      cout << Tag() << "Writing quals to '" << qualb_fn << "'." << endl;
      WriteAll(quals, qualb_fn);
    }


  }

 
  cout << Tag() << "END" << endl;
}
