///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CleanCorrectedReads

#include "MainTools.h"
#include "reporting/PerfStat.h"

#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/QualNibbleVec.h"

#include "util/ReadTracker.h"

#include "kmers/KmerSpectrumCore.h"

#include "kmers/naif_kmer/NaifKmerizer.h"
#include "kmers/naif_kmer/KernelUniqueFinder.h"
#include "kmers/naif_kmer/KernelKmerSpectralizer.h"

#include "paths/PolymorphismRemoveCore.h"





static inline 
String Tag(String S = "CCR") { return Date() + " (" + S + "): "; } 













void heuristics_estimate(const u_int ESTIMATED_COVERAGE,
                         const unsigned ploidy, 
                         const BaseVecVec & bases,
                         const KmerSpectrum & kspec,
                         int & MAX_KMER_FREQ_TO_MARK)
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

  // Set kmer cutoff halfway between 1 and the spectrum minimum (log scale)
  MAX_KMER_FREQ_TO_MARK = round(sqrt(kf_lo));
  cout << Tag() << "Setting MAX_KMER_FREQ_TO_MARK to " << MAX_KMER_FREQ_TO_MARK << endl;
}









void report_on_traced_reads(const BaseVecVec & bases,
                            const vec<BaseVec> & traced_reads_orig,
                            const vec<Bool> & remove_read, 
                            const vec<int64_t> & trace_ids) 
{
  if (trace_ids.nonempty()) {
    cout << endl << "STATUS OF TRACED READS:" << endl;
    for (size_t i = 0; i < trace_ids.size(); i++) {
      size_t id = trace_ids[i];
      cout << ">" << id << "_orig (" << (remove_read[id] ? "DISCARD" : "KEEP") << ")" << endl;
      cout << traced_reads_orig[i].ToString() << endl;
      cout << ">" << id << "_edit (" << (remove_read[id] ? "DISCARD" : "KEEP") << ")" << endl; 
      cout << bases[id].ToString() << endl;    
    }
    cout << endl;
  }
}








template <class BOOL_t>
void find_reads_with_rare_kmers_parallel(const size_t       K, 
                                         const size_t       kf_rare_max, // kmer frequency
                                         const BaseVecVec & bases, 
                                         vec<BOOL_t>      * remove_p,
                                         KmerSpectrum     * kspec_p,
                                         const size_t       VERBOSITY,
                                         const size_t       NUM_THREADS,
                                         const size_t       mem_mean_ceil = 0)
{
  const size_t n_reads = bases.size();
  if (VERBOSITY > 0) cout << Tag() << "Finding reads with unique kmers." << endl;
    
  // ---- Find reads with unique kmers

  ForceAssertLe(K, 124u);
  if (K <= 29) {
    UniqueFinder<Kmer29, BOOL_t> unique_finder(K, kf_rare_max, bases, 
                                               remove_p, kspec_p, NUM_THREADS);
    naif_kmerize(&unique_finder, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 60) {
    UniqueFinder<Kmer60, BOOL_t> unique_finder(K, kf_rare_max, bases, 
                                               remove_p, kspec_p, NUM_THREADS);
    naif_kmerize(&unique_finder, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 124) {
    UniqueFinder<Kmer124, BOOL_t> unique_finder(K, kf_rare_max, bases, 
                                               remove_p, kspec_p, NUM_THREADS);
    naif_kmerize(&unique_finder, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }

    
  if (VERBOSITY > 0) cout << Tag() << "Done finding reads with unique kmers." << endl;


  size_t n_reads_low_kf = 0;
  for (size_t i = 0; i < n_reads; ++i) 
    if ((*remove_p)[i]) n_reads_low_kf++;

  if (VERBOSITY > 0) cout << Tag() << "n_reads_low_kf = " << n_reads_low_kf << endl;

}
















void remove_reads(const vec<Bool>  & remove_read,
                  BaseVecVec       * p_bases,
                  QualNibbleVecVec * p_quals,
                  PairsManager     * p_pairs,
                  vec<size_t>      * p_i_reads_in,
                  const bool         SAVE_IN_READ_IDS,
                  const unsigned     VERBOSITY)
{
  if (VERBOSITY > 0) cout << Tag() << "Removing unique reads." << endl;
  const size_t n_reads_in = p_bases->size();
  const size_t n_pairs_in = p_pairs->nPairs();

  const size_t n_reads_paired_in = 2 * n_pairs_in;
  const size_t n_reads_unpaired_in = n_reads_in - n_reads_paired_in;


  // Tally the number of reads being removed from each part of each pair.

  size_t n_A_reads_removed = 0;
  size_t n_B_reads_removed = 0;
  size_t n_pairs_removed = 0;
  for (size_t ip = 0; ip < n_pairs_in; ++ip) {
    const bool remove_A = remove_read[p_pairs->ID1(ip)];
    const bool remove_B = remove_read[p_pairs->ID2(ip)];
    
    if ( remove_A && !remove_B) n_A_reads_removed++;
    if (!remove_A &&  remove_B) n_B_reads_removed++;
    if ( remove_A &&  remove_B) n_pairs_removed++;
  }
  

  // ---- Remove reads from dataset

  // Shift the elements of reads and lane_index into their new locations.
  // Along the way, create a mapping from old read ID to new read ID.
  
  size_t n_reads_out = 0;
  
  for (size_t i_read_in = 0; i_read_in < n_reads_in; ++i_read_in) {
    if (!remove_read[i_read_in]) {

      if (SAVE_IN_READ_IDS)
        (*p_i_reads_in).push_back(i_read_in);

      if (i_read_in != n_reads_out) {
        (*p_bases).SwapElements(n_reads_out, i_read_in);
        (*p_quals).SwapElements(n_reads_out, i_read_in);
      }
      n_reads_out++;
    }
  }

  // Resize the data structures.

  (*p_bases).resize(n_reads_out);
  (*p_quals).resize(n_reads_out);


  // update the pairs info.

  const bool allow_pair_severance = true;
  p_pairs->removeReads(remove_read, allow_pair_severance);
  
  ForceAssertEq(n_reads_out, (*p_bases).size());
  ForceAssertEq(n_reads_out, (*p_quals).size());
  ForceAssertEq(n_reads_out, p_pairs->nReads());
  


  
  // Add up the tallies.

  const size_t n_pairs_out = p_pairs->nPairs();

  const size_t n_reads_removed = n_reads_in - n_reads_out;
  const size_t n_reads_paired_removed = (n_A_reads_removed + 
                                         n_B_reads_removed + 
                                         2 * n_pairs_removed);
  
  const size_t n_reads_unpaired_removed = n_reads_removed - n_reads_paired_removed;
  
  const size_t n_reads_unpaired_gained = n_A_reads_removed + n_B_reads_removed;
  const size_t n_reads_unpaired_out = (n_reads_unpaired_in - 
                                       n_reads_unpaired_removed +
                                       n_reads_unpaired_gained);
  // Report to cout.

  if (VERBOSITY > 0) {
    cout << endl;
    cout << Tag() << "REPORT" << endl;

    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_in)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 << " %)" 
         << "  reads in." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(2 * n_pairs_in)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * 2 * n_pairs_in / n_reads_in << " %)" 
         << "  paired reads in." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_unpaired_in)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_reads_unpaired_in / n_reads_in << " %)" 
         << "  unpaired reads in." << endl;

    cout << Tag() << endl;

    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_out)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_reads_out / n_reads_in << " %)" 
         << "  reads out." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(2 * n_pairs_out)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * 2 * n_pairs_out / n_reads_in << " %)" 
         << "  paired reads out." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_unpaired_out)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_reads_unpaired_out / n_reads_in << " %)" 
         << "  unpaired reads out." << endl;

    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_unpaired_gained)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_reads_unpaired_gained / n_reads_in << " %)" 
         << "  unpaired reads gained from broken pairs." << endl;

    cout << Tag() << endl;


    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_removed)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_reads_removed / n_reads_in << " %)" 
         << "  reads removed." << endl;

    cout << Tag() 
         << setw(14) << ToStringAddCommas(2 * n_pairs_removed)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * 2 * n_pairs_removed / n_reads_in << " %)" 
         << "  A and B reads removed." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_A_reads_removed)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_A_reads_removed / n_reads_in << " %)" 
         << "  A reads removed where B read is kept." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_B_reads_removed)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_B_reads_removed / n_reads_in << " %)" 
         << "  B reads removed where A read is kept." << endl;
    cout << Tag() 
         << setw(14) << ToStringAddCommas(n_reads_unpaired_removed)
         << " (" << setw(5) << fixed << setprecision(1) << 100.0 * n_reads_unpaired_removed / n_reads_in << " %)" 
         << "  unpaired reads removed." << endl;

    
    PerfStat::log() << std::fixed << std::setprecision(1)
                    << PerfStat("frac_reads_removed", 
                                "% of reads removed because of low frequency kmers", 
                                100.0 * (float(n_reads_in) - float(n_reads_out)) / float(n_reads_in));
  }



  // Sanity checks: make sure the output is not degenerate.
  if (n_reads_in > 0 && n_reads_out == 0)
    FatalErr("FindErrors removed all reads!");
  if (n_pairs_in > 0 && n_pairs_out == 0)
    FatalErr("FindErrors removed all read pairings!");
  
}




















int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  
  // File locations:
  
  CommandArgument_String_Doc(HEAD_IN, "Looks for 'HEAD_IN.{fastb,qualb}'.");
  CommandArgument_String_Doc(HEAD_OUT, "Generates 'HEAD_OUT.{fastb,keep,log}'.");
  CommandArgument_String_OrDefault_Doc(PLOIDY_FILE, "", 
     "If file exists, reads ploidy.");
  

  // Computational performance:

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0).");
  CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB, 0);

  // General control:

  CommandArgument_Int_OrDefault_Doc(K, 24, 
     "Kmer size for error correction.");
  CommandArgument_Int_OrDefault_Doc(MAX_KMER_FREQ_TO_MARK, 2, 
     "Reads with minimum kmer frequencies at and below this will be marked and discarded in _corr.");
  CommandArgument_String_OrDefault_Doc(HEAD_KMER_SPECTRUM_IN, "",
     "If set, looks for a kmer spectrum in '<HEAD_KMER_SPECTRUM_IN>.<K_PC>.kspec'" 
     " and, based on its features, automatically adjust MAX_KMER_FREQUENCY_TO_MARK.");
  CommandArgument_UnsignedInt_OrDefault_Doc(ESTIMATED_COVERAGE, 0, 
     "If positive, use this as estimated coverage for adjusting heuristics instead of kmer spectrum.");

  CommandArgument_Bool_OrDefault_Doc(HAPLOIDIFY, False, 
     "Optionally remove polymorphisms after error correction.");
  CommandArgument_UnsignedInt_OrDefault_Doc(K_HAP, 25, 
     "Kmer size for polymorphism removal.");


  // Logging control:

  CommandArgument_UnsignedInt_OrDefault_Doc(VERBOSITY, 1, "Verbosity of logging.");

  CommandArgument_String_OrDefault_Doc(TRACE_IDS, "", "filter processing: if "
     "specified (in ParseIntSet format), only process stacks intersecting these "
                                       "read ids - note that this affects output!");


  // Arguments for removed_marked_reads()
  CommandArgument_Bool_OrDefault_Doc(DELETE, True, 
     "actually delete reads, rather than mark them for deletion");
  CommandArgument_Bool_OrDefault_Doc(SAVE_IN_READ_IDS, False,
    "for each OUT read_ID save the IN read_ID to a text file");
  CommandArgument_Bool_OrDefault_Doc(TRACK_READS, True,
    "use ReadTracker to ouput a mapping of old reads to new");



  EndCommandArguments;






  // Thread control and memory limits

  NUM_THREADS = configNumThreads(NUM_THREADS);
  SetMaxMemory(MAX_MEMORY_GB << 30);

  
    

  // ---- Getting ploidy
  
  const unsigned ploidy = (PLOIDY_FILE != "" && IsRegularFile(PLOIDY_FILE)) ?
    StringOfFile(PLOIDY_FILE, 1).Int() : 1;
  
  cout << Tag() << " K       = " << setw(12) << K << endl;
  cout << Tag() << " ploidy  = " << setw(12) << ploidy << endl;


  // ---- Adjust heuristics from kmer spectrum or estimated coverage
  
  if (HEAD_KMER_SPECTRUM_IN != "" || ESTIMATED_COVERAGE > 0) {
    
    // get kmer spectrum before error correction
    const unsigned K_PC = (K & 1u) ? K : K + 1; 

    KmerSpectrum kspec(K_PC);
    kspec.from_text_file(HEAD_KMER_SPECTRUM_IN, VERBOSITY);

    void heuristics_estimate(const u_int ESTIMATED_COVERAGE,
                             const unsigned ploidy, 
                             const BaseVecVec & bases,
                             const KmerSpectrum & kspec,
                             int & MAX_KMER_FREQ_TO_MARK);
  }



  // ---- Load reads

  const String fastb_fn = HEAD_IN + ".fastb";
  const String qualb_fn = HEAD_IN + ".qualb";

  cout << Tag() << "Reading bases from '" << fastb_fn << "'." << endl;
  BaseVecVec bases(fastb_fn);;

  cout << Tag() << "Reading quals from '" << qualb_fn << "'." << endl;
  QualNibbleVecVec quals;
  LoadQualNibbleVec(qualb_fn, &quals);

  ForceAssertEq(bases.size(), quals.size());

  const size_t n_reads = bases.size();
  cout << Tag() << setw(12) << n_reads << "  reads read." << endl;


  // ---- Read tracing

  vec<int64_t> trace_ids;
  ParseLongLongSet(TRACE_IDS, trace_ids, true, 0, 
                   MastervecFileObjectCount(HEAD_IN + ".fastb"));

  vec<BaseVec> traced_reads_orig;
  for (size_t i = 0; i < trace_ids.size(); i++)
    traced_reads_orig.push_back(bases[trace_ids[i]]);




  // ==== Find reads with unique kmers
 
  vec<Bool> bad_read(n_reads, false);
    
  if (true) {

    KmerSpectrum kspec(K); // kmer spectrum of input reads

    find_reads_with_rare_kmers_parallel(K, MAX_KMER_FREQ_TO_MARK, bases, 
                                        &bad_read, &kspec,
                                        VERBOSITY, NUM_THREADS);
    
    kspec.to_text_file(HEAD_IN, VERBOSITY);
    
    if (TRACE_IDS != "")
      report_on_traced_reads(bases, traced_reads_orig, bad_read, trace_ids);
  }







  // ==== Delete bad reads and optionally take out ambiguities



  if (DELETE || HAPLOIDIFY) {

    // ---- Remove reads adjusting also the paring info

    cout << Tag() << "Loading pairing info from '" << HEAD_IN << ".pairs'." << endl;
    PairsManager pairs(HEAD_IN + ".pairs");
    cout << Tag() << "Loaded " << pairs.nPairs() << " pairs." << endl;
    vec<size_t> i_reads_in;
    remove_reads(bad_read, &bases, &quals, &pairs, &i_reads_in,
                 SAVE_IN_READ_IDS, VERBOSITY);
    
    
    // ---- Remove polymorphisms

    if (HAPLOIDIFY) {
      cout << Tag() << "Removing ambiguities/polymorphisms from reads at K=" << K_HAP << "." << endl;

      // ---- Find polymorphisms/ambiguities

      Polymorphisms polys(K_HAP);
      KmerSpectrum kspec_hap(K_HAP);

      polymorphisms_find_parallel(K_HAP, bases, & polys, & kspec_hap, 
                                  VERBOSITY, NUM_THREADS);
      polys.print_stats();
      
      kspec_hap.to_text_file(HEAD_OUT + ".diploid", VERBOSITY);

      
      // ---- Save polymorphisms/ambiguities to fasta files
      {
        const String fn_poly_stats = HEAD_OUT + ".k" + ToString(K_HAP) + ".poly_stats";
        cout << Tag() << "Saving polymorphism statistics to '" << fn_poly_stats << "'." << endl;
        polys.write_stats(fn_poly_stats); 
        
        const String head_poly_in  = HEAD_OUT + ".poly.k" + ToString(K_HAP) + ".in";
        const String head_poly_out = HEAD_OUT + ".poly.k" + ToString(K_HAP) + ".out";

        cout << Tag() << "Saving kept ambiguity branches to '" << head_poly_in << ".fasta'." << endl;
        cout << Tag() << "Saving discarded ambiguity branches to '" << head_poly_out << ".fasta'." << endl;
        polys.write_fastas(head_poly_in, head_poly_out);
      }

      // ----- Remove polymorphisms/ambiguities from corrected reads 

      cout << Tag() << "Removing polymorphisms from reads." << endl;
      const size_t n_fixed = 
        polymorphisms_remove_parallel(polys, &bases, &quals, VERBOSITY, NUM_THREADS);
      
      cout << Tag() << setw(12) << n_fixed << "  reads changed by polymorphism removal." << endl;
      

    }


    // ---- Write reads and pairing info to disk

    cout << Tag() << "Writing output files." << endl;
    
    if (true) {
      const String fastb_fn = HEAD_OUT + ".fastb";
      const String qualb_fn = HEAD_OUT + ".qualb";
      const String pairs_fn = HEAD_OUT + ".pairs";
      cout << Tag() << "Writing bases to '" << fastb_fn << "'." << endl;
      bases.WriteAll(fastb_fn);
      cout << Tag() << "Writing quals to '" << qualb_fn << "'." << endl;
      WriteAll(quals, qualb_fn);    
      cout << Tag() << "Writing pairs to '" << pairs_fn << "'." << endl;
      pairs.Write(pairs_fn);
    }


    // ---- Compute final kmer spectrum
    if (true) {
      const unsigned K_odd = (K & 1) ? K : K + 1;
      
      cout << Tag() << "Computing kmer spectrum of corrected reads at K=" << K_odd << "." << endl;
      KmerSpectrum kspec(K_odd);
      kmer_spectrum_compute(bases, & kspec, 
                            VERBOSITY, NUM_THREADS);
      kspec.to_text_file(HEAD_OUT, VERBOSITY);
    }

    
    // ---- Some optional stuff
    
    if (SAVE_IN_READ_IDS) {
      ofstream outfs;
      outfs.open((HEAD_OUT + ".in_read_IDs").c_str());
      outfs << "# out_read_ID in_read_ID" << endl;
      
      const size_t n = bases.size();
      for (size_t i = 0; i < n; i++)
        outfs << i << " " << i_reads_in[i] << endl;
      outfs.close();
    }
    
    if (TRACK_READS) {
      ReadTracker rt;
      rt.AddReadSet(HEAD_IN, bad_read);
      rt.Dump(HEAD_OUT);
    }
    
        
  }





  cout << Tag() << "END" << endl;
}

















