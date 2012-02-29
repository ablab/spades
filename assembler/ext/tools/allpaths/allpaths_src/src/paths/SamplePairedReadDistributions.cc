///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: SamplePairedReadDistributions

   Samples the invariant size distributions of paired reads.
   This is done by aligning a random sample of pairs to a reference.

   The reference can be a set of unibases.

   NOTES:
   
     "Invariant size" means the size that is invariant under read trimming.
     Because some sheared jump reads are, at times, flipped prior to being fed 
     to this module, the parameter FLIP has been added to account for this. 


     FRAGMENT READS:
 
        FLIP=False     invariant size is POSITIVE = full fragment size


     SHEARED JUMPS:   

        FLIP=True      invariant size is NEGATIVE = -|separation|


     NON-SHEARED JUMPS:  

        FLIP=False     invariant size is POSITIVE = full insert size


   INPUT FILES:
     reference.fastb
     HEAD_READS.{fastb,pairs}

   OUTPUT FILES:
     HEAD_READS.*

*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#include <omp.h>

#include "MainTools.h"
#include "system/SysConf.h"  // processorsOnline()
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "reporting/PerfStat.h"


#include "PairsManager.h"
#include "lookup/FirstLookupFinderECJ.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTabBuilder.h"
#include "random/Shuffle.h"
#include "paths/UnibaseUtils.h"

#include "math/IntDistribution.h"
#include "math/IntFrequencies.h"

static inline 
String Tag(String S = "SPRD") { return Date() + " (" + S + "): "; } 




String pct_blank_str(const float x, const unsigned precision = 1)
{
  return (x < 0.001) ? "" : ToString(100.0 * x, precision) + "%";
}


String cov_blank_str(const float x, const unsigned precision = 1)
{
  return (abs(x) < 0.1) ? "" : ToString(x, precision) + "x";
}


String genome_size_str(const float G)
{
  if (G >= 1e9) return ToString(G/1e9, 1) + " Gb";
  if (G >= 1e8) return ToString(G/1e6, 0) + " Mb";
  if (G >= 1e7) return ToString(G/1e6, 0) + " Mb";
  if (G >= 1e6) return ToString(G/1e6, 1) + " Mb";
  if (G >= 1e5) return ToString(G/1e3, 0) + " kb";
  if (G >= 1e4) return ToString(G/1e3, 0) + " kb";
  if (G >= 1e3) return ToString(G/1e3, 1) + " kb";
  return               ToString(G,     1) + "  b";
}



void print_cov_line(const IntFunction<double> & cov)
{
  cout << "  " << setw(7) << cov_blank_str(cov.sum_below(0));
  cout << "  " << setw(7) << cov_blank_str(cov.sum_in(    0,   499));
  cout << " "  << setw(7) << cov_blank_str(cov.sum_in(  500,  1999));
  cout << " "  << setw(7) << cov_blank_str(cov.sum_in( 2000,  7999));
  cout << " "  << setw(7) << cov_blank_str(cov.sum_in( 8000, 31999));
  cout << " "  << setw(7) << cov_blank_str(cov.sum_above(32000));
  cout << "  " << setw(7) << cov_blank_str(cov.sum());
}



    
void print_cov_line2(const IntFunction<double> & cov)
{
  cout << " "  << setw(8) << cov_blank_str(cov.sum_above(    0));
  cout << " "  << setw(8) << cov_blank_str(cov.sum_above(  300));
  cout << " "  << setw(8) << cov_blank_str(cov.sum_above( 1000));
  cout << " "  << setw(8) << cov_blank_str(cov.sum_above( 3000));
  cout << " "  << setw(8) << cov_blank_str(cov.sum_above(10000));
  cout << " "  << setw(8) << cov_blank_str(cov.sum_above(30000));
}


void print_cov_line3(const IntFunction<double> & cov)
{
  cout << " "  << setw(8) << ToString(cov.sum_above(    0), 1);
  cout << " "  << setw(8) << ToString(cov.sum_above( 1000), 1);
  cout << " "  << setw(8) << ToString(cov.sum_above( 2000), 1);
  cout << " "  << setw(8) << ToString(cov.sum_above( 3000), 1);
  cout << " "  << setw(8) << ToString(cov.sum_above( 4000), 1);
  cout << " "  << setw(8) << ToString(cov.sum_above( 5000), 1);
}

void print_n_gap_links_line(const IntFunction<double> & cov,
                            const IntFunction<double> & n_links,
                            const vec<size_t> & gaps)
{
  for (size_t i = 0; i < gaps.size(); i++) {
    const size_t gap = gaps[i];
    const size_t n = cov.sum_above(gap) - gap * n_links.sum_above(gap) + 0.5;
    cout << " "  << setw(5) << (n > 0 ? ToString(n) : "     ");
  }
}
    





class Library
{
  String name;
  size_t n0;
  size_t n_sampled;
  size_t n_hits;
  size_t n_ilogical;

};





















// ---- this filters out outliers with few hits from the -inf and +inf ends
void distribution_trim_limits(const IntFrequencies & hits,
                              int * p_x_min,
                              int * p_x_max)
{
  if (hits.size() > 0) {
    const int x_h_min = hits.x_min();
    const int x_h_max = hits.x_max();
    size_t freq_max = 0;
    size_t n_hits = 0;
    for (int x = x_h_min; x <= x_h_max; x++) {
      const size_t freq = hits.freq(x);
      if (freq_max < freq) 
        freq_max = freq;
      n_hits += freq;
    }
    
    const int diameter = n_hits / freq_max;
    const int dens_min = 20;
    
    *p_x_min = x_h_max;
    {
      int dens = 0;
      for (int x = x_h_min; x <= x_h_max && *p_x_min == x_h_max; x++) {
        dens += hits.freq(x);
        if (x - diameter >= x_h_min)
          dens -= hits.freq(x - diameter);
        if (dens >= dens_min)
          *p_x_min = x - diameter;
      }
    }
    if (*p_x_min < x_h_min) *p_x_min = x_h_min;
    
    
    *p_x_max = x_h_min;
    {
      int dens = 0;
      for (int x = x_h_max; x >= x_h_min && *p_x_max == x_h_min; x--) {
        dens += hits.freq(x);
        if (x + diameter <= x_h_max)
          dens -= hits.freq(x + diameter);
        if (dens >= dens_min)
          *p_x_max = x + diameter;
      }
    }
    if (*p_x_max > x_h_max) *p_x_max = x_h_max;
    
    //ForceAssertGt(*p_x_max, *p_x_min);
  }
}





IntDistribution distribution_compute(const IntFrequencies & inv_sz_freq,
                                     const IntFrequencies & n_inserts,
                                     const int x_min,
                                     const int x_max)
{
  IntFunction<double> counts;
  for (int x = x_min; x <= x_max; x++)
    counts[x] = double(inv_sz_freq[x]) / double(n_inserts[x]); 

  // the constructor of IntDistribution will convert and NORMALIZE the IntFunction<double>
  return counts;
}


IntFunction<double> physical_coverage_compute(const IntFrequencies & inv_sz_freq, 
                                              const IntFrequencies & n_inserts,
                                              const double frac_sampled,
                                              const int x_min,
                                              const int x_max)
{
  IntFunction<double> cov;
  for (int x = x_min; x <= x_max; x++) 
    cov[x] = double(inv_sz_freq[x]) * double(x) / (double(n_inserts[x]) * frac_sampled);
  
  return cov;
}


IntFunction<double> number_links_compute(const IntFrequencies & inv_sz_freq, 
                                         const IntFrequencies & n_inserts,
                                         const double frac_sampled,
                                         const int x_min,
                                         const int x_max)
{
  IntFunction<double> n_links;
  for (int x = x_min; x <= x_max; x++) 
    n_links[x] = double(inv_sz_freq[x]) / (double(n_inserts[x]) * frac_sampled);
  
  return n_links;  // just the raw counts scaled up from the sampled counts
}




// boxcar, if you're wondering, is just a rectangle function:
//                     _________________
// ___________________|                 |_______________
//              -radius        0        radius
//
void function_convolute_boxcar(IntFunction<double> * p_func,
                               const size_t radius)
{
  if (radius > 0) {
    const double diameter = 2 * radius + 1;
    const int x0 = p_func->x_min() - radius;
    const int x1 = p_func->x_max() + radius;
    
    const IntFunction<double> & func_in = *p_func;
    IntFunction<double> func_out(x0, x1);
    
    for (int x = x0; x <= x1; x++) {
      double mean = 0.0;
      const int y0 = x - radius;
      const int y1 = x + radius;

      for (int y = y0; y <= y1; y++) 
        mean += func_in[y];

      mean /= diameter;
      
      func_out[x] = mean;
    }
    
    *p_func = func_out;
  }
}



IntDistribution distribution_smooth(const IntDistribution & dist_orig,
                                    const IntFrequencies & n_inserts)
{
  const int x_prob_max = dist_orig.x_prob_max();
  const double rad_conv = 100.0 * 100.0 / n_inserts.freq(x_prob_max);
  //cout << Tag() << fixed << setw(16) << setprecision(1) 
  //     << rad_conv << " smoothing radius." << endl;
        
        
  IntFunction<double> func = dist_orig.probs(); 
  if (rad_conv >= 4) {
    // ---- convolute a few times with a rectangle function
    //      it's almost like convoluting with a gaussian
    function_convolute_boxcar(& func, int(0.25 * rad_conv)); 
    function_convolute_boxcar(& func, int(0.25 * rad_conv)); 
    function_convolute_boxcar(& func, int(0.25 * rad_conv)); 
    function_convolute_boxcar(& func, int(0.25 * rad_conv)); 
  }
  else {
    function_convolute_boxcar(& func, int(0.5 + rad_conv)); 
  }

  // ---- make sure that there are no p(x) = 0.0 in range [x_min, x_max]

  const int x_min = func.x_min();
  const int x_max = func.x_max();
  double f_small = func.f_max();
  vec<int> x_null;
  for (int x = x_min; x <= x_max; x++) {
    const double f = func[x];
    if (f == 0.0) 
      x_null.push_back(x);
    else if (f < f_small)
      f_small = f;    
  }
    
  for (size_t i = 0; i != x_null.size(); i++)
    func[x_null[i]] = f_small;

  // here I rely on the IntDistribution constructor that accepts an IntFunction
  return func;
}







void build_unibase_reference(const String & UNIBASES,
                             const unsigned UNIBASES_K,
                             const String & ref_head,
                             const double TARGET_UNIBASE_COVERAGE,
                             const size_t MIN_UNIBASE_LENGTH,
                             const bool NEW_LOOKUP_TABLE)
{
  
  const String ref_fn = ref_head + ".fastb";

  if (! IsRegularFile(ref_fn) || NEW_LOOKUP_TABLE) {
    cout << Tag() << "Loading unibases." << endl;
    BaseVecVec unibases(UNIBASES);
    const size_t n_unibases = unibases.size();
    cout << Tag() << n_unibases << " loaded." << endl;      
    vec<int> toRc;
    cout << Tag() << "Tagging rc copies." << endl;
    UnibaseInvolution(unibases, toRc, UNIBASES_K);
      
    cout << Tag() << "Building unibase reference file." << endl;

    const size_t ub_len_min = (MIN_UNIBASE_LENGTH < 2 * UNIBASES_K) ?
      2 * UNIBASES_K : MIN_UNIBASE_LENGTH;

    size_t ub_len_total = 0; 
    size_t ub_len_usable = 0; 
    vec<size_t> i_ub_sorted(n_unibases, vec<size_t>::IDENTITY);
    vec<size_t> ub_lens(n_unibases);
    vec<bool> ub_keep(n_unibases, true);

    for (size_t i_ub = 0; i_ub < n_unibases; i_ub++) {
      const size_t ub_len = unibases[i_ub].size();
      ub_lens[i_ub] = ub_len;
      ub_len_total += ub_len;
      if (ub_len > ub_len_min) ub_len_usable += ub_len;
      if (ub_keep[i_ub])
        ub_keep[toRc[i_ub]] = false; //don't use rc, we will not worry about palindromes
    }
    
    vec<size_t> ub_lens_sorted(ub_lens);
      
    ReverseSortSync(ub_lens_sorted, i_ub_sorted); 
      
    const size_t ub_len_total_approx = (ub_len_total - n_unibases * UNIBASES_K) / 2;
    const size_t ub_len_total_target = round(ub_len_total_approx * TARGET_UNIBASE_COVERAGE);
      
    cout << Tag() << setw(14) << ub_lens_sorted.front() << "  length of largest unibase." << endl;
    cout << Tag() << setw(14) << ub_len_total           << "  total unibases length." << endl;
    cout << Tag() << setw(14) << ub_len_usable          << "  total usable unibases length." << endl;
    cout << Tag() << setw(14) << ub_len_total_target    << "  target total unibases length (" 
         << ToString(100.0 * TARGET_UNIBASE_COVERAGE) << " %)." << endl;
    

    size_t ub_len_sum = 0;
    BaseVecVec ub_selected;
    for (size_t i = 0; i < n_unibases && ub_len_sum < ub_len_total_target; i++) {
      const size_t i_ub = i_ub_sorted[i];
      const size_t ub_len = unibases[i_ub].size();
      if (ub_keep[i_ub] && ub_len > ub_len_min) {
        ub_selected.push_back(unibases[i_ub]);
        ub_len_sum += ub_len;
      }
    }

    if (ub_selected.size() > 0) {
      cout << Tag() << setw(14) << ub_selected.size() << "  selected unibases." << endl;
      cout << Tag() << setw(14) << ub_selected.back().size() << "  length of shortest selected unibase." << endl;
      cout << Tag() << setw(14) << ub_selected.front().size() << "  length of largest selected unibase." << endl;
      cout << Tag() << setw(14) << ub_lens_sorted.front() << "  length of largest unibase." << endl;
      cout << Tag() << setw(14) << ub_len_sum << "  total bases covered." << endl;
    }
    else {
      FatalErr("ERROR: Couldn't find any usable unibases. Aborting."); 
    }
      
    ub_selected.WriteAll(ref_fn);
  }
}




IntFrequencies n_inserts_in_unibases(const BaseVecVec & ref, 
                                     const int read_sz, 
                                     const int MAX_SIZE)
{
  IntFrequencies n_inserts;
  const size_t n_bvs = ref.size();

  for (size_t i_bv = 0; i_bv != n_bvs; i_bv++) {
    const size_t bv_sz = ref[i_bv].size();
    for (int inv_sz = -MAX_SIZE; inv_sz <= MAX_SIZE; inv_sz++) {
      
      const size_t min_insert_sz = read_sz + abs(inv_sz - read_sz);
      if (min_insert_sz <= bv_sz)
        n_inserts[inv_sz] += bv_sz - min_insert_sz + 1;
    }
  }
  return n_inserts;
}






bool found_pair_alignment(const PairsManager & pairs, 
			  const BaseVecVec & reads, 
			  const QualNibbleVecVec & quals,
			  const size_t i_pair, 
                          const BaseVecVec & ref,
			  const FirstLookupFinderECJ & lfinder, 
			  const size_t K_lookup,
			  int * p_inv_sz,
			  bool * p_ilogical,
			  const bool FLIP, 
			  const int TRIM)
{
  size_t    id1 = pairs.ID1(i_pair);
  size_t    id2 = pairs.ID2(i_pair);
  BaseVec read1 = reads[id1];
  BaseVec read2 = reads[id2];

  list<first_look_align> hits1;
  list<first_look_align> hits2;
  if (FLIP) {
    read1.ReverseComplement();
    read2.ReverseComplement();
  }
  if (TRIM > 0) {
    read1.SetToSubOf(read1, TRIM, -1);
    read2.SetToSubOf(read2, TRIM, -1);
  }
  
  if (quals.size() > 0) {
    QualNibbleVec qual1 = quals[id1];
    QualNibbleVec qual2 = quals[id2];

    if (FLIP) {
      qual1.ReverseMe();
      qual2.ReverseMe();
    }

    lfinder.getAlignments(read1, qual1, id1, &hits1);
    lfinder.getAlignments(read2, qual2, id2, &hits2);
  }
  else {
    lfinder.getAlignments(read1, id1, &hits1);
    lfinder.getAlignments(read2, id2, &hits2);
  }


  if (hits1.size() != 1 || hits2.size() != 1)
    return false;

  const first_look_align & hit1 = hits1.front();
  const first_look_align & hit2 = hits2.front();
  
  if (hit1.target_loc.getContig() != hit2.target_loc.getContig())
    return false;  // not the same contig
  
  *p_ilogical = false;  
  if (hit1.is_FW() == hit2.is_FW()) { // throw away ilogical orientations
    *p_ilogical = true;
    return false;
  }
  
  const int s1 = hit1.get_start_on_target(read1.isize(), K_lookup);
  const int s2 = hit2.get_start_on_target(read2.isize(), K_lookup);
  
  *p_inv_sz = (hit1.is_FW()) ? 
    (s2 + read2.size()) - s1 :
    (s1 + read1.size()) - s2;
 
  // ---- HACK so that inariant sizes are positive for sheared jumps
  if (FLIP) *p_inv_sz = -*p_inv_sz;

  // ---- check alignment
  if (false) {
    const size_t i_tig = hit1.target_loc.getContig();
    BaseVec ref_bv;
    BaseVec read_f;
    BaseVec read_b;
    if (hit1.is_FW()) {    
      ref_bv.SetToSubOf(ref[i_tig], s1, *p_inv_sz);
      read_f = read1;
      read_b = read2;
    }
    else {
      ref_bv.SetToSubOf(ref[i_tig], s2, *p_inv_sz); 
      read_f = read2;
      read_b = read1;
    }
    read_b.ReverseComplement();

    cout << "1fw = " << hit1.is_FW() << endl;
    cout << "id1 = " << id1 << endl;
    cout << "id2 = " << id2 << endl;
    cout << "s1 = " << s1 << " " << endl;
    cout << "s2 = " << s2 << " " << endl;
    cout << "tig = " << hit1.target_loc.getContig() << endl;
    exit(1);
  }

  return true;
}




IntFrequencies invariant_sizes_from_aligments(const PairsManager         & pairs,
                                              const BaseVecVec           & reads,
                                              const QualNibbleVecVec     & quals,
                                              const BaseVecVec           & ref,
                                              const FirstLookupFinderECJ & lfinder,
                                              const size_t                 K_lookup,
                                              const size_t               & i_lib,
                                              const vec<size_t>          & is_pairs,
                                              const size_t                 MAX_SIZE,
                                              const size_t                 TARGET_SAMPLE_SIZE,
					      size_t                     * p_n_pairs_sampled,
					      size_t                     * p_n_hits,
					      size_t                     * p_n_ilogical,
                                              const bool                   FLIP,
                                              const int                    TRIM,
                                              const size_t                 NUM_THREADS)
{
  // ---- Estimate number of samples
  const size_t n_pairs = is_pairs.size();
  const size_t n_pairs_test = (n_pairs > 100) ? 100 : n_pairs;

  size_t n_pairs_align = 0;
  int inv_sz = 0;
  for (size_t ii_pair = 0; ii_pair < n_pairs_test; ii_pair ++) {
    
    const size_t i_pair = is_pairs[ii_pair];
    ForceAssertEq(i_lib, (size_t)pairs.libraryID(i_pair));
    
    bool ilogical = false;
    if (found_pair_alignment(pairs, reads, quals, i_pair, 
                             ref, lfinder, K_lookup,
                             & inv_sz, & ilogical,
                             FLIP, TRIM))
      if (unsigned(abs(inv_sz)) <= MAX_SIZE)
        n_pairs_align ++;
  }

  //cout << endl;
  //cout << "i_lib= " << i_lib << "   n_pairs= " << n_pairs << endl;
  //cout << "i_lib= " << i_lib << "   n_pairs_align= " << n_pairs_align << endl;

  *p_n_pairs_sampled = (n_pairs_align > 0) ? 
    TARGET_SAMPLE_SIZE * n_pairs_test / n_pairs_align :
    n_pairs;
  
  //cout << "i_lib= " << i_lib << "   n_pairs_sample= " << n_pairs_sample << endl;
  if (*p_n_pairs_sampled > n_pairs) 
    *p_n_pairs_sampled = n_pairs;
  //cout << "i_lib= " << i_lib << "   n_pairs_sample= " << n_pairs_sample << endl;


  // ---- Sample reads in parallel

  vec< vec<unsigned> > inv_szss(NUM_THREADS);  // one for each thread

  #pragma omp parallel for 
  for (size_t i_thread = 0; i_thread < NUM_THREADS; i_thread++) {
    vec<unsigned> & inv_szs = inv_szss[i_thread];

    const size_t ii0_pair = (*p_n_pairs_sampled *  i_thread     ) / NUM_THREADS;
    const size_t ii1_pair = (*p_n_pairs_sampled * (i_thread + 1)) / NUM_THREADS;

    //cout << "i_thread= " << i_thread 
    //     << "   ii0_pair= " << ii0_pair 
    //     << "   ii1_pair= " << ii1_pair << endl;
    
    for (size_t ii_pair = ii0_pair; ii_pair != ii1_pair; ii_pair++) {
      const size_t i_pair = is_pairs[ii_pair];

      //cout << "GREP " << i_pair;
      int inv_sz = 0;
      bool ilogical = false;
      if (found_pair_alignment(pairs, reads, quals, i_pair, 
			       ref, lfinder, K_lookup,
			       & inv_sz, & ilogical,
			       FLIP, TRIM)) {
        if (unsigned(abs(inv_sz)) <= MAX_SIZE) {
          inv_szs.push_back(inv_sz);
          #pragma omp critical
          if (ilogical) (*p_n_ilogical)++;
        }
        //cout << " " << inv_sz << endl;
      }
      //else cout << " " << 0 << " not" << endl;
    }
  }

  // ---- bring the results from all threads together
  IntFrequencies inv_sz_freq;
  for (size_t it = 0; it != NUM_THREADS; it++)
    for (size_t i = 0; i != inv_szss[it].size(); i++) {
      inv_sz_freq[inv_szss[it][i]]++;
      (*p_n_hits)++;
    }

  return inv_sz_freq;
}




int main(int argc, char *argv[]) 
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String_OrDefault_Doc
    (ROOT, "", "If specified, all paths start from here.");
  CommandArgument_String_OrDefault_Doc
    (HEAD_REF, "", "Looks for HEAD_REF.fastb.");
  CommandArgument_String_OrDefault_Doc
    (UNIBASES, "", "Looks for UNIBASES as a reference.");
  CommandArgument_UnsignedInt_OrDefault_Doc
    (UNIBASES_K, 96, "Needed by involution algorithm.");
  CommandArgument_String_Doc
    (HEAD_READS, "Looks for HEAD_READS.{fastb,qualb,pairs/pairto}.");
  CommandArgument_Bool_OrDefault_Doc
    (WRITE, True, "Write distributions into HEAD_READS.<lib_name>.distrib files.");
  CommandArgument_String_OrDefault_Doc
    (OUT_SUFFIX,"", "If specified, output to HEAD_READS.OUT_SUFFIX.... instead.");

  CommandArgument_Bool_OrDefault_Doc
    (BUILD_LOOKUP_TABLE, False, "Build lookup table.");
  CommandArgument_Bool_OrDefault_Doc
    (KEEP_LOOKUP_TABLE, False, "Keep lookup table if generated internally.");

  CommandArgument_Bool_Doc
    (FLIP, "True: high quality at read end; False: high quality at read start.");
  CommandArgument_Int_OrDefault_Doc
    (TRIM, 0, "Trim that many bases from the beginning of a read (just before alignment).");



  CommandArgument_UnsignedInt_OrDefault_Doc
    (TARGET_SAMPLE_SIZE, 50000, "Sample size per library.");
  CommandArgument_UnsignedInt_OrDefault_Doc
    (MIN_SAMPLE_SIZE, 1000, "Minimum sample size per library.");
  CommandArgument_Double_OrDefault_Doc
    (TARGET_UNIBASE_COVERAGE, .2, "Unibases covering this fraction of total are enough.");
  CommandArgument_UnsignedInt_OrDefault_Doc
    (MAX_SIZE, 100000, "Maximum aligned pair invariant size considered.");
  CommandArgument_UnsignedInt_OrDefault_Doc
    (MIN_UNIBASE_LENGTH, 0, "Only include unibases of this length or longer for alignments.");
  CommandArgument_Int_OrDefault_Doc
    (RANDOM_SEED, 133333, "Seed value for the random generator.");

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");

  CommandArgument_Bool_OrDefault(VERBOSE, False);

  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  cout << Tag() << "Using " << NUM_THREADS << " threads." << endl;
  
  ForceAssertGe(TARGET_SAMPLE_SIZE, MIN_SAMPLE_SIZE);
  ForceAssertGt(TARGET_SAMPLE_SIZE, 0ul);

  ForceAssert(!HEAD_REF.empty() || !UNIBASES.empty());
  if (!UNIBASES.empty())
    ForceAssertNe(UNIBASES_K, 0u);

  if ( ROOT != "" )
  {    if ( HEAD_REF.nonempty( ) ) HEAD_REF = ROOT + "/" + HEAD_REF;
       if ( UNIBASES.nonempty( ) ) UNIBASES = ROOT + "/" + UNIBASES;
       HEAD_READS = ROOT + "/" + HEAD_READS;    }

  const size_t K_lookup = 12;
  
  if (UNIBASES.empty() && HEAD_REF.empty())
    FatalErr("ERROR: You must specify either UNIBASES or HEAD_REF.");

  if (!UNIBASES.empty() && !HEAD_REF.empty())
    FatalErr("ERROR: You can't specify both UNIBASES and HEAD_REF.");
  
  const String fastb_fn = HEAD_READS + ".fastb";
  const String qualb_fn = HEAD_READS + ".qualb";
  const String pairs_fn = HEAD_READS + ".pairs";
  const String pairto_fn = HEAD_READS + ".pairto";

  String ref_head;
  
  if (!HEAD_REF.empty()) {      // ---- Use a reference file
    cout << Tag() << "Using reference provided." << endl;
    ref_head = HEAD_REF;
  }
  else {       // ---- Use unibases as reference
    cout << Tag() << "Using unibases as reference." << endl;
    ref_head = UNIBASES + ".SPRD";
    build_unibase_reference(UNIBASES, UNIBASES_K, ref_head,
                            TARGET_UNIBASE_COVERAGE, MIN_UNIBASE_LENGTH, BUILD_LOOKUP_TABLE);
  }

  const String ref_fn = ref_head + ".fastb";
  const String lookup_fn = ref_head + ".lookuptab";


  const String head_out = (OUT_SUFFIX == "") ? HEAD_READS : HEAD_READS + "." + OUT_SUFFIX;
  const String dir_ascii = head_out + ".libs.ascii";
  Mkdir777(dir_ascii);


  // ---- Build lookup table if needed
  bool built_new_lookup = false;
  if (BUILD_LOOKUP_TABLE) {
    cout << Tag() << "Creating lookup table." << endl;
    LookupTabBuilder lookup(K_lookup);
    lookup.addFastb(ref_fn.c_str());
    cout << Tag() << "Writing lookup table to '" << lookup_fn << "'." << endl;
    lookup.write(lookup_fn.c_str());
    built_new_lookup = true;
  } else if (!IsRegularFile(lookup_fn))
    FatalErr("Unable to find lookup table: " + lookup_fn);


  // ---- Load reference
  cout << Tag() << "Loading reference '" << ref_fn << "'." << endl;
  const BaseVecVec ref(ref_fn);

  // ---- Load lookup table
  cout << Tag() << "Loading lookup table '" << lookup_fn << "'." << endl;
  LookupTab lookup_tab(lookup_fn.c_str());

  // ---- Set alignment filtering parameters
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


  // ---- Load pairing info
  PairsManager pairs;
  {
    if (IsRegularFile(pairs_fn)) {
      cout << Tag() << "Loading pair info '" << pairs_fn << "'." << endl;
      pairs.Read(pairs_fn);
    }
    else if (IsRegularFile(pairto_fn)) {
      cout << Tag() << "Loading pair info '" << pairto_fn << "'." << endl;
      const size_t n_reads = MastervecFileObjectCount(fastb_fn);
      pairs.ReadFromPairtoFile(pairto_fn, n_reads);
    }
    else {
      FatalErr("ERROR: Can't find pairs file.");
    }
  }
  const size_t n_pairs = pairs.nPairs();
  const size_t n_libs = pairs.nLibraries();
  cout << Tag() << "Found " << n_pairs << " pairs in " << n_libs << " libraries." << endl;

  
  // ---- Shuffling pair indices
  cout << Tag() << "Selecting pairs in each library." << endl;
  vec< vec<size_t> > i_pairs_lib(n_libs);
  for (size_t i_pair = 0; i_pair < n_pairs; i_pair++)
    i_pairs_lib[pairs.libraryID(i_pair)].push_back(i_pair);

  srand(RANDOM_SEED);
  for (size_t i_lib = 0; i_lib < n_libs; i_lib++)
    random_shuffle(i_pairs_lib[i_lib].begin(), i_pairs_lib[i_lib].end());


  // ---- Load reads and quals
  cout << Tag() << "Reading reads from '" << fastb_fn << "'." << endl;
  BaseVecVec reads(fastb_fn);
  const size_t n_reads = reads.size();
  
  QualNibbleVecVec quals;
  if (IsRegularFile(qualb_fn)) {
    cout << Tag() << "Reading quals from '" << qualb_fn << "'." << endl;
    LoadQualNibbleVec(qualb_fn, &quals);
    ForceAssertEq(reads.size(), quals.size());
    for (size_t i = 0; i != n_reads; i++)
      ForceAssertEq(reads[i].size(), quals[i].size());
  }
  else {
    cout << Tag() << "Not using quals.  '" << qualb_fn << "' not found." << endl;
  }


  // ---- Compute density of possible hits
  const size_t read_sz = (reads[0].size() + reads[1].size()) / 2;
  const IntFrequencies n_inserts = n_inserts_in_unibases(ref, read_sz, MAX_SIZE);
  n_inserts.to_text_file(ref_head + ".SPRD.sizes");

  const int inv_sz_max = n_inserts.x_max(); // the maximum possible invariant size


  vec<IntDistribution> inv_sz_smooth_dist(n_libs);
  vec<IntFunction<double> > cov_phys(n_libs);
  vec<IntFunction<double> > n_links(n_libs);

  // ---- Computing invariant sizes for libraries
  cout << Tag() << "Computing invariant size distributions for " << n_libs << " libraries." << endl;
  cout << Tag() << setw(16) << inv_sz_max << " max unibase size." << endl;
  cout << Tag() << setw(16) << TARGET_SAMPLE_SIZE << " target sample size for each lib." << endl;;
  cout << Tag() << endl;

  cout << Tag() << endl;
  cout << Tag() << " id         library name   inv size orig        pairs   sampled   hits !logic    inv size < 0   frac    inv size > 0   frac" << endl;
  cout << Tag() << " -- -------------------- --------------- ------------ --------- ------ ------ --------------- ------ --------------- ------" << endl;


  for (size_t i_lib = 0; i_lib < n_libs; i_lib++) {
    const String head = (dir_ascii + "/lib" + ToString(i_lib) + ".inv_sz");
    const size_t n_pairs_lib = i_pairs_lib[i_lib].size();
    const size_t i0_pair = i_pairs_lib[i_lib][0];
    const size_t id1 = pairs.ID1(i0_pair);
    const size_t id2 = pairs.ID2(i0_pair);
    const int inv_sz_mean = (FLIP) ? 
      -pairs.getLibrarySep(i_lib) : 
      pairs.getLibrarySep(i_lib) + reads[id1].size() + reads[id2].size();
    const int inv_sz_sd   = pairs.getLibrarySD(i_lib);

    bool short_unibases = false;
    bool enough_pairs = true;
    bool enough_samples = true;

    size_t n_ilogical = 0;
    size_t n_pairs_sampled = 0;
    size_t n_hits = 0;
    
    float pct_neg = 0; 

    if (inv_sz_mean > 0.8 * inv_sz_max) {
      short_unibases = true;
      inv_sz_smooth_dist[i_lib] = IntDistribution::gaussian(inv_sz_mean, inv_sz_sd, 4.0);
      cov_phys[i_lib] = IntFunction<double>(0, 1, 0);
    }
    else if (n_pairs_lib < MIN_SAMPLE_SIZE) {
      enough_pairs = false;
      inv_sz_smooth_dist[i_lib] = IntDistribution::gaussian(inv_sz_mean, inv_sz_sd, 4.0);
      cov_phys[i_lib] = IntFunction<double>(0, 1, 0);
    }
    else {
      
      // ---- compute alignment hits frequencies
      
      const IntFrequencies inv_sz_freq = 
        invariant_sizes_from_aligments(pairs, reads, quals, 
                                       ref, lfinder, K_lookup,
                                       i_lib, i_pairs_lib[i_lib], 
                                       MAX_SIZE, TARGET_SAMPLE_SIZE,
				       & n_pairs_sampled, & n_hits, & n_ilogical,
                                       FLIP, TRIM, NUM_THREADS);

      const float frac_sampled = float(n_pairs_sampled) / float(n_pairs_lib);

      // ---- find the lower and upper limits of the distribution
      
      int x_min = 1;
      int x_max = 0;

      if (n_hits > 0)
        distribution_trim_limits(inv_sz_freq, &x_min, &x_max);


      if (n_hits < MIN_SAMPLE_SIZE || x_max < x_min) {
        enough_samples = false;
        inv_sz_smooth_dist[i_lib] = IntDistribution::gaussian(inv_sz_mean, inv_sz_sd, 4.0);
        cov_phys[i_lib] = IntFunction<double>(0, 1, 0);
      }
      else {

        // ---- compute distribution from hits frequency and unibase weigths

        const IntDistribution inv_sz_dist = distribution_compute(inv_sz_freq, n_inserts,
                                                                 x_min, x_max);
        
        // ---- smooth out probability distribution

        inv_sz_smooth_dist[i_lib] = distribution_smooth(inv_sz_dist, inv_sz_freq);

        pct_neg = 0.5 + 100.0 * inv_sz_dist.prob_lt(0);
        
        
        // ---- coverage function

        cov_phys[i_lib] = physical_coverage_compute(inv_sz_freq, n_inserts, frac_sampled,
                                                    x_min, x_max);

        n_links[i_lib] = number_links_compute(inv_sz_freq, n_inserts, frac_sampled,
                                              x_min, x_max);
        
        // ---- output invariant size frequencies and 'raw' distributions

        inv_sz_freq.to_text_file(head);
        inv_sz_dist.to_text_file(head + ".raw");

      }

    }


    // ---- print report line
      
    cout << Tag();
    cout << " "     << setw(2)  << i_lib;
    cout << " "     << setw(20) << pairs.getLibraryName(i_lib);
    cout << " "     << setw(6)  << inv_sz_mean;
    cout << " +/- " << setw(4)  << inv_sz_sd;
    cout << " "     << setw(12) << n_pairs_lib;

    if (short_unibases) {
      cout << " **** Mean inv_sz > 0.8 x max unibase size => assume gaussian.";
    }
    else if (! enough_pairs) {
      cout << " **** not enough pairs to sample (< " << MIN_SAMPLE_SIZE << ") => assume gaussian.";
    }
    else {
      
      cout << " " << setw(9) << n_pairs_sampled;
      cout << " " << setw(6) << n_hits;
      const int ilogical_pct = 0.5 + float(100 * n_ilogical)/float(n_pairs_sampled + n_ilogical);
      cout << " " << setw(5) << ilogical_pct << "%";
   
      if (! enough_samples) {
        cout << " **** Poorly sampled distribution => assume gaussian.";
      }
      else {

        // ---- split negative and positive distributions
        
        IntDistribution dist_neg;
        IntDistribution dist_pos;
        inv_sz_smooth_dist[i_lib].split(&dist_neg, &dist_pos);
        
        if (dist_neg) cout << " " << setw(6) << int(dist_neg.mean()) 
                           << " +/- " << setw(4) << int(sqrt(dist_neg.variance()));
        else          cout << " " << setw(6) << "-" 
                           << "     " << setw(4) << "-";
        cout << " " << setw(5) << fixed << setprecision(1) << pct_neg << "%";

        if (dist_pos) cout << " " << setw(6) << int(dist_pos.mean()) 
                           << " +/- " << setw(4) << int(sqrt(dist_pos.variance()));
        else          cout << " " << setw(6) << "-" 
                           << "     " << setw(4) << "-";
        cout << " " << setw(5) << fixed << setprecision(1) << 100.0 - pct_neg << "%";
      }
    }
    cout << endl;
        
    // ---- output 'smooth' distribution

    inv_sz_smooth_dist[i_lib].to_text_file(head + ".smooth");
  }
  cout << Tag() << endl;


  // ---- output 'smooth' distributions in binary format

  if (true) {

   const String fn = head_out + ".distribs";
    cout << Tag() << "Writing binary distributions to '" << fn << "'." << endl;
    BinaryWriter::writeFile(fn.c_str(), inv_sz_smooth_dist);
    
    // To read:
    //   vec<IntDistribution> dists;
    //   BinaryReader::readFile(fn.c_str(), &dists);
  }



  // ---- PERFSTAT BLOCK for distribution


  // ---- PERFSTAT BLOCK for coverage

  IntFunction<double> cov_phys_total;
  IntFunction<double> n_links_total;
  size_t n_pairs_lib_total = 0;
  for (size_t i_lib = 0; i_lib < n_libs; i_lib++) {
    cov_phys_total    += cov_phys[i_lib];
    n_links_total     += n_links[i_lib];
    n_pairs_lib_total += i_pairs_lib[i_lib].size();
  }




  PerfStat::log() << PerfStatBlockStart("Libraries statistics tables");








  cout << endl;
  cout << "Table 1: library names, number of pairs (N), original (L0) and new sizes (L)" << endl;
  cout << endl;
  
  //      "........10........20........30........40........50........60........70........80"
  cout << "--------------------------------------------------------------------------" << endl;
  cout << " id          library name  num pairs N    orig size L0       new size L" << endl;
  cout << "--- --------------------- ------------ ----------------- -----------------" << endl;
  for (size_t i_lib = 0; i_lib < n_libs; i_lib++) {
    const size_t n_pairs_lib = i_pairs_lib[i_lib].size();
    const size_t i0_pair = i_pairs_lib[i_lib][0];
    const size_t id1 = pairs.ID1(i0_pair);
    const size_t id2 = pairs.ID2(i0_pair);
    const int inv_sz_mean = (FLIP) ? 
      pairs.getLibrarySep(i_lib) : 
      pairs.getLibrarySep(i_lib) + reads[id1].size() + reads[id2].size();
    const int inv_sz_sd   = pairs.getLibrarySD(i_lib);

    const IntDistribution & distrib = inv_sz_smooth_dist[i_lib];
    const IntFunction<double> & cov = cov_phys[i_lib];

    cout << " "     << setw(2)  << i_lib;
    cout << " "     << setw(21) << pairs.getLibraryName(i_lib);
    cout << " "     << setw(12) << n_pairs_lib;

    cout << " "     << setw(7)  << inv_sz_mean;
    cout << " +/- " << setw(5)  << inv_sz_sd;

    cout << " "     << setw(7)  << int(distrib.mean());
    cout << " +/- " << setw(5)  << int(sqrt(distrib.variance()));
    cout << endl;
  }
  if (n_libs > 1) {
    cout << endl;
    cout << "tot                 total";
    cout << " "     << setw(12) << n_pairs_lib_total << endl;
  }
  cout << "--------------------------------------------------------------------------" << endl;
  cout << endl << endl;









  cout << "Table 2: fraction of reads in each length interval" << endl;
  cout << endl;

  //      "........10........20........30........40........50........60........70........80"
  cout << "---------------------------------------------------------------------------" << endl;
  cout << " id   <L>    L < 0    0-500  500-1k   1k-2k   2k-4k   4k-8k  8k-16k    >16k" << endl;
  cout << "--- -----  -------  ------- ------- ------- ------- ------- ------- -------" << endl;
  for (size_t i_lib = 0; i_lib < n_libs; i_lib++) {
    const size_t n_pairs_lib = i_pairs_lib[i_lib].size();
    const size_t i0_pair = i_pairs_lib[i_lib][0];
    const size_t id1 = pairs.ID1(i0_pair);
    const size_t id2 = pairs.ID2(i0_pair);
    const int inv_sz_mean = (FLIP) ? 
      pairs.getLibrarySep(i_lib) : 
      pairs.getLibrarySep(i_lib) + reads[id1].size() + reads[id2].size();
    const int inv_sz_sd   = pairs.getLibrarySD(i_lib);

    const IntDistribution & distrib = inv_sz_smooth_dist[i_lib];

    cout << setw(3)  << i_lib;
    cout << " "     << setw(5)  << int(distrib.mean());

    cout << "  " << setw(7) << pct_blank_str(distrib.prob_lt(0));
    cout << "  " << setw(7) << pct_blank_str(distrib.prob_in(    0,   499));
    cout << " "  << setw(7) << pct_blank_str(distrib.prob_in(  500,   999));
    cout << " "  << setw(7) << pct_blank_str(distrib.prob_in( 1000,  1999));
    cout << " "  << setw(7) << pct_blank_str(distrib.prob_in( 2000,  3999));
    cout << " "  << setw(7) << pct_blank_str(distrib.prob_in( 4000,  7999));
    cout << " "  << setw(7) << pct_blank_str(distrib.prob_in( 8000, 15999));
    cout << " "  << setw(7) << pct_blank_str(distrib.prob_ge(16000));
    
    cout << endl;
  }
  cout << "---------------------------------------------------------------------------" << endl;







  if (false) { 
    cout << endl << endl;
    
    cout << "Table 3: physical coverage (C) in each length interval" << endl;
    cout << endl;

    //      "........10........20........30........40........50........60........70........80"
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << " id  <L>    L < 0    0-500   500-2k  2k-8k   8k-32k  > 32k    total    G=NL/C" << endl;
    cout << "--- -----  -------  ------- ------- ------- ------- -------  -------  -------" << endl;
    for (size_t i_lib = 0; i_lib < n_libs; i_lib++) {
      const size_t n_pairs_lib = i_pairs_lib[i_lib].size();
      const size_t i0_pair = i_pairs_lib[i_lib][0];
      const size_t id1 = pairs.ID1(i0_pair);
      const size_t id2 = pairs.ID2(i0_pair);
      const int inv_sz_mean = (FLIP) ? 
        pairs.getLibrarySep(i_lib) : 
        pairs.getLibrarySep(i_lib) + reads[id1].size() + reads[id2].size();
      const int inv_sz_sd   = pairs.getLibrarySD(i_lib);

      const IntFunction<double> & cov = cov_phys[i_lib];
      const IntDistribution & distrib = inv_sz_smooth_dist[i_lib];
    
      const double cov_total = cov.sum();
      const size_t G = (cov_total > 0 ? float(n_pairs_lib * distrib.mean()) / cov_total : 0);


      cout << setw(3)  << i_lib;
      cout << " "  << setw(5) << int(distrib.mean());
      print_cov_line(cov);
      cout << "  " << setw(7) << (G > 0 ? genome_size_str(G) : "n/a  ");
      cout << endl;
    }
    if (n_libs > 1) {
      cout << endl;
      cout << "tot      ";
      print_cov_line(cov_phys_total);    
      cout << endl;
    }
    cout << "-----------------------------------------------------------------------------" << endl;


    cout << endl;

  }
  



  cout << endl << endl;
  
  cout << "Table 3: number of bridging links over a specific gap size" << endl;
  cout << endl;

  //      "........10........20........30........40........50........60........70........80"

  cout << "--------------------------------------------------------------------" << endl;
  cout << " id   <L> <= 0     0    1k    2k    3k    4k    6k    8k   12k   16k" << endl;
  cout << "--- ----- ---- ----- ----- ----- ----- ----- ----- ----- ----- -----" << endl;

  vec<size_t> gaps(9);
  gaps[0] =     0;  gaps[1] =  1000;  gaps[2] =  2000;
  gaps[3] =  3000;  gaps[4] =  4000;  gaps[5] =  6000;  
  gaps[6] =  8000;  gaps[7] = 12000;  gaps[8] = 16000;
  for (size_t i_lib = 0; i_lib < n_libs; i_lib++) {
    const size_t n_pairs_lib = i_pairs_lib[i_lib].size();
    const size_t i0_pair = i_pairs_lib[i_lib][0];
    const size_t id1 = pairs.ID1(i0_pair);
    const size_t id2 = pairs.ID2(i0_pair);
    const int inv_sz_mean = (FLIP) ? 
      pairs.getLibrarySep(i_lib) : 
      pairs.getLibrarySep(i_lib) + reads[id1].size() + reads[id2].size();
    const int inv_sz_sd   = pairs.getLibrarySD(i_lib);

    const IntDistribution & distrib = inv_sz_smooth_dist[i_lib];
    const IntFunction<double> & cov = cov_phys[i_lib];
    const IntFunction<double> & nlinks = n_links[i_lib];
    
    //cout <<         setw(15) << pairs.getLibraryName(i_lib);
    cout << setw(3) << i_lib;
    cout << " "  << setw( 5) << int(distrib.mean());
    cout << " "  << setw( 4) << pct_blank_str(distrib.prob_lt(0), 0);
  
    print_n_gap_links_line(cov, nlinks, gaps);
    cout << endl;
  }
  cout << "tot           ";
  print_n_gap_links_line(cov_phys_total, n_links_total, gaps);
  cout << endl;
  cout << "--------------------------------------------------------------------" << endl;

  PerfStat::log() << PerfStatBlockStop();



  // ---- clean up

  
  if (built_new_lookup && !KEEP_LOOKUP_TABLE)
    Remove(lookup_fn);
  cout << Tag() << "Done!" << endl;

} // main() ends here.








