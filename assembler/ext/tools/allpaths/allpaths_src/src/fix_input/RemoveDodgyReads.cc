///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* RemoveDodgyReads
 *
 *
 * This module identifies 'dodgy' reads (that's a Britishism) in an input
 * dataset and removes them.  It is designed to be run at the very beginning
 * of the pipeline - in particular, it should be before error correction, since
 * these problematic reads can damage the error-correction algorithm.
 *
 *
 * 1. POLY-A'S
 * Controllable via REMOVE_POLY_A
 *
 * Poly-A's are reads of which more than 90% of the bases are A (reading forward
 * only).   These reads are suspect because an Illumina machine sometimes
 * reports large numbers of A's when it has a speck of dust in it, preventing
 * proper reading.
 *
 *
 * 2. DUPLICATE MOLECULES
 * Controllable via REMOVE_DUPLICATES
 *
 * Sometimes an Illumina machine reads the same insert repeatedly and creates
 * several copies of the same read or read pair, which are identical except for
 * sequencing errors.  We identify these duplicate molecules and keep only the one 
 * with the highest quality score.
 *
 *
 * As of Feb. 2010, this module removes about 5% of reads from E.coli datasets
 * (taking 5 minutes to run) and 15% from Microstick datasets (taking 5 minutes
 * per 10MB.)
 *
 *
 *
 * Josh Burton    2010-02 original author
 * Filipe Ribeiro 2010-10 remove_duplicate_reads()
 * 
 ******************************************************************************/

// MakeDepend: library PTHREAD
#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "util/ReadTracker.h"



# define K_DODGY 16


static inline 
String Tag(String S = "RDR") { return Date() + " (" + S + "): "; } 



size_t find_poly_A(const Bool RC, 
                   const double MAX_FRAC_A, 
                   const BaseVecVec & bases, 
                   vec<Bool> * p_is_poly_A)
{
  size_t n_poly_A = 0;
  const size_t n_reads = bases.size();
  const char poly_base = RC ? 'T' : 'A';
  const float max_pct_A = 100.0 * MAX_FRAC_A; // percent
  cout << Tag() << "Cheching reads for poly-A" << endl;
    
  pair< float, unsigned char > homopol_pct;
  for (size_t i = 0; i < n_reads; i++) {
    homopol_pct = bases[i].HomopolPercent();
    if (homopol_pct.first >= max_pct_A &&
        as_base(homopol_pct.second) == poly_base) {
      //PRINT3(i, homopol_pct.first, as_base(homopol_pct.second));
      
      n_poly_A++;
      (*p_is_poly_A)[i] = true;
    }
  }
  
  cout << Tag() << setw(12) << n_poly_A << " poly-" << poly_base << " reads found." << (RC ? "  (Used poly-T because RC=True)" : "") << endl;
  return n_poly_A;
}











// -------------------------------------------------------------------------
// PairKey class
//
// Stores a pair id (m_pid) and an associated  RC-invariant key (m_key).
//
// The 'm_key' is built from the first 16 bases (highest quality) of 
// both reads in the pair.
//
// Variables of this type are put in an array and sorted according to the 
// 'm_key' (see operator<()).  If two pairs have the same 'm_key' they are 
// assumed to be the same (yes, this is a strech, but hey).
// -------------------------------------------------------------------------

class PairKey
{
private:
  uint64_t m_pid;
  uint64_t m_key[2];
  
  static const size_t n_bases = K_DODGY;
  
public:

  // ---- Constructor ----
  PairKey(const size_t pid,
          const BaseVec & bv0,
          const BaseVec & bv1,
	  const size_t step,
          const bool from_end = false)
    : m_pid(pid)
  {
    uint64_t key0_fw = 0;
    uint64_t key0_rc = 0;
    uint64_t key1_fw = 0;
    uint64_t key1_rc = 0;
    

    if (from_end) {
      size_t b0_hi_id = bv0.size() - 1;
      size_t b1_hi_id = bv1.size() - 1;
      
      for (size_t i = 0; i != n_bases; i++) {
        const size_t b_id = step * i;
        const size_t b0 = bv0[b0_hi_id - b_id];
        const size_t b1 = bv1[b1_hi_id - b_id];

        const size_t shift = i << 1;  // i * 2

        key0_fw |= b0 << shift;
        key1_fw |= (3ul - b1) << shift;

        key0_rc |= b1 << shift;
        key1_rc |= (3ul - b0) << shift;
      }
    }
    else {
      for (size_t i = 0; i != n_bases; i++) {
        const size_t b_id = step * i;
        const size_t b0 = bv0[b_id];
        const size_t b1 = bv1[b_id];

        const size_t shift = i << 1;  // i * 2

        key0_fw |= b0 << shift;
        key1_fw |= (3ul - b1) << shift;

        key0_rc |= b1 << shift;
        key1_rc |= (3ul - b0) << shift;
      }
    }

    // keep canonical version
    if (key0_fw < key0_rc || 
        (key0_fw == key0_rc && key1_fw < key1_rc)) { 
      m_key[0] = key0_fw;
      m_key[1] = key1_fw;
    }
    else {
      m_key[0] = key0_rc;
      m_key[1] = key1_rc;
    }
  }
    
    
  uint64_t get_pair_id() const { return m_pid; }

  
  friend bool operator==(const PairKey & a, const PairKey & b)
  { 
    return (a.m_key[0] == b.m_key[0] && 
            a.m_key[1] == b.m_key[1]); 
  }

  friend bool operator<(const PairKey & a, const PairKey & b)
  { 
    if (a.m_key[0] < b.m_key[0]) return true;
    if (a.m_key[0] > b.m_key[0]) return false;
    return (a.m_key[1] < b.m_key[1]);
  }


};








void find_duplicate_pairs(const PairsManager & pairs, 
                          const BaseVecVec   & bases, 
                          const QualVecVec   & quals, 
                          vec<Bool>          * p_read_is_duplicate,
                          const Bool           FROM_END,
                          unsigned             nThreads)
{
  // ---- find minimum read lenght
  const size_t n_reads = bases.size();
  const size_t n_pairs = pairs.nPairs();
  const size_t n_libs = pairs.nLibraries();

  cout << Tag() << "Finding duplicate pairs." << endl;
  cout << Tag() << setw(12) << n_reads << "  reads." << endl;
  cout << Tag() << setw(12) << 2 * n_pairs << "  paired reads." << endl;
  cout << Tag() << setw(12) << fixed << setprecision(1) 
       << (200.0 * n_pairs) / n_reads << "  % reads are paired." << endl;
    
  vec<PairKey> pair_keys;
  pair_keys.reserve(n_pairs);

  vec<size_t> n_duplicates_counts(32, 0); // count dups in exponential bins, e.g. 32 means 2**32
  size_t n_duplicates_total = 0;
  size_t n_non_duplicates_total = 0;

  typedef vec<PairKey>::iterator Itr;
  InPlaceParallelSorter<Itr,Comparator<PairKey> > ipps(nThreads);

  // ---- cycle through all the libraries
  for (size_t lib_id = 0; lib_id != n_libs; lib_id++) {
    cout << Tag() << "---- Finding duplicates in library " 
         << pairs.getLibraryNames()[lib_id] << " (" << lib_id << "/" << n_libs << ")"
         << endl; 


    // ---- find minimum read length for this library
    unsigned len_lo = 1ul << 20;
    for (size_t p_id = 0; p_id < n_pairs; p_id++) {
      if ((unsigned)pairs.libraryID(p_id) == lib_id) {
        const unsigned len1 = bases[pairs.ID1(p_id)].size();
        const unsigned len2 = bases[pairs.ID2(p_id)].size();
        
        if (len1 < len_lo) len_lo = len1;
        if (len2 < len_lo) len_lo = len2;
      }
    }

    const size_t step = 1; // this is good enough
    
    cout << Tag() << setw(12) << len_lo << " shortest read length in library." << endl; 

    const unsigned len_min = step * K_DODGY;
    cout << Tag() << setw(12) << len_min << " minimum acceptable read length." << endl; 

    // ---- find pair keys for this library
    pair_keys.clear();
    for (size_t p_id = 0; p_id < n_pairs; p_id++) {
      if ((unsigned)pairs.libraryID(p_id) == lib_id) {
        const BaseVec & bv1 = bases[pairs.ID1(p_id)];
        const BaseVec & bv2 = bases[pairs.ID2(p_id)];
        if (bv1.size() >= len_min ||
            bv2.size() >= len_min)    // only look at reads big enough
          
          pair_keys.push_back(PairKey(p_id, bv1, bv2, step, FROM_END));
      }
    }
    
    // ---- sort pair keys 
    //cout << Tag() << "Sorting pair keys." << endl;
    ipps.sort(pair_keys.begin(), pair_keys.end());
    
    
    //cout << Tag() << "Finding duplicate pairs." << endl;
    
    // ---- find duplicate pairs
    size_t n_duplicates = 0;
    size_t n_non_duplicates = 0;
    typedef vec<PairKey>::const_iterator PKIter;
    for (PKIter it0 = pair_keys.begin(), it1 = it0; it0 != pair_keys.end(); it0 = it1) {
      n_non_duplicates++;
      
      // ---- serch for consecutive pair key matches
      it1++;
      while (it1 != pair_keys.end() && *it0 == *it1)
        it1++;
      
      const size_t n_reds = it1 - it0;
      {
        size_t n = 1;
        size_t i = 0;
        while (n < n_reds) { n <<= 1; i++; }
        n_duplicates_counts[i]++;
      }
      
      // ---- found more than 1 matching pair keys 
      if (n_reds >= 2) {
        n_duplicates += n_reds - 1;
        
        PKIter it_max = it0;    // iterator of highest quality score
        size_t q_sum_max = 0;   // highest quality score
        for (PKIter it = it0; it != it1; it++) {
          const size_t p_id = it->get_pair_id();
          const QualVec & qv1 = quals[pairs.ID1(p_id)];
          const QualVec & qv2 = quals[pairs.ID2(p_id)];
          const size_t q_sum = Sum(qv1) + Sum(qv2);
          if (q_sum > q_sum_max) {
            q_sum_max = q_sum;
            it_max = it;
          }
        }
        
        // ---- label duplicates if qual sum is not the maximum
        for (PKIter it = it0; it != it1; it++) {
          if (it != it_max) {
            const size_t p_id = it->get_pair_id();
            (*p_read_is_duplicate)[pairs.ID1(p_id)] = 
              (*p_read_is_duplicate)[pairs.ID2(p_id)] = true;
          }
        }
        
        
        // ---- debug: output reads
        if (false && n_reds >= 100) {
          const size_t sz = bases[0].size();
          cout << setw(12) << n_duplicates << "  redundancies found." << endl;
          for (PKIter it = it0; it != it1; it++) 
            bases[pairs.ID1(it->get_pair_id())].PrintBases(cout,0, sz, false, sz);
          cout << endl;
          for (PKIter it = it0; it != it1; it++) 
            bases[pairs.ID2(it->get_pair_id())].PrintBases(cout,0, sz, false, sz);
        }
      }
      
    }

    cout << Tag() << setw(12) << pair_keys.size()  << "  pairs in library." << endl;
    cout << Tag() << setw(12) << n_duplicates << "  redundant pairs found." << endl;
    cout << Tag() << setw(12) << fixed << setprecision(1)  
         << (100.0 * n_duplicates) / pair_keys.size() << "  % pairs are duplicates." << endl;

    n_duplicates_total     += n_duplicates;
    n_non_duplicates_total += n_non_duplicates;

  }


  // ---- report
  {
    cout << Tag() << "Totals:" << endl;
    cout << Tag() << setw(12) << n_pairs  << "  pairs." << endl;
    cout << Tag() << setw(12) << n_duplicates_total << "  redundant pairs found." << endl;
    cout << Tag() << setw(12) << fixed << setprecision(1)  
         << (100.0 * n_duplicates_total) / n_pairs << "  % pairs are duplicates." << endl;
    
    size_t i_max = 0;
    for (size_t i = 0; i != n_duplicates_counts.size(); i++) 
      if (n_duplicates_counts[i] != 0) 
        i_max = i + 1;

    for (size_t i = 0; i != i_max; i++) 
      cout << setw(12) << (1ul << i)  << " " << n_duplicates_counts[i] << endl;
    
    cout << endl;
    
  }
    
}

































int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String_Doc(IN_HEAD,
    "Input is in <RUN>/<IN_HEAD>.{fastb,qualb,pairs}");
  CommandArgument_String_Doc(OUT_HEAD,
    "Output goes to <RUN>/<OUT_HEAD>.{fastb,qualb,pairs}");

  // Arguments to do read-removing steps.
  CommandArgument_Bool_OrDefault_Doc(REMOVE_POLY_A, True,
    "Remove reads with over 90% As; these are an artifact of Illumina sequencing failures.");
  CommandArgument_Bool_OrDefault_Doc(REMOVE_DUPLICATES, True,
    "Determine if multiple reads appear to come from the same molecule. "
    "If so, merge their quality scores into one, and remove the others.");
  
  // Parallelization.

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  
  // POLY_A arguments.
  CommandArgument_Double_OrDefault_Doc(MAX_FRAC_A, 0.9,
    "In REMOVE_POLY_A, the maximum allowable fraction of A's present in pairs of reads."); 
  CommandArgument_Bool_OrDefault_Doc(RC, False,
     "True if the reads have been RC'ed."); 
 
  
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);  
  
  // Load reads and pairing info.
  const String fn_bases = IN_HEAD + ".fastb";
  cout << Tag() << "Loading bases from '" << fn_bases << "'."<< endl;
  BaseVecVec bases(fn_bases);
  const size_t n_reads = bases.size();
  

  // Check to see if any reads are Poly-A.
  // A read is defined as "Poly-A" if more than MAX_FRAC_A (default: 90%) of its
  // bases are A's (or T's, if RC=True).
  // These reads are created erroneously in Illumina when a speck of dust gets
  // in the machine and prevents proper reading of the bases.

  vec<Bool> is_poly_A(n_reads, false);
  const size_t n_poly_A = (REMOVE_POLY_A) ? 
    find_poly_A(RC, MAX_FRAC_A, bases, &is_poly_A) : 0;


  // Load quals.
  cout << Tag() << "Loading quals" << endl;
  QualVecVec quals(IN_HEAD + ".qualb");
  ForceAssertEq(n_reads, quals.size());


  // Load pairing information
  cout << Tag() << "Loading pairings" << endl;
  PairsManager pairs(IN_HEAD + ".pairs");
  pairs.makeCache(); // this prevents segfaults in parallelized processing later
  
  if (true) { // sanity check
    size_t n_pairs = pairs.nPairs();
    for (size_t p_id = 0; p_id < n_pairs; p_id++) {
      uint64_t r1_id = pairs.ID1(p_id);
      ForceAssertLt(r1_id, n_reads);
      uint64_t r2_id = pairs.ID2(p_id);
      ForceAssertLt(r2_id, n_reads);
    }
  }
  


  
  vec<Bool> is_duplicate(n_reads, false);

  if (REMOVE_DUPLICATES) {
    find_duplicate_pairs(pairs, bases, quals, &is_duplicate, RC, NUM_THREADS);
  }



  
  
  // Merge the various problem flags into a single 'dodgy' flag.
  cout << Tag() << "Synthesizing problems" << endl;
  vec<Bool> dodgy(n_reads, false);
  for (size_t i = 0; i < n_reads; i++) {
    const BaseVec & bv = bases[i];
    if (is_poly_A[i] || 
        is_duplicate[i] || 
        bv.size() < K_DODGY)
      dodgy[i] = true;
  }
  size_t n_dodgy_initial = Sum(dodgy);
  
  // For all reads marked as dodgy, mark their partners as well.
  cout << Tag() << "Marking partners of dodgy reads as dodgy" << endl;
  for (size_t i = 0; i < pairs.nPairs(); i++) {
    int64_t ID1 = pairs.ID1(i);
    int64_t ID2 = pairs.ID2(i);
    if ( dodgy[ID1] && !dodgy[ID2]) dodgy[ID2] = true;
    if (!dodgy[ID1] &&  dodgy[ID2]) dodgy[ID1] = true;
  }
  
  
  // Remove dodgy reads from the bases, quals, and PairsManager.
  bases.EraseIf(dodgy);
  quals.EraseIf(dodgy);
  pairs.removeReads(dodgy, false);
  
  
  
  // Save memory by cleaning up some data structures that are no longer needed.
  size_t n_duplicate = Sum(is_duplicate);
  Destroy(is_poly_A);
  Destroy(is_duplicate);
  
  
  
  // Write output files.
  if (true) {
    cout << Tag() << "Writing output files" << endl;
    // BinaryWrite3(IN_HEAD + ".is_dodgy", dodgy); // this file is never used; it's for reference only
    bases.WriteAll(OUT_HEAD + ".fastb");
    quals.WriteAll(OUT_HEAD + ".qualb");
    pairs.Write(OUT_HEAD + ".pairs");
    
    // Write a ReadTracker.
    ReadTracker rt;
    rt.AddReadSet(IN_HEAD, dodgy);
    rt.Dump(OUT_HEAD); // creates the file <OUT_HEAD>.readtrack
  }
  
  
  // Write a report to cout.
  size_t n_dodgy = Sum(dodgy);
  
  String pct_poly_A        = ToString(100.0 * n_poly_A        / n_reads, 3);
  String pct_duplicate     = ToString(100.0 * n_duplicate     / n_reads, 3);
  String pct_dodgy_initial = ToString(100.0 * n_dodgy_initial / n_reads, 3);
  String pct_dodgy         = ToString(100.0 * n_dodgy         / n_reads, 3);
  
  
  cout << endl;
  cout << Tag() << "REPORT" << endl;
  cout << "Out of " << n_reads << " reads..." << endl;
  const char poly_base = RC ? 'T' : 'A';
  cout << "\t" << n_poly_A << " (" << pct_poly_A
       << "%) were marked as poly-" << poly_base << " (REMOVE_POLY_A="
       << ToStringBool(REMOVE_POLY_A) << ")" << (RC ? " (RC=True)" : "")
       << endl;
  cout << "\t" << n_duplicate << " (" << pct_duplicate
       << "%) were marked as duplicates (REMOVE_DUPLICATES="
       << ToStringBool(REMOVE_DUPLICATES) << ")" << endl;
  cout << "In sum, " << n_dodgy_initial << " (" << pct_dodgy_initial
       << "%) were marked as dodgy for one or more reasons." << endl;
  cout << "Then these reads partners' were marked, bringing the number up to "
       << n_dodgy << " (" << pct_dodgy << "%)." << endl;
  cout << "These marked reads were removed, leaving "
       << bases.size() << " output reads." << endl;
  cout << endl;
  
  
  cout << Tag() << "Done!" << endl;
  return 0;
}
