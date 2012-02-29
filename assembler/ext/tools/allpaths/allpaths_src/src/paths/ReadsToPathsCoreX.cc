/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <map>

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "FeudalMimic.h"
#include "kmers/KmerRecord.h"
#include "kmers/SortKmers.h"
#include "pairwise_aligners/BalancedMutmerGraph.h"
#include "paths/KmerPath.h"
#include "paths/MakeAlignsPathsParallelX.h"
#include "paths/ReadsToPathsCoreX.h"

static inline 
String Tag(String S = "RTPCX") { return Date() + " (" + S + "): "; } 

// Returns true if the K-mer at BaseVec b, in the range [offset, offset+K),
// is palindromic. 

inline bool is_palindrome(const BaseVec & bv, 
                          const size_t K, 
                          const size_t offset)
{
  if (K % 2) return false; // kmers cannot be palindromic if K is odd
  
  //  ForceAssertLt(offset, bv.size());
  const size_t rev_offset = offset + K - 1;
  //  ForceAssertLt(rev_offset, bv.size());
  
  const size_t half_K = K >> 1;
  for (size_t i = 0; i < half_K; i++)
    if ((bv[offset + i] ^ bv[rev_offset - i]) != 3) 
      return false;
  return true;
}




// MutmerHit: this class takes 128 bits (16 bytes) distributed as follows:
// 
// first 64 bits:
//
// read1_id   bits  0 .. 38  (39)
// rc         bits 39        (1)
// len        bits 40 .. 47  (8) (low order bits)
// pos1       bits 48 .. 63  (16)
//
// second 64 bits:
//
// read2_id   bits  0 .. 38  (39)
// unused     bits 39        (1)
// len        bits 40 .. 47  (8) (high order bits)
// pos2       bits 48 .. 63  (16)
//
//
// NOTE: this class can store read ids up to 512 G (39 bits)
//       that should be enough for a while. 

class MutmerHit
{
private:
  union 
  {
    uint64_t _id[2];
    uint16_t _pos[8];
    uint8_t  _len[16];
  };
    
public:
  MutmerHit() { memset(reinterpret_cast<char*>(_len), 0, sizeof(MutmerHit)); }

  void set(uint64_t id1, 
           uint64_t id2, 
           bool rc,
           uint16_t pos1, 
           uint16_t pos2, 
           uint16_t len)
  {
    set_id1(id1);
    set_id2(id2);
    set_rc2(rc);
    set_pos1(pos1);
    set_pos2(pos2);
    set_len(len);
  }

  inline uint64_t get_id1()  const { return _id[0] & 0x0000007fffffffff; } // 39 bits
  inline uint64_t get_id2()  const { return _id[1] & 0x0000007fffffffff; } // 39 bits
  inline uint16_t get_pos1() const { return _pos[3]; } // 16 bits;
  inline uint16_t get_pos2() const { return _pos[7]; } // 16 bits;
  inline uint16_t get_len()  const { return _len[5] | (_len[13] << 8); } // 8 + 8 bits
  inline bool     get_rc2()  const { return _len[4] >> 7; }  // 1 bit

  inline void set_id1 (const uint64_t id1)  { _id[0] &= 0xffffff8000000000; _id[0] |= id1; } // 39 bits
  inline void set_id2 (const uint64_t id2)  { _id[1] &= 0xffffff8000000000; _id[1] |= id2; } // 39 bits
  inline void set_pos1(const uint16_t pos1) { _pos[3] = pos1; } // 16 bits;
  inline void set_pos2(const uint16_t pos2) { _pos[7] = pos2; } // 16 bits;
  inline void set_len (const uint16_t len)  
  { _len[5] = len & 0x00ff;
    _len[13] = len >> 8; }  // 8 + 8 bits;
  inline void set_rc2 (const bool rc2)  { if (rc2) _len[4] |= 0x80; else _len[4] &= 0x7f; } // 1 bit




  // Each mutmer hit (mhit) represents a kmer equivalence: it claims that
  // the bases of read #id1, starting at pos1, are equal to (or rc of)
  // the bases of read #id2, starting at pos2, for a length of len bases.
  // If these asserts fail, the bases are not actually equal.
  void AssertValid(const BaseVecVec & bases) const 
  {
    const BaseVec & bv1 = bases[get_id1()];
    const BaseVec & bv2 = bases[get_id2()];
    const size_t len  = get_len();
    const size_t pos1 = get_pos1();
  
    if (get_rc2()) {   // RC
      const size_t pos2 = bv2.size() - get_pos2() - 1;
      for (size_t i = 0; i != len; i++)
        ForceAssertEq(bv1[pos1 + i], 3 - bv2[pos2 - i]);
    }
    else {             // FW
      const size_t pos2 = get_pos2();
      for (size_t i = 0; i != len; i++)
        ForceAssertEq(bv1[pos1 + i], bv2[pos2 + i]);
    }
  }


  friend bool operator<(const MutmerHit & h1, 
                        const MutmerHit & h2) 
  {
    if (h1.get_id2() < h2.get_id2()) return true;
    if (h1.get_id2() > h2.get_id2()) return false;
    return h1.get_len() > h2.get_len();    
  }
  
  friend std::istream& operator>>(std::istream & is, MutmerHit & mhit) 
  { is.read(reinterpret_cast<char*>(&mhit), sizeof(MutmerHit)); return is; }

  friend std::ostream& operator<<(std::ostream & os, const MutmerHit & mhit)
  { os.write(reinterpret_cast<const char*>(&mhit), sizeof(MutmerHit)); return os; }

};














template<int I>
void mutmer_graph_to_mutmer_hits(const BaseVecVec & sbvv,
                                 BMG<I> & bmg,      // can't be const, need to be sorted
                                 vec<MutmerHit> * mhits,
                                 const bool VERBOSE = true)
{
  size_t n_bvs = sbvv.size();
  // please don't comment out these outputs. 
  // they help when things go wrong.
  if (VERBOSE) {
    cout << Tag() << "(MG2MH): Start." << endl;
    cout << Tag() << "(MG2MH): bvvs.size() = " << n_bvs << endl;
    cout << Tag() << "(MG2MH): bmg.size() = " << bmg.size() << endl;
    cout << Tag() << "(MG2MH): bmg.SizeOf() = " << bmg.SizeOf()/(1u<<20) << " MB" << endl;
  }

  if (bmg.size() == 0) return;

  // Generate mutmers.
  size_t n_hits = 0;
  for (size_t id1 = 0; id1 < n_bvs; id1++)
    n_hits += bmg[id1].size();

  if (VERBOSE) {
    cout << Tag() << "(MG2MH): n_hits = " << n_hits << endl;
    cout << Tag() << "(MG2MH): allocating " << ((n_hits * sizeof(MutmerHit)) >> 20) << " MB to store MutmerHits." << endl;
  }

  mhits->resize(n_hits);

  size_t lt = 0;
  size_t eq = 0;
  size_t mt = 0;
  

  size_t i_hit = 0;
  for (size_t id1 = 0; id1 < n_bvs; id1++) {

    BMG_node<I> & mmids = bmg[id1];

    // sorting mutmer_id over id and rc
    sort(mmids.begin(), mmids.end());

    const size_t n_mmids = mmids.size();
    for (size_t im0 = 0, im1 = 1; im0 < n_mmids; im0 = im1++) {
      
      // find continuous block of == mmids ('==' means same (id,rc))
      while (im1 < n_mmids && mmids[im0] == mmids[im1]) 
        im1++;

      const bool    RC = mmids[im0].Rc();
      const size_t id2 = mmids[im0].UnpackReadId();

      for (size_t im = im0; im != im1; im++) {
        MutmerHit & h = (*mhits)[i_hit++];

        uint16_t pos1, pos2, len;
        mmids[im].Unpack(pos1, pos2, len);
        if (id1 < id2) h.set(id1, id2, RC, pos1, pos2, len);
        else if (!RC)  h.set(id2, id1, RC, pos2, pos1, len);
        else {
          size_t nb1 = sbvv[id1].size();
          size_t nb2 = sbvv[id2].size();
          h.set(id2, id1, RC, nb2 - pos2 - len, nb1 - pos1 - len, len);    
        }

	if (id1 < id2) lt++;
	if (id1 == id2) eq++;
	if (id1 > id2) mt++;

      }
    }

  }
  
  // cout << "lt = " << lt << endl;
  // cout << "eq = " << eq << endl;
  // cout << "mt = " << mt << endl;


  if (VERBOSE) cout << Tag() << "(MG2MH): Done." << endl;

}






bool mutmer_hits_on_file(const String CHECKPOINT_HEAD)
{
  return IsRegularFile(CHECKPOINT_HEAD + ".mhits");
}

void mutmer_hits_to_file(const vec<MutmerHit> & mhits, const String CHECKPOINT_HEAD) 
{
  ofstream os;
  os.open((CHECKPOINT_HEAD + ".mhits").c_str());
  uint64_t n_mhits = mhits.size();
  os.write(reinterpret_cast<const char*>(&(n_mhits)), sizeof(uint64_t));
  for (size_t i = 0; i != n_mhits; i++)
    os << mhits[i];
  os.close();
}


void mutmer_hits_from_file(vec<MutmerHit> & mhits, const String CHECKPOINT_HEAD)
{
  String chkpt_fn = CHECKPOINT_HEAD + ".mhits";
 
  ifstream is;
  is.open((CHECKPOINT_HEAD + ".mhits").c_str());
  uint64_t n_mhits;
  is.read(reinterpret_cast<char*>(&n_mhits), sizeof(uint64_t));
  mhits.resize(n_mhits);
  for (size_t i = 0; i != n_mhits; i++)
    is >> mhits[i];
  is.close();
}














void base_vec_vec_to_mutmer_hits(const size_t          K, 
                                 BaseVecVec          & bvv, 
                                 BitVecVec           * chosen,
                                 vec<MutmerHit>      * mhits,
                                 const String        & PARCEL_DIR,
                                 const size_t          N_THREADS,
                                 const String          CHECKPOINT_HEAD,
                                 // please leave the default for VERBOSE at 'true'. 
                                 // when CommonPather calls this routine it uses 
                                 // a lot of memory and it is somewhat unstable.
                                 // verbosity allow faster diagnostics.
                                 const bool            VERBOSE)
{
  const size_t n_bvs = bvv.size();



  // ---- whole bunch of sanity checks

  ForceAssertLt(n_bvs, (1ul << 39));  // mutmer_hits have only 39 bits for ids



  // Determine whether this dataset will require a big BalancedMutmerGraph
  // (e.g., one with I=2 rather than I=1.)
  bool big_BMG = false;
  for (size_t id = 0; id < n_bvs; id++) {
    if (bvv[id].size() >= 1024) big_BMG = true;
    ForceAssertLe(bvv[id].size(), static_cast<unsigned>(USHRT_MAX)); // ReadsToPathsCoreY should ensure that this is true
  }
  if (n_bvs > INT_MAX) big_BMG = true;

  if (bvv.sumSizes() > static_cast<size_t>(halfway_kmer))
    FatalErr("Five byte limit on k-mer indices exceeded.");
 



  // ---- need sorting?  sorting/not sorting changes results...
  bool must_sort = true;
  size_t len0 = bvv[0].size();
  for (size_t id = 1; id < n_bvs && !must_sort; id++)
    if (bvv[id].size() != len0)
      must_sort = true;
  




  // sort BaseVec's in descending order by size, if necessary.  

  vec<size_t> id_sid;       // maps sorted indices 'sid' to unsorted 'id'
  BaseVecVec & sbvv = bvv;  // 'sbvv': sorted base vectors
    

  if (must_sort) {
    id_sid.resize(n_bvs);
    {
      vec< pair<int, size_t> > lengths(n_bvs);
      for (size_t id = 0; id < n_bvs; id++)
        lengths[id] = make_pair(bvv[id].size(), id);
      ReverseSort(lengths);
      for (size_t sid = 0; sid < n_bvs; sid++)  // running over sorted ids
        id_sid[sid] = lengths[sid].second;
    }
    
    BaseVecVec sbvv_tmp(n_bvs);
    for (size_t sid = 0; sid < n_bvs; sid++) {  // running over sorted ids
      sbvv_tmp[sid] = bvv[id_sid[sid]];  
    }
    sbvv = sbvv_tmp;  // note: 'sbvv' is just a reference to 'bvv'
  }








      
  if (0 &&    // for now, this never happens. Also, must save 'chosen'
      CHECKPOINT_HEAD != "" && 
      mutmer_hits_on_file(CHECKPOINT_HEAD)) {
    
    if (VERBOSE) cout << Tag("BVV2MH") << "CHECKPOINT: retrieving previously computed 'mhits'." << endl;
    mutmer_hits_from_file(*mhits, CHECKPOINT_HEAD);
    if (VERBOSE) cout << Tag("BVV2MH") << "CHECKPOINT: retrieved " << mhits->size() << " mhits." << endl;
    
  } 
  else {
    if (VERBOSE) cout << Tag("BVV2MH") << "Computing 'mhits' and 'chosen'." << endl;
    Mimic(sbvv, *chosen);
    
    
    // Call MakeAligns with the goal of creating a BMG (Balanced Mutmer Graph)
    // The BMG has a parameter I, indicating internal storage space
    // If the the reads are long (>1024bp), we must use I=2;
    // otherwise, we can use I=1, which is more efficient
    // This is the big runtime/memory bottleneck of pathing!
    
#    define CALL_MAKE_ALIGNS(I, __K)                                \
      {                                                             \
         BMG<I> bmg;                                                \
         MakeAlignsPathsParallelX<I>                                \
           (__K, sbvv, chosen, &bmg,                                \
            PARCEL_DIR, N_THREADS, CHECKPOINT_HEAD);                \
                                                                    \
         mutmer_graph_to_mutmer_hits<I>(sbvv, bmg, mhits, VERBOSE); \
       }
     
#    define CASE(_K)                               \
      {                                            \
        if (big_BMG)  { CALL_MAKE_ALIGNS(2, _K); } \
        else          { CALL_MAKE_ALIGNS(1, _K); } \
      }
     

    // ---- Call MakeAlignsPathsParallelX(...) and mutmer_graph_to_mutmer_hits(...)
    //      See #define's above.
    DISPATCH_ON_K(K, CASE);
    
    if (CHECKPOINT_HEAD != "") {
      //if (VERBOSE) cout << Tag("BVV2MH") << "CHECKPOINT: saving 'mhits'." << endl;
      //mutmer_hits_to_file(*mhits, CHECKPOINT_HEAD);
      //if (VERBOSE) cout << Tag("BVV2MH") << "CHECKPOINT: saved " << mhits->size() << " mhits." << endl;
    }     
  }
  





  const size_t n_mhits = mhits->size();



  
  
  // Re-order the reads, undoing the sorting by length that we did earlier.

  if (must_sort) {
    if (VERBOSE)
      cout << Tag("BVV2MH") << "Recover original read order" << endl;

    for (size_t ih = 0; ih < n_mhits; ih++) {
      (*mhits)[ih].set_id1(id_sid[(*mhits)[ih].get_id1()]);
      (*mhits)[ih].set_id2(id_sid[(*mhits)[ih].get_id2()]);    
    }
    
    // 'chosen' was computed as a function of 'sid' (sorted base vectors)
    // now we want 'chosen' as a function of 'id' (unsorted base vectors)
    BitVecVec chosen_tmp(n_bvs);
    BaseVecVec bvv_tmp(n_bvs);
    for (size_t sid = 0; sid < n_bvs; sid++) {  // running over sorted ids
      size_t id = id_sid[sid];
      chosen_tmp[id] = (*chosen)[sid];
      bvv_tmp[id]    = sbvv[sid];      // remember: sbvv is just a ref to bvv
    }

    // return 'chosen' and 'bvv' to its original order
    (*chosen) = chosen_tmp;    
    bvv = bvv_tmp;  
  }






  // Sanity check.  Make sure all mutmer hits make sense.  
  // Each mutmer hit (mhit) represents a kmer equivalence: it claims that
  // the bases of read #id1, starting at pos1, are equal to (or rc of)
  // the bases of read #id2, starting at pos2, for a length of len bases.
  // If these asserts fail, the bases are not actually equal.
  if (VERBOSE) cout << Tag("BVV2MH") << "Validate mhits" << endl;
  for (size_t ih = 0; ih != n_mhits; ih++)
    (*mhits)[ih].AssertValid(bvv);   







  
  // Sort the mhits over 'id2', and then from larger to smaller 'len'
  if (VERBOSE) cout << Tag("BVV2MH") << "Sort mhits" << endl;
  Sort(*mhits);
  

}









void mutmer_hits_to_paths(const size_t           K, 
                          const BaseVecVec     & bvv, 
                          const BitVecVec      & chosen,
                          const vec<MutmerHit> & mhits, 
                          vecKmerPath          * paths,
                          // please leave the default for VERBOSE at 'true'. 
                          // when CommonPather calls this routine it uses 
                          // a lot of memory and it is somewhat unstable.
                          // verbosity allow faster diagnostics.
                          const bool            VERBOSE = true)
{
  if (VERBOSE) cout << Tag("MH2P") << "Start." << endl;

  const size_t n_bvs = bvv.size();
  const size_t n_mhits = mhits.size();


  // Find which canonical kmers (i.e., kmers marked as chosen) are
  // palindromic.

  if (VERBOSE) cout << Tag("MH2P") << "Finding palindromes." << endl;

  BitVecVec chosen_palindrome;
  Mimic(bvv, chosen_palindrome);
  for (size_t id = 0; id < n_bvs; id++) {
    const size_t nb = bvv[id].size();
    const size_t nk = (nb >= K) ? nb - K + 1 : 0;
    for (size_t ik = 0; ik < nk; ik++)
      if (chosen[id][ik])
	if (is_palindrome(bvv[id], K, ik))
	  chosen_palindrome[id].Set(ik, true);
  }



  // Number positions on the reads.
  if (VERBOSE) cout << Tag("MH2P") << "Finding starting kmer numbers for each read." << endl;
  vec<size_t> knums0(n_bvs);
  {
    size_t knum = 0;
    for (size_t id = 0; id < n_bvs; id++) {
      knums0[id] = knum;
      knum += bvv[id].size();    
    }
    ForceAssertEq(knum, bvv.sumSizes());
  }


 



  // Build read paths for the reads themselves.
  if (VERBOSE)
    cout << Tag("MH2P") << "Building read paths from mutmer hits." << endl;
  
  paths->reserve(n_bvs);
     
  // Set up tracking of palindromes.
  size_t n_palindromes = 0;
  map<size_t, size_t> knum_palindrome;

  size_t n_segments = 0;

  // Loop over all reads, filling a 'knum' vector for each 
  // that indicates what the kmer numbers on that read will be.
  // cout << "n_bvs = " << n_bvs << endl;
  // cout << "n_mhits = " << n_mhits << endl;

  // buffer variables cleared each time
  vec<longlong> knums;
  KmerPath path;

  size_t ih0 = 0;  // starting hit index 

  for (size_t id2 = 0; id2 < n_bvs; id2++) {

    // find first mhit for which id2 is different (mhits is sorted over id2)
    size_t ih1 = ih0;  // ending hit index
    while (ih1 < n_mhits && mhits[ih1].get_id2() == id2) 
      ih1++;

    if (VERBOSE) {
      if ((n_bvs >= 1e6 && id2 % (n_bvs / 20) == 0) ||
          n_bvs < 20) {
        cout << "id2 = " << id2 << " (" << (100.0 * id2) / n_bvs << " %)  ih01 = " << ih0 << " " << ih1 << endl;
      }
    }     


    size_t nb2 = bvv[id2].size();
    size_t nk2 = (nb2 >= K) ? nb2 - K + 1 : 0;  // number of kmers in bv2

    knums.clear();
    knums.resize(nk2, -1);  // -1 means unset; checked at the end

    // Find canonical kmers on this read (there may be any number, including 0)
    // and mark them with a kmer number corresponding to their location.
    for (size_t ik = 0; ik < nk2; ik++)
      if (chosen[id2][ik]) 
        knums[ik] = halfway_kmer + knums0[id2] + ik + 1;
       

    // Look through all the mutmer hits involving this read, in search of
    // correspondences between a kmer on this read and its canonical version
    // on other reads.  For all such correspondences found, mark the location
    //  on this read with the same kmer_id used by the canonical kmer on the
    // other read.

    for (size_t ih = ih0; ih < ih1; ih++) {

      const size_t id1  = mhits[ih].get_id1();
      const size_t pos1 = mhits[ih].get_pos1();
      const size_t pos2 = mhits[ih].get_pos2();
      const size_t len  = mhits[ih].get_len();
      const bool   rc2  = mhits[ih].get_rc2();

      if (rc2) {   // RC
        for (size_t ik = 0; ik < len - K + 1; ik++) {

          const size_t ik1 = pos1 + ik;
          const size_t ik2 = nb2 - K - pos2 - ik;

          if (chosen[id1][ik1] && knums[ik2] < 0) {

            if (chosen_palindrome[id1][ik1])
              knums[ik2] = halfway_kmer + knums0[id1] + ik1 + 1;
            else
              knums[ik2] = halfway_kmer - knums0[id1] - ik1;
          }
        }
      }
      else {       // FW
        for (size_t ik = 0; ik < len - K + 1; ik++) {

          const size_t ik1 = pos1 + ik;
          const size_t ik2 = pos2 + ik;

          if (chosen[id1][ik1] && knums[ik2] < 0) 
            knums[ik2] = halfway_kmer + knums0[id1] + ik1 + 1;
        }
      }

    }




    // Check for palindromes.
    for (size_t ik = 0; ik < nk2; ik++) {
      if (is_palindrome(bvv[id2], K, ik)) {
        size_t knum = knum_palindrome[knums[ik]];
        if (knum == 0) { // not defined
          ForceAssertLt(n_palindromes, (size_t)max_palindromes);
          knum = first_palindrome + n_palindromes;
          knum_palindrome[knums[ik]] = knum;
          n_palindromes++;    
        }
        knums[ik] = knum;
      }
    }

    // Build the new KmerPath, segment by segment.
    path.clear();
    for (size_t ik = 0, jk = 1; ik < nk2; ik = jk++) {

      while (jk < nk2 && knums[jk] == knums[jk - 1] + 1) 
        jk++;

      path.AddSegment(knums[ik], knums[jk-1]);
      //path.AddSegment(knums[ik], knums[ik] + jk - ik - 1);

      n_segments++;
    }
    paths->push_back(path);

       
    // Sanity check.  If this assert fails, the knums vector was not built correctly.
    // Most likely, there are missing mutmer hits, and some base-locations are
    // not marked with kmer IDs as they should have been.
    for (size_t ik = 0; ik < nk2; ik++) {
      if (knums[ik] < 0) PRINT3_TO(cout, id2, ik, nb2);
      ForceAssert(knums[ik] >= 0);
    }


    // update the starting hit index (ih0) to the previous ending index (ih1)
    ih0 = ih1;    
  }
  
   if (VERBOSE) cout << Tag("MH2P") << "Done." << endl;
    
}













void ReadsToPathsCoreX(const size_t     K, 
                       BaseVecVec     & bvv, // not const, change in place
                       vecKmerPath    * paths, 
                       const String   & PARCEL_DIR,
                       const size_t     N_THREADS,
                       const String     CHECKPOINT_HEAD,
                       const bool       VERBOSE)
{
  if (bvv.size() > 0) {
    if (VERBOSE) cout << Tag("RTPCX") << "Start." << endl;
    ForceAssertSupportedK(K);
    
 
    if (VERBOSE) cout << Tag("RTPCX") << "Calling base_vec_vec_to_mutmer_hits()." << endl;
    vec<MutmerHit> mhits;
    BitVecVec chosen;
    base_vec_vec_to_mutmer_hits(K, bvv, & chosen, & mhits, 
                                PARCEL_DIR, N_THREADS, CHECKPOINT_HEAD, VERBOSE);
    
    
    
    
    if (VERBOSE) cout << Tag("RTPCX") << "Calling mutmer_hits_to_paths()." << endl;
    mutmer_hits_to_paths(K, bvv, chosen, mhits, paths, VERBOSE);
    


    if (VERBOSE) cout << Tag("RTPCX") << "Done." << endl;
  }
}
















// ReadsToPathsCoreX will fail if it faces the prospect of a KmerPathInterval
// with length greater than 65535.
// ReadsToPathsCoreY is a wrapper which solves this:
// it breaks long kmer paths, call ReadsToPathsCoreX, 
// and puts broken things back together.

void ReadsToPathsCoreY(BaseVecVec      & bvv, // not const, change in place
                       const size_t      K, 
                       vecKmerPath     & paths, 
                       const String    & PARCEL_DIR,
                       const size_t      N_THREADS,
                       const String      CHECKPOINT_HEAD,
                       const bool        VERBOSE)
{

  // Are there any long reads?
  // "Long": greater than USHRT_MAX, the maximum unsigned short, which is 65535.
  bool long_read = false;
  for(size_t id = 0; id < bvv.size() && !long_read; id++)
    long_read = (bvv[id].size() > USHRT_MAX);
 
  // If not, Y=X.
  if(! long_read) {
    ReadsToPathsCoreX(K, bvv, &paths, 
                      PARCEL_DIR, N_THREADS, CHECKPOINT_HEAD, VERBOSE);
  }
  else {

    // If so, we need to break up the reads.
    size_t n_reads = 0;
    size_t reads_rawsize = 0;
    vec<size_t> contig;
    BaseVecVec sheared_reads;
    BaseVec bv;

    for (int pass = 1; pass <= 2; pass++) {
      if (pass == 2) {
        sheared_reads.Reserve(reads_rawsize, n_reads);
        contig.reserve(n_reads);
      }
      for (size_t id = 0; id < bvv.size(); id++) {
        const size_t n_pos = bvv[id].size(); 
        for (size_t pos = 0; pos < n_pos; pos++) {
	  size_t len = n_pos - pos;
	  if (len > USHRT_MAX) len = USHRT_MAX;
	  
          if (pass == 1) {
            n_reads++;
            reads_rawsize += (len + 15) / 16;
          }
          else {
            bv.SetToSubOf(bvv[id], pos, len);
            sheared_reads.push_back(bv);
            contig.push_back(id);
          }
	  
          if (pos + len == n_pos) break;
          pos += len - K;
        }
      }
    }

    // Now the BaseVecVec sheared_reads holds the sheared version of bvv.
    vecKmerPath sheared_paths;
    ReadsToPathsCoreX(K, sheared_reads, & sheared_paths, 
                      PARCEL_DIR, N_THREADS, CHECKPOINT_HEAD, VERBOSE);

    // Now release the memory that held the sheared reads...
    Destroy(sheared_reads);

    // Finally, reassemble.  We have to do this gently,
    //  so that nothing ends up self-owned.
    paths.clear().reserve(bvv.size());

    KmerPath one_path;
    size_t i_shear = 0;

    for(size_t i_path = 0; i_path < bvv.size(); i_path++) {
      one_path.Clear();
      while(i_shear < contig.size() && 
             contig[i_shear] == i_path) {
        one_path.Append(sheared_paths[i_shear]);
        i_shear++;
      }
      paths.push_back(one_path);
    }

  }
}














// Wrapper functions, which also create paths_rc and/or pathsdb.

void ReadsToPathsCoreY(BaseVecVec          & bvv, // not const, change in place
                       const size_t          K,
                       vecKmerPath         & paths, 
                       vecKmerPath         & paths_rc, 
                       vec<tagged_rpint>   & pathsdb,
                       const String        & PARCEL_DIR,
                       const size_t          N_THREADS,
                       const String          CHECKPOINT_HEAD,
                       const bool            VERBOSE)
{
  ReadsToPathsCoreY(bvv, K, paths, 
                    PARCEL_DIR, N_THREADS, CHECKPOINT_HEAD, VERBOSE);
  paths_rc = paths;
  for (size_t i = 0; i < paths_rc.size(); i++)
    paths_rc[i].Reverse();
  CreateDatabase(paths, paths_rc, pathsdb);    
}
  


void ReadsToPathsCoreY(BaseVecVec            & bvv, // not const, change in place
                       const size_t            K,
                       vecKmerPath           & paths, 
                       vecKmerPath           & paths_rc, 
                       vec<big_tagged_rpint> & pathsdb,
                       const String          & PARCEL_DIR,
                       const size_t            N_THREADS,
                       const String            CHECKPOINT_HEAD,
                       const bool              VERBOSE)
{ 
  ReadsToPathsCoreY(bvv, K, paths, 
                    PARCEL_DIR, N_THREADS, CHECKPOINT_HEAD, VERBOSE);
  paths_rc = paths;
  for (size_t i = 0; i < paths_rc.size(); i++)
    paths_rc[i].Reverse();
  CreateDatabase(paths, paths_rc, pathsdb);    
}










