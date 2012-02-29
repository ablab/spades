///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Filipe Ribeiro - Nov 2011 - <crdhelp@broadinstitute.org>
//

#ifndef PATHS__POLYMORPHISM_REMOVE_CORE_H
#define PATHS__POLYMORPHISM_REMOVE_CORE_H

#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/QualNibbleVec.h"
#include "PairsManager.h"

#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/KmerFreqAffixesMap.h"



typedef Kmer29                     Kmer_t;
typedef KmerFreqAffixes<Kmer_t>    KmerRec_t;






class PolymorphismStats
{
  unsigned sz_min_lo;
  unsigned sz_min_hi;
  unsigned sz_max_lo;
  unsigned sz_max_hi;
  map<pair<unsigned, unsigned>, unsigned> sz_alts;

public:
  PolymorphismStats() : sz_min_lo(0), sz_min_hi(0), sz_max_lo(0), sz_max_hi(0) {}

  void add(const unsigned sz_a, const unsigned sz_b) 
  {
    const unsigned sz_min = sz_a < sz_b ? sz_a : sz_b;
    const unsigned sz_max = sz_a < sz_b ? sz_b : sz_a;
    sz_alts[make_pair(sz_min, sz_max)]++;
    if (sz_min < sz_min_lo || sz_min_lo == 0) sz_min_lo = sz_min; 
    if (sz_min > sz_min_hi) sz_min_hi = sz_min; 
    if (sz_max < sz_max_lo || sz_max_lo == 0) sz_max_lo = sz_max; 
    if (sz_max > sz_max_hi) sz_max_hi = sz_max; 
  }


  void to_text_file(const String & fn) const 
  {
    ofstream os;
    os.open(fn.c_str());
    os << "#   sz_min     sz_max       freq" << endl; 
    for (unsigned sz_min = sz_min_lo; sz_min < sz_min_hi; sz_min++) {
      for (unsigned sz_max = sz_max_lo; sz_max < sz_max_hi; sz_max++) {

	map<pair<unsigned, unsigned>, unsigned>::const_iterator it = 
	  sz_alts.find(make_pair(sz_min, sz_max));

	if (it != sz_alts.end()) {
	  const unsigned n_alts = it->second;
	  os << " " << setw(10) << sz_min
	     << " " << setw(10) << sz_max 
	     << " " << setw(10) << n_alts  
	     << endl;
	}
      }
    }
    os << endl;
    os.close();
  }

  void print(const unsigned K) const 
  {
    const unsigned sz_snp = 2 * K - 1;
    unsigned n_indels_simple = 0;
    unsigned n_indels_complex = 0;
    unsigned n_snps = 0;
    unsigned n_mnps = 0;
    
    for (unsigned sz_min = sz_min_lo; sz_min < sz_min_hi; sz_min++) {
      for (unsigned sz_max = sz_max_lo; sz_max < sz_max_hi; sz_max++) {

	map<pair<unsigned, unsigned>, unsigned>::const_iterator it = 
	  sz_alts.find(make_pair(sz_min, sz_max));

	if (it != sz_alts.end()) {
	  const unsigned n_alts = it->second;

          if (sz_min == sz_max) {
            if (sz_min == sz_snp) n_snps += n_alts;
	    else                  n_mnps += n_alts;
          }
          else {
            if (sz_min < sz_snp)  n_indels_simple  += n_alts;
            else                  n_indels_complex += n_alts;
          }         
        }
      }
    }

    cout << "  n_indels_simple  = " << setw(10) << n_indels_simple << endl;
    cout << "  n_indels_complex = " << setw(10) << n_indels_complex << endl;
    cout << "  n_snps           = " << setw(10) << n_snps << endl;
    cout << "  n_mnps           = " << setw(10) << n_mnps << endl;
    cout << "  n_polys_total    = " << setw(10) 
         << (n_indels_simple + n_indels_complex + n_snps + n_mnps) << endl;
    
  }


}; 






// -------------------- Class KmerIPoly -------------------
//
//  Associates a kmer with a polymorphism 
//  In this context, a 'polymorphism' is just a set of sequential kmers, or a BaseVec
//  The polymorphism is identified by its index and it is stored externally 
//
template<class KMER_t> 
class KmerIPoly : public KMER_t
{
  int64_t _i_poly  : 31; // poly index
  int64_t _is_fw   :  1; 
  int64_t _ik      : 16; // the index of the kmer inside the polymorphism
  int64_t _nk      : 16; // the number of total kmers in the polymorphism
  
public:
  KmerIPoly(const unsigned K = 0) : 
    KMER_t(K), _i_poly(0), _is_fw(0), _ik(0), _nk(0) {}
  KmerIPoly(const KMER_t & kmer, 
            const unsigned i_poly, 
            const bool is_fw,
            const unsigned ik,
            const unsigned nk) : 
    KMER_t(kmer), _i_poly(i_poly), _is_fw(is_fw), _ik(ik), _nk(nk) {}
  
  size_t i_poly()       const { return _i_poly; }
  size_t n_kmers()      const { return _nk; }
  size_t i_kmer()       const { return _ik; }
  bool   is_fw()        const { return _is_fw; }
  bool   is_end_point() const { return _ik == 0 || _ik == _nk - 1; }
  
  void print() const
  {
    cout << " i_poly= " << setw(8) << i_poly() << (is_fw() ? " fw" : " rc")
         << " ik= "     << setw(3) << i_kmer() 
         << " nk= "     << setw(3) << n_kmers();
  }
};
TRIVIALLY_SERIALIZABLE(KmerIPoly<Kmer_t>);










/*

class Alleles 
{
  BaseVecVec                  _bvv;  // Alleles in BaseVec format
  vec<float>                  _kfv;  // mean kmer freq in allele
  KmerMap<KmerIPoly<Kmer_t> > _poly; // map of allele kmers to the poly index in _bvv

public:
  void add(const BaseVec & bv,
            const float kf,
            BaseVecVec * bvv_p,
            vec<float> * kfv_p,
            vec<KmerIPoly<Kmer_t> > * i_poly_vec_p);
  

// TODO later.




};
*/




// ------------------- class Polymorphisms --------------------
//
//  Stores all the polymorphisms at a specific K found in a data set
//  Stores polys as BaseVecVec
//  Stores a mapping of kmers to the corresponding polys
//


class Polymorphisms
{
  unsigned                    _K;
  BaseVecVec                  _bvv_a;  // A polymorphisms in BaseVec format
  BaseVecVec                  _bvv_b;  // B polymorphisms in BaseVec format
  vec<float>                  _kfv_a;  // mean kmer freq in A polys
  vec<float>                  _kfv_b;  // mean kmer freq in B polys
  KmerMap<KmerIPoly<Kmer_t> > _poly_a; // map of A poly kmers to the poly index in _bvv_a 
  KmerMap<KmerIPoly<Kmer_t> > _poly_b; // map of B poly kmers to the poly index in _bvv_b 

  vec<bool>                  _a_is_strong;

  PolymorphismStats          _stats;

public:
  Polymorphisms(const size_t K) : _K(K) {}
  
private:
  void _add_single(const BaseVec & bv,
                   const float kf,
                   BaseVecVec * bvv_p,
                   vec<float> * kfv_p,
                   vec<KmerIPoly<Kmer_t> > * i_poly_vec_p);

  bool _add_poly(const BaseVec & bv_a,
                 const BaseVec & bv_b,
		 const float kf_a,
		 const float kf_b,
                 vec<KmerIPoly<Kmer_t> > * poly_a_vec_p,
                 vec<KmerIPoly<Kmer_t> > * poly_b_vec_p);

  bool _kmer_frequency_valid(const size_t kf,
                             const size_t kf_min,
                             const size_t kf_max)
  {
    return (kf_min == 0 || kf >= kf_min) && (kf_max == 0 || kf <= kf_max);
  }



public:
  unsigned K() const { return _K; }
  unsigned size() const { return _bvv_a.size(); }

  KmerIPoly<Kmer_t> kmer_poly_a(const Kmer_t & kmer) const { return _poly_a(kmer); }
  KmerIPoly<Kmer_t> kmer_poly_b(const Kmer_t & kmer) const { return _poly_b(kmer); }
  const BaseVec & base_vec_a(const unsigned i_poly) const { return _bvv_a[i_poly]; }
  const BaseVec & base_vec_b(const unsigned i_poly) const { return _bvv_b[i_poly]; }
  float kmer_freq_a(const unsigned i_poly) const { return _kfv_a[i_poly]; }
  float kmer_freq_b(const unsigned i_poly) const { return _kfv_b[i_poly]; }
  
  bool a_is_strong(const size_t i_poly) const { return _a_is_strong[i_poly]; }
  bool a_is_weak(const size_t i_poly)   const { return !_a_is_strong[i_poly]; }
  bool b_is_strong(const size_t i_poly) const { return !_a_is_strong[i_poly]; }
  bool b_is_weak(const size_t i_poly)   const { return _a_is_strong[i_poly]; }


  KmerIPoly<Kmer_t> kmer_poly_strong(const Kmer_t & kmer) const 
  { 
    KmerIPoly<Kmer_t> kp_undef;

    const KmerIPoly<Kmer_t> & kp_a = _poly_a(kmer); // only defined if kmer is in poly_a
    if (kp_a.is_valid_kmer()) 
      return (a_is_strong(kp_a.i_poly())) ? kp_a : kp_undef;

    const KmerIPoly<Kmer_t> & kp_b = _poly_b(kmer); // only defined if kmer is in poly_b
    if (kp_b.is_valid_kmer()) 
      return (b_is_strong(kp_b.i_poly())) ? kp_b : kp_undef;
    
    return kp_undef;
  }

  KmerIPoly<Kmer_t> kmer_poly_weak(const Kmer_t & kmer) const 
  { 
    KmerIPoly<Kmer_t> kp_undef;

    const KmerIPoly<Kmer_t> & kp_a = _poly_a(kmer); // only defined if kmer is in poly_a
    if (kp_a.is_valid_kmer()) 
      return (a_is_weak(kp_a.i_poly())) ? kp_a : kp_undef;

    const KmerIPoly<Kmer_t> & kp_b = _poly_b(kmer); // only defined if kmer is in poly_b
    if (kp_b.is_valid_kmer()) 
      return (b_is_weak(kp_b.i_poly())) ? kp_b : kp_undef;
    
    return kp_undef;
  }

  const BaseVec & base_vec_strong(const size_t i_poly) const 
  { return a_is_strong(i_poly) ? _bvv_a[i_poly] : _bvv_b[i_poly]; }
  const BaseVec & base_vec_weak(const size_t i_poly)   const 
  { return a_is_weak(i_poly)   ? _bvv_a[i_poly] : _bvv_b[i_poly]; }


  // ---- kmap not a const because follow_kmers() needs to tag the visited.
  void build(KmerMap<KmerRec_t> & kmap,
             const size_t kf_min = 0,
             const size_t kf_max = 0);


  void write_fastas(const String & HEAD_STRONG, 
                    const String & HEAD_WEAK) const
  {
    const size_t n = size();

    ofstream os_strong((HEAD_STRONG + ".fasta").c_str());
    ofstream os_weak((HEAD_WEAK + ".fasta").c_str());
    for (size_t i = 0; i < n; i++) {
      const BaseVec & bv_strong = a_is_strong(i) ? _bvv_a[i] : _bvv_b[i];
      const BaseVec & bv_weak   = a_is_weak(i)   ? _bvv_a[i] : _bvv_b[i];
      os_strong << ">ambiguity_" << i << endl << bv_strong.ToString() << endl;
      os_weak   << ">ambiguity_" << i << endl << bv_weak.ToString() << endl;
    }
    os_strong.close();
    os_weak.close();
  }


  void write(const String & HEAD) const
  {
    _bvv_a.WriteAll(HEAD + ".A.fastb"); 
    _bvv_b.WriteAll(HEAD + ".B.fastb"); 
    _poly_a.write_binary(HEAD + ".kmer.poly.A.map");
    _poly_b.write_binary(HEAD + ".kmer.poly.B.map");
  }

  void read(const String & HEAD)
  {
    _bvv_a.ReadAll(HEAD + ".A.fastb"); 
    _bvv_b.ReadAll(HEAD + ".B.fastb"); 
    _poly_a.read_binary(HEAD + ".kmer.poly.A.map");
    _poly_b.read_binary(HEAD + ".kmer.poly.B.map");
  }

  void print_stats() const { _stats.print(_K); }
  void write_stats(const String & fn) const { _stats.to_text_file(fn); }

  void print(const unsigned i_poly) const 
  {
    cout << "a: i_poly= " << setw(6) << i_poly 
	 << " nb= " << _bvv_a[i_poly].size() 
	 << " " << hieroglyphs(_bvv_a[i_poly]) 
	 << " kf_mean= " << _kfv_a[i_poly] 
	 << endl;
    cout << "b: i_poly= " << setw(6) << i_poly 
	 << " nb= " << _bvv_b[i_poly].size() 
	 << " " << hieroglyphs(_bvv_b[i_poly]) 
	 << " kf_mean= " << _kfv_b[i_poly] 
	 << endl;
  }

};







void polymorphisms_find_parallel(const unsigned      K,
                                 const BaseVecVec  & bvv, 
                                 Polymorphisms     * polys_p,
                                 KmerSpectrum      * kspec_p,
                                 const unsigned      VERBOSITY,
                                 const unsigned      NUM_THREADS,
                                 const size_t        mean_mem_ceil = 0);

void polymorphisms_correlate(const PairsManager  & pairs, 
                             const BaseVecVec    & bvv,
                             Polymorphisms       * polys_p,
                             const unsigned        VERBOSITY);
  


size_t polymorphisms_remove_parallel(const Polymorphisms  & polys,
                                     BaseVecVec           * bvv_p,
                                     QualNibbleVecVec     * qvv_p,
                                     const unsigned         VERBOSITY,
                                     const unsigned         NUM_THREADS);





#endif
