///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef KMERS__NAIF_KMER__KMER_FREQ_AFFIXES_MAP_H
#define KMERS__NAIF_KMER__KMER_FREQ_AFFIXES_MAP_H


#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelKmerStorer.h" 
#include "kmers/naif_kmer/KmerMap.h" 



static inline 
String TagKFAM(String S = "KFAM") { return Date() + " (" + S + "): "; } 

class FreqAffixes
{
  uint64_t _freq : 56;
  uint64_t _pre  :  4;  // 1 bit per possible base
  uint64_t _suf  :  4;  // 1 bit per possible base

  unsigned _n_bases(const unsigned bit4) const 
  { 
    return ((bit4 & 1u) + ((bit4 >> 1) & 1u) + ((bit4 >> 2) & 1u) + ((bit4 >> 3) & 1u)); 
  }

  unsigned _base(const unsigned bit4, const unsigned i) const
  {
    unsigned n = 0;
    if (bit4 & 1u) n++;   if (n > i) return 0u;
    if (bit4 & 2u) n++;   if (n > i) return 1u;
    if (bit4 & 4u) n++;   if (n > i) return 2u;
    if (bit4 & 8u) n++;   if (n > i) return 3u;
    return 4u;
  }

  String _str_affixes(const unsigned aff, const bool fw) const 
  {
    String s = "[" + (((aff & 1u) ? hieroglyph(fw ? 0 : 3) : " ") +
		      ((aff & 2u) ? hieroglyph(fw ? 1 : 2) : " ") +
		      ((aff & 4u) ? hieroglyph(fw ? 2 : 1) : " ") +
		      ((aff & 8u) ? hieroglyph(fw ? 3 : 0) : " ")) + "]";
    return s;
  }



public:
  FreqAffixes() : _freq(0), _pre(0), _suf(0) {}

  void set_freq  (const size_t freq)    { _freq = freq; }
  void set_prefix(const unsigned base)  { _pre |= (1u << (base & 3)); }
  void set_suffix(const unsigned base)  { _suf |= (1u << (base & 3)); }


  size_t freq() const { return _freq; }

  bool has_prefix(const unsigned base) const { return _pre & (1u << (base & 3)); }
  bool has_suffix(const unsigned base) const { return _suf & (1u << (base & 3)); }

  size_t n_prefixes(const bool fw = true) const 
  { return (fw ? _n_bases(_pre) : _n_bases(_suf)); } 

  size_t n_suffixes(const bool fw = true) const
  { return (fw ? _n_bases(_suf) : _n_bases(_pre)); }

  size_t n_affixes() const 
  { return _n_bases(_pre) + _n_bases(_suf); }


  unsigned prefix(const unsigned i, const bool fw = true) const
  {
    return (fw ? _base(_pre, i) : 3u ^ _base(_suf, i));
  }

  unsigned suffix(const unsigned i, const bool fw = true) const
  {
    return (fw ? _base(_suf, i) : 3u ^ _base(_pre, i));
  }

  String str_prefixes(const bool fw = true) const 
  {
    return _str_affixes(fw ? _pre : _suf, fw);
  }

  String str_suffixes(const bool fw = true) const 
  {
    return _str_affixes(fw ? _suf : _pre, fw);
  }
};





template<class KMER_t>
class KmerFreqAffixes : public KMER_t, public FreqAffixes
{
public:
  KmerFreqAffixes(const KMER_t & kmer) : KMER_t(kmer), FreqAffixes() {}
  explicit KmerFreqAffixes(const unsigned K = 0) : KMER_t(K), FreqAffixes() {}
 
  friend
  bool operator < (const KmerFreqAffixes & a, const KmerFreqAffixes & b)
  { 
    return (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b));
  }
 
};





TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer29>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer60>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer124>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer248>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer504>);













template<class MAP_t>
class KmerAffixesMapNavigator
{
  typedef typename MAP_t::rec_type  KmerRec_t;
  typedef typename MAP_t::kmer_type Kmer_t;

  const MAP_t & _kmap;
  Kmer_t _kmer_fw;
  Kmer_t _kmer_rc;
  bool _fw;
  KmerRec_t _kmer_rec;

  KmerAffixesMapNavigator(const MAP_t & kmap,
			  const Kmer_t & kmer) : 
    _kmap(kmap),
    _kmer_fw(kmer),
    _kmer_rc(reverse_complement(kmer)),
    _fw(_kmer_fw < _kmer_rc),
    _kmer_rec(kmap(_fw ? _kmer_fw : _kmer_rc))
  {}

  KmerRec_t & kmer_rec() { return _kmer_rec; }

  unsigned n_prefixes() const { return _kmer_rec.n_prefixes(_fw); }
  unsigned n_suffixes() const { return _kmer_rec.n_suffixes(_fw); }

  unsigned prefix(unsigned i_pre) const { return _kmer_rec.prefix(i_pre, _fw); }
  unsigned suffix(unsigned i_suf) const { return _kmer_rec.suffix(i_suf, _fw); }

  void next_by_prefix(const unsigned i_pre)
  {
    ForceAssertLt(i_pre, n_prefixes());

    unsigned base_fw = prefix(i_pre);
    unsigned base_rc = 3u ^ base_fw;
    _kmer_fw.push_left (base_fw);
    _kmer_rc.push_right(base_rc);
    _fw = (_kmer_fw < _kmer_rc);
    _kmer_rec = _kmap(_fw ? _kmer_fw : _kmer_rc);

    ForceAssert(_kmer_rec);
  }

  void next_by_suffix(const unsigned i_suf)
  {
    ForceAssertLt(i_suf, n_suffixes());

    unsigned base_fw = suffix(i_suf);
    unsigned base_rc = 3u ^ base_fw;
    _kmer_fw.push_right(base_fw);
    _kmer_rc.push_left (base_rc);
    _fw = (_kmer_fw < _kmer_rc);
    _kmer_rec = _kmap(_fw ? _kmer_fw : _kmer_rc);

    ForceAssert(_kmer_rec);
  }

};



































template<class KMER_REC_t>
void kmer_freq_affixes_map_build_parallel(const size_t K,
                                          const BaseVecVec & bvv,
                                          const Validator & validator_kf,
                                          const double hash_table_ratio,
                                          KmerMap<KMER_REC_t> * kmap_p,
                                          const size_t verbosity,
                                          const size_t num_threads,
                                          const size_t mean_mem_ceil = 0)
{
  // ---- build kmer vector with frequencies and affixes
  
  vec<KMER_REC_t> kvec;
  if (true) {
    KernelKmerAffixesStorer<KMER_REC_t> storer(bvv, K, &kvec, &validator_kf);
    naif_kmerize(&storer, num_threads, verbosity, mean_mem_ceil);
  }

  if (verbosity > 0)
    cout << TagKFAM() << setw(14) << kvec.size() << " records found." << endl;

  
  // ---- sort kvec by highest kmer frequency 
  //      since we are building a chain hash we want to add first the 
  //      high frequency kmers so that, when we recall them (which will happen
  //      often) they'll come up first.
  
  if (verbosity > 0) cout << TagKFAM() << "Sorting records." << endl; 
  sort(kvec.begin(), kvec.end(), &(kmer_freq_gt<KMER_REC_t>));

  // ---- convert from vec<kmer> to hash table
  
  if (verbosity > 0) cout << TagKFAM() << "Building " << K << "-mer hash table." << endl;
  kmap_p->from_kmer_vec(kvec, hash_table_ratio, verbosity);
  
  if (verbosity > 0) cout << TagKFAM() << "Done building " << K << "-mer hash table." << endl;
  
  if (verbosity > 0) 
    kmer_affixes_map_freq_table_print(*kmap_p);

  //kmer_affixes_map_verify(*kmap_p);

}








template <class KMER_REC_t>
void kmer_freq_affixes_map_build_parallel(const size_t K,
                                          const BaseVecVec & bvv,
                                          const double hash_table_ratio,
                                          KmerMap<KMER_REC_t> * kmap_p,
                                          const size_t verbosity,
                                          const size_t num_threads,
                                          const size_t mean_mem_ceil = 0)
{
  Validator validator_kf;
  kmer_freq_affixes_map_build_parallel(K, bvv, validator_kf, hash_table_ratio, 
                                       kmap_p,
                                       verbosity, num_threads, mean_mem_ceil);
}








template <class KMER_REC_t>
void kmer_spectrum_from_kmer_freq_map(const KmerMap<KMER_REC_t> & kmap,
                                      KmerSpectrum * kspec_p)
{
  const size_t nh = kmap.size_hash();
  for (size_t ih = 0; ih != nh; ih++) {
    const KMER_REC_t & krec = kmap[ih];
    if (krec.is_valid_kmer()) {
      const size_t freq = krec.freq();
      if (kspec_p->size() <= freq)
        kspec_p->resize(freq + 1, 0);
      (*kspec_p)[freq]++;
    }
  }
}

  



// ---- print a table of the kmer frequency regarding number of affixes

template<class KMER_REC_t>
void kmer_affixes_map_freq_table_print(const KmerMap<KMER_REC_t> & kmap)
{
  const size_t nh = kmap.size_hash();
  vec<vec<size_t> > n_kmers(5, vec<size_t>(5, 0));
  vec<vec<size_t> > kf(5, vec<size_t>(5, 0));
    
  size_t n_kmers_total = 0;

  for (size_t ih = 0; ih < nh; ih++) {
    const KMER_REC_t & krec = kmap[ih];
    if (krec.is_valid_kmer()) { 
      unsigned n_pre = krec.n_prefixes();
      unsigned n_suf = krec.n_suffixes();
      n_kmers[n_pre][n_suf]++;
      n_kmers_total++;
      kf[n_pre][n_suf] += krec.freq();
    }
  }
  // ---- output table of flows
  cout << TagKFAM() << endl;
    
  for (size_t n_pre = 0; n_pre <= 4; n_pre++) {
    cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_pre << ")=  " 
	 << setw(10) << n_kmers[n_pre][n_pre]
	 << " (mean_freq= " << setw(6) << kf[n_pre][n_pre] / (n_kmers[n_pre][n_pre] + 1) << ")"
	 << endl;
    for (size_t n_suf = n_pre + 1; n_suf <= 4; n_suf++) {
      size_t nk = (n_kmers[n_pre][n_suf] + n_kmers[n_suf][n_pre]);
      size_t kf2 = (kf[n_pre][n_suf] + kf[n_suf][n_pre]);
      cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_suf << ")=  " 
	   << setw(10) << nk
	   << " (mean_freq= " << setw(6) << kf2 / (nk + 1) << ")"
	   << endl;
    }
    cout << TagKFAM() << endl;
  }
  cout << TagKFAM() << "n_kmers(tot)=  " << setw(10) << n_kmers_total << endl;
  cout << TagKFAM() << endl;
}




template<class KMER_REC_t>
void kmer_affixes_vec_freq_table_print(const vec<KMER_REC_t> & kvec)
{
  const size_t nh = kvec.size();
  vec<vec<size_t> > n_kmers(5, vec<size_t>(5, 0));
  vec<vec<size_t> > kf(5, vec<size_t>(5, 0));
    
  size_t n_kmers_total = 0;

  for (size_t ih = 0; ih < nh; ih++) {
    const KMER_REC_t & krec = kvec[ih];
    unsigned n_pre = krec.n_prefixes();
    unsigned n_suf = krec.n_suffixes();
    n_kmers[n_pre][n_suf]++;
    n_kmers_total++;
    kf[n_pre][n_suf] += krec.freq();
 }
  // ---- output table of flows
  cout << TagKFAM() << endl;
    
  for (size_t n_pre = 0; n_pre <= 4; n_pre++) {
    cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_pre << ")=  " 
	 << setw(10) << n_kmers[n_pre][n_pre]
	 << " (mean_freq= " << setw(6) << kf[n_pre][n_pre] / (n_kmers[n_pre][n_pre] + 1) << ")"
	 << endl;
    for (size_t n_suf = n_pre + 1; n_suf <= 4; n_suf++) {
      size_t nk = (n_kmers[n_pre][n_suf] + n_kmers[n_suf][n_pre]);
      size_t kf2 = (kf[n_pre][n_suf] + kf[n_suf][n_pre]);
      cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_suf << ")=  " 
	   << setw(10) << nk
	   << " (mean_freq= " << setw(6) << kf2 / (nk + 1) << ")"
	   << endl;
    }
    cout << TagKFAM() << endl;
  }
  cout << TagKFAM() << "n_kmers(tot)=  " << setw(10) << n_kmers_total << endl;
  cout << TagKFAM() << endl;
}








// ---- Makes sure that the affixes make sense



template<class KMER_REC_t>
void kmer_affixes_map_verify(const KmerMap<KMER_REC_t> & kmap)
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;

  const size_t nh = kmap.size_hash();
  size_t n_good = 0;
  size_t n_total = 0;
  for (size_t ih = 0; ih != nh; ih++) {
    dots_pct(ih, nh);

    const KMER_REC_t & krec0 = kmap[ih];
    if (krec0) { 
      const Kmer_t kmer0FW(krec0);
      const Kmer_t kmer0RC(reverse_complement(krec0));
      unsigned n_pres = krec0.n_prefixes();
      unsigned n_sufs = krec0.n_suffixes();

      for (unsigned i = 0; i != n_pres; i++) {
	n_total++;
	Kmer_t kmerFW = kmer0FW;
	Kmer_t kmerRC = kmer0RC;
	kmerFW.push_left(krec0.prefix(i));
	kmerRC.push_right(3u ^ krec0.prefix(i));
	Kmer_t & kmer = (kmerFW < kmerRC) ? kmerFW : kmerRC;
	if (kmap(kmer))
	  n_good++;
      }

      for (unsigned i = 0; i != n_sufs; i++) {
	n_total++;
	Kmer_t kmerFW = kmer0FW;
	Kmer_t kmerRC = kmer0RC;
	kmerFW.push_right(krec0.suffix(i));
	kmerRC.push_left(3u ^ krec0.suffix(i));
	Kmer_t & kmer = (kmerFW < kmerRC) ? kmerFW : kmerRC;
	if (kmap(kmer))
	  n_good++;
      }
      
    }
  }      
    
  cout << "n_good= " << n_good << " n_total= " << n_total << endl;
}





#endif
