///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#ifndef KMERS__KMER_SPECTRA_H
#define KMERS__KMER_SPECTRA_H


class KmerSpectrum : public vec<size_t>
{
  size_t _K;

  mutable size_t _kf_min1;
  mutable size_t _kf_min2;
  mutable size_t _kf_min3;
  mutable size_t _kf_max1;
  mutable size_t _kf_max2;
  mutable size_t _kf_hi;

  mutable double _kf_ave_uniq;
    
  mutable size_t _nk_total;
  mutable size_t _nk_bad_low_kf;
  mutable size_t _nk_good_uniq;
  mutable size_t _nk_good_snp;
  mutable size_t _nk_good_rep;
  mutable size_t _nk_bad_high_kf;
  
  mutable size_t _ndk_total;
  mutable size_t _ndk_bad_low_kf;
  mutable size_t _ndk_good_snp;
  mutable size_t _ndk_good_uniq;
  mutable size_t _ndk_good_rep;
  mutable size_t _ndk_bad_high_kf;

  static const size_t _coverage_min = 7;
  mutable float  _coverage;
  mutable size_t _genome_size_unique;
  mutable size_t _genome_size_repetitive;
  mutable size_t _genome_size;

  mutable double _stddev_bias;

  mutable size_t _d_SNP;

  
  String tag(String S = "KS") const { return Date() + " (" + S + "): "; } 

public:

  KmerSpectrum(const size_t K, const size_t n = 0, const size_t val = 0)
    : vec<size_t>(n, val),
      _K(K)
  {}

  friend KmerSpectrum operator+(const KmerSpectrum & ksa, 
                                const KmerSpectrum & ksb)
  {
    KmerSpectrum ks = ksa;
    if (ksb.size() > ksa.size()) ks.resize(ksb.size(), 0);
    for (size_t i = 0; i != ksb.size(); i++)
      ks[i] += ksb[i];
    return ks;
  }

  friend KmerSpectrum operator-(const KmerSpectrum & ksa, 
                                const KmerSpectrum & ksb)
  {
    KmerSpectrum ks = ksa;
    if (ksb.size() > ksa.size()) ks.resize(ksb.size(), 0);
    for (size_t i = 0; i != ksb.size(); i++)
      ks[i] -= ksb[i];
    return ks;
  }


  void analyze(const unsigned ploidy, 
               const unsigned read_len,
               const size_t kf_min1_arg = 10,
               const unsigned verbosity = 0) const;

  size_t sum() const;

  // accessors
  size_t K()       const { return _K; }

  void increment(const size_t kf, const size_t n = 1) 
  {
    if (size() <= kf) resize(kf + 1, 0);
    (*this)[kf] += n;
  }


  size_t kf_min1() const { return _kf_min1; }
  size_t kf_min2() const { return _kf_min2; }
  size_t kf_min3() const { return _kf_min3; }
  size_t kf_max1() const { return _kf_max1; }
  size_t kf_max2() const { return _kf_max2; }
  size_t kf_hi()   const { return _kf_hi; }

  size_t genome_size_unique() const { return _genome_size_unique; }
  size_t genome_size_repetitive() const { return _genome_size_repetitive; }
  size_t genome_size() const { return _genome_size; }
  float  coverage() const { return _coverage; }

  float bias_stddev() const { return _stddev_bias; }
  size_t d_SNP() const { return _d_SNP; }

  float fraction_of_error_kmers() const { return float(_nk_bad_low_kf) / float(_nk_total); }


  String head_K(const String & head) const 
  {
    return head + "." + ToString(_K) + "mer";
  }

  void to_text_file(const String & head, const unsigned verbosity = 1) const;
  void from_text_file(const String & head, const unsigned verbosity = 1);

};










struct KmerBiSpectrum 
{
  KmerSpectrum AB;
  KmerSpectrum A_in;
  KmerSpectrum A_out;
  KmerSpectrum B_in;
  KmerSpectrum B_out;
  
  KmerBiSpectrum(const unsigned K)
  : AB(K), A_in(K), A_out(K), B_in(K), B_out(K) {}

  unsigned K() const { return AB.K(); }

};





class KmerSpectraAffixes : public vec<KmerSpectrum>
{
  unsigned i_ks(const unsigned n_pre, 
                const unsigned n_suf) const
  {
    return (n_pre < n_suf) ?
      n_pre * (9 - n_pre) / 2 + n_suf :
      n_suf * (9 - n_suf) / 2 + n_pre;
  }
  
public:
  KmerSpectraAffixes(const unsigned K) : vec<KmerSpectrum>(15, KmerSpectrum(K)) {}

  KmerSpectrum & operator()(const unsigned n_pre, 
                            const unsigned n_suf)
  { return (*this)[i_ks(n_pre, n_suf)];  }

  const KmerSpectrum & operator()(const unsigned n_pre, 
                                  const unsigned n_suf) const
  { return (*this)[i_ks(n_pre, n_suf)];  }
  
  unsigned K() const { return (*this)[0].K(); }

  KmerSpectrum all() const 
  {
    KmerSpectrum ks(K());
    for (vec<KmerSpectrum>::const_iterator it = begin(); it != end(); it++) {

      const size_t nkf = it->size();
 
      if (ks.size() < nkf) 
        ks.resize(nkf, 0);

      for (size_t kf = 0; kf != nkf; kf++)
        ks[kf] += (*it)[kf];
    }
    return ks;
  }

};









#define __NQ__  51


class KmerQualitySpectra : public vec< vec<size_t> >
{
private:
  const size_t _K;

public:
  KmerQualitySpectra(const size_t K)
    : vec< vec<size_t> >(__NQ__),
      _K(K)
  {}

  unsigned K() const { return _K; }

  void to_text_file(const String & head) const;
};













#endif
