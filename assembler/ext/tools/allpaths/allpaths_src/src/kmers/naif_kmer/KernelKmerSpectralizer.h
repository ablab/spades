///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Filipe Ribeiro - Feb 2011 - <crdhelp@broadinstitute.org>
//   



#ifndef KMERS__NAIF_KMER__KERNEL_KMER_SPECTRALIZER_H
#define KMERS__NAIF_KMER__KERNEL_KMER_SPECTRALIZER_H

#include "Basevector.h"
#include "Qualvector.h"

#include "system/LockedData.h"

#include "kmers/KmerSpectra.h"
#include "kmers/naif_kmer/Kmers.h"







// -------- KERNEL_t class --------

template<class KMER_t>
class KernelKmerSpectralizer
{
private:
  const BaseVecVec & _bvv;
  const size_t       _K;

  KmerSpectrum     * _p_kspec;     // final results go here; only merge() adds to it
  KmerSpectrum       _kspec_tmp;   // only used in clones
  LockedData         _lock;        // lock for merge()


public:
  typedef KMER_t      rec_type;

  KernelKmerSpectralizer(const BaseVecVec & bvv, 
                         const size_t       K,
                         KmerSpectrum     * p_kspec) 
    : _bvv(bvv), 
      _K(K), 
      _p_kspec(p_kspec), 
      _kspec_tmp(_K),
      _lock() 
  {}

  // copy constructor for temporary kernels
  explicit KernelKmerSpectralizer(const KernelKmerSpectralizer & that)
    : _bvv(that._bvv),
      _K(that._K),
      _p_kspec(0),
      _kspec_tmp(_K),
      _lock()
  {
  }

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 
  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_parcels, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_parcels;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcels.in_one_parcel(kmer))
        parcels.add(kmer);
      kmer_cur.next();          
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & kmer_recs,  // in this case, unused
                 const size_t i0,
                 const size_t i1)
  {
    const size_t kf = i1 - i0;
    if (_kspec_tmp.size() <= kf) _kspec_tmp.resize(kf + 1, 0);
    _kspec_tmp[kf]++;
  }
  

  
  // interface function needed by naif_kmerize()
  void merge(const KernelKmerSpectralizer & kernel_tmp,
	     const size_t i_parcel)
  {
    const size_t nkf = kernel_tmp._kspec_tmp.size();

    Locker lock(_lock);
    
    if (_p_kspec->size() < nkf) 
      _p_kspec->resize(nkf, 0);
    
    for (size_t kf = 0; kf != nkf; kf++)
      (*_p_kspec)[kf] += kernel_tmp._kspec_tmp[kf];
  }

};






// -------- KERNEL_t class --------

  
template<class KMER_t>
class KernelKmerAffixesSpectralizer
{
private:
  const BaseVecVec   & _bvv;
  const size_t         _K;
  const size_t         _Ks;

  KmerSpectraAffixes * _p_kspecs;     // final results go here; only merge() adds to it
  KmerSpectraAffixes   _kspecs_tmp;   // only used in clones
  LockedData           _lock;         // lock for merge()

  const Validator * _p_validator;
  
  
public:
  typedef KmerBVLoc<KMER_t>       rec_type;

  KernelKmerAffixesSpectralizer(const BaseVecVec        & bvv, 
                                const size_t              K,
                                KmerSpectraAffixes      * p_kspecs,
                                const Validator         * p_validator = 0) 
    : _bvv(bvv), 
      _K(K), 
      _Ks(K - 2),
      _p_kspecs(p_kspecs), 
      _kspecs_tmp(_K),
      _lock(),
      _p_validator(p_validator) 
  {
    ForceAssert(K % 2 == 1);
  }

  // copy constructor for temporary kernels
  explicit KernelKmerAffixesSpectralizer(const KernelKmerAffixesSpectralizer & that)
    : _bvv(that._bvv),
      _K(that._K),
      _Ks(that._Ks),
      _p_kspecs(0),
      _kspecs_tmp(_K),
      _lock(),
      _p_validator(that._p_validator)
  {
  }

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 
  // interface function needed by naif_kmerize()
  // Here we look for (k-2)mers and store their locations
  void parse_base_vec(ParcelBuffer<rec_type> * p_parcels, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_parcels;
    SubKmers<BaseVec, rec_type> kmer_cur(_Ks, _bvv[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcels.in_one_parcel(kmer)) {
        kmer.set_ibv(ibv);
        kmer.set_ib(kmer_cur.index_start());
        kmer.set_fw(kmer_cur.is_canonical_fw());
        kmer.set_palindrome(kmer_cur.is_palindrome());
        parcels.add(kmer);
      }
      kmer_cur.next();          
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs, 
                 const size_t i0k,
                 const size_t i1k)
  {
    // palidromes are not implemented right yet.
    
    // compute frequencies for various kmers around central (k-2)mer
    vec<vec<size_t> > kf_pre_suf(4, vec<size_t>(4, 0));
    vec<vec<size_t> > kf_suf_suf(4, vec<size_t>(4, 0));
    vec<vec<size_t> > kf_pre_pre(4, vec<size_t>(4, 0));

    for (size_t ik = i0k; ik < i1k; ik++) {
      const rec_type & krec    = krecs[ik];
      const BaseVec  & bv      = _bvv[krec.ibv()];
      const unsigned   nb      = bv.size();
      const unsigned   ib_kmer = krec.ib();

      if (krec.is_fw() || krec.is_palindrome()) {         // kmer has FW orientation in read

        if (ib_kmer >= 2) {
          const unsigned pre1 = bv[ib_kmer - 1];
          const unsigned pre2 = bv[ib_kmer - 2];
          kf_pre_pre[pre2][pre1]++;
        }
        if (ib_kmer >= 1 && ib_kmer + _Ks < nb) {
          const unsigned pre1 = bv[ib_kmer - 1];
          const unsigned suf1 = bv[ib_kmer + _Ks];
          kf_pre_suf[pre1][suf1]++;
        }
        if (ib_kmer + _Ks + 1 < nb) {
          const unsigned suf1 = bv[ib_kmer + _Ks];
          const unsigned suf2 = bv[ib_kmer + _Ks + 1];
          kf_suf_suf[suf1][suf2]++;
        }
      }         

      if (krec.is_rc() || krec.is_palindrome()) {         // kmer has RC orientation in read

        if (ib_kmer >= 2) {
          const unsigned suf1 = 3u ^ bv[ib_kmer - 1];
          const unsigned suf2 = 3u ^ bv[ib_kmer - 2];
          kf_suf_suf[suf1][suf2]++;
        }
        if (ib_kmer >= 1 && ib_kmer + _Ks < nb) {
          const unsigned suf1 = 3u ^ bv[ib_kmer - 1];
          const unsigned pre1 = 3u ^ bv[ib_kmer + _Ks];
          kf_pre_suf[pre1][suf1]++;
        }
        if (ib_kmer + _Ks + 1 < nb) {
          const unsigned pre1 = 3u ^ bv[ib_kmer + _Ks];
          const unsigned pre2 = 3u ^ bv[ib_kmer + _Ks + 1];
          kf_pre_pre[pre2][pre1]++;
        }
      }         

    }

    // summarize kmer frequencies based on number of affixes

    for (unsigned pre = 0; pre < 4; pre++) {
      for (unsigned suf = 0; suf < 4; suf++) {

        const size_t kf = kf_pre_suf[pre][suf];

        if (!_p_validator || (*_p_validator)(kf)) {
          
          unsigned n_pre = 0;
          for (unsigned pre2 = 0; pre2 < 4; pre2++)
            if (!_p_validator || (*_p_validator)(kf_pre_pre[pre2][pre]))
              n_pre++;

          unsigned n_suf = 0;
          for (unsigned suf2 = 0; suf2 < 4; suf2++)
            if (!_p_validator || (*_p_validator)(kf_suf_suf[suf][suf2]))
              n_suf++;
          
          // pick kmer spectrum corresponding to n_pre and n_suf
          // (invariant under n_pre <-> n_suf)
          KmerSpectrum & kspec_tmp = _kspecs_tmp(n_pre, n_suf);

          if (kspec_tmp.size() <= kf) kspec_tmp.resize(kf + 1, 0);
          kspec_tmp[kf]++;

        }

      }
    }



  }
  

  
  // interface function needed by naif_kmerize()
  void merge(const KernelKmerAffixesSpectralizer & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
 
    for (unsigned n_pre = 0; n_pre <= 4; n_pre++) {
      for (unsigned n_suf = n_pre; n_suf <= 4; n_suf++) {
    
        KmerSpectrum & kspec = (*_p_kspecs)(n_pre, n_suf);
          const KmerSpectrum & kspec_tmp = kernel_tmp._kspecs_tmp(n_pre, n_suf);
          const size_t nkf = kspec_tmp.size();

          if (kspec.size() < nkf) 
            kspec.resize(nkf, 0);
    
          for (size_t kf = 0; kf != nkf; kf++)
            kspec[kf] += kspec_tmp[kf];
      }
    }

  }

};






// -------- KERNEL_t class --------


template<class KMER_t>
class KernelKmerBiSpectralizer
{
private:
  
  const BaseVecVec & _bvv;      // all bv's, first from set A followed by set B
  const size_t       _nbvA;       // number of bv's in set A
  const size_t       _K;

  KmerBiSpectrum   * _p_kbspec;    // final results go here; only merge() adds to it
  KmerBiSpectrum     _kbspec_tmp;  // only used in clones
  LockedData         _lock;        // lock for merge()


public:
  typedef KMER_t      rec_type;

  KernelKmerBiSpectralizer(const BaseVecVec & bvv, 
                           const size_t       nbvA, 
                           const size_t       K,
                           KmerBiSpectrum   * p_kbspec) 
    : _bvv(bvv), 
      _nbvA(nbvA), 
      _K(K), 
      _p_kbspec(p_kbspec), 
      _kbspec_tmp(_K),
      _lock() 
  {}

  // copy constructor for temporary kernels
  explicit KernelKmerBiSpectralizer(const KernelKmerBiSpectralizer & that)
    : _bvv(that._bvv),
      _nbvA(that._nbvA), 
      _K(that._K),
      _p_kbspec(0),
      _kbspec_tmp(_K),
      _lock()
  {
  }

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 

  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_parcels, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_parcels;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
    const unsigned in_A = (ibv < _nbvA); 
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      kmer.tag_set(in_A); 
      if (parcels.in_one_parcel(kmer))
        parcels.add(kmer);
      kmer_cur.next();          
    }
  }

  

  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & parcel,  // in this case, unused
                 const size_t i0,
                 const size_t i1)
  {
    size_t kfA = 0;
    size_t kfB = 0;
    for (size_t i = i0; i != i1; i++)
      if (parcel[i].tag()) kfA++;
      else                 kfB++;
    
    if      (kfB == 0)  _kbspec_tmp.A_out.increment(kfA);
    else if (kfA == 0)  _kbspec_tmp.B_out.increment(kfB);
    else {
      _kbspec_tmp.A_in.increment(kfA);
      _kbspec_tmp.B_in.increment(kfB);
    }
    
    _kbspec_tmp.AB.increment(kfA + kfB);

  }
  

  
  // interface function needed by naif_kmerize()
  void merge(const KernelKmerBiSpectralizer & kks,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    {
      const KmerSpectrum & ks_tmp = kks._kbspec_tmp.AB;
      KmerSpectrum       & ks     = (*_p_kbspec).AB;
      const size_t nkf  = ks_tmp.size();
      if (ks.size() < nkf) ks.resize(nkf, 0);
      for (size_t kf = 0; kf != nkf; kf++)
        ks[kf] += ks_tmp[kf];
    }
    {
      const KmerSpectrum & ks_tmp = kks._kbspec_tmp.A_in;
      KmerSpectrum       & ks     = (*_p_kbspec).A_in;
      const size_t nkf  = ks_tmp.size();
      if (ks.size() < nkf) ks.resize(nkf, 0);
      for (size_t kf = 0; kf != nkf; kf++)
        ks[kf] += ks_tmp[kf];
    }
    {
      const KmerSpectrum & ks_tmp = kks._kbspec_tmp.A_out;
      KmerSpectrum       & ks     = (*_p_kbspec).A_out;
      const size_t nkf  = ks_tmp.size();
      if (ks.size() < nkf) ks.resize(nkf, 0);
      for (size_t kf = 0; kf != nkf; kf++)
        ks[kf] += ks_tmp[kf];
    }
    {
      const KmerSpectrum & ks_tmp = kks._kbspec_tmp.B_in;
      KmerSpectrum       & ks     = (*_p_kbspec).B_in;
      const size_t nkf  = ks_tmp.size();
      if (ks.size() < nkf) ks.resize(nkf, 0);
      for (size_t kf = 0; kf != nkf; kf++)
        ks[kf] += ks_tmp[kf];
    }
    {
      const KmerSpectrum & ks_tmp = kks._kbspec_tmp.B_out;
      KmerSpectrum       & ks     = (*_p_kbspec).B_out;
      const size_t nkf  = ks_tmp.size();
      if (ks.size() < nkf) ks.resize(nkf, 0);
      for (size_t kf = 0; kf != nkf; kf++)
        ks[kf] += ks_tmp[kf];
    }
  }

};




// -------- KERNEL_t class --------


template<class KMER_t>
class KernelKmerQualitySpectralizer
{
private:
  const BaseVecVec   & _bvv;
  const QualVecVec   & _qvv;
  const size_t         _K;
  
  KmerQualitySpectra * _p_kqspec;   // final results go here; only merge() adds to it
  KmerQualitySpectra   _kqspec_tmp; // only used in temps
  LockedData           _lock;       // lock for merge()
 
public:
  typedef KmerKmerQual<KMER_t> rec_type;
  
  KernelKmerQualitySpectralizer(const BaseVecVec   & bvv, 
                                const QualVecVec   & qvv,
                                const size_t         K,
                                KmerQualitySpectra * p_kqspec)
    :  _bvv(bvv), 
       _qvv(qvv), 
       _K(K),
       _p_kqspec(p_kqspec),
       _kqspec_tmp(_K),
       _lock()
  {}
  
  explicit KernelKmerQualitySpectralizer(const KernelKmerQualitySpectralizer & that)
    : _bvv(that._bvv),
      _qvv(that._qvv),
      _K(that._K),
      _p_kqspec(0),
      _kqspec_tmp(_K),
      _lock()
  {}


  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
  
  



  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_parcels, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_parcels;
    const BaseVec & bv = _bvv[ibv];
    const QualVec & qv = _qvv[ibv];
    const unsigned nb = bv.size();
    rec_type kmerFW(_K);
    rec_type kmerRC(_K);
    vec<uint8_t> q(_K);
      
    // ---- load the 1st K - 1 bases
    for (unsigned ib = 0; ib != _K - 1; ib++) {
      const uint64_t base = bv[ib];
      kmerFW.push_right(base);
      kmerRC.push_left(3ul ^ base); // the complement of base
      q[ib] = qv[ib];
    }
    
    // ---- build kmers by adding one base at a time
    for (unsigned ib = _K - 1; ib != nb; ib++) {
      const uint64_t base = bv[ib];
      kmerFW.push_right(base);
      kmerRC.push_left(3ul ^ base);  // the complement of base
      
      // ---- replace the quality for base left off from previous kmer
      //      think about it... it works! 
      //      (because we only care about the minimum, not the position of the qual)
      q[ib % _K] = qv[ib]; 
      
      rec_type & kmerq = (kmerFW < kmerRC) ? kmerFW : kmerRC; 
      if (parcels.in_one_parcel(kmerq)) {
	
        // ---- find the MINIMUM base quality in the kmer
        size_t qmin = q[0]; 
        for (size_t iq = 1; iq != _K; iq++)
          if (qmin < q[iq]) qmin = q[iq];
        kmerq.qual(qmin);
	
        parcels.add(kmerq);
      }
      //cout << this->back().hieroglyphs() << " " << kmerFW.hieroglyphs() << endl;
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & kmer_recs,
                 const size_t ik0,
                 const size_t ik1)
  {
    // ---- frequency is easy
    const size_t kf = ik1 - ik0;
    
    // ---- the kmer quality is the MEDIAN of all the kmer qualities
    vec<unsigned> qs(kf);
    for (size_t ik = ik0; ik != ik1; ik++)
      qs[ik - ik0] = kmer_recs[ik].qual();
    sort(qs.begin(), qs.end());

    ForceAssert(kf > 0u);
    const unsigned kq = (kf & 1u) ? qs[kf / 2] : (qs[kf / 2 - 1] + qs[kf / 2]) / 2;
    ForceAssert(kq < 41u);

    if (_kqspec_tmp[0].size()  <= kf) _kqspec_tmp[0].resize(kf + 1, 0);
    if (_kqspec_tmp[kq].size() <= kf) _kqspec_tmp[kq].resize(kf + 1, 0);
    _kqspec_tmp[0][kf]++;
    _kqspec_tmp[kq][kf]++;
  }



  // interface function needed by naif_kmerize()
  void merge(const KernelKmerQualitySpectralizer & kqs,
	     const size_t i_parcel)  // i_parcel for debug purposes only
  {

    Locker lock(_lock);

    const size_t nkq = kqs._kqspec_tmp.size();
    size_t nkf_max = 0;
    for (size_t kq = 0; kq < nkq; kq++) {
      const size_t nkf = kqs._kqspec_tmp[kq].size();
      if (nkf_max < nkf) nkf_max = nkf;
    }
    
    for (size_t kq = 0; kq < nkq; kq++) {
      if ((*_p_kqspec)[kq].size() < nkf_max)
        (*_p_kqspec)[kq].resize(nkf_max, 0);
      
      const size_t nkf = kqs._kqspec_tmp[kq].size();
      for (size_t kf = 0; kf != nkf; kf++) 
	(*_p_kqspec)[kq][kf] += kqs._kqspec_tmp[kq][kf];
    }    
  }
  


};






#endif
