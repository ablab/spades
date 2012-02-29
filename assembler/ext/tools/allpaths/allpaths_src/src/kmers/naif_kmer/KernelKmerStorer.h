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



#ifndef KMERS__NAIF_KMER__KERNEL_KMER_STORER_H
#define KMERS__NAIF_KMER__KERNEL_KMER_STORER_H


#include "kmers/naif_kmer/Kmers.h"





// -------- KERNEL_t class --------


template<class ELEM_t>
class KernelKmerStorer
{
private:
  const BaseVecVec   & _bvv;
  const size_t         _K;

  vec<ELEM_t>          _kvec_tmp; // only used in temporaries
  vec<ELEM_t>        * _p_kvec;     // final results go here; only merge() adds to it
  LockedData           _lock;     // lock for merge()

  const Validator    * _p_validator;
  
public:
  typedef typename ELEM_t::kmer_type  rec_type;

  KernelKmerStorer(const BaseVecVec & bvv, 
                   const size_t       K,
                   vec<ELEM_t>      * p_kvec,
                   const Validator  * p_validator = 0) 
    : _bvv(bvv), 
      _K(K), 
      _kvec_tmp(), 
      _p_kvec(p_kvec), 
      _lock(),
      _p_validator(p_validator)
  {}

  // copy constructor for temporary kernels
  explicit KernelKmerStorer(const KernelKmerStorer<ELEM_t> & that)
    : _bvv(that._bvv),
      _K(that._K),
      _kvec_tmp(),
      _p_kvec(0),
      _lock(),
      _p_validator(that._p_validator)
  {}

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 

  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcel_buf = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcel_buf.in_one_parcel(kmer))
        parcel_buf.add(kmer);
      kmer_cur.next();          
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs,  
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t kf = i1k - i0k;

    if (!_p_validator || (*_p_validator)(kf)) {
      ELEM_t krec(krecs[i0k]);
      krec.set_freq(kf);
      _kvec_tmp.push_back(krec);
    }
  }
 

  // interface function needed by naif_kmerize()
  void merge(const KernelKmerStorer<ELEM_t> & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    _p_kvec->insert(_p_kvec->end(), kernel_tmp._kvec_tmp.begin(), kernel_tmp._kvec_tmp.end());
  }

};













// -------- KERNEL_t class --------



template<class ELEM_t>
class KernelKmerAffixesStorer
{
private:
  const BaseVecVec   & _bvv;
  const size_t         _K;
  const size_t         _Ks;

  vec<ELEM_t>        * _p_kvec;   // final results go here; only merge() adds to it
  vec<ELEM_t>          _kvec_tmp; // only used in temporaries
  LockedData           _lock;     // lock for merge()

  const Validator    * _p_validator;
  
public:
  typedef KmerBVLoc<typename ELEM_t::kmer_type>  rec_type;

  KernelKmerAffixesStorer(const BaseVecVec & bvv, 
                          const size_t       K,
                          vec<ELEM_t>      * p_kvec,
                          const Validator  * p_validator = 0) 
    : _bvv(bvv), 
      _K(K), 
      _Ks(K - 2),
      _p_kvec(p_kvec), 
      _kvec_tmp(), 
      _lock(),
      _p_validator(p_validator)
  { 
    ForceAssert(K % 2 == 1);
  }

  // copy constructor for temporary kernels
  explicit KernelKmerAffixesStorer(const KernelKmerAffixesStorer<ELEM_t> & that)
    : _bvv(that._bvv),
      _K(that._K),
      _Ks(that._Ks),
      _p_kvec(0),
      _kvec_tmp(),
      _lock(),
      _p_validator(that._p_validator)
  {}

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

    // summarize kmers and store them with number of affixes and kmer frequencies

    const rec_type & ks_mer_fw = krecs[i0k];
    const rec_type   ks_mer_rc = reverse_complement(ks_mer_fw);


    for (unsigned pre = 0; pre < 4; pre++) {
      for (unsigned suf = 0; suf < 4; suf++) {

        const size_t kf = kf_pre_suf[pre][suf];
        
        if (!_p_validator || (*_p_validator)(kf)) { // kmer frequency is valid
          
          // build kmer with affix info
          ELEM_t kmerFW(ks_mer_fw);
          kmerFW.add_left (pre);
          kmerFW.add_right(suf);

          ELEM_t kmerRC(ks_mer_rc);
          kmerRC.add_left (3u ^ suf);
          kmerRC.add_right(3u ^ pre);
          
          if (kmerFW < kmerRC) {
            kmerFW.set_freq(kf);
            for (unsigned pre2 = 0; pre2 < 4; pre2++)
              if (!_p_validator || (*_p_validator)(kf_pre_pre[pre2][pre]))
                kmerFW.set_prefix(pre2);

            for (unsigned suf2 = 0; suf2 < 4; suf2++)
              if (!_p_validator || (*_p_validator)(kf_suf_suf[suf][suf2]))
                kmerFW.set_suffix(suf2);

            _kvec_tmp.push_back(kmerFW);

          }
          else {
            kmerRC.set_freq(kf);
            for (unsigned pre2 = 0; pre2 < 4; pre2++)
              if (!_p_validator || (*_p_validator)(kf_pre_pre[pre2][pre]))
                kmerRC.set_suffix(3u ^ pre2);

            for (unsigned suf2 = 0; suf2 < 4; suf2++)
              if (!_p_validator || (*_p_validator)(kf_suf_suf[suf][suf2]))
                kmerRC.set_prefix(3u ^ suf2);

            _kvec_tmp.push_back(kmerRC);
          }

        }
      }
    }
    
  }
    

  // interface function needed by naif_kmerize()
  void merge(const KernelKmerAffixesStorer<ELEM_t> & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    _p_kvec->insert(_p_kvec->end(), kernel_tmp._kvec_tmp.begin(), kernel_tmp._kvec_tmp.end());
  }

};















#endif
