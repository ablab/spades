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


#ifndef KMERS__NAIF_KMER__KERNEL_UNIQUE_FINDER_H
#define KMERS__NAIF_KMER__KERNEL_UNIQUE_FINDER_H

#include "kmers/naif_kmer/LockedBlocks.h"
#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/Kmers.h"



// -------- KERNEL_t class --------





template<class KMER_t, class BOOL_t>
class UniqueFinder
{
  const size_t             _n_threads;
  const size_t             _n_blocks;
  LockedBlocks             _blocks;
  LockedBlocks           * _p_blocks;
  vec< vec<size_t> >       _i_read_unique;

  const size_t             _K;
  const size_t             _how_unique;
  const BaseVecVec       & _bases;
 
  vec<BOOL_t>            & _is_unique;


  KmerSpectrum           * _p_kspec;     // final results go here; only merge() adds to it
  KmerSpectrum             _kspec_tmp;   // only used in clones
  LockedData               _lock;        // lock for merge()


public:

  typedef KmerBVLoc<KMER_t>        rec_type;

  UniqueFinder(const size_t       K,
               const size_t       how_unique,
	       const BaseVecVec & bases,
	       vec<BOOL_t> *      p_is_unique,
               KmerSpectrum     * p_kspec,
               const size_t       n_threads)
    : _n_threads(n_threads),
      _n_blocks(2 * n_threads),
      _blocks(_n_blocks),
      _p_blocks(&_blocks),
      _i_read_unique(_n_blocks),
      _K(K),
      _how_unique(how_unique),
      _bases(bases),
      _is_unique(*p_is_unique), 
      _p_kspec(p_kspec), 
      _kspec_tmp(_K),
      _lock()
  {
    ForceAssertEq(_bases.size(), _is_unique.size());
    const size_t n = bases.size();
    for (size_t i = 0; i < n; i++)
      _is_unique[i] = false;
  }
    

  ~UniqueFinder()
  {}
      

  
  // copy constructor for temporary kernels
  explicit UniqueFinder(const UniqueFinder & that)
    : _n_threads(that._n_threads),
      _n_blocks(that._n_blocks),
      _blocks(),                   // no local blocks
      _p_blocks(that._p_blocks),   // get the blocks from 'that'
      _i_read_unique(_n_blocks),
      _K(that._K),
      _how_unique(that._how_unique),
      _bases(that._bases),
      _is_unique(that._is_unique),
      _p_kspec(0),
      _kspec_tmp(_K),
      _lock()
  {}
 
  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bases; }
  
  

  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcel_buf = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bases[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcel_buf.in_one_parcel(kmer)) {
        kmer.set_ibv(ibv);          // just need the read id (bv index) for uniqueness
        parcel_buf.add(kmer);
      }
      kmer_cur.next();          
    }
  } 




  void flush_uniques(const size_t max_buf_size)
  {
    bool verbose = false;//(max_buf_size == 0);

    std::set<size_t> blocks_todo;
    for (size_t i = 0; i < _n_blocks; i++)
      if (_i_read_unique[i].size() > max_buf_size)  // buffers too big? flush them!
        blocks_todo.insert(i);

    
    
    while (blocks_todo.size() > 0) {

      // ---- lock a subset of the data to update

      const size_t i_blk = _p_blocks->lock_some_block(blocks_todo);

      

      // ---- set is_unique[i_read] = true

      vec<size_t> & uniqs = _i_read_unique[i_blk];
      const size_t n_uniqs = uniqs.size();
      for (size_t i = 0; i < n_uniqs; i++)
	_is_unique[uniqs[i]] = true;
      uniqs.clear();



      // ---- unlock block and remove it from to-do list

      _p_blocks->unlock_block(i_blk);
      blocks_todo.erase(i_blk);

    } // while (blocks_todo.size() > 0) 

  } // void flush_uniques(...)



  // _is_unique[] is a vec<BOOL_t> which could be bit aligned
  // so we must make sure that every 64 entries are in the same block
  // hence the '>> 6' (i.e. '/64')
  size_t i_block(const size_t i, const size_t n) 
  { 
    return ((i >> 6) * _n_blocks) / ((n + 63) >> 6);  
  }



  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & parcel, 
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t kfreq = i1k - i0k;

    if (_kspec_tmp.size() <= kfreq) _kspec_tmp.resize(kfreq + 1, 0);
    _kspec_tmp[kfreq]++;

    if (kfreq <= _how_unique) {
      const size_t n_reads = _bases.size();
      for (size_t ik = i0k; ik != i1k; ik++) {
	const size_t i_read = parcel[ik].ibv();
	const size_t i_blk  = i_block(i_read, n_reads);
	_i_read_unique[i_blk].push_back(i_read);
      }
      flush_uniques(5000);
    }
  }
  
 










  // ---- this merges all the recommendations and confirmations
  // interface function needed by naif_kmerize()
  void merge(UniqueFinder & kernel_tmp,
	     const size_t i_parcel)
  {
    // final flush of uniques
    kernel_tmp.flush_uniques(0);


    // merge the local kmer spectrum into the global
    const size_t nkf = kernel_tmp._kspec_tmp.size();

    Locker lock(_lock);
    
    if (_p_kspec->size() < nkf) 
      _p_kspec->resize(nkf, 0);
    
    for (size_t kf = 0; kf != nkf; kf++)
      (*_p_kspec)[kf] += kernel_tmp._kspec_tmp[kf];
  }






};  // class UniqueFinder
  







#endif
