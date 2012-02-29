///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Filipe Ribeiro - Aug 2011 - <crdhelp@broadinstitute.org>
//


#ifndef KMERS__NAIF_KMER__KERNEL_PRE_CORRECTOR_H
#define KMERS__NAIF_KMER__KERNEL_PRE_CORRECTOR_H

#include <set>
#include <iostream>
#include <sstream>

#include "system/LockedData.h"

#include "kmers/naif_kmer/LockedBlocks.h"

#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/Kmers.h"




class BaseColumn 
{
public:
  const PC_Params &     _params;
  vec< vec<unsigned> >  nq; 
  unsigned              q_sum[4];
  unsigned              n_q_high[4];
  unsigned              base_winner;
  bool                  fix_base[4];

  BaseColumn(const PC_Params & params) : _params(params) , nq(4, vec<unsigned>(42, 0)) { reset(); }

  void reset() 
  {
    base_winner = 0;
    q_sum   [0] = q_sum   [1] = q_sum   [2] = q_sum   [3] = 0;
    n_q_high[0] = n_q_high[1] = n_q_high[2] = n_q_high[3] = 0;
    fix_base[0] = fix_base[1] = fix_base[2] = fix_base[3] = false;
  }

  void add_base_qual(const unsigned base, const unsigned qual) 
  {
    nq[base][qual]++;
    q_sum[base] += qual;
    if (qual >= _params.q_high) n_q_high[base]++;
  }
    
  void process()
  {
    // find winner base
    base_winner = 0;
    if (q_sum[1] > q_sum[base_winner]) base_winner = 1; 
    if (q_sum[2] > q_sum[base_winner]) base_winner = 2; 
    if (q_sum[3] > q_sum[base_winner]) base_winner = 3; 
    
    const unsigned q_sum_winner = q_sum[base_winner];
    const unsigned q_sum_max_to_lose = q_sum_winner * _params.ratio_loser;

    if (q_sum_winner >= _params.q_sum_min_to_win) {  // winner is strong to win
      
      for (unsigned base = 0; base < 4; base++)
        if (base           != base_winner &&
            n_q_high[base] <= _params.n_q_high_max_to_lose &&
            q_sum[base]     < q_sum_max_to_lose)
          fix_base[base] = true;
            
    } // if (q_sum_winner >= _params.min_q_sum_to_win) 
    
  } // void process()


  void print() 
  {
    cout << "  q_sum: ";
    for (size_t base = 0; base != 4; base++) cout << " " << setw(5) << q_sum[base];
    cout << "  n_q_high: ";
    for (size_t base = 0; base != 4; base++) cout << " " << setw(2) << n_q_high[base];
    cout << "  Win: " << base_winner;
    cout << "  fix_base: ";
    for (size_t base = 0; base != 4; base++) cout << unsigned(fix_base[base]);
    cout << endl;
  }
  



};  // class BaseColumn








// -------- KERNEL_t class --------

template<class KMER_t, class QVV_t>
class PreCorrector
{
  const size_t             _n_threads;
  const size_t             _n_blocks;
  LockedBlocks             _blocks;
  LockedBlocks           * _p_blocks;
  vec< vec<BVLocBase> >    _recommends;
  
  const PC_Params        & _params;
  
  const size_t             _K;

  const BaseVecVec       & _bases;
  const QVV_t            & _quals;
  BaseVecVec             & _bases_new;
  
  KmerSpectrum           * _p_kspec;     // final results go here; only merge() adds to it
  KmerSpectrum             _kspec_tmp;   // only used in clones
  LockedData               _lock;        // lock for merge()
  
  
public:
  
  typedef KmerBVLocBase<KMER_t>       rec_type;
  
  PreCorrector(const PC_Params      & params,
               const size_t           K,
               const BaseVecVec     & bases,
               const QVV_t          & quals,
               BaseVecVec           * p_bases_new,
               KmerSpectrum         * p_kspec,
               const size_t           n_threads) :
    _n_threads(n_threads),
    _n_blocks(2 * n_threads),
    _blocks(_n_blocks),
    _p_blocks(&_blocks),
    _recommends(_n_blocks),
    _params(params),
    _K(K),
    _bases(bases),
    _quals(quals),
    _bases_new(*p_bases_new),
    _p_kspec(p_kspec), 
    _kspec_tmp(_K),
    _lock()
  {
    ForceAssert(_K & 1);  // K must be odd
  }
  
  ~PreCorrector()
  {}
      

  
  // copy constructor for temporary kernels
  explicit PreCorrector(const PreCorrector & that) :
    _n_threads(that._n_threads),
    _n_blocks(that._n_blocks),
    _blocks(),                   // no local blocks
    _p_blocks(that._p_blocks),   // get the blocks from 'that'
    _recommends(_n_blocks),
    _params(that._params),
    _K(that._K),
    _bases(that._bases),
    _quals(that._quals),
    _bases_new(that._bases_new),
    _p_kspec(0),
    _kspec_tmp(_K),
    _lock()
  {
    //cout << "clone: n_blocks = " << _n_blocks << endl; 
    //cout << "clone: blocks.size() = " << _p_blocks->size() << endl; 
  }
  
  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }
  
  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bases; }
  
  

  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bases[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      //if (p_buf->in_one_parcel(kmer) && 
      //    !kmer_cur.is_palindrome()) {  // skip palindromes
      if (p_buf->in_one_parcel(kmer)) { // deal with palindromes in summarize()
        kmer.i_read     = ibv;
        kmer.i_base     = kmer_cur.index_start();
        kmer.rc         = kmer_cur.is_canonical_rc();
        kmer.palindrome = kmer_cur.is_palindrome();
        p_buf->add(kmer);
      }
      kmer_cur.next();          
    }
  } 



  void flush_recommends(const size_t max_buf_size, 
                        const size_t parcel_size) // parcel_size is for debug
  {
    bool debug = false;//(max_buf_size == 0);

    std::set<size_t> blocks_todo;
    for (size_t i = 0; i < _n_blocks; i++)
      if (_recommends[i].size() > max_buf_size) // buffers too big? flush them! 
        blocks_todo.insert(i);

    
    
    while (blocks_todo.size() > 0) {

      // ---- lock a subset of the data to update

      const size_t iblk = _p_blocks->lock_some_block(blocks_todo);

 
      // ---- apply recommendations to '_bases_new'

      vec<BVLocBase> & recoms = _recommends[iblk];
      const size_t n_recoms = recoms.size();
      for (size_t i = 0; i < n_recoms; i++) { 
        const size_t i_read = recoms[i].i_read;
        const size_t i_base = recoms[i].i_base;

        if (_bases[i_read][i_base] != recoms[i].base) {
          _bases_new[i_read].Set(i_base, recoms[i].base);
        }
      }
      recoms.clear();


      // ---- unlock block and remove it from to-do list

      _p_blocks->unlock_block(iblk);
      blocks_todo.erase(iblk);

    } // while (blocks_todo.size() > 0) 

  } // void flush_recommends(...)
  



  size_t i_block(const size_t i, const size_t n) { return (i * _n_blocks) / n;  }
  

  bool stack_debug_print(const vec<rec_type> & krecs, 
                         const size_t i0k,
                         const size_t i1k,
                         const BaseColumn & col)
  {

    bool print_stack = false;
    for (size_t ik = i0k; ik != i1k && !print_stack; ik++)
      if (_params.i_reads_print.count(krecs[ik].i_read))
        print_stack = true;

    if (print_stack) {
      stringstream oss;

      if (false) {
	vec< pair<unsigned, size_t> > i_base_i_krec;
	for (size_t ik = i0k; ik != i1k; ik++) {
	  const rec_type & krec = krecs[ik];
	  const size_t     ir   = krec.i_read;
	  const unsigned   ib   = krec.i_base;
	  const BaseVec  & bv   = _bases[ir];
	  const unsigned   nb   = bv.size();
	  const unsigned   ibc  = (krec.rc) ? (nb - ib - _K) : ib; // canonical ib
	  i_base_i_krec.push_back(make_pair(ibc, ik));
	}
	sort(i_base_i_krec.begin(), i_base_i_krec.end());
	
	const size_t nk = i1k - i0k;
	for (size_t iik = 0; iik != nk; iik++) {
	  const size_t     ibc  = i_base_i_krec[iik].first;
	  const size_t     ik   = i_base_i_krec[iik].second;
	  const rec_type & krec = krecs[ik];
	  const size_t     ir   = krec.i_read;
	  const unsigned   ib   = krec.i_base;
	  const BaseVec  & bv   = _bases[ir];
	  const unsigned   nb   = bv.size();
	  oss << setw(12) << ir << (krec.rc ? " rc " : " fw ")
	      << (_params.i_reads_print.count(ir) ? "* " : "  ");
	  for (unsigned i = 0; i < nb - _K - ibc; i++) 
	    oss << " ";
	  for (unsigned i = 0; i < ibc; i++) 
	    oss << hieroglyph(krec.rc ? 3u ^ bv[nb - i - 1] : bv[i]);
	  oss << "  ";
	  for (unsigned i = ibc; i < ibc + _K / 2; i++) 
	    oss << hieroglyph(krec.rc ? 3u ^ bv[nb - i - 1] : bv[i]);
	  oss << "  " << hieroglyph(krec.rc ? 3u ^ bv[nb - (ibc + _K / 2) - 1] : bv[ibc + _K / 2])
	      << "  ";
	  for (unsigned i = ibc + _K / 2 + 1; i < ibc + _K; i++) 
	    oss << hieroglyph(krec.rc ? 3u ^ bv[nb - i - 1] : bv[i]);
	  oss << "  ";
	  for (unsigned i = ibc + _K; i < nb; i++) 
	    oss << hieroglyph(krec.rc ? 3u ^ bv[nb - i - 1] : bv[i]);
	  oss << endl;
	}
	
	oss << endl;
	
	for (unsigned base = 0; base < 4; base++)
	  oss << "base " << hieroglyph(base) 
	      << ((col.base_winner == base) ? " winner" : "       ")
	      << " fix= " << col.fix_base[base]
	      << " q_sum= " << setw(8) << col.q_sum[base]
	      << " n_q_hi= " << setw(4) << col.n_q_high[base]
	      << endl;
	oss << endl;
      }
      else {
	oss << "quals: " << endl;
	for (unsigned q = 0; q < 42; q++) {
	  oss << "quals: " << setw(3) << q;
	  for (unsigned b = 0; b < 4; b++)
	    oss << " " << setw(10) << col.nq[b][q];
	  oss << endl;
	}
	oss << "quals: " << endl;
      }
      
      cout << oss.str() << endl;
    }
    
    return print_stack;
  }



  // ---- this looks for errors and issues recommendations
  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs, 
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t n_rows = i1k - i0k;

    const unsigned ibk_center = _K / 2; // the base at the kmer center
    const bool is_flank_palindrome = krecs[i0k].palindrome;


    // ---- kmer spectrum step

    if (true) {
      vec<size_t> freq(4, 0);
      for (size_t ik = i0k; ik != i1k; ik++) {
        const rec_type & krec  = krecs[ik];
        const unsigned   ib    = krec.i_base + ibk_center;
        unsigned         base  = _bases[krec.i_read][ib];
        if (krec.rc) base = 3u ^ base;
        if (is_flank_palindrome) {
          if      (base == 3) base = 0;
          else if (base == 2) base = 1;
        }
        freq[base]++;
      }
      for (unsigned base = 0; base != 4; base++) {
        const unsigned kf = freq[base];
        if (kf > 0) {
          if (_kspec_tmp.size() <= kf) _kspec_tmp.resize(kf + 1, 0);
          _kspec_tmp[kf]++;
        }
      }
    }


    // ---- center column error correction step

    if (n_rows >= _params.n_kmers_min &&
        ! is_flank_palindrome) { 

      // ---- build center base column statistics
      BaseColumn col(_params);
      
      for (size_t ik = i0k; ik != i1k; ik++) {
        const rec_type & krec = krecs[ik];
        const unsigned   ib   = krec.i_base + ibk_center;
        const unsigned   base = _bases[krec.i_read][ib];
        const unsigned   qual = _quals[krec.i_read][ib];
        col.add_base_qual(((krec.rc) ? 3u ^ base : base), qual); 
      }

      // ---- process base column statistics
      col.process();

      // ---- debugging stacks
      if (_params.i_reads_print.size() > 0) 
        stack_debug_print(krecs, i0k, i1k, col);


      // ---- buffer recommendations for these krecs
      
      const size_t n_reads = _bases.size();

      for (size_t ik = i0k; ik != i1k; ik++) {
        rec_type        rloc = krecs[ik];
        const BaseVec & bv   = _bases[rloc.i_read];
        const size_t    ib   = rloc.i_base + ibk_center;
        const size_t    iblk = i_block(rloc.i_read, n_reads);
        const unsigned  base = (rloc.rc) ? 3u ^ bv[ib] : bv[ib];
        if (col.fix_base[base]) {
          rloc.i_base = ib;
          rloc.base = (rloc.rc) ? 3u ^ col.base_winner : col.base_winner;
          _recommends[iblk].push_back(rloc);
        }

      } 
      
      // ---- if too many buffered recommends, flush them
      flush_recommends(5000, krecs.size());

    }
    //if (krecs[i0k].palindrome) 
    //  cout << "FOUND PALINDROME!!!!" << endl;


  }
  




  // ---- this issues a final flush to merge all recommendations and confirmations
  // interface function needed by naif_kmerize()
  void merge(PreCorrector & kernel_tmp,
	     const size_t i_parcel)
  {
    // ---- merge corrections

    kernel_tmp.flush_recommends(0, 0);

    // ---- merge kmer spectrum

    const size_t nkf = kernel_tmp._kspec_tmp.size();

    Locker lock(_lock);
    
    if (_p_kspec->size() < nkf) 
      _p_kspec->resize(nkf, 0);
    
    for (size_t kf = 0; kf != nkf; kf++)
      (*_p_kspec)[kf] += kernel_tmp._kspec_tmp[kf];
  }


};  // class PreCorrector
  







#endif
