///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Filipe Ribeiro - Feb 2011 - <crdhelp@broadinstitute.org>
//


#ifndef KMERS__NAIF_KMER__KERNEL_ERROR_FINDER_H
#define KMERS__NAIF_KMER__KERNEL_ERROR_FINDER_H

#include "kmers/naif_kmer/LockedBlocks.h"
#include "kmers/naif_kmer/Kmers.h"








class ColumnQualityScores 
{
public:
  const EF_Params & _params;
  unsigned QSsum[4], nQ[4], QBcount[4];
  unsigned highest_base;
  bool recommended[4];
  bool confirmed[4];

  ColumnQualityScores(const EF_Params & params)
    : _params(params) 
  {}


  void reset() 
  {
    highest_base = 0;
    QSsum      [0] = QSsum      [1] = QSsum      [2] = QSsum      [3] = 0;
    nQ         [0] = nQ         [1] = nQ         [2] = nQ         [3] = 0;
    QBcount    [0] = QBcount    [1] = QBcount    [2] = QBcount    [3] = 0;
    confirmed  [0] = confirmed  [1] = confirmed  [2] = confirmed  [3] = false;
    recommended[0] = recommended[1] = recommended[2] = recommended[3] = false;
  }

  void add_base_qual(const unsigned base, const unsigned qual) 
  {
    QSsum[base] += qual;
    nQ[base]++;
    if (qual >= unsigned(_params.min_Q_to_support)) 
      QBcount[base]++;
  }
    
  void process(const bool under_kmer, bool * branched)
  {
    const unsigned col_size = nQ[0] + nQ[1] + nQ[2] + nQ[3];
    if (col_size >= unsigned(_params.min_basestack_depth)) {

      highest_base = 0;
      if (QSsum[1] > QSsum[highest_base]) highest_base = 1; 
      if (QSsum[2] > QSsum[highest_base]) highest_base = 2; 
      if (QSsum[3] > QSsum[highest_base]) highest_base = 3; 
	
      const unsigned highest_QSsum = QSsum[highest_base];
	
      if (highest_QSsum >= unsigned(_params.min_QSsum_to_win)) {

	bool supported[4] = { false, false, false, false };
	bool found_main_support = false;
	bool found_alt_support  = false;
	  
	// Experimental: if min_qss_to_support is set, use it instead of
	// n_bases_to_support.  --bruce
	if (_params.min_QSsum_to_support > 0) {
	  for (unsigned base = 0; base < 4; base++) {
	    if (QSsum[base] >= unsigned(_params.min_QSsum_to_support)) {
	      supported[base] = true;
	      if (base == highest_base) found_main_support = true;
	      else                      found_alt_support  = true;
	    }
	  }
	} 
	else {
	  for (unsigned base = 0; base < 4; base++) {
	    if (QBcount[base] >= unsigned(_params.min_n_bases_to_support)) {
	      supported[base] = true;
	      if (base == highest_base) found_main_support = true;
	      else                      found_alt_support  = true;
	    }
	  }
	}
	  
	*branched = found_alt_support;
	  
	// ---------------- recommendations ---------------------
	// Make recommendations for replacing bases with other bases.
	// Criteria for recommendations are as follows...
	// Each base-location has an associated quality-score-sum ("QSsum").
	// We replace an A with a C iff the following criteria are ALL true:
	// -- QSsum(C) is the highest of the four QSsum's
	// -- QSsum(A) / QSsum(C) is lower than max_QSsum_ratio_to_correct
	// -- A is not supported [see definition of "supported" above]

	double QSsum_ratio[4];
	double full_max_QSsum_ratio = 0;
	for (unsigned base = 0; base < 4; base++) {
	  QSsum_ratio[base] = double(QSsum[base]) / double(highest_QSsum);
	  if (base != highest_base && full_max_QSsum_ratio < QSsum_ratio[base]) 
	    full_max_QSsum_ratio = QSsum_ratio[base];
	}
	  
	double max_QSsum_ratio = 0;
	if (full_max_QSsum_ratio <= _params.max_QSsum_ratio_to_correct2) {
	  for (unsigned base = 0; base < 4; base++) {
	    if (base != highest_base && !supported[base] && QSsum[base] != 0) {
	      if (max_QSsum_ratio < QSsum_ratio[base]) 
		max_QSsum_ratio = QSsum_ratio[base];
		
	      // means that QSsum[base] is much smaller than highest_QSsum
              if (QSsum_ratio[base] <= _params.max_QSsum_ratio_to_correct)
                recommended[base] = true;
	    }
	  }
	}
	

	// ---------------- confirmations ---------------------
	// Auto-confirm.  A large enough quality score sum gets automatic confirmation.
	// See next section for 'standard' confirmation.
	  
	if (_params.auto_confirm >= 0)
	  for (unsigned base = 0; base < 4; base++)
	    if (QSsum[base] >= unsigned(_params.auto_confirm))
	      confirmed[base] = true;
	


	// Mark some base-locations as confirmed.  We do this (say, with C) iff:
	// -- ALL bases in this column are C, or are recommended for change to C
	//    (this requires all bases other than C to be non-supported)
	// -- QSsum(C) is at least min_QSsum_to_confirm (if set), or
	//    QBcount(C) is at least min_n_bases_to_confirm if QSsum_to_confirm is not set
	// -- C is supported, or this column is NOT part of the key kmer
	// A confirmed base-location will not be error-corrected even if it is
	// marked in another stack as in need of correction.
	  
	if (!confirmed[highest_base]            && 
	    !found_alt_support                  && 
	    (!under_kmer || found_main_support) &&
	    full_max_QSsum_ratio       <= _params.max_QSsum_ratio_to_confirm &&
	    max_QSsum_ratio            <= _params.max_QSsum_ratio_to_correct &&
	    ((0                         < _params.min_QSsum_to_confirm && 
              highest_QSsum            >= unsigned(_params.min_QSsum_to_confirm)) ||
             QBcount[highest_base]     >= unsigned(_params.min_n_bases_to_confirm))) {

	  confirmed[highest_base] = true;

	}

      } // if highest_QSsum >= _params.min_QSsum_to_win) 

    } // if (col_size >= _MIN_COLUMN_SIZE)

  } // void recommendations_and_confirmations()

  void print() 
  {
    cout << "  QSsum: ";
    for (size_t base = 0; base != 4; base++) cout << " " << setw(5) << QSsum[base];
    cout << "  nQ: ";
    for (size_t base = 0; base != 4; base++) cout << " " << setw(2) << nQ[base];
    cout << "  QBcount: ";
    for (size_t base = 0; base != 4; base++) cout << " " << setw(2) << QBcount[base];
    cout << "  Win: " << highest_base;
    cout << "  recom: ";
    for (size_t base = 0; base != 4; base++) cout << unsigned(recommended[base]);
    cout << "  conf: ";
    for (size_t base = 0; base != 4; base++) cout << unsigned(confirmed[base]);
    cout << endl;
  }


};  // class ColumnQualityScores








// -------- SUMMARIZER_t class --------

template<class KMER_t, class QVV_t>
class ErrorFinder
{
  const size_t             _n_threads;
  const size_t             _n_blocks;
  LockedBlocks             _blocks;
  LockedBlocks           * _p_blocks;
  vec< vec<BVLocBase> >    _recommends;
  vec< vec<BVLocBase> >    _confirms;

  const EF_Params        & _params;
  
  const size_t             _K;

  const BaseVecVec       & _bases;
  const QVV_t            & _quals;
  BaseVecVec             & _bases_new;
  BitVecVec              & _bases_locked;

  const size_t             _nb_max;  // the largest base vec size 
  const size_t             _n_cols;  // n_cols = K + 2 (nb_max - K)

  vec<ColumnQualityScores> _qs;
  
  size_t test_recoms_counts;
  size_t test_confs_counts;
  
  typedef typename QVV_t::value_type QV_t;

public:

  typedef KmerBVLocBase<KMER_t>       rec_type;

  ErrorFinder(const EF_Params        & params,
              const size_t             K,
              const BaseVecVec       & bases,
              const QVV_t            & quals,
              BaseVecVec             * p_bases_new,
              BitVecVec              * p_bases_locked,
              const size_t             n_threads) 
    : _n_threads(n_threads),
      _n_blocks(2 * n_threads),
      _blocks(_n_blocks),
      _p_blocks(&_blocks),
      _recommends(_n_blocks),
      _confirms(_n_blocks),
      _params(params),
      _K(K),
      _bases(bases),
      _quals(quals),
      _bases_new(*p_bases_new),
      _bases_locked(*p_bases_locked),
      _nb_max(bases.MaxSize()),
      _n_cols(2 * _nb_max - _K),
      _qs(_n_cols, ColumnQualityScores(params))
  {
    test_recoms_counts = test_confs_counts = 0;
  }

  ~ErrorFinder()
  {
    if (0) {
      cout << "test_confs_counts = " << test_confs_counts << endl;
      cout << "test_recoms_counts = " << test_recoms_counts << endl;
    }
  }
      

  
  // copy constructor for temporary kernels
  explicit ErrorFinder(const ErrorFinder & that)
    : _n_threads(that._n_threads),
      _n_blocks(that._n_blocks),
      _blocks(),                   // no local blocks
      _p_blocks(that._p_blocks),   // get the blocks from 'that'
      _recommends(_n_blocks),
      _confirms(_n_blocks),
      _params(that._params),
      _K(that._K),
      _bases(that._bases),
      _quals(that._quals),
      _bases_new(that._bases_new),
      _bases_locked(that._bases_locked),
      _nb_max(that._nb_max),
      _n_cols(that._n_cols),
      _qs(that._qs)
  {
    //cout << "clone: qs.size = " << _qs.size() << endl; 
    //cout << "clone: n_blocks = " << _n_blocks << endl; 
    //cout << "clone: blocks.size() = " << _p_blocks->size() << endl; 
  }
 
  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bases; }
  
  

  String todo_str(const size_t n, const std::set<size_t> & todo) const
  {
    String s = "";
    for (size_t i = 0; i != n; i++) 
      s += (todo.count(i)) ? "|" : "_";
    return s;
  }
    


  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bases[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (p_buf->in_one_parcel(kmer)) {
        kmer.i_read     = ibv;
        kmer.i_base     = kmer_cur.index_start();
        kmer.rc         = kmer_cur.is_canonical_rc();
        kmer.palindrome = kmer_cur.is_palindrome();
        p_buf->add(kmer);
      }
      kmer_cur.next();          
    }
  } 



  void flush_recommends_confirms(const size_t max_buf_size, 
                                 const size_t parcel_size) // parcel_size is for debug
  {
    bool debug = false;//(max_buf_size == 0);

    std::set<size_t> blocks_todo;
    for (size_t i = 0; i < _n_blocks; i++)
      if (_recommends[i].size() > max_buf_size || 
          _confirms[i].size() > max_buf_size)  // buffers too big? flush them!
        blocks_todo.insert(i);

    
    
    while (blocks_todo.size() > 0) {


      if (debug)
        cout << Date() << ":[" << setw(10) << parcel_size << "] need lock for           " 
             << " " << todo_str(_p_blocks->size(), blocks_todo)
             <<  " blocks: " << _p_blocks->str() << endl; 

      // ---- lock a subset of the data to update

      const size_t iblk = _p_blocks->lock_some_block(blocks_todo);


      if (debug)
        cout << Date() << ":[" << setw(10) << parcel_size << "] locked block " << setw(2) << iblk << " from    " 
             << " " << todo_str(_p_blocks->size(), blocks_todo)
             <<  " blocks: " << _p_blocks->str() << endl; 


      

      // ---- apply confirmations to '_bases_locked'

      vec<BVLocBase> & confs = _confirms[iblk];
      const size_t n_confs = confs.size();
      for (size_t i = 0; i < n_confs; i++) {
        const size_t i_read = confs[i].i_read;
        const size_t i_base = confs[i].i_base;
        _bases_locked[i_read].Set(i_base, true);

        // revert previous corrections, if any
        _bases_new   [i_read].Set(i_base, _bases[i_read][i_base]); 
      }
      confs.clear();




      // ---- apply recommendations to '_bases_new'

      vec<BVLocBase> & recoms = _recommends[iblk];
      const size_t n_recoms = recoms.size();
      for (size_t i = 0; i < n_recoms; i++) { 
        const size_t i_read = recoms[i].i_read;
        const size_t i_base = recoms[i].i_base;

        if (!_bases_locked[i_read][i_base] && _bases[i_read][i_base] != recoms[i].base) {
          const unsigned char base_new   = _bases_new[i_read][i_base]; 
          const unsigned char base_orig  = _bases    [i_read][i_base]; 
          const unsigned char base_newer = recoms[i].base; 
          

          if (base_new == base_orig) { // no previous correction
            _bases_new[i_read].Set(i_base, base_newer);
          }
          else if (base_new != base_newer) { 
            // two different recommended corrections 
            //    => ignore corrections and set 'bases_locked'
            // NOTE: let base_new and base_orig stay different so that 
            //       at the end we can count the number of conflicts
            _bases_locked[i_read].Set(i_base, true);
          }
        }
      }
      recoms.clear();



      if (debug)
        cout << Date() << ":[" << setw(10) << parcel_size << "] unlocking block " << setw(2) << iblk << " from " 
             << " " << todo_str(_p_blocks->size(), blocks_todo)
             <<  " blocks: " << _p_blocks->str() << endl; 

      // ---- unlock block and remove it from to-do list

      _p_blocks->unlock_block(iblk);
      blocks_todo.erase(iblk);


      if (debug)
        cout << Date() << ":[" << setw(10) << parcel_size << "] unlocked block  " << setw(2) << iblk << " from " 
             << " " << todo_str(_p_blocks->size(), blocks_todo)
             <<  " blocks: " << _p_blocks->str() << endl; 

    } // while (blocks_todo.size() > 0) 

  } // void flush_recommends_confirms(const size_t max_buf_size)




  size_t i_block(const size_t i, const size_t n) { return (i * _n_blocks) / n;  }

  unsigned i0_col_fw(const unsigned ibk) { return _nb_max - _K - ibk; }
  unsigned i0_col_rc(const unsigned ibk) { return _nb_max - 1 + ibk; }
  
 


  void build_column_stats(const vec<rec_type> & parcel, 
			  const size_t i0k,
			  const size_t i1k)
  {
    const size_t i0k_col = _nb_max - _K;
    const size_t i1k_col = _nb_max;

    // ---- reset all _qs[i_col]
    for (size_t i_col = 0; i_col != _n_cols; i_col++)
      _qs[i_col].reset(); 
    
    // ---- go through all rows in stack and build _qs[i_col]
    for (size_t ik = i0k; ik != i1k; ik++) {
      const rec_type & kloc = parcel[ik];
      const BaseVec  & bv   = _bases[kloc.i_read];
      const QV_t     & qv   = _quals[kloc.i_read];
      const size_t     nb   = bv.size();
      const size_t     ibk  = kloc.i_base;
      
      if (kloc.rc) {
	for (size_t ib = 0, i_col = i0_col_rc(ibk); ib != nb; ib++, i_col--)
	  //if (i_col < i0k_col || i_col >= i1k_col)
          _qs[i_col].add_base_qual(3u - bv[ib], qv[ib]); // complement of base
      }
      else {
	for (size_t ib = 0, i_col = i0_col_fw(ibk); ib != nb; ib++, i_col++)
	  //if (i_col < i0k_col || i_col >= i1k_col)
          _qs[i_col].add_base_qual(bv[ib], qv[ib]);
      }
    }
  }



  // ---- obtain all recomendations and confirmations
  void process_column_stats()
  {
    // ---- obtain confirmations on the kmer (or not)
    if (!_params.skip_under) {
      bool branched = false;
      const bool under_kmer = true;
      const size_t i0_col = _nb_max - _K;
      const size_t i1_col = _nb_max;
      for (size_t i_col = i0_col; i_col < i1_col; i_col++)
	_qs[i_col].process(under_kmer, &branched);
    }

    // ---- obtain recommendations and confirmations to the left of the kmer
    {
      bool branched = false;
      const bool under_kmer = false;
      const size_t j0_col = 0;
      const size_t j1_col = _nb_max - _K;
      const size_t i0_col = _nb_max - _K - 1;
      for (size_t j_col = j0_col; j_col < j1_col; j_col++)
	if (!branched || _params.do_branches)
	  _qs[i0_col - j_col].process(under_kmer, &branched);
          
    }

    // ---- obtain recommendations and confirmations to the right of the kmer
    {
      bool branched = false;
      const bool under_kmer = false;
      const size_t i0_col = _nb_max;
      const size_t i1_col = _n_cols;
      for (size_t i_col = i0_col; i_col < i1_col; i_col++) 
	if (!branched || _params.do_branches)
	  _qs[i_col].process(under_kmer, &branched);
    }
  }



  void confirm_and_recommend(const vec<rec_type> & parcel, 
			     const size_t i0k,
			     const size_t i1k)
  {
    // ---- buffer recommendations and confirmations for this parcel/stack

    const size_t n_reads = _bases.size();
    for (size_t ik = i0k; ik != i1k; ik++) {
      BVLocBase       rloc = parcel[ik];
      const BaseVec & bv   = _bases[rloc.i_read];
      const size_t    nb   = bv.size();
      const size_t    ibk  = rloc.i_base;
      const size_t    iblk = i_block(rloc.i_read, n_reads);

      if (rloc.rc) {
	for (size_t ib = 0, i_col = i0_col_rc(ibk); ib != nb; ib++, i_col--) {
	  const unsigned base = 3u - bv[ib];
	  if (_qs[i_col].confirmed[base]) {
	    rloc.i_base = ib;
	    _confirms[iblk].push_back(rloc);
	  } 
	  if (_qs[i_col].recommended[base]) {
	    rloc.i_base = ib;
	    rloc.base = 3u - _qs[i_col].highest_base;
	    _recommends[iblk].push_back(rloc);
	  }
	}
      }
      else {
	for (size_t ib = 0, i_col = i0_col_fw(ibk); ib != nb; ib++, i_col++) {
	  const unsigned base = bv[ib];
	  if (_qs[i_col].confirmed[base]) {
	    rloc.i_base = ib;
	    _confirms[iblk].push_back(rloc);
	  } 
	  if (_qs[i_col].recommended[base]) {
	    rloc.i_base = ib;
	    rloc.base = _qs[i_col].highest_base;
	    _recommends[iblk].push_back(rloc);
	  }
	}
      }
    } // for (size_t ik = i0k; ik != i1k; ik++)

  }


  bool stack_debug_print(const vec<rec_type> & krecs, 
                         const size_t i0k,
                         const size_t i1k)
  {

    bool print_stack = false;
    for (size_t ik = i0k; ik != i1k && !print_stack; ik++)
      if (_params.i_reads_print.count(krecs[ik].i_read))
        print_stack = true;

    if (print_stack) {
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
      stringstream oss;
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
        oss << "   ";
        for (unsigned i = ibc; i < ibc + _K; i++) 
          oss << hieroglyph(krec.rc ? 3u ^ bv[nb - i - 1] : bv[i]);
        oss << "   ";
        for (unsigned i = ibc + _K; i < nb; i++) 
          oss << hieroglyph(krec.rc ? 3u ^ bv[nb - i - 1] : bv[i]);
        oss << endl;
      }
      
      oss << endl;
      /*
      for (unsigned base = 0; base < 4; base++)
        oss << "base " << hieroglyph(base) 
            << ((col.base_winner == base) ? " winner" : "       ")
            << " fix= " << col.fix_base[base]
            << " q_sum= " << setw(8) << col.q_sum[base]
            << " n_q_hi= " << setw(4) << col.n_q_high[base]
            << endl;
      oss << endl;
      */
      
      cout << oss.str() << endl;
    }
    /*
    for (size_t i_col = 0; i_col != _n_cols; i_col++) {
      cout << parcel[ikt].i_read << "[ " << parcel[ikt].i_base << " " << ((parcel[ikt].rc) ? "RC ]:" : "FW ]:");
      _qs[i_col].print();
    }
    */
    return print_stack;
  }



  // ---- this looks for errors and issues recommendations
  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & parcel, 
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t n_rows = i1k - i0k;
    
    if (n_rows >= unsigned(_params.min_readstack_depth) &&
        (_params.max_readstack_depth == 0 || n_rows <= unsigned(_params.max_readstack_depth)) &&
        n_rows >= unsigned(_params.min_basestack_depth) &&
	(!parcel[i0k].palindrome || _params.do_palindromes)) {

      build_column_stats(parcel, i0k, i1k);
      
      process_column_stats();

      confirm_and_recommend(parcel, i0k, i1k);

      if (_params.i_reads_print.size() > 0)
        stack_debug_print(parcel, i0k, i1k);
  
      flush_recommends_confirms(5000, parcel.size());
    }
    if (false && parcel[i0k].palindrome) {
      cout << "FOUND PALINDROME!!!!" << endl;
    }

  }
  




  // ---- this issues a final flush to merge all recommendations and confirmations
  // interface function needed by naif_kmerize()
  void merge(ErrorFinder & tmp_corrector,
	     const size_t i_parcel)
  {
    tmp_corrector.flush_recommends_confirms(0, 0);
  }


};  // class ErrorFinder
  







#endif
