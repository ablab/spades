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

#include "paths/PolymorphismRemoveCore.h"


static inline 
String Tag(String S = "PRC") { return Date() + " (" + S + "): "; } 













// ---- Starts from a kmer with several suffixes and follows one path of 1-1 kmers

// ---- kmap not a const because follow_kmers() needs to tag the visited kmers

void follow_kmers(KmerMap<KmerRec_t> & kmap,
		  const Kmer_t & kmer0,
		  const unsigned isuf,
		  BaseVec * bv_p,
		  float * mean_freq_p)
{
  bool verbose = false;
  const unsigned K = kmer0.K();


  // ---- krec is the walking rec; the original one has 2 suffixes; we'll pick one.
  Kmer_t kmerFW(kmer0);
  Kmer_t kmerRC(reverse_complement(kmerFW));
  bool fw = (kmerFW < kmerRC);

  bv_p->clear();
  *mean_freq_p = 0;
  size_t nk = 0;

  for (unsigned i = 1; i < K; i++) 
    bv_p->push_back(kmer0[i]); 


  // ---- get 0th krec
  KmerRec_t & krec0 = kmap(fw ? kmerFW : kmerRC);

  String patch = (isuf == 0) ? "a  " : "b  ";
  if (verbose) {
    cout << patch << hieroglyphs(kmerFW) 
         << " exists= " << (krec0.is_valid_kmer() ? 1 : 0)
         << (fw ? " fw" : " rc")
         << " pre= " << krec0.str_prefixes() 
         << " suf= " << krec0.str_suffixes() << endl;
  }

  // ---- pick a suffix 
  unsigned baseFW = krec0.suffix(isuf, fw);
  unsigned baseRC = 3u ^ baseFW;

  bool done = false;
  do {
    
    // ---- build next kmer
    kmerFW.push_right(baseFW);
    kmerRC.push_left (baseRC);
    fw = (kmerFW < kmerRC);

    // ---- get next krec
    KmerRec_t & krec = kmap(fw ? kmerFW : kmerRC);
    ForceAssert(krec.is_valid_kmer());

    if (!krec.tag() &&     // kmer not visited yet
        krec.n_prefixes() == 1 && 
        krec.n_suffixes() == 1) {

      // ---- add base to bv
      bv_p->push_back(baseFW); 
      *mean_freq_p += krec.freq();
      nk++;
      
      krec.tag_set(); // mark this 1-1 kmer as visited

      baseFW = krec.suffix(0, fw);
      baseRC = 3u ^ baseFW;

      if (verbose) {
        patch += " ";
        cout << patch << hieroglyphs(kmerFW)
	     << (fw ? " fw" : " rc")
             << " pre= " << krec.str_prefixes() 
             << " suf= " << krec.str_suffixes() << endl;
      }
    }
    else {
      done = true;
    }
  } while (!done);
           
  if (nk > 0) 
    *mean_freq_p /= float(nk);

  if (verbose)
    cout << endl;

}







void print_alt_debug(const unsigned K,
                     const String & type, 
                     const unsigned n_alts_total,
                     const BaseVec & bv_a, 
                     const BaseVec & bv_b)
{
  
  cout << "---------------------- " << type << " " << n_alts_total << endl;
  
  for (unsigned i = 0; i < K; i++) cout << hieroglyph(bv_a[i]);
  cout << "   ";
  for (unsigned i = 0; i < K; i++) cout << hieroglyph(bv_a[bv_a.size() - K - 1 + i]);
  cout << endl;
  
  for (unsigned i = 0; i < K; i++) cout << hieroglyph(bv_b[i]);
  cout << "   ";
  for (unsigned i = 0; i < K; i++) cout << hieroglyph(bv_b[bv_b.size() - K - 1 + i]);
  cout << endl;
}
























void Polymorphisms::_add_single(const BaseVec           & bv,
                                const float               kf,
                                BaseVecVec              * bvv_p,
                                vec<float>              * kfv_p,
                                vec<KmerIPoly<Kmer_t> > * i_poly_vec_p)
{
  const size_t i_poly = bvv_p->size(); 
  bvv_p->push_back(bv); 
  kfv_p->push_back(kf);
  
  SubKmers<BaseVec, Kmer_t> kmers_bv(_K, bv);
  const size_t nk = kmers_bv.n_kmers();
  for (unsigned ik = 0; ik != nk; ik++, kmers_bv.next())
    i_poly_vec_p->push_back(KmerIPoly<Kmer_t>(kmers_bv.canonical(), 
                                              i_poly, 
                                              kmers_bv.is_canonical_fw(), 
                                              ik, nk));
}

bool Polymorphisms::_add_poly(const BaseVec           & bv_a,
                              const BaseVec           & bv_b,
                              const float               kf_a,
                              const float               kf_b,
                              vec<KmerIPoly<Kmer_t> > * poly_a_vec_p,
                              vec<KmerIPoly<Kmer_t> > * poly_b_vec_p)
{
  // ---- polymorphism only if last K-1 bases are identical
  
  const unsigned nb_a = bv_a.size();
  const unsigned nb_b = bv_b.size();
  
  bool is_poly = false;
  
  if (nb_a >= _K + 2 && nb_b >= _K + 2)
    is_poly = true;
  
  if (is_poly) {
    
    // ---- verify that final K-1 bases are identical

    for (unsigned i = 1; i < _K; i++)
      if (bv_a[nb_a - i] != bv_b[nb_b - i]) 
        is_poly = false;
    
    if (is_poly) {     // ---- add A and B poly entries 

      _add_single(bv_a, kf_a, & _bvv_a, & _kfv_a, poly_a_vec_p);
      _add_single(bv_b, kf_b, & _bvv_b, & _kfv_b, poly_b_vec_p);
      
      _a_is_strong.push_back(true); // ---- pick A as strong by default

      _stats.add(bv_a.size(), bv_b.size());
    }
  }
  
  return is_poly;
}




// ---- kmap not a const because follow_kmers() needs to tag the visited.

void Polymorphisms::build(KmerMap<KmerRec_t> & kmap,
                          const size_t         kf_min,
                          const size_t         kf_max)
{
  bool verbose = false;
    
  vec<KmerIPoly<Kmer_t> > poly_a_vec;
  vec<KmerIPoly<Kmer_t> > poly_b_vec;
    
  BaseVec bv_a;
  BaseVec bv_b;
  float kf_a;
  float kf_b;
    
  cout << Tag() << "Looking for polymorphisms." << endl;
  const size_t nh = kmap.size_hash();
  for (size_t ih = 0; ih != nh; dots_pct(ih++, nh)) {
      
    const KmerRec_t & krec0 = kmap[ih];
    if (krec0.is_valid_kmer() && 
        (krec0.n_suffixes() == 2 ||
         krec0.n_prefixes() == 2)) {
        
      // ---- follow kmers forward

      if (krec0.n_suffixes() == 2) {
          
        Kmer_t kmer0 = krec0;
        follow_kmers(kmap, kmer0, 0, &bv_a, &kf_a);
        follow_kmers(kmap, kmer0, 1, &bv_b, &kf_b);
          
        if (_kmer_frequency_valid(kf_a + kf_b, kf_min, kf_max))
          _add_poly(bv_a, bv_b, kf_a, kf_b, & poly_a_vec, & poly_b_vec);
      }

      // ---- follow kmers backward
        
      if (krec0.n_prefixes() == 2) { 
          
        Kmer_t kmer0 = reverse_complement(krec0);
        follow_kmers(kmap, kmer0, 0, &bv_a, &kf_a);
        follow_kmers(kmap, kmer0, 1, &bv_b, &kf_b);
          
        if (_kmer_frequency_valid(kf_a + kf_b, kf_min, kf_max))
          _add_poly(bv_a, bv_b, kf_a, kf_b, & poly_a_vec, & poly_b_vec);
      }
    }
  }
    
  // ---- Building maps from vectors

  cout << Tag() << "Building A allele kmer map." << endl;
  _poly_a.from_kmer_vec(poly_a_vec, 1.5);
  cout << Tag() << "mapA.size= " << _poly_a.num_recs() << endl;

  cout << Tag() << "Building B allele kmer map." << endl;
  _poly_b.from_kmer_vec(poly_b_vec, 1.5);
  cout << Tag() << "mapB.size= " << _poly_b.num_recs() << endl;

}












void polymorphisms_find_parallel(const unsigned K,
                                 const BaseVecVec & bvv, 
                                 Polymorphisms * polys_p,
                                 KmerSpectrum * kspec_p,
                                 const unsigned VERBOSITY,
                                 const unsigned NUM_THREADS,
                                 const size_t mean_mem_ceil)
{
  ForceAssertEq(K, kspec_p->K());
  
  if (!(K & 1u)) {
    cout << Tag() << "K = " << K << " is not odd. Changing K to " << (K + 1) << "." << endl;
    exit(1);
  }
  const size_t nbv = bvv.size();
  size_t len_read = 0;
  for (size_t ibv = 0; ibv < nbv; ibv++) 
    len_read += bvv[ibv].size();
  len_read /= nbv;


  {
    // ---- build the kmer hash map

    cout << Tag() << "Building " << K << "-kmer frequency/affixes hash table." << endl;
    Validator validator_kf(1, 0);
    const double hash_table_ratio = 1.5;
    KmerMap<KmerRec_t> kmap;
    
    kmer_freq_affixes_map_build_parallel(K, bvv, validator_kf, hash_table_ratio, 
                                         & kmap,
                                         VERBOSITY, NUM_THREADS, mean_mem_ceil);


    // ---- get kmer spectrum and determine upper limit on CN=1 kmer frequencies
    
    kmer_spectrum_from_kmer_freq_map(kmap, kspec_p);


    const unsigned ploidy = 2;
    kspec_p->analyze(ploidy, len_read);
    ForceAssertGt(kspec_p->genome_size(), 0u);
    const size_t kf_min = 0;
    const size_t kf_max = 1.5 * kspec_p->kf_max2(); 
    cout << Tag() << setw(12) << kspec_p->kf_max2() << "  kmer frequency of mode." << endl;
    cout << Tag() << setw(12) << kf_max << "  kmer frequency upper cutoff for CN=1." << endl;
  

    // ---- find ambiguities
    
    cout << Tag() << "Finding polymorphisms." << endl;
    polys_p->build(kmap, kf_min, kf_max);
    

  }
}




















// ---- returns the approximate kmer index in the alternative path

size_t i_kmer_alt(const BaseVec& bv_alt,
		  const unsigned K, 
		  const unsigned ik,
		  const unsigned nk) 
{
  return (bv_alt.size() - K) * ik / (nk - 1);
}




void base_vec_insert_alternative(const unsigned            K,
                                 const unsigned            i0k_bv0,
                                 const unsigned            i1k_bv0,
                                 const unsigned            nb_bv0,
                                 const bool                is_fw_kbv0,
                                 const KmerIPoly<Kmer_t> & k0_del,
                                 const BaseVec           & bv_ins,
                                 BaseVec                 * bv_p)
{
  bool verbose = false;

  const unsigned nk_bv_del = i1k_bv0 - i0k_bv0 + 1;
  const unsigned nk_del    = k0_del.n_kmers();

  if (verbose) {
    cout << "base_vec_insert_alternative():" << endl;
    cout << " i0k_bv0= " << i0k_bv0 << endl;
    cout << " i1k_bv0= " << i1k_bv0 << endl;
    cout << " nb_bv0= " << nb_bv0 << endl;
    BaseVec bv_ins_rc = bv_ins;
    bv_ins_rc.ReverseComplement();
    cout << "bv_ins_fw = " << hieroglyphs(bv_ins) << endl;
    cout << "bv_ins_rc = " << hieroglyphs(bv_ins_rc) << endl;
    cout << " nk_del= " << nk_del << endl;
    cout << " nk_bv_del= " << nk_bv_del << endl;
    cout << "before:     " << hieroglyphs(*bv_p) << endl;
  }

  // ---- if the last kmer to delete is the last kmer in the read 
  //      we will need to add the last K-1 bases of the alternative
  const bool not_end = (i1k_bv0 < nb_bv0 - K); 

  // ---- the ambiguity is FW-FW or RC-RC

  if (is_fw_kbv0 == k0_del.is_fw()) {  
    const unsigned i0k_del = k0_del.i_kmer();
    const unsigned i1k_del = i0k_del + nk_bv_del - 1;

    if (verbose) {
      cout << " fw-fw" << endl;
      cout << " i0k_del= " << i0k_del << endl;
      cout << " i1k_del= " << i1k_del << endl;
      cout << " bv_ins.size()= " << bv_ins.size() << endl;
    }
    if (i1k_del >= nk_del) cout << "oopsie! fw" << endl;
    else {
      const unsigned i0b_ins = i_kmer_alt(bv_ins, K, i0k_del, nk_del);
      const unsigned i1b_ins = i_kmer_alt(bv_ins, K, i1k_del, nk_del) + 
        ((not_end) ? 0 : K - 1);        
      const unsigned nb_ins = i1b_ins + 1 - i0b_ins;

      if (verbose) {
        cout << " i0b_ins= " << i0b_ins << endl;
        cout << " i1b_ins= " << i1b_ins << endl;
        cout << " nb_ins=  " << nb_ins << endl;
        cout << " not_end= " << not_end << endl;
      }
      for (unsigned i = 0; i < nb_ins; i++) {
        bv_p->push_back(bv_ins[i0b_ins + i]);
      }
      
    }
  } 
  
  // ---- the ambiguity is FW-RC or RC-FW
  
  else {

    const unsigned i1k_del = k0_del.i_kmer();
    const unsigned i0k_del = i1k_del - nk_bv_del + 1;
    
    if (verbose) {
      cout << " fw-rc" << endl;
      cout << " i0k_del= " << i0k_del << endl;
      cout << " i1k_del= " << i1k_del << endl;
      cout << " bv_ins.size()= " << bv_ins.size() << endl;
    }
    if (i0k_del >= nk_del) cout << "oopsie! rc" << endl;
    else {
      const unsigned i0b_ins = i_kmer_alt(bv_ins, K, i0k_del, nk_del)  + 
        ((not_end) ? K - 1 : 0);
      const unsigned i1b_ins = i_kmer_alt(bv_ins, K, i1k_del, nk_del) + K - 1;
      const unsigned nb_ins = i1b_ins + 1 - i0b_ins;
      
      if (verbose) {
        cout << " i0b_ins= " << i0b_ins << endl;
        cout << " i1b_ins= " << i1b_ins << endl;
        cout << " nb_ins=  " << nb_ins << endl;
        cout << " not_end= " << not_end << endl;
      }
      for (unsigned i = 0; i < nb_ins; i++) {
        bv_p->push_back(3u ^ bv_ins[i1b_ins - i]);
      }
    }
  }
  if (verbose)
    cout << "after:     " << hieroglyphs(*bv_p) << endl;

}















  
// ---- Takes a table of ambiguity conversions and fixes a read

bool base_vec_fix(const Polymorphisms & polys,
                  BaseVec             * bv_p,
                  QualNibbleVec       * qv_p = 0)
{
  const unsigned K = polys.K();
  const unsigned nb_bv0 = bv_p->size();

  bool fixed = false; // the return value, to be updated

  if (nb_bv0 >= K) {
    BaseVec       bv_new;
    QualNibbleVec qv_new;
    
    // ---- compute mean read quality 
    //      qual vecs get updated too, and using the mean is the easiest way, 
    //      albeit quite meaningless
    
    unsigned q_mean = 0;
    if (qv_p) {
      for (unsigned i = 0; i < nb_bv0; i++) q_mean += (*qv_p)[i];
      q_mean /= nb_bv0;
    }
    
    // ---- cycle through all the kmers in the read
   
    const unsigned nk_bv0 = nb_bv0 - K + 1;
    bool is_good = true;
    bool to_correct = false;
    bool fill_end = true;
    
    // ---- cycle through all bv0 kmers
    
    unsigned ik_bv0 = 0;
    SubKmers<BaseVec, Kmer_t> kmer_bv0(K, *bv_p);
    while (kmer_bv0.not_done()) {

      // ---- get the kmer record with deletion information
      KmerIPoly<Kmer_t> k0_del = polys.kmer_poly_weak(kmer_bv0.canonical());

      if (k0_del.is_valid_kmer()) {  // there is an alternative! delete and replace
        to_correct = true;

        const unsigned i0k_bv0 = ik_bv0; // index of the first kmer to delete
        const bool is_fw_kbv0 = kmer_bv0.is_canonical_fw(); // fw info on first kmer
        
        KmerIPoly<Kmer_t> k_del = k0_del;
        const unsigned  i_poly0 = k0_del.i_poly();
        unsigned        i_poly  = i_poly0;

        // ---- find extension of kmer region to be deleted

        bool done = false;
        while (!done) {
          kmer_bv0.next();
          ik_bv0++;
          if (kmer_bv0.not_done()) {
            k_del = polys.kmer_poly_weak(kmer_bv0.canonical());

            if (k_del.is_valid_kmer())  i_poly = k_del.i_poly();
            else                        done   = true;
          } 
          else 
            done = true;
        }
        const unsigned i1k_bv0 = ik_bv0 - 1; // index of the last kmer to delete


	// ---- if all is well insert alternative into basevector

	if (i1k_bv0 >= i0k_bv0 && i_poly == i_poly0) {
          
          base_vec_insert_alternative(K, i0k_bv0, i1k_bv0, nb_bv0, is_fw_kbv0,
                                      k0_del, polys.base_vec_strong(i_poly), & bv_new);


          // -- add read mean quality score to fill unknown qualities
          //    as mentioned above this is quite meaningless, but very easy

          const unsigned nb_ins = bv_new.size() - qv_new.size(); 
          for (unsigned i = 0; i < nb_ins; i++) 
            qv_new.push_back(q_mean);

          if (i1k_bv0 == nk_bv0 - 1)
            fill_end = false;
        }
        else {
          is_good = false;
        }
      }
      
      if (ik_bv0 < nk_bv0) {
        bv_new.push_back((*bv_p)[ik_bv0]);
        if (qv_p)
          qv_new.push_back((*qv_p)[ik_bv0]);
      }
     
      ik_bv0++;
      kmer_bv0.next();
    } 


    // ---- only update base vector if all went well
    
    if (to_correct && is_good) {

      // ---- fill end if necessary

      if (fill_end) {
        for (unsigned i = 1; i < K; i++) {
          const unsigned ib = nb_bv0 - K + i;
          bv_new.push_back((*bv_p)[ib]);
          if (qv_p)
            qv_new.push_back((*qv_p)[ib]);
        }
      }

      // ---- update the base vector

      *bv_p = bv_new;
      if (qv_p) {
        ForceAssertEq(bv_new.size(), qv_new.size());
        *qv_p = qv_new;
      }

      fixed = true;
    }



  }

  return fixed;
  
}
  







class PolyRemoverProc
{
  const Polymorphisms & _polys;
  BaseVecVec          & _bvv;
  QualNibbleVecVec    & _qvv;
  size_t              & _n_fixed;
  size_t              & _jv;
  LockedData          & _lock;
  const size_t          _n_threads;
  
public:
  PolyRemoverProc(const Polymorphisms & polys,
                  BaseVecVec          * bvv_p,
                  QualNibbleVecVec    * qvv_p,
                  size_t              * n_fixed_p,
                  size_t              * jv_p,
                  LockedData          & lock,
                  const unsigned        n_threads)
    : _polys(polys),
      _bvv(*bvv_p),
      _qvv(*qvv_p),
      _n_fixed(*n_fixed_p),
      _jv(*jv_p),
      _lock(lock),      
      _n_threads(n_threads)
  {}
  
  void operator()(const size_t i_thread) 
  {
    const size_t nv = _bvv.size();
    const size_t i0v = ( i_thread      * nv) / _n_threads;
    const size_t i1v = ((i_thread + 1) * nv) / _n_threads;

    // ---- Remove polymorphisms from the reads
    BaseVec bv;
    QualNibbleVec qv;
    for (size_t iv = i0v; iv < i1v; iv++) {
      bv = _bvv[iv];
      qv = _qvv[iv];

      const bool fixed = base_vec_fix(_polys, & bv, & qv);

      Locker locker(_lock);
      if (fixed) {
        _n_fixed++;
        _bvv[iv] = bv;
        _qvv[iv] = qv;
      }

      dots_pct(_jv++, nv);
    }

    /*
    Locker lock(_lock);
    cout << " i_thread= " << i_thread 
         << " nv= " << nv 
         << " i0v= " << i0v 
         << " i1v= " << i1v 
         << " n_threads= " << _n_threads 
         << " _n_fixed= " << _n_fixed 
         << endl;
    */
  }
  

};




size_t polymorphisms_remove_parallel(const Polymorphisms & polys,
                                     BaseVecVec          * bvv_p,
                                     QualNibbleVecVec    * qvv_p,
                                     const unsigned        VERBOSITY,
                                     const unsigned        NUM_THREADS)
{
  size_t n_fixed = 0;
  size_t jv = 0;
  LockedData lock;
  {
    PolyRemoverProc remover(polys, bvv_p, qvv_p, &n_fixed, &jv, lock, NUM_THREADS);
    Worklist<size_t, PolyRemoverProc> worklist(remover, NUM_THREADS - 1);
    
    for (size_t i = 0; i != NUM_THREADS; i++) 
      worklist.add(i);
  }
  return n_fixed;
}


