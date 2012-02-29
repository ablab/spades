///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#ifndef KMERS__NAIF_KMER__KMER_MAP_H
#define KMERS__NAIF_KMER__KMER_MAP_H

#include "feudal/BinaryStream.h"

/*
template<class HASH_TABLE_t, class VAL_t>
kmer_hash_table_build_parallel(const unsigned K,
			       const BaseVecVec & bvv,
			       const VAL_t & validator,
			       const double size_ratio,
			       HASH_TABLE_t * hash_p,
			       const unsigned n_threads,
			       const unsigned verbosity)
{
  typedef typename HASH_TABLE_t::rec_type KmerRec_t;
  vec<KmerRec_t> kvec;
  KernelKmerStorer<KmerRec_t, VAL_t > storer(bvv, K, &kvec, &validator);
  naif_kmerize(&storer, n_threads, verbosity);

  sort(kvec.begin(), kvec.end(), &(kmer_freq_gt<KmerRec_t>));
  hash_p->from_kmer_vec(kvec, size_ratio);
}
*/



// ---- this hash table is a vec<REC_t> so it more than duplicates the size of the data.

 // KmerMap implemented as a chain hash table

template<class REC_t>
class KmerMap
{
  vec<REC_t>     _hash;
  size_t         _n_rec;
  
  size_t _ih_next(size_t & ih, size_t & inc) const
  { 
    ih += ++inc; // (inc + 5) * (++inc);
    while (ih >= _hash.size()) ih -= _hash.size();
    return ih;
  }

  template<class KMER_t>
  size_t _ih0_from_kmer(const KMER_t & kmer) const 
  { 
    return kmer.hash_64bits() % _hash.size(); 
  }
public:

  typedef REC_t rec_type;

  KmerMap() : _hash(), _n_rec(0) {}

  KmerMap(const vec<REC_t> & kvec, 
          const float ratio = 1.5) : 
    _hash(),
    _n_rec(0)
  { from_kmer_vec(kvec, ratio); }


  // ---- kmer record index if kmer record exists; otherwise, index of empty record

  template<class KMER_t>
  size_t ih_of_kmer(const KMER_t & kmer) const
  {
    size_t ih = _ih0_from_kmer(kmer);
    REC_t rec = _hash[ih];
    size_t inc = 1;
    while (rec.is_valid_kmer() && !kmer.match(rec)) 
      rec = _hash[_ih_next(ih, inc)];
  
    return ih;
  }


  // ---- build the hash from a kmer record vector

  void from_kmer_vec(const vec<REC_t> & kvec, 
                     const float ratio = 1.5,
		     const unsigned verbosity = 1)
  {
    ForceAssertGt(ratio, 1.0);
    _n_rec = kvec.size();

    _hash.clear();
    _hash.resize(ratio * _n_rec, REC_t(0));
    if (verbosity) cout << "nh= " << _hash.size() << endl;

    // ---- set up the hash for each value in _kvec
    for (size_t i = 0; i != _n_rec; i++) {
      const REC_t & rec = kvec[i];
      size_t ih = _ih0_from_kmer(rec);
      size_t inc = 1;
      while (_hash[ih].is_valid_kmer())
        _ih_next(ih, inc);
      _hash[ih] = rec;

      if (verbosity) dots_pct(i, _n_rec);
    }
  }

  void write_binary(const String & fn) const
  {
    BinaryWriter::writeFile(fn.c_str(), _hash);
  }
  
  size_t read_binary(const String & fn) 
  {
    BinaryReader::readFile(fn.c_str(), &_hash);
    return _hash.size();
  }
  
 

  template<class KMER_t>
  REC_t   operator()(const KMER_t & kmer) const { return _hash[ih_of_kmer(kmer)]; }

  template<class KMER_t>
  REC_t & operator()(const KMER_t & kmer)       { return _hash[ih_of_kmer(kmer)]; }

  REC_t   operator[](const size_t i_h) const { return _hash[i_h]; }
  REC_t & operator[](const size_t i_h)       { return _hash[i_h]; }

  size_t  size_hash()                  const { return _hash.size(); }
  size_t  num_recs()                   const { return _n_rec; }
  float   ratio()                      const { return float(size_hash()) / float(num_recs()); }


  
  // ---- outputs the frequencies of number of steps to find a record 
  void report() const
  {
    vec<size_t> freqs(1, 0);
    size_t not_found = 0;
    const size_t n_h = _hash.size();
    for (size_t i_h = 0; i_h < n_h; i_h++) {
      const REC_t & rec0 = _hash[i_h];
      if (rec0) {
        size_t j_h = _ih0_from_kmer(rec0);
        REC_t rec = _hash[j_h];
        size_t inc = 1;
        while (rec && !rec0.match(rec))
          rec = _hash[_ih_next(j_h, inc)];
        
        
        
        if (_hash[j_h]) {
          if (inc + 1 > freqs.size()) freqs.resize(inc + 1, 0);
          freqs[inc]++;
        }
        else {
          not_found++;
        }
      }
      else {
        freqs[0]++;
      }
    }
    for (size_t i = 0; i != freqs.size(); i++) 
      cout << setw(10) << i << " " << setw(10) << freqs[i] << endl;
    cout << setw(10) << not_found << " not found." <<endl;
  }


};

















#endif

