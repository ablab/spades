///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef KMERS__KMER_PARCELS__KMER_PARCELS_STATISTICS_H
#define KMERS__KMER_PARCELS__KMER_PARCELS_STATISTICS_H



// --------------------------------------
// MapOfCounters
// --------------------------------------
class MapOfCounters : public map<size_t, size_t>
{
public:
  MapOfCounters() : map<size_t, size_t>() {}

  template <class T>
  MapOfCounters(const vec<T> & v) 
    : map<size_t, size_t>()
  {
    for (typename vec<T>::const_iterator it = v.begin(); it != v.end(); it++)
      (*this)[*it]++;
  }


  size_t Total() const 
  {
    size_t total = 0;
    for (MapOfCounters::const_iterator it = (*this).begin(); 
         it != (*this).end(); it++)
      total += it->second;
    return total;
  }

  void Write(const String & filename) const
  {
    double total = static_cast<double>(Total());

    ofstream ofs;
    ofs.open(filename.c_str());

    ofs << "# 1:x 2:y 3:y/total 4:cum(y) 5:cum(y)/total" << endl;

    size_t cum = 0;
    for (MapOfCounters::const_iterator it = (*this).begin(); 
       it != (*this).end(); it++) {
      cum += it->second;
      ofs << " " << it->first
          << " " << it->second
          << " " << setprecision(10) << static_cast<double>(it->second) / total
          << " " << cum
          << " " << setprecision(10) << static_cast<double>(cum) / total
          << endl;
    }
    ofs.close();
  }


  MapOfCounters & operator += (const MapOfCounters & m) 
  {
    for (MapOfCounters::const_iterator it = m.begin(); it != m.end(); it++)
      (*this)[(*it).first] += (*it).second;
    return *this;
  } 

};




// --------------------------------------
// KmerFrequencyCounters
// --------------------------------------
class KmerFrequencyCounters : public MapOfCounters 
{
public:
  void Write(const String & filename) const
  {
    // ---- counting kmers and distinct kmers
    const size_t total_num_distinct_kmers = (*this).Total();
  
    size_t total_num_kmers = 0;
    for (KmerFrequencyCounters::const_iterator it = (*this).begin(); 
         it != (*this).end(); it++) {
      
      const size_t & kmer_frequency = it->first;
      const size_t & num_distinct_kmers = it->second;
      
      total_num_kmers += kmer_frequency * num_distinct_kmers;
    }
    
    // ---- writing kmer frequency stats to disk
    ofstream outfs;
    outfs.open(filename.c_str());
    outfs << ("#"
              " 1:kmer_frequency"
              " 2:num_distinct_kmers"
              " 3:2/total_num_distinct_kmers"
              " 4:cummulative(2)"
              " 5:4/total_num_distinct_kmers"
              " 6:1*2"
              " 7:6/total_num_kmers"
              " 8:cummulative(6)"
              " 9:8/total_num_kmers") 
          << endl;


    size_t cum_num_distinct_kmers = 0;
    size_t cum_num_kmers = 0;
    for (KmerFrequencyCounters::const_iterator it = (*this).begin(); 
         it != (*this).end(); it++) {
    
      const size_t & kmer_frequency = it->first;
      const size_t & num_distinct_kmers = it->second;

      cum_num_distinct_kmers += num_distinct_kmers;
      const size_t num_kmers = kmer_frequency * num_distinct_kmers;
      cum_num_kmers += num_kmers;
    
      outfs << " " << kmer_frequency 
            << " " << num_distinct_kmers
            << " " << setprecision(10) << static_cast<double>(num_distinct_kmers) / static_cast<double>(total_num_distinct_kmers) 
            << " " << cum_num_distinct_kmers 
            << " " << setprecision(10) << static_cast<double>(cum_num_distinct_kmers) / static_cast<double>(total_num_distinct_kmers) 
            << " " << num_kmers
            << " " << setprecision(10) << static_cast<double>(num_kmers) / static_cast<double>(total_num_kmers) 
            << " " << cum_num_kmers 
            << " " << setprecision(10) << static_cast<double>(cum_num_kmers) / static_cast<double>(total_num_kmers) 
            << endl;
    }
  
    outfs.close();
  }
};



class GenomeSizeEstimator
{
  size_t _kf_mode;

  size_t _nk_total;
  size_t _nk_bad_low_kf;
  size_t _nk_good_uniq;
  size_t _nk_good_rep;
  size_t _nk_bad_high_kf;

  size_t _ndk_total;
  size_t _ndk_bad_low_kf;
  size_t _ndk_good_uniq;
  size_t _ndk_good_rep;
  size_t _ndk_bad_high_kf;
  
public:
  GenomeSizeEstimator(const KmerFrequencyCounters & kfc)
  {
    // ---- First convert map to a set of convenient vectors.

    const size_t nkf = kfc.rbegin()->first + 2;
    cout << "nkf = " << nkf << endl;

    vec<size_t> nk(nkf, 0);     // number of kmers at freq kf
    vec<size_t> ndk(nkf, 0);    // number of distinct kmers at freq kf
    vec<size_t> cnk(nkf, 0);    // cumulative number of kmers at freq kf
    vec<size_t> cndk(nkf, 0);   // cumulative number of distinct kmers at freq kf

    typedef KmerFrequencyCounters::const_iterator FreqIter;
    for (FreqIter it = kfc.begin(); it != kfc.end(); it++) {
      const size_t kf = it->first;
      ndk[kf] = it->second;
      nk [kf] = kf * ndk[kf];
    }
    for (size_t kf = 1; kf != nkf; kf++) {
      cndk[kf] = cndk[kf - 1] + 0.5 * (ndk[kf - 1] + ndk[kf]);
      cnk [kf] = cnk [kf - 1] + 0.5 * (nk [kf - 1] + nk [kf]);
    }


    // ---- Separate kmer spectrum in 4 regions based on the kf:
    //    1      ... kf_lo    : bad kmers with low frequency
    //    kf_lo  ... kf_hi    : good kmers CN = 1
    //    kf_hi  ... kf_max   : godd kmers CN > 1 (repetitive
    //    kf_max ... inf      : bad kmers with high frequency
 

    // ---- Find first guess for low kmer frequency
  
    size_t kf_lo = 2;
    while (kf_lo + 1 != nkf && nk[kf_lo + 1] < nk[kf_lo]) {
      kf_lo++;
    }

    //cout << "kf_lo = " << kf_lo << endl;

    // ---- Find the kmer frequency mode above low kf 
    
    _kf_mode = kf_lo;
    for (size_t kf = kf_lo + 1; kf != nkf; kf++)
      if (nk[kf] > nk[_kf_mode])
        _kf_mode = kf;

    //cout << "kf_mode = " << _kf_mode << endl;

    // ---- Refine the low kf 
    for (size_t kf = _kf_mode - 1; kf != 0; kf--)
      if (nk[kf] < nk[kf_lo])
        kf_lo = kf;

    //cout << "kf_lo = " << kf_lo << endl;


    // ---- Define the high kf above the mode 
    size_t kf_hi = _kf_mode + (_kf_mode - kf_lo);

    //cout << "kf_hi = " << kf_hi << endl;
  
  
    // ---- Define maximum kf above which we neglect data.
    //      At times there are a lot of very high frequency kmers that
    //      interfere with the genome size estimate. Those kmers 
    //      need to be thrown away.
    //
    //      Just a rough estimate for kf_max is good enough.
    //      Assume
    //  
    //          ndk[kf >= 2]  ~=  4 * ndk[2] / kf^2
    //         
    //      then
    // 
    //          ndk[kf_max] = 1   =>  kf_max = sqrt(4 * ndk[2])
    //
    //      Because  'kf' in the data is equivalent to  
    //      'kf * coverage' in the genome, 
    //      we need to multiply kf_max by the coverage freq. 
    //

    size_t kf_max = _kf_mode * sqrt(4 * ndk[2 * _kf_mode] * _kf_mode);
    if (kf_max >= nkf) 
      kf_max = nkf - 1;
    
    //cout << "kf_max = " << kf_max << endl;

    _nk_total       = cnk[nkf - 1];
    _nk_bad_low_kf  = cnk[kf_lo];
    _nk_good_uniq   = cnk[kf_hi]  - cnk[kf_lo];
    _nk_good_rep    = cnk[kf_max] - cnk[kf_hi];
    _nk_bad_high_kf = _nk_total   - cnk[kf_max];
     
    //cout << "nk_total = " << _nk_total << endl;
    //cout << "nk_bad_low_kf = " << _nk_bad_low_kf << endl;
    //cout << "nk_good_uniq = " << _nk_good_uniq << endl;
    //cout << "nk_good_rep = " << _nk_good_rep << endl;
    //cout << "nk_bad_high_kf = " << _nk_bad_high_kf << endl;



    _ndk_total       = cndk[nkf - 1];
    _ndk_bad_low_kf  = cndk[kf_lo];
    _ndk_good_uniq   = cndk[kf_hi]  - cndk[kf_lo];
    _ndk_good_rep    = cndk[kf_max] - cndk[kf_hi];
    _ndk_bad_high_kf = _ndk_total   - cndk[kf_max];

    //cout << "ndk_total = " << _ndk_total << endl;
    //cout << "ndk_bad_low_kf = " << _ndk_bad_low_kf << endl;
    //cout << "ndk_good_uniq = " << _ndk_good_uniq << endl;
    //cout << "ndk_good_rep = " << _ndk_good_rep << endl;
    //cout << "ndk_bad_high_kf = " << _ndk_bad_high_kf << endl;
  }



  size_t genome_size_unique()
  {
    return _ndk_good_uniq;
  }

  size_t genome_size_repetitive()
  {
    float kf_ave_uniq = float(_nk_good_uniq) / float(_ndk_good_uniq);
    return size_t(0.5 + float(_nk_good_rep) / kf_ave_uniq);
  }

  size_t genome_size_total()
  {
    return genome_size_unique() + genome_size_repetitive();  
  }
 
};



// --------------------------------------
// ComparativeKmerFrequencyCounters
// --------------------------------------
class ComparativeKmerFrequencyCounters
{
public:
  KmerFrequencyCounters both1, both2, only1, only2;

  ComparativeKmerFrequencyCounters() {}

  ComparativeKmerFrequencyCounters & operator += (const ComparativeKmerFrequencyCounters & c) 
  {
    this->both1 += c.both1;
    this->both2 += c.both2;
    this->only1 += c.only1;
    this->only2 += c.only2;
    return *this;
  } 
};




// --------------------------------------
// KmerFrequencyStatistics
// --------------------------------------

class KmerFrequencyStatistics
{
private: 
  LockedData mutexes_min[512];
  LockedData mutexes_max[512];

  vec<KmerFrequencyCounters> kmer_freq_counters;
  vec<vec<KmerFrequencyCounters> > kmer_gc_freq_counters;
  
  size_t m_n_gc;

  bool collapsed;

public:
  vec<uint32_t> reads_min_kmer_freq;
  vec<uint32_t> reads_max_kmer_freq;

  KmerFrequencyStatistics(const size_t n_reads, 
                          const size_t n_blocks,
                          const size_t K = 0) 
  {
    reads_min_kmer_freq.resize(n_reads, 1 << 30);
    reads_max_kmer_freq.resize(n_reads, 0);

    kmer_freq_counters.resize(n_blocks);
    kmer_gc_freq_counters.resize(n_blocks);

    m_n_gc = 0;
    if (K != 0) {
      m_n_gc = K + 1;
      for (size_t ib = 0; ib != n_blocks; ib++)
        kmer_gc_freq_counters[ib].resize(m_n_gc);
    }
  }


  void UpdateMinMax(const size_t read_ID, 
                    const size_t kmer_freq)
  {
    if (true) {
      Locker lock(mutexes_max[read_ID & 0x1FF]);
      if (reads_max_kmer_freq[read_ID] < kmer_freq) 
        reads_max_kmer_freq[read_ID] = kmer_freq;
    }
    if (true) {
      Locker lock(mutexes_min[read_ID % 0x1FF]);
      if (reads_min_kmer_freq[read_ID] > kmer_freq)
        reads_min_kmer_freq[read_ID] = kmer_freq;
    }
  }

  void UpdateFreqs(const size_t block_ID, 
                   const size_t kmer_freq,
                   const size_t i_gc = 0)
  {
    kmer_freq_counters[block_ID][kmer_freq]++;
    if (m_n_gc) 
      kmer_gc_freq_counters[block_ID][i_gc][kmer_freq]++;
    collapsed = false;
  }


  void Update(const vec<const KmerLoc*> & kmer_locs_p,
              const size_t block_ID, 
              const size_t i_gc = 0)
  {
    size_t kmer_freq = kmer_locs_p.size();
    
    for (vec<const KmerLoc*>::const_iterator it = kmer_locs_p.begin(); it != kmer_locs_p.end(); ++it) 
      UpdateMinMax((*it)->GetReadID(), kmer_freq);
    UpdateFreqs(block_ID, kmer_freq, i_gc);
  }


  

  KmerFrequencyCounters & AllKmerFrequencyCounters()
  {
    Collapse();
    return kmer_freq_counters[0];
  }

  
  vec<KmerFrequencyCounters> & AllKmerGCFrequencyCounters()
  {
    Collapse();
    return kmer_gc_freq_counters[0];
  }





private:
  void Collapse() 
  {
    if (!collapsed) {
      // accumulate all freqs on kmer_freq_counters[0]
      
      for (size_t ib = 1; ib != kmer_freq_counters.size(); ++ib) {
        kmer_freq_counters[0] += kmer_freq_counters[ib];       
        kmer_freq_counters[ib].clear();
        
        for (size_t i_gc = 0; i_gc != m_n_gc; i_gc++) {
          kmer_gc_freq_counters[0][i_gc] += kmer_gc_freq_counters[ib][i_gc];       
          kmer_gc_freq_counters[ib][i_gc].clear();
        }
      }
      collapsed = true;
    }  
  }

};





#endif
