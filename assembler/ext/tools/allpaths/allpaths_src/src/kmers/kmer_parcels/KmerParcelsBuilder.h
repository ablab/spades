///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ParcelsBuilder
//  
//    class KmerParcelLock<K> 
//       : public vec< KmerRecord<K> >, private LockedData
//
//          void PushKmerRecords(buf);
//          void SortKmerParcel();
//          * -- Writing parcels is the main time bottle neck -- *
//          void WriteSortedKmerParcel(thread_ID, parcel_ID, parcels_store, stats); 
//
//
//    class KmerParcelVec<K>(parcel0_ID, parcel1_ID, n_parcels, n_blocks) 
//       : public vec<KmerParcelLock<K> >, private LockedData
//
//          void ParseReadKmersForParcelIDs(bases, read1_IDs, block_ID);
//          bool WaitingForBlocks();
//          bool RunNextTask(thread_ID, bases, read1_IDs, parcels_store, stats);
//
//
//    class KmerParcelVecVec<K>(n_parcels, n_threads) 
//       : public vec< KmerParcelVec<K> >
//
//          void RunTasks(thread_ID, bases, read1_IDs, parcels_store, stats);
//
//
//    class ParcelProcessor<K>(bases, read1_IDs, parcel_vec_vec, parcels_store, stats);
//
//          void operator()(thread_ID)
//
//
//    class KmerParcelsBuilder(bases, parcels_store, n_threads);
//
//          void   ComputeNumParcels();
//          size_t GetNumParcels() const { return _n_parcels; }
//          void   SetVerbose(verbose);
//          void   SetCanonicalize(canonicalize);
//          void   Build(n_parcels);
//
//
//    class KmerParcelsCoupledBuilder(bases1, bases2, parcels1_store, parcels2_store)
//
//          void   SetVerbose(verbose);
//          void   Build(n_parcels);
//
//

#ifndef KMER_PARCELS_BUILDER_H
#define KMER_PARCELS_BUILDER_H

#include "kmers/KmerParcels.h" // KmerParcelsStore

// ----------------------------------------
//  KmerParcel<K>
// ----------------------------------------

// This class is used by the parcel builder to add kmer records to the parcel
// in a thread safe way.
// 
// When done, there are methods to sort and write to a KmerParcelStore.


template<size_t K>
class KmerParcel : public vec< KmerRecord<K> >, private LockedData
{
private:
  NaifTimer _timer_sort;
  NaifTimer _timer_write;

public:
  KmerParcel() : _timer_sort(), _timer_write() {} 
  ~KmerParcel() {}

  // ---- copy constructor
  //      specified so that vec<> can be copied without copying LockedData
  KmerParcel(KmerParcel<K> const & that) 
    : vec<KmerRecord<K> >(that), 
      LockedData()
  { }

  // ---- = operator
  //      specified so that vec<> can be copied without copying LockedData
  KmerParcel & operator=(KmerParcel<K> const & that) 
  {
    vec< KmerRecord<K> >::operator=(that);
    return *this;
  }


  void PushKmerRecords(const vec< KmerRecord<K> > & buf);
  void SortKmerParcel();

  void WriteSortedKmerParcel(const size_t thread_ID, 
                             const size_t parcel_ID, 
                             KmerParcelsStore & parcels_store,
                             KmerFrequencyStatistics & stats);

};  //  KmerParcel<K>





// ----------------------------------------
//  KmerParcelVec<K>
// ----------------------------------------


// This class stores the set of KmerParcel's that are processed together.  


template<size_t K>
class KmerParcelVec : public vec<KmerParcel<K> >, private LockedData
{
private:

  size_t _parcel0_ID;  // id of the first parcel on the set
  size_t _parcel1_ID;  // id of the parcel after the last one
  size_t _n_parcels; // total number of parcels 

  // the read data is divided into blocks when
  // searching for the subset of parcel IDs
  // each read block is processed by its own thread 
  size_t _n_blocks;  
  size_t _n_blocks_started; 
  size_t _n_blocks_done;

  // after the parcel IDs are computed, then each 
  // thread sorts and writes each parcel.
  size_t _n_sub_parcels; // parcel_ID1 - parcel_ID0
  size_t _n_sub_parcels_started;
  size_t _n_sub_parcels_done;

  //LockedData _state_mutex;

public:
  KmerParcelVec() {}
  ~KmerParcelVec() {}

  // the copy constructor and operator= are critical!  
  // LockedData can't be copied but vec<> can.
  KmerParcelVec(KmerParcelVec<K> const & that) 
    : vec< KmerParcel<K> >(that), LockedData()
  {}
    
  KmerParcelVec & operator=(KmerParcelVec<K> const & that)
  {
    vec< KmerParcel<K> >::operator=(that);
    return *this;
  }

  KmerParcelVec(const size_t parcel0_ID,
                const size_t parcel1_ID,
                const size_t n_parcels,
                const size_t n_blocks) 
  { Init(parcel0_ID, parcel1_ID, n_parcels, n_blocks); }


  void Init(const size_t parcel0_ID,
            const size_t parcel1_ID,
            const size_t n_parcels,
            const size_t n_blocks);

private:
  void ParseReadKmersForParcelIDs(const BaseVecVec & bases,
                                  const vec<size_t> & read1_IDs,
                                  const size_t block_ID);
public:
  bool WaitingForBlocks();

  bool RunNextTask(const size_t thread_ID,
                   const BaseVecVec & bases, 
                   const vec<size_t> & read1_IDs,
                   KmerParcelsStore & parcels_store,
                   KmerFrequencyStatistics & stats);

};   //  KmerParcelVec<K>







// ----------------------------------------
//  KmerParcelVecVec<K>
// ----------------------------------------


// This class stores all the sets of KmerParcel.

template<size_t K>
class KmerParcelVecVec : public vec< KmerParcelVec<K> >
{
public:

  KmerParcelVecVec(const size_t n_parcels,
                   const size_t n_threads);

  void RunTasks(const size_t thread_ID,
                const BaseVecVec & bases, 
                const vec<size_t> & read1_IDs,
                KmerParcelsStore & parcels_store,
                KmerFrequencyStatistics & stats);

};







// ----------------------------------------
//  ParcelProcessor<K>
// ----------------------------------------

template<size_t K>
class ParcelProcessor 
{
private:
  // Inputs
  const BaseVecVec & _bases;
  const vec<size_t> & _read1_IDs;

  // Outputs
  KmerParcelVecVec<K>     * _parcel_vec_vec;
  KmerParcelsStore        & _parcels_store;
  KmerFrequencyStatistics & _stats;

public:
  ParcelProcessor(const BaseVecVec & bases,
                  const vec<size_t> & read1_IDs,
                  KmerParcelVecVec<K> * parcel_vec_vec,
                  KmerParcelsStore & parcels_store,
                  KmerFrequencyStatistics & stats)
    : _bases(bases),
      _read1_IDs(read1_IDs),
      _parcel_vec_vec(parcel_vec_vec),
      _parcels_store(parcels_store),
      _stats(stats)
  {
    //cout << "sizeof _parcel_vec_vec = " << (*_parcel_vec_vec).size() << endl;
  }

  void operator() (const size_t thread_ID) 
  {
    //cout << "thread_ID = " << thread_ID << endl;
    //cout << "sizeof _parcel_vec_vec = " << (*_parcel_vec_vec).size() << endl;

    _parcel_vec_vec->RunTasks(thread_ID,
                              _bases, 
                              _read1_IDs,
                              _parcels_store,
                              _stats);
  }
};













// ----------------------------------------
// KmerParcels    Builder
// ----------------------------------------

class KmerParcelsBuilder
{
private:
  const BaseVecVec &      _bases;
  KmerParcelsStore &      _parcels_store;
  const size_t            _n_threads;
  KmerFrequencyStatistics _stats;

  size_t _K;
  size_t _n_parcels;
  size_t _n_kmers_total;

  vec<size_t> _read1_IDs;

  size_t _n_bytes_max; 
  bool _verbose;
  bool _canonicalize;

  bool _write_gc_stats;

  ofstream _log;


public:
  // ---- constructor ----

  KmerParcelsBuilder(const BaseVecVec & bases, 
                     KmerParcelsStore & parcels_store, 
                     const size_t n_threads) 
    : _bases(bases),
      _parcels_store(parcels_store),
      _n_threads(n_threads),
      _stats(bases.size(), n_threads, parcels_store.GetK())
  { 
    _K = parcels_store.GetK();
    _n_parcels = 0;
    
    _parcels_store.SetNumReads(bases.size());

    // ---- compute total number of kmers in 'bases'
    _n_kmers_total = 0;
    for (size_t read_ID = 0; read_ID != _bases.size(); read_ID++) {
      const size_t n_bases = _bases[read_ID].size();
      if (n_bases >= _K) 
        _n_kmers_total += n_bases - _K + 1;  // add the number of kmers on read
    }
    
    _n_bytes_max = 64ul << 30; // 64 GB
    _verbose = false;
    _write_gc_stats = false;
    _canonicalize = true;
  }

  // ---- destructor ----
  ~KmerParcelsBuilder() {}


private:
  size_t EstimateTotalMemory();
  void ComputeRead1IDs();

  template<size_t K>
  void BuildTemplate();
  
public:
  void ComputeNumParcels();
  void SetVerbose(const bool verbose = true) { _verbose = verbose; }
  void SetWriteGCStats(const bool write_gc_stats = true) { _write_gc_stats = write_gc_stats; }
  void SetCanonicalize(const bool canonicalize = true) { _canonicalize = canonicalize; }

  size_t GetNumParcels() const { return _n_parcels; }

  void Build(const size_t n_parcels = 0);

  const KmerFrequencyCounters & GetKmerFrequencyCounters()
  { return _stats.AllKmerFrequencyCounters(); }

};





// ----------------------------------------
// KmerParcels    Coupled    Builder
// ----------------------------------------

class KmerParcelsCoupledBuilder
{
private:
  KmerParcelsBuilder _builder1;
  KmerParcelsBuilder _builder2;

public:
  // ---- constructor ----

  KmerParcelsCoupledBuilder(const BaseVecVec & bases1, 
                            const BaseVecVec & bases2, 
                            KmerParcelsStore & parcels1, 
                            KmerParcelsStore & parcels2, 
                            const size_t n_threads) 
    : _builder1(bases1, parcels1, n_threads),
      _builder2(bases2, parcels2, n_threads)

  { 
    ForceAssertEq(parcels1.GetK(), parcels2.GetK());
  }

public:
  void SetVerbose(const bool verbose = true) 
  { 
    _builder1.SetVerbose(verbose);
    _builder2.SetVerbose(verbose);
  }

  void Build(size_t n_parcels = 0) 
  {
    if (n_parcels == 0) {

      _builder1.ComputeNumParcels();
      _builder2.ComputeNumParcels();
      
      n_parcels  = Max(_builder1.GetNumParcels(), 
                       _builder2.GetNumParcels());
    }

    _builder1.Build(n_parcels);
    _builder2.Build(n_parcels);
  }
};


#endif
