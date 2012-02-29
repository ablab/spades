///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
  Kmer Parcels

  2009-04 Josh Burton (initial developer, algorithm design)
  2009-07 Filipe Ribeiro (paralelization, optimization)
  2010-03 Filipe Ribeiro (objectification, parcelization in memory)
  
  Description:  

    A tool for handling HUGE numbers of kmers.  Enables fast, low-memory,
    parallelizable processing of all of the kmers in a BaseVecVec, and
    catalogues each kmer with a list of its appearances.  
 
  Usage:

      // ---- build on disk
      KmerParcelsDiskStore parcels_disk(K, parcel_head);  
      KmerParcelsBuilder   parcels_disk_builder(bases, parcels_disk, n_threads);
      
      parcels_disk_builder.Build(n_parcels); 
      
      

      // ---- build in memory
      KmerParcelsMemStore parcels_mem(K);                // in memory
      KmerParcelsBuilder  parcels_mem_builder(bases, parcels_mem, n_threads);
      
      parcels_mem_builder.Build(n_parcels);



      // ---- use parcels (either in memory or on disk)
      n_parcels = parcels_store.GetNumParcels();

      for (size_t parcel_ID = 0; parcel_ID != n_parcels; parcel_ID++) {
        KmerParcelReader reader(parcels_store, parcel_ID);

        while (reader.GetNextKmerBatch()) {
          const KmerBatch & batch = reader.CurrentKmerBatch();
          BaseVec kmer = batch.GetKmerBaseVec();
          size_t kmer_freq = batch.GetKmerFreq();
          vec<KmerLoc> kmer_locs = batch;  // KmerBatch is a public vec<KmerLoc>
          ...
        }
      }


  Terminology:

    BaseVec: a read

    BaseVecVec: a vector of reads

    KmerLoc: the position of a kmer in the set of reads: (read_ID, pos)
             where if pos >= 0 the kmer is in the FW direction, 
             and if pos < 0 the kmer is in the RC direction at the base position ~pos.

    KmerBatch: a canonicalized kmer + vec<KmerLoc>

    KmerRecord<K>: a kmer + its KmerLoc
                   note: templatization is use for performance reasons. 
                         sorting a vector of templatized kmers is 2 to 3 times faster 
                         than sorting a vector of non-templatized BaseVec kmers.

    Parcel: a set of KmerBatch written to disk or stored in memory as a vec<KmerBatch>



   Algorithm:

     1. Number of parcels: n_parcels

        Either provided by user or computed based on memory usage restrictions. 
         
        If all the parcels do not fit within the memory limit, then the algorithm is 
        run over several passes, each dealing with a parcel subset.

     2. Kmers:

        Divide the input BaseVecVec into kmers and assign to each kmer a 
        parcel_ID in the range [0..n_parcels).

        Keep only the kmers with parcel_IDs within the parcel subset of the current pass.
        Store the kmers and their locations as a vec< KmerRecord<K> > for each parcel.     

     3. Sort:

        Sort the vec< KmerRecord<K> > over the kmers to bring identical kmers together.
        Each block of identical kmers and their respective locations constitutes a KmerBatch.
             
  
     4. Repeat:
       
        Go back to step 2 for the next pass (i.e., the next parcel subset)

*/



#ifndef KMER_PARCELS_H
#define KMER_PARCELS_H

#include "Basevector.h"            // BaseVec, BaseVecVec
#include "feudal/QualNibbleVec.h"  // QualNibbleVec

#include "Vec.h"              // vec
#include "VecUtilities.h"     // ReverseSortSync
#include "kmers/KmerShape.h"  // DISPATCH_ON_K


#include "system/System.h"    // ForceAssert, etc.

#include "system/Worklist.h"  // worklist
#include "system/LockedData.h"

#include <map>




double NumSigmaBinomial(size_t n, size_t m, double p = 0.5);

void PrintKmerBaseVec(const BaseVec & bv, 
                      const size_t color = 3, ostream& out = cout);

void PrintKmerBaseVec(const BaseVec & bv, const QualNibbleVec & qv, 
                      const size_t color, ostream& out);


#include "kmers/kmer_parcels/KmerParcelsClasses.h" 
// NaifTimer
// NaifBuffer 
// KmerLoc
// KmerRecord<K> 
// KmerBatch



#include "kmers/kmer_parcels/KmerParcelsStatistics.h"
// MapOfCounters
// KmerFrequencyCounters
// ComparativeKmerFreqeucnyCounters
// KmerFrequencyStatistics



#include "kmers/kmer_parcels/KmerParcelAccessor.h" 
// KmerParcel Accessor
//            ReaderImpl
//            WriterImpl
//            Files
//            DiskReaderImpl
//            DiskWriterImpl
//            MemReaderImpl
//            MemWriterImpl








// ---------------------------------------------
//  KmerValidator
// ---------------------------------------------
// example:
// {
//   Validator min10(10);
//   
//   if (min10(x)) {  // same as if (x < 10)
//     ...
//   }
// }
//
// - useful as a parameter to a function
// - expand upon need to include other criteria, e.g. max
//
class Validator
{
private:
  bool   _all;   // all are valid 
  bool   _none;  // all are invalid
  size_t _min;   // only values greater or equal than _min are valid

public: 
  Validator(const size_t & min = 0) 
    : _all(false), _none(false), _min(min)
  {}

  void SetMinimum(const size_t min) 
  { _min = min; _all = _none = false; }
  
  
  void SetNone(const bool none = true)
  {
    _none = none;
    if (_all && _none) _all = false;
  }

  void SetAll(const bool all = true)
  {
    _all = all;
    if (_all && _none) _none = false;
  }

  bool operator()(const size_t val) 
  {
    if (_all)  return true;
    if (_none) return false;
    return (val >= _min);
  }
};












// ---------------------------------------------
// KmerParcels    Store
// ---------------------------------------------
//
//  This is a semi-abstract base class that:
//
//   - is agnostic regarding how the parcels are stored.
//   - provides the interface between a user and the kmer parcels 
//
class KmerParcelsStore
{
private:
  size_t    _K;
  Validator _kmer_validator;


public:
  // ---- constructors
  //KmerParcelsStore();
  KmerParcelsStore(const size_t K) : _K(K), _kmer_validator() {}

  // ---- destructors
  virtual ~KmerParcelsStore() {}


  // ---- pure virtual
  virtual size_t GetNumParcels() const = 0;
  virtual void   SetNumParcels(const size_t n) = 0;

  virtual size_t GetNumReads() const = 0;
  virtual void   SetNumReads(const size_t n) = 0;

  virtual KmerParcelWriterImpl * GetKmerParcelWriter(const size_t parcel_ID) = 0;
  virtual KmerParcelReaderImpl * GetKmerParcelReader(const size_t parcel_ID) const = 0;
  
  
  virtual void Open() = 0;  // to signal that we're about to create parcels
  virtual void Close() = 0; // to signal that we're done creating parcels

  // ---- virtual
  virtual void Log(String const & mess) {}
  virtual String GetDirectoryName() const { return String(""); }


  // ---- misc
  inline size_t GetK() const { return _K; }

  vec<size_t> GetParcelsIDsSizeSorted() const;
  size_t      GetTotalNumKmerBatches() const;
  
  Validator & KmerValidator() { return _kmer_validator; }
  bool KmerValidator(const size_t kf) { return _kmer_validator(kf); }

};  // KmerParcels    Store













// ---------------------------------------------
// KmerParcels    Disk    Store
// ---------------------------------------------
//
//  Derived class that manages parcels stored on disk
//
class KmerParcelsDiskStore : public KmerParcelsStore
{
 private:
  String   _dn;
  String   _n_reads_fn;
  String   _n_parcels_fn;
  String   _log_fn;
  size_t   _n_reads;
  size_t   _n_parcels;
  bool     _keep_on_disk;
  ofstream _log;
  

public:
  // ---- constructors
  KmerParcelsDiskStore(const size_t K, 
                       const String & head_fastb)
    : KmerParcelsStore(K)
  {
    _dn = head_fastb + "." + ToString(GetK()) + "merParcels"; 
    _n_parcels_fn = _dn + "/num_parcels";
    _n_reads_fn = _dn + "/num_reads";
    _log_fn = _dn + "/log";
    _n_reads = 0;
    _n_parcels = 0;
    _keep_on_disk = false;

    if (IsRegularFile(_n_parcels_fn)) {    
      _keep_on_disk = true;

      ifstream ifs(_n_parcels_fn.c_str());
      ifs >> _n_parcels;
      ifs.close();
      
      if (IsRegularFile(_n_reads_fn)) {    
        ifstream ifs(_n_reads_fn.c_str());
        ifs >> _n_reads;
        ifs.close();
      }
    }
  }

  // ---- no copy constructor nor operator=
  //      if one of the copies got destructed and the files were erased 
  //      the copy would be meaningless
  KmerParcelsDiskStore(const KmerParcelsDiskStore & p);
  KmerParcelsDiskStore & operator=(const KmerParcelsDiskStore & p);



  // ---- destructor
  ~KmerParcelsDiskStore()
  {
    if (!_keep_on_disk) {
      /*
      size_t K = GetK();
      for (size_t parcel_ID = 0; parcel_ID != _n_parcels; parcel_ID++) {
        KmerParcelFiles files(_dn, parcel_ID);
        files.Remove();
      }
      */
      System("rm -fr " + _dn);
    }
  }

  // ---- pure virtual
  size_t GetNumParcels() const { return _n_parcels; }

  void SetNumParcels(const size_t n) { _n_parcels = n; }

  size_t GetNumReads() const { return _n_reads; }
  void SetNumReads(const size_t n) { _n_reads = n; }

  

  KmerParcelWriterImpl * GetKmerParcelWriter(const size_t parcel_ID)
  {
    ForceAssertLt(parcel_ID, _n_parcels);
    return new KmerParcelDiskWriterImpl(GetK(), parcel_ID, _dn, _log);
  }    

  KmerParcelReaderImpl * GetKmerParcelReader(const size_t parcel_ID) const
  {
    ForceAssertLt(parcel_ID, _n_parcels);
    //cout << "KmerParcelDiskStore.GetKmerParcelReader(): parcel_ID = " << parcel_ID << endl;
    return new KmerParcelDiskReaderImpl(GetK(), parcel_ID, _dn);
  }    

  
  void Open()  // writing: creates directory and opens log
  {
    if (!IsDirectory(_dn)) Mkpath(_dn);
    _log.open((_log_fn).c_str());
  }

  void Close() // writing: saves n_parcels and closes log
  {
    if (_keep_on_disk) {
      ofstream ofs(_n_parcels_fn.c_str());
      ofs << _n_parcels << endl;
      ofs.close();

      if (_n_reads != 0) {
        ofstream ofs(_n_reads_fn.c_str());
        ofs << _n_reads << endl;
        ofs.close();
      }
    }
    _log.close();
  }


  void Log(String const & mess) { _log << mess << endl; }

  // ---- virtual
  inline String GetDirectoryName() const { return _dn; }


  // ---- Other
  inline void SetKeepOnDisk(const bool keep) { _keep_on_disk = keep; }

  inline bool ExistOnDisk()     // reading: checks if parcels exist
  { return (IsRegularFile(_n_parcels_fn)); }

}; // KmerParcels    Disk    Store




// ---------------------------------------------
// KmerParcels    Mem    Store
// ---------------------------------------------
//
//  Derived class that manages parcels stored in memory
//
class KmerParcelsMemStore : public KmerParcelsStore
{
 private:
  vec< vec<KmerBatch> > _parcels;
  size_t _n_reads;


public:
  // ---- constructors
  KmerParcelsMemStore(const size_t K) : KmerParcelsStore(K), _n_reads(0) {}
  // ---- destructor 
  ~KmerParcelsMemStore() {}
  
  // ---- pure virtual
  size_t GetNumParcels() const { return _parcels.size(); }

  void SetNumParcels(const size_t n) 
  {
    ForceAssertGt(n, _parcels.size());
    _parcels.resize(n);
  }

  size_t GetNumReads() const { return _n_reads; }
  void   SetNumReads(const size_t n) { _n_reads = n; }




  KmerParcelWriterImpl * GetKmerParcelWriter(const size_t parcel_ID)
  {
    size_t n_parcels = GetNumParcels();
    ForceAssertLt(parcel_ID, n_parcels);
    return new KmerParcelMemWriterImpl(GetK(), parcel_ID, _parcels[parcel_ID]);

  }    




  KmerParcelReaderImpl * GetKmerParcelReader(const size_t parcel_ID) const
  {
    size_t n_parcels = GetNumParcels();
    ForceAssertLt(parcel_ID, n_parcels);
    return new KmerParcelMemReaderImpl(GetK(), parcel_ID, _parcels[parcel_ID]);
  }    


  void Open() {}
  void Close() {}

};  // KmerParcels    Mem    Store 















/*
  The two classes below, KmerParcelReader and KmerParcelWriter are a bit tricky to
  understand.  I (Filipe Ribeiro, 4/2010) am not entirely sure this is the best solution.
  
  A KmerParcelReader/Writer must provide access to one parcel alone, without knowing 
  where that parcel is stored (disk or memory).
  
  It is obtained from a KmerParcelsStore like this:
  
  { 
    KmerParcelsStoreDisk parcels(K, "frag_reads_orig");
    size_t n_parcels = parcels.GetNumParcels();
    for (size_t parcel_ID = 0; parcel_ID != n_parcels; parcel_ID++) {
      KmerParcelReader parcel_reader(parcels, parcel_ID);
      ...
    } 
  }

  So the concrete object 'parcels' needs to provide an abstract KmerParcelReader for
  each of its parcels. 

  However, it can only provide either a derived KmerParcelDiskReader/MemReader or 
  a pointer to an abstract KmerParcelReader.  You can't return a reference to a local variable.

  Returning a concrete KmerParcelDiskReader is no good because we want to be 
  agnostic with respect to where the parcels are stored.

  Returning a pointer to a KmerParcelReader means that there is a call to 'new' by 
  the KmerParcelsDiskStore to create a new KmerParcelDiskReader, which in turn means 
  that we must remember to call 'delete' on that pointer to KmerParcelReader.

  So, the following two classes, KmerParcelReader and KmerParcelWriter, are wrappers 
  whose destructors make sure 'delete' is called on the pointers returned by KmerParcelsStore.

  Uff.
  
  See 'KmerParcelAccessor.h' for details on the following classes:
  
    class KmerParcelAccessor

    class KmerParcelReaderImpl : public KmerParcelAccessor
    class KmerParcelWriterImpl : public KmerParcelAccessor
  
    class KmerParcelDiskReaderImpl : public KmerParcelReaderImpl
    class KmerParcelDiskWriterImpl : public KmerParcelWriterImpl

    class KmerParcelMemReaderImpl : public KmerParcelReaderImpl
    class KmerParcelMemWriterImpl : public KmerParcelWriterImpl

*/


// ---------------------------------------------
// KmerParcel    Reader
// ---------------------------------------------
class KmerParcelReader
{
private: 
  const KmerParcelReaderImpl * _reader;

public:
  KmerParcelReader() : _reader(0) {}


  KmerParcelReader(const KmerParcelsStore & parcels, 
                   const size_t parcel_ID) 
    : _reader(parcels.GetKmerParcelReader(parcel_ID)) 
  {
    //cout << "KmerParcelReader constructor: parcelID = " << parcel_ID << endl; 
    //cout << "KmerParcelReader constructor: reader->parcelID = " << GetParcelID()  << endl; 
  }

  // the 'delete' matches the 'new' in parcels.GetKmerParcelReader
  ~KmerParcelReader() { 
    //cout << "KmerParcelReader destructor" << endl;
    delete _reader; } 


  void Init(const KmerParcelsStore & parcels, 
            const size_t parcel_ID) 
  {
    _reader = parcels.GetKmerParcelReader(parcel_ID);
  }

  inline size_t GetK() const { return _reader->GetK(); }
  inline size_t GetParcelID() const { return _reader->GetParcelID(); }
  inline double GetTimer() const { return _reader->GetTimer(); }

  inline size_t GetNumKmerBatches() const 
  { return _reader->GetNumKmerBatches(); }

  inline size_t GetNumKmerBatchesRead() const 
  { return _reader->GetNumKmerBatchesRead(); }
  
  inline size_t GetKmerBatchID() const
  { return _reader->GetKmerBatchID(); }

  inline KmerBatch const & CurrentKmerBatch() const 
  { return _reader->CurrentKmerBatch(); }

  inline vec<KmerLoc> const & CurrentKmerLocs() const 
  { return _reader->CurrentKmerBatch(); }

  inline bool GetNextKmerBatch() const 
  { return _reader->GetNextKmerBatch(); }

  inline bool GetNextKmer(NaifBuffer * kmer, size_t * freq) const
  { return _reader->GetNextKmer(kmer, freq); }

};













// ---------------------------------------------
// KmerParcel    Writer
// ---------------------------------------------

class KmerParcelWriter
{
private: 
  KmerParcelWriterImpl * _writer;

public:
  KmerParcelWriter(KmerParcelsStore & parcels, 
                   const size_t parcel_ID) 
    : _writer(parcels.GetKmerParcelWriter(parcel_ID)) 
  {}

  // the 'delete' matches a 'new' in parcels.GetKmerParcelWriter 
  ~KmerParcelWriter() { delete _writer; } 

  inline size_t GetNumKmerBatches() const 
  { return _writer->GetNumKmerBatches(); }


  inline void PutKmerKmerLocs(const NaifBuffer & kmer, 
                              const vec<const KmerLoc*> & kmer_locs)
  { _writer->PutKmerKmerLocs(kmer, kmer_locs); }

  inline void PutKmerKmerLocs(const NaifBuffer & kmer, 
                              const vec<KmerLoc> & kmer_locs)
  { _writer->PutKmerKmerLocs(kmer, kmer_locs); }
  
  inline void PutKmerBatch(const KmerBatch & batch)
  { _writer->PutKmerBatch(batch); }

};






#include "kmers/kmer_parcels/KmerParcelsBuilder.h" 




















#endif






