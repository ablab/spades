///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// -------------------------------------------------------------------------------
// KmerParcel    Accessor
// -------------------------------------------------------------------------------
class KmerParcelAccessor
{
private:
  size_t       _K;
  size_t       _parcel_ID;

  NaifBuffer   _kmer_buf;
  
  mutable NaifTimer _timer;
  
public:
  // ---- constructor
  KmerParcelAccessor(const size_t K, 
                     const size_t parcel_ID) 
    : _K(K), 
      _parcel_ID(parcel_ID),
      _kmer_buf(2*K),
      _timer() 
  {}

  // ---- destructor
  virtual ~KmerParcelAccessor() 
  { 
    //cout << "KPAccessor destructor" << endl; 
  }
  
protected:
  
  inline NaifBuffer & KmerBuf() { return _kmer_buf; }
  inline const NaifBuffer & KmerBuf() const { return _kmer_buf; }

  inline void TimerStart() const { _timer.Start(); } 
  inline void TimerStop() const { _timer.Stop(); } 

public:
  // --- accessors
  inline size_t GetK()           const { return _K; }
  inline size_t GetParcelID()    const { return _parcel_ID; }
  inline double GetTimer()       const { return _timer.GetTimer(); }
  
  virtual size_t GetNumKmerBatches() const = 0;
};






// -------------------------------------------------------------------------------
// KmerParcel    Reader    Impl
// -------------------------------------------------------------------------------
class KmerParcelReaderImpl : public KmerParcelAccessor
{
public:
  KmerParcelReaderImpl(const size_t K, 
                       const size_t parcel_ID) 
    : KmerParcelAccessor(K, parcel_ID) {}

  virtual size_t            GetNumKmerBatchesRead() const = 0;

  virtual size_t            GetKmerBatchID() const = 0;

  virtual KmerBatch const & CurrentKmerBatch() const = 0;

  virtual bool              GetNextKmerBatch() const = 0;

  virtual bool              GetNextKmer(NaifBuffer * kmer, size_t * freq) const = 0;
};






// -------------------------------------------------------------------------------
// KmerParcel    Writer    Impl
// -------------------------------------------------------------------------------
class KmerParcelWriterImpl : public KmerParcelAccessor
{
public:
  KmerParcelWriterImpl(const size_t K, 
                       const size_t parcel_ID) 
    : KmerParcelAccessor(K, parcel_ID) {}

  virtual void PutKmerKmerLocs(const NaifBuffer & kmer, 
                               const vec<const KmerLoc*> & kmer_locs) = 0;

  virtual void PutKmerKmerLocs(const NaifBuffer & kmer, 
                               const vec<KmerLoc> & kmer_locs) = 0;

  virtual void PutKmerBatch(const KmerBatch & batch) = 0;
};






// -------------------------------------------------------------------------------
// KmerParcel    Files  
// -------------------------------------------------------------------------------

class KmerParcelFiles
{
private:
  String _kkf_fn;  // file name for storage of kmer + kmer_freq 
  String _kl_fn;   // file name for storage of vec<kmer_loc>
  String _size_fn; // file name for storage of number of kmer batches  
  
public:
  // Constructor
  KmerParcelFiles(const String & dn, const size_t parcel_ID) 
  {
    const String base = dn + "/" + ToString(parcel_ID);

    _kkf_fn  = base + ".v3.kkf";
    _kl_fn   = base + ".v3.kl";
    _size_fn = base + ".v3.size";
  }
  
  String KKFFileName() const { return _kkf_fn; }
  String KLFileName() const { return _kl_fn; }
  String SizeFileName() const { return _size_fn; }


  void Remove()
  {
    ::Remove(KKFFileName());
    ::Remove(KLFileName());
    ::Remove(SizeFileName());
  }

};








// -------------------------------------------------------------------------------
// KmerParcel    Disk    ReaderImpl
// -------------------------------------------------------------------------------

class KmerParcelDiskReaderImpl : public KmerParcelReaderImpl
{
private:
  const KmerParcelFiles _files;
  mutable KmerBatch     _batch;
  mutable uint32_t      _kmer_freq;
  
  mutable size_t        _n_batches_read; // batches read so far
  mutable size_t        _n_batches;   // total number of KmerBatches

  mutable ifstream      _is_kl;      // input stream kmerloc
  mutable ifstream      _is_kkf;     // input stream kmerkmerfreq

  
public:
  // Constructor
  KmerParcelDiskReaderImpl(const size_t K, 
                           const size_t parcel_ID,
                           const String & dn) 
    : KmerParcelReaderImpl(K, parcel_ID),
      _files(dn, parcel_ID),
      _batch(K)
  {
    ifstream is_size(_files.SizeFileName().c_str());
    is_size >> _n_batches;
    is_size.close();
    _n_batches_read = 0;

    _is_kl.open(_files.KLFileName().c_str()); 
    _is_kkf.open(_files.KKFFileName().c_str()); 
  }
  
  ~KmerParcelDiskReaderImpl()
  {
    //cout << "KPDiskReaderImpl destructor" << endl;  
    _is_kl.close();
    _is_kkf.close();
  }
  

private:
  inline void BinaryReadKmer() const
  { 
    _batch.Kmer().BinaryStreamRead(_is_kkf); 
  }

  inline void BinaryReadKmerFreq() const
  {
    _is_kkf.read(reinterpret_cast<char *>(&_kmer_freq), sizeof(uint32_t));
  }

public:
  // Accessors
  inline size_t GetNumKmerBatches() const 
  { return _n_batches; }

  inline size_t GetNumKmerBatchesRead() const 
  { return _n_batches_read; }

  inline size_t GetKmerBatchID() const 
  { return _n_batches_read - 1; }

  inline KmerBatch const & CurrentKmerBatch() const 
  { return _batch; }

  // return false if no batches left
  bool GetNextKmerBatch() const
  {
    if (_n_batches_read == _n_batches)
      return false; 

    TimerStart();

    BinaryReadKmer();
    BinaryReadKmerFreq();

    _batch.resize(_kmer_freq);
    
    for (size_t i = 0; i != _kmer_freq; i++) 
      _batch[i].BinaryStreamRead(_is_kl);

    _n_batches_read++;
   
    TimerStop();

    return true;
  }


  // return false if no batches left
  bool GetNextKmer(NaifBuffer * kmer, size_t * freq) const
  {
    if (_n_batches_read == _n_batches)
      return false; 

    TimerStart();

    kmer->BinaryStreamRead(_is_kkf);

    BinaryReadKmerFreq();
    *freq = _kmer_freq;

    _batch.resize(0);

    _n_batches_read++;
   
    TimerStop();

    return true;
  }

};








// -------------------------------------------------------------------------------
// KmerParcel    Disk    WriterImpl
// -------------------------------------------------------------------------------
class KmerParcelDiskWriterImpl : public KmerParcelWriterImpl
{
private:
  const KmerParcelFiles _files;
  ofstream   & _log;

  size_t       _n_batches_written;   // total number of KmerBatches written

  ofstream     _os_kl;      // output stream kmerloc
  ofstream     _os_kkf;     // output stream kmerkmerfreq

public:
  // Constructor
  KmerParcelDiskWriterImpl(const size_t K,
                           const size_t parcel_ID,
                           const String & dn,
                           ofstream & log)
    : KmerParcelWriterImpl(K, parcel_ID),
      _files(dn, parcel_ID),
      _log(log),
      _n_batches_written(0)
  {
    _os_kl.open(_files.KLFileName().c_str()); 
    _os_kkf.open(_files.KKFFileName().c_str()); 
  }


  // Destructor
  ~KmerParcelDiskWriterImpl()
  {
    //cout << "KPDiskWriterImpl destructor" << endl;  

    _os_kkf.close();
    _os_kl.close();
      
    //cout << _files.SizeFileName() << endl;

    ofstream os_size(_files.SizeFileName().c_str());
    os_size << _n_batches_written << "\n";
    os_size.close();

    _log << (Date() + " (KPDW " + ToString(GetParcelID()) + 
             "): wrote " + ToString(_n_batches_written) + 
             " kmer batches") << endl;
  }
  

private:
  inline
  void BinaryWriteKmer(const NaifBuffer & kmer_buf)
  {
    kmer_buf.BinaryStreamWrite(_os_kkf);
  }

  inline 
  void BinaryWriteKmerFreq(const size_t kmer_freq)
  {
    uint32_t kf = kmer_freq;
    _os_kkf.write(reinterpret_cast<char *>(&kf), sizeof(uint32_t));
  }



public:

  size_t GetNumKmerBatches() const { return _n_batches_written; }
  

  void PutKmerKmerLocs(const NaifBuffer & kmer, 
                       const vec<const KmerLoc*> & kmer_locs)
  {
    TimerStart();
    size_t kmer_freq = kmer_locs.size();
    
    BinaryWriteKmer(kmer);
    BinaryWriteKmerFreq(kmer_freq);
    
    typedef vec<const KmerLoc*>::const_iterator KmerLocIter;
    for (KmerLocIter it = kmer_locs.begin(); it != kmer_locs.end(); it++)
      (*(*it)).BinaryStreamWrite(_os_kl);
    
    _n_batches_written++;
    TimerStop();
  }

  void PutKmerKmerLocs(const NaifBuffer & kmer, 
                       const vec<KmerLoc> & kmer_locs)
  {
    TimerStart();
    size_t kmer_freq = kmer_locs.size();

    BinaryWriteKmer(kmer);
    BinaryWriteKmerFreq(kmer_freq);
    
    typedef vec<KmerLoc>::const_iterator KmerLocIter;
    for (KmerLocIter it = kmer_locs.begin(); it != kmer_locs.end(); it++)
      (*it).BinaryStreamWrite(_os_kl);
    
    _n_batches_written++;
    TimerStop();
  }
  
  void PutKmerBatch(const KmerBatch & batch)
  {
    PutKmerKmerLocs(batch.Kmer(), batch);
  }
};









// -------------------------------------------------------------------------------
// KmerParcel    Mem    ReaderImpl
// -------------------------------------------------------------------------------
class KmerParcelMemReaderImpl : public KmerParcelReaderImpl
{
private:
  const vec<KmerBatch> & _parcel;
  mutable size_t _n_batches_read;
  
public:
  // Constructor
  KmerParcelMemReaderImpl(const size_t K, 
                          const size_t parcel_ID,
                          const vec<KmerBatch> & parcel) 
    : KmerParcelReaderImpl(K, parcel_ID),
      _parcel(parcel)
  {
    _n_batches_read = 0;
  }

public:
  inline size_t GetNumKmerBatches() const 
  { return _parcel.size(); }

  inline size_t GetNumKmerBatchesRead() const 
  { return _n_batches_read; }

  inline size_t GetKmerBatchID() const 
  { return _n_batches_read - 1; }

  inline KmerBatch const & CurrentKmerBatch() const 
  { return _parcel[GetKmerBatchID()]; }

  bool GetNextKmerBatch() const 
  {
    if (_n_batches_read == GetNumKmerBatches())
      return false; 
      
    _n_batches_read++;
    return true;
  }


  // return false if no batches left
  bool GetNextKmer(NaifBuffer * kmer, size_t * freq) const
  {
    if (_n_batches_read == GetNumKmerBatches())
      return false; 

    const KmerBatch & batch = CurrentKmerBatch();
    *kmer = batch.Kmer();

    *freq = batch.GetKmerFreq();

    _n_batches_read++;
   
    return true;
  }



};





// -------------------------------------------------------------------------------
// KmerParcel   Mem    WriterImpl
// -------------------------------------------------------------------------------
class KmerParcelMemWriterImpl : public KmerParcelWriterImpl
{
private: 
  vec<KmerBatch> & _parcel;

public:
  // Constructor
  KmerParcelMemWriterImpl(const size_t K, 
                          const size_t parcel_ID,
                          vec<KmerBatch> & parcel) 
    : KmerParcelWriterImpl(K, parcel_ID),
      _parcel(parcel)      
  {}
  
  
  inline size_t GetNumKmerBatches() const { return _parcel.size(); }
  
  
  void PutKmerBatch(const KmerBatch & batch)
  {
    TimerStart();
    _parcel.push_back(batch);
    TimerStop();
  }
  

  void PutKmerKmerLocs(const NaifBuffer & kmer, 
                       const vec<const KmerLoc*> & kmer_locs_p)
  {
    TimerStart();

    size_t kmer_freq = kmer_locs_p.size();
    vec<KmerLoc> kmer_locs(kmer_freq);

    for (size_t i = 0; i != kmer_freq; i++)
      kmer_locs[i] = *(kmer_locs_p[i]);

    _parcel.push_back(KmerBatch(kmer, kmer_locs));

    TimerStop();
  }



  void PutKmerKmerLocs(const NaifBuffer & kmer, 
                       const vec<KmerLoc> & kmer_locs)
  {
    TimerStart();
    _parcel.push_back(KmerBatch(kmer, kmer_locs));
    TimerStop();
  
  }
  
};



















