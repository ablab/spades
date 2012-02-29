///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ---------------------------------
// NaifTimer
// ---------------------------------

class NaifTimer 
{
private:
  double _t;

public:
  NaifTimer(const double t0 = 0) : _t(t0) {}
  NaifTimer(const NaifTimer & t) : _t(t._t) {}
  NaifTimer & operator=(const NaifTimer & t) { _t = t._t; return *this; }

  inline void Start() { _t -= WallClockTime(); }
  inline void Stop()  { _t += WallClockTime(); }
  inline double GetTimer() const { return _t; }
};



// ---------------------------------
// NaifBuffer
// ---------------------------------

class NaifBuffer
{
private:
  size_t           _n_bits;
  size_t           _n_bytes;
  uint8_t        * _buf;
  const uint8_t  * _cbuf;
  bool             _owned;

  inline void SetSizes() 
  {
    _n_bytes = (_n_bits + 7) / 8;
  } 

public:
  // ---- constructor 
  NaifBuffer(const size_t n_bits, uint8_t * const buf) 
    : _n_bits(n_bits), 
      _n_bytes((n_bits + 7) / 8), 
      _buf(buf), 
      _cbuf(buf),
      _owned(false) {}

  // ---- constructor from a pointer to const
  NaifBuffer(const size_t n_bits, uint8_t const * cbuf) 
    : _n_bits(n_bits), 
      _n_bytes((n_bits + 7) / 8), 
      _buf(0),     // _buf is not const
      _cbuf(cbuf),
      _owned(false) {}

  // ---- constructor
  NaifBuffer(const size_t n_bits) 
    : _n_bits(n_bits), 
      _n_bytes((n_bits + 7) / 8)
  {
    _cbuf = _buf = new uint8_t[_n_bytes];
    memset(_buf,0,_n_bytes);
    _owned = true;
  }

  // ---- copy constructor
  NaifBuffer(NaifBuffer const & b) 
  {
    _n_bits = b._n_bits;
    _n_bytes = b._n_bytes;
    _buf = new uint8_t[_n_bytes]; 
    _cbuf = _buf; 
    _owned = true;
    memcpy(_buf, b._cbuf, _n_bytes);
  }

  // ---- destructor
  ~NaifBuffer() { if (_owned) delete[] _buf; }
  
  // ---- =
  NaifBuffer & operator=(const NaifBuffer & b)
  {
    ForceAssertEq(_n_bits, b._n_bits);
    memcpy(_buf, b._cbuf, _n_bytes);
    return *this;
  } 

  // ---- accessors
  inline size_t NumBytes() const { return _n_bytes; }
  inline size_t NumBits()  const { return _n_bits; }

  inline       uint8_t * UInt8Pointer()       { return _buf; }
  inline const uint8_t * UInt8Pointer() const { return _cbuf; }
  
  inline char * CharPointer() 
  { return reinterpret_cast<char *>(_buf); }

  inline const char * CharPointer() const 
  { return reinterpret_cast<char const *>(_cbuf);  }


  // ---- cmp
  inline int Compare(const NaifBuffer & b) const
  {
    ForceAssertEq(_n_bits, b._n_bits);
    return memcmp(_cbuf, b._cbuf, _n_bytes);
  } 
  
  inline void Reset(const uint8_t val = 0) 
  {
    memset(_buf, val, _n_bytes);
  }

  void BinaryStreamRead(istream & is) 
  { 
    ForceAssertGt(_n_bytes, 0u);
    ForceAssert(_buf != 0);
    is.read(reinterpret_cast<char*>(_buf), _n_bytes);  
  }

  void BinaryStreamWrite(ostream & os) const 
  { 
    ForceAssertGt(_n_bytes, 0u);
    os.write(reinterpret_cast<char const *>(_cbuf), _n_bytes);
  }
};


// ---- <
inline bool operator<(const NaifBuffer & a, const NaifBuffer & b)
{  return (a.Compare(b) < 0); } 

// ---- >
inline bool operator>(const NaifBuffer & a, const NaifBuffer & b)
{  return (a.Compare(b) > 0); } 
  
// ---- ==
inline bool operator==(const NaifBuffer & a, const NaifBuffer & b)
{  return (a.Compare(b) == 0); } 






/*
class NaifKmer : public NaifBuffer
{
  NaifKmer(const BaseVec & read, const size_t K, const size_t offset)
    : NaifBuffer(K << 1)
  {
    Reset();
    uint8_t buf = UInt8Pointer();
    for (size_t ib = 0; ib != K; ib++) {
      
      
      size_t i = ib >> 2;  //   ib / 4
      size_t j = ib & 3;   //   ib % 4
      buf[i] 
    }
  }
 





};


*/



// ----------------------------------------------------------------------------------
//
// KmerLoc
//
//   A KmerLoc stores the address of a kmer in a BaseVecVec in the form of 
//   a read_ID, a position, and the kmer orientation (forward or reverse complement).
//
//   In this implementation, the kmer orientation is associated with the sign of 
//   the position:
//
//   FW: _pos > 0   and   position = _pos
//   RC: _pos < 0   and   position = ~_pos = -_pos - 1
//
// ----------------------------------------------------------------------------------
class KmerLoc 
{
private:
  uint64_t _read_ID;
  int32_t  _pos;
  
public:
  // constructors
  KmerLoc() {}
  KmerLoc(const size_t read_ID, const int32_t base_ID, const bool RC = false) 
    : _read_ID(read_ID) 
  {
    _pos = ((RC) ? ~base_ID : base_ID); 
  }


  // destructor
  ~KmerLoc() {}

  // accessors
  inline uint64_t GetReadID()      const { return _read_ID; }
  inline int32_t  GetSignedPos()   const { return _pos; }
  inline size_t   GetUnsignedPos() const { return (_pos < 0) ? ~_pos : _pos; } // ~x = -x-1
  inline size_t   GetBaseID()      const { return (_pos < 0) ? ~_pos : _pos; } // ~x = -x-1
  inline bool     IsRC()           const { return (_pos < 0); }

  // mutators
  inline KmerLoc & SetReadID(const size_t read_ID) { _read_ID = read_ID; return *this; }
  inline KmerLoc & SetSignedPos(const int32_t pos)   { _pos = pos; return *this; }

  inline KmerLoc GetRCKmerLoc() const { return KmerLoc(_read_ID, ~_pos); }
  inline KmerLoc GetRC() const { return GetRCKmerLoc(); } // backward compatibility


  void BinaryStreamRead(istream & is) 
  {
    is.read(reinterpret_cast<char *>(&_read_ID), sizeof(uint64_t));
    is.read(reinterpret_cast<char *>(&_pos), sizeof(int32_t));
  }

  void BinaryStreamWrite(ostream & os) const 
  {
    os.write(reinterpret_cast<const char *>(&_read_ID), sizeof(uint64_t));
    os.write(reinterpret_cast<const char *>(&_pos), sizeof(int32_t));
  }
  
  // output
  void Print(ostream& out, bool alt = false ) const
  {
    if ( alt ) {
      out << "[" << _read_ID
	  << " " << ( IsRC( ) ? "+" : "-" )
	  << " " << GetUnsignedPos( )
	  << "]";
    }
    else {
      out << setw(12) << _read_ID << " " 
	  << setw(12) << GetUnsignedPos() << " " 
	  << (IsRC() ? "(RC)" : "(FW)") << endl;
    }
  }
  
  inline static bool KmerLocEQ(const KmerLoc & a, const KmerLoc & b)
  { 
    return (a._read_ID == b._read_ID && a._pos == b._pos);
  }

  // for sorting purposes
  inline static bool ReadIDSignedPosLT(const KmerLoc & a, const KmerLoc & b)
  { 
    if (a._read_ID < b._read_ID) return true;
    if (a._read_ID > b._read_ID) return false;
    return (a._pos < b._pos);
  }

  inline static bool ReadIDUnsignedPosLT(const KmerLoc & a, const KmerLoc & b)
  { 
    if (a._read_ID < b._read_ID) return true;
    if (a._read_ID > b._read_ID) return false;
    return (a.GetUnsignedPos() < b.GetUnsignedPos());
  }

  inline friend bool operator == (const KmerLoc & a, const KmerLoc & b)
  { 
    return (a._read_ID == b._read_ID && a._pos == b._pos);
  }

};





// ----------------------------------------------------------------------------------
//
// KmerRecord 
//
//   A KmerRecord stores a kmer, its read_ID, and its signed position.
//
//   Originally designed to replace kmer_record<K,2> which stores read_IDs as int32
//   as opposed to int64 (which is the future! - as of Sep 2009)
//
// ----------------------------------------------------------------------------------

template<int K> 
class KmerRecord : public KmerLoc 
{
private:
  
  // static data
  static const size_t _n_bytes = (K + 3) / 4;

  // object data
  uint8_t _kmer[_n_bytes];

public:
  KmerRecord(const BaseVec & kmer_bv) { SetKmerBaseVec(kmer_bv); }

  KmerRecord(const BaseVec & kmer_bv, const KmerLoc & loc) 
    : KmerLoc(loc)
  {
    SetKmerBaseVec(kmer_bv);
  }

  KmerRecord(const BaseVec & kmer_bv, 
             const int64_t read_ID, 
             const int32_t base_ID, 
             const bool RC) 
    : KmerLoc(read_ID, base_ID, RC)
  {
    SetKmerBaseVec(kmer_bv);
  }
  

  KmerRecord() {}

  
  inline void SetKmerBaseVec(const BaseVec & kmer_bv) 
  { kmer_bv.extractBaseBits(_kmer, _n_bytes); }

  inline void GetKmerBaseVec(BaseVec * kmer_bv) const 
  { kmer_bv->assignBaseBits(K, _kmer); }
  
  inline size_t KmerNumBytes() const 
  { return _n_bytes; }

  inline NaifBuffer Kmer()
  { return NaifBuffer(2*K, reinterpret_cast<uint8_t *>(_kmer)); }

  inline const NaifBuffer Kmer() const 
  { return NaifBuffer(2*K, reinterpret_cast<uint8_t const *>(_kmer)); }
  
  size_t GetKmerNumGC() const 
  {
    BaseVec bv;
    bv.assignBaseBits(K, _kmer);
    size_t n_gc = 0;
    for (size_t i = 0; i != K; i++)
      if (bv[i] == 1 || bv[i] == 2) n_gc++;
    return n_gc;
  }


  friend bool operator == (const KmerRecord & a, const KmerRecord & b)
  { return (a._loc == b._loc && 0 == memcmp(a._kmer, b._kmer, _n_bytes)); }

  static bool KmersEQ(const KmerRecord & a, const KmerRecord & b)
  { return (memcmp(a._kmer, b._kmer, _n_bytes) == 0); }
  
  static bool KmersLT(const KmerRecord & a, const KmerRecord & b)
  { return (memcmp(a._kmer, b._kmer, _n_bytes) < 0); }

  static bool KmersGT(const KmerRecord & a, const KmerRecord & b)
  { return  (memcmp(a._kmer, b._kmer, _n_bytes) > 0); }

  static bool KmerReadIDSignedPosLT(const KmerRecord & a, const KmerRecord & b)
  {
    int cmp = memcmp(a._kmer, b._kmer, _n_bytes);
    if (cmp < 0) return true;
    if (cmp > 0) return false;
    return KmerLoc::ReadIDSignedPosLT(a, b);
  }




  inline void PrintKmer(std::ostream& out = cout) 
  { 
    out << "[";
    for (size_t i = 0; i != _n_bytes; i++)
      out << setw(16) << hex << static_cast<uint16_t>(_kmer[i]) << " ";
    out << "]" << dec; 
  }



};







// ----------------------------------------------------------------------------------
//
// KmerBatch
//
//   A KmerBatch stores a kmer and vector of all all its occurences (vec<KmerLoc>).
//
//   The KmerBatch is not templatized and only stores the kmer once.
//
// ----------------------------------------------------------------------------------


class KmerBatch : public vec<KmerLoc> 
{
private:
  NaifBuffer _kmer;
  
public:
  KmerBatch(const size_t K) : vec<KmerLoc>(), _kmer(2*K) {} 

  template<size_t K>
  KmerBatch(typename vec<KmerRecord<K> >::const_iterator const & kr_begin, 
            typename vec<KmerRecord<K> >::const_iterator const & kr_end) 
    : vec<KmerLoc>(), _kmer(2*K)
  {
    _kmer = (*kr_begin).Kmer();
    typedef vec<KmerLoc>::const_iterator KLIter;
    for (KLIter it = kr_begin; it != kr_end; it++)
      (*this).push_back(*it);
  }


  KmerBatch(const NaifBuffer & kmer, 
            const vec<KmerLoc> & kmer_locs)
    : vec<KmerLoc>(kmer_locs), _kmer(kmer) {}

  ~KmerBatch() {}


  inline size_t GetK() const { return (_kmer.NumBits() >> 1); } // K = n_bits / 2
  inline size_t GetKmerFreq()   const { return size(); }

  inline       NaifBuffer & Kmer()       { return _kmer; }
  inline const NaifBuffer & Kmer() const { return _kmer; }

  
  inline void SetKmerBaseVec(const BaseVec & kmer_bv) 
  { kmer_bv.extractBaseBits(_kmer.CharPointer(), _kmer.NumBytes()); }

  inline void GetKmerBaseVec(BaseVec * kmer_bv) const 
  { kmer_bv->assignBaseBits(GetK(), _kmer.CharPointer()); }




  size_t GetKmerFreqRC() const 
  { 
    size_t kf = 0;
    typedef vec<KmerLoc>::const_iterator KLVIter;
    for (KLVIter it = begin(); it != end(); it++)
      if ((*it).IsRC()) 
        kf++;
    return kf; 
  }

  size_t GetKmerFreqMinFWRC() const 
  {
    const size_t kfRC = GetKmerFreqRC();
    const size_t kfFW = GetKmerFreq() - kfRC;
    return (kfFW < kfRC) ? kfFW : kfRC;
  }



  size_t GetKmerFreqPosLT(const size_t mid_pos) const 
  { 
    size_t kf = 0;
    typedef vec<KmerLoc>::const_iterator KLVIter;
    for (KLVIter it = begin(); it != end(); it++)
      if ((*it).GetUnsignedPos() < mid_pos) 
        kf++;
    return kf; 
  }

  bool IsPalindrome() const  // requires some processing...
  {
    BaseVec kmer_bv;
    kmer_bv.assignBaseBits(GetK(), _kmer.CharPointer());
    return (kmer_bv.Canonicalize() == BaseVec::palindrome);
  }
    


  void PrintBinnedPos(ostream& out, const size_t read_len = 101) const
  {
    const size_t n_bins = 9;  // last bin is for overflow
    vec<size_t> pos_fw(n_bins);
    vec<size_t> pos_rc(n_bins);
    const size_t n_pos = read_len - GetK() + 1;

    for (size_t i = 0; i != size(); i++) {
      size_t pos = (*this)[i].GetUnsignedPos();

      size_t i_bin = (pos < n_pos) ? pos * (n_bins - 1) / n_pos : n_bins - 1;
      
      if ((*this)[i].IsRC())
        pos_rc[i_bin]++;
      else
        pos_fw[i_bin]++;
    }

    cout << " bin          FW          RC" << endl;
    for (size_t i = 0; i != n_bins; i++) 
      cout << setw(4) << i << "  " 
           << setw(10) << pos_fw[i] << "  " 
           << setw(10) << pos_rc[i] << endl; 
  }


  void PrintKmerLocs(ostream& out, const size_t max_freq = 0)
  {
    if (max_freq) {
      size_t n = (max_freq < size()) ? max_freq : size();
      for (size_t i = 0; i < n; i++)
        (*this)[i].Print(out);
      
      cout << ((max_freq < size()) ? "..." : "") << endl;
    }
  }
  
  void PrintKmerStats(ostream& out, const size_t read_len = 101) const
  {
    const size_t f = GetKmerFreq();

    const size_t nsig_fr = 0.5 + NumSigmaBinomial(f, GetKmerFreqRC());

    const size_t n_pos = read_len - GetK() + 1;
    size_t f_beg = 0;
    bool good = true;

    for (size_t i = 0; i != size(); i++) {
      size_t pos = (*this)[i].GetUnsignedPos();
      if (pos >= n_pos) good = false;
      if (pos < n_pos/2) f_beg++;
    }
    
    const size_t nsig_be = 0.5 + NumSigmaBinomial(f, f_beg);
    

    out << GetK() << "-mer [ " 
        << GetKmerFreq() << " kmers ["
        << GetKmerFreqMinFWRC() << " FWRC] | "
        << nsig_fr << " sig_FWRC";

    if (good) out << " | " << nsig_be << " sig_BE";
    out << " ]";
  }

  // output 
  void Print(ostream& out, const size_t read_len = 101) const
  {
    BaseVec kmer_bv;
    kmer_bv.assignBaseBits(GetK(), _kmer.CharPointer());

    PrintKmerBaseVec(kmer_bv, 3, out);
    out << "   ";
    PrintKmerStats(out, read_len);
  }



  void Print(ostream& out, const QualNibbleVec & qv, const size_t read_len = 101) const
  {
    BaseVec kmer_bv;
    kmer_bv.assignBaseBits(GetK(), _kmer.CharPointer());

    PrintKmerBaseVec(kmer_bv, qv, 3, out);
    out << "   ";
    PrintKmerStats(out, read_len);
  }

};


