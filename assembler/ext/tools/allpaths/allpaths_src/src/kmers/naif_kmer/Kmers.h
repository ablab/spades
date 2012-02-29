///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#ifndef KMERS__NAIF_KMER__KMERS_H
#define KMERS__NAIF_KMER__KMERS_H



const uint64_t hash_mixers[256] = {
  0xd5efa3673459f0b0, 0xe2bfb7185ba0dd69, 0xfa1a07dc2338c504, 0x98c417061852c314,
  0xeda1d46041d5c177, 0x1baab4f318bade9e, 0x4a5fb2b27f0deeb7, 0x6f1a17677caa2063,
  0xd715b252b9e4e1f1, 0xfb65375470d5f231, 0xe170d79fb92aa052, 0x3f45f7465e45b072,
  0xafac0820a39db6ca, 0x786b59d83ebc1190, 0x2c84e821d2b08a12, 0xe3102b6996451f94,
  0x12287d82a70522f0, 0xb2e61310d977a810, 0x1cc6abca4d22fab9, 0x2bec3eccb783463f,
  0xc0d8dd82c71e8c17, 0xb7881effb5126b6a, 0xa5196ec0bd22ac0c, 0x9ca9e51d439a098c,
  0xa9eead017747ef36, 0x382db9a06e26c439, 0xa0e5803c67269f54, 0x8045dd3fdb51cbbb,
  0x5c05113303540508, 0x01c9295e3205c67d, 0x5f86010ba2c409de, 0xd404a3c7765e6cac,
  0x060d061ffe14768c, 0x9e66199a2b73ba1a, 0x293202100157517b, 0x5640f582888bda66,
  0x8b19e5206e462186, 0xe14dc529da3d335a, 0x137af2c03d201901, 0x0e9821eddcbc3193,
  0x07c83b66be116cfc, 0x1295f0d30226bb67, 0x64ba5fa4c6225cc7, 0xcf5e2dbbd8596e02,
  0x6438f0760590adb9, 0x0ea4b5d52bf908d2, 0x970900db5a5e9f2d, 0x8383c553cc1cab25,
  0x32dabca6250e9acd, 0xf5731c6089abc80f, 0xc499189577332512, 0x9757ff50f14cf491,
  0x297befa247eed20a, 0xf49886186cf0f3f6, 0xe88e219999f8425c, 0x6dcfeeffb4d7a680,
  0xf75e86ead5726c89, 0x4e68e89565a2be44, 0xe8abcdc03629b674, 0x9f3e06e4fbe45e52,
  0x1848945188caa324, 0xffafd7e58bd60818, 0xcd765477ffaf18c7, 0xb1a34b351fc67b07,
  0x89caf47f662b82f9, 0x4cd261086f936a76, 0xe79e1141a8155546, 0x1604545d796e2dc7,
  0x9e62135f56f283fd, 0xf56ed15b4a85f474, 0x97c03af511b0cd99, 0x3a6deafa9de36e54,
  0xbc0d08e54026d814, 0x42e7c6c95f99575f, 0x3a3716c59dc18174, 0xd650ab4d9e7c6c98,
  0x51d3fad57141e500, 0xa8714699cfeb367d, 0x0bdb5d988638834a, 0x947a5003c4619ab3,
  0x6d8554a62ebc6670, 0x4944111aec8163e2, 0x853316c02eb9c69e, 0x434909d77040e3de,
  0xe0265d252de7f2d6, 0xce2a2411a15f1406, 0x2f3d0623a82e81f3, 0xf17e2a36ba4ed7f1,
  0x962bccf0ba3e0fea, 0x3a95753789c14646, 0xdd0c2aba89a10b53, 0x11f605c4ab111bea,
  0x5cbeb4fb18cd5727, 0xf6e37cb9ab4fbd61, 0x74e48d12d107c4cb, 0x8889a6d70bd66a6d,
  0xeb47b50fe976aa4d, 0x9a978fb76c4487fb, 0x19b35dce94d118ec, 0xcd844ffa60fe243f,
  0xf538557086d0b1ed, 0x0e61b473a678cec7, 0xe2d0723375cb38e6, 0x796c217b76bb1a5a,
  0xf6068cf442264326, 0xcc68a08ef5c59757, 0xef4836772cb30e21, 0xe386bae15e5287b0,
  0xa1943f2d0e879b90, 0x96fddc36a59357d6, 0x3cb6fcc9979b4974, 0xaffda49518db205b,
  0x5a8cd5139046665a, 0x11e5213ba03d1dd8, 0xfa6221fda27045df, 0x2b0020d390c04310,
  0x4764c2f3b1eb9976, 0x584916e39f374955, 0xf18cef775d391168, 0x7ca142887bdee32a,
  0x78a81eec4ce48a1c, 0x34d80c66ccd9f529, 0x988a19b170314998, 0x99f5cfd9acbdf927,
  0xe9158771de0e1147, 0xa9cd836ff0dfef95, 0x2e0e8eb8819207fc, 0x261bf49a5b348329,
  0x06bf59cafdcfc036, 0xebb28aa1739c42bc, 0x059f94b0cddbdea5, 0x8896c4d82e9f1575,
  0x61ad3892d7b924ed, 0xc0dee81042284d27, 0xaf36384cbbb3e4a9, 0x01a78f503606faf7,
  0x71b7e23aa9d31eb2, 0xaf0513c4e40c6a42, 0x0d720b5684a5d2a0, 0x308efff564dce5bd,
  0x2299578568174290, 0x7f53023ae1252add, 0xbef5332d245d2f27, 0x52e0126e9d292f78,
  0xfa6b520afcbf4fd7, 0x68f1bee2ecc71dab, 0x78c3c742a0278b60, 0x7f4ed495ca53a9fe,
  0x3c0e0ab6425eb49a, 0xa0edd6a2296327c6, 0x10c07f9a60d1cf70, 0x41a4f8f911e4ffca,
  0x5c315047c0f8222f, 0x25d292513e2d6d29, 0x00a0b34fc9329700, 0x48c6395da9e3a47a,
  0x7cf9edd26b8131b3, 0x5c82013de688c935, 0x1af1a70f6bf39dbd, 0x4aee893776ae5550,
  0x697d613d17861883, 0x7b13d3a61938d72d, 0x05092b9b436c3bd9, 0xdf9418322c7629b2,
  0xb78de1c52cdf6cc4, 0xcdc50540e39c8ab9, 0xc9fb8c49169ef287, 0xd19724ac52ba33ab,
  0x4b90b67a17abf9dd, 0x994257b59b735867, 0x83f0cd83c6ac0882, 0x57a99e383266b067,
  0x7758dcbe5b00a4f8, 0xbcd79d2117eeee25, 0xaa733b46984b3687, 0xb418981e6979c8bb,
  0x523ffecabf006184, 0x3e4dd29439137bc9, 0x0d9147a4ad1e63d6, 0xf8618bd6994fe19a,
  0xb5f7ba29f94e35b5, 0x98820a9422a88b89, 0xe6ecb942c8ef72b5, 0x0c3d6b9166f97ca1,
  0x80cd363b5b0b4f03, 0xfcff2798363c6a83, 0xb31129db2f2b1ced, 0xa30b84b0cf3ca88c,
  0x5b590ab483d326a7, 0xb824668e05fc2136, 0xb6abd4bbdd2ae04b, 0x7dc6324aa22401be,
  0x44d5761ad567b621, 0x8ea2b9553061fc07, 0x2c7eb544a61ffb20, 0x32eb5fabd43d42be,
  0x9694e44a67ffbbb5, 0x3a5af14452fca3be, 0x1a28f56dbdb178c2, 0xb902dacf35e964b4,
  0x0452c4f38b690ee8, 0x23e6b9a208d0c606, 0xb31eac3d0e9f4e09, 0x1bac76eb0b604ff0,
  0x5c58b21975710508, 0x1556622b37174df0, 0x0269ec5552fc85d2, 0x165ffee5f24a1c64,
  0x6fb1e6949a44e5a2, 0xe9ed7992c15a2c70, 0x4d0f226443d27d80, 0xcb48f7d77c0fe225,
  0x97e7048e4fbc6f0a, 0x2de837f91c31b0df, 0x402ec3b05e3f3176, 0x6fe83390e64619ee,
  0xc4fc2c09c8c06ad6, 0x3893b97a1c247177, 0x873e4c94ac60999e, 0x9465331264d28d9c,
  0xdbab61572a3f5263, 0x8a4308a2cc6bc3da, 0xd7f28efcb89f15e5, 0x9db9591583afd8b3,
  0x711543a07762034f, 0x1eeaf8f00f9ec58c, 0xdab550ea5228eabc, 0x261eec83407078d7,
  0xf5ba8246c9ec77dd, 0xd0960ac516cc191d, 0x06fc6422d742874c, 0xd856ce6c027170ca,
  0x615bbd914aa60fa0, 0x26aa4c6c821ff0b3, 0xda96bc5a1cfcc103, 0x2d729c9fb0efee9d,
  0x72fb1a325a390721, 0xb8e17faa481dd40d, 0x58d2e6518036d52f, 0x3e76d4b97911f3e9,
  0xf2dbd442a7caae2e, 0x7084366e445d965e, 0xf0e423fada826d70, 0x4718f3162f225f04,
  0xcfc5772c218b935b, 0xfdb8d41df8077140, 0x68a727db472b1a6c, 0xf47368c35131f6ea,
  0xbff6f430a06f187e, 0xca18c810d3eced2f, 0x2020138db082b04c, 0x6fa3fe05a7065bb7,
  0xbbf608e106911734, 0xe1c2981358a65a09, 0x4e49c63e52acb33d, 0xd7e404d189c6872e };






/*
class AnyKmer
{
public:
  bool     is_valid_kmer () const = 0;
  unsigned K             () const = 0;
  unsigned operator []   (const unsigned i) const = 0;
  void     set           (const unsigned i, const unsigned base) = 0;
  void     push_right    (const unsigned base) = 0;
  void     push_left     (const unsigned base) = 0;
  void     add_right     (const unsigned base) = 0;
  void     add_left      (const unsigned base) = 0;
  void     del_right     () = 0;
  void     del_left      () = 0;

  size_t   hash_64bits   () const = 0;
  bool     match         (const AnyKmer & a) const = 0;


  friend bool operator <  (const AnyKmer & a, const AnyKmer & b) = 0;
  friend bool operator <= (const AnyKmer & a, const AnyKmer & b) = 0;
  friend bool operator >  (const AnyKmer & a, const AnyKmer & b) = 0;
  friend bool operator >= (const AnyKmer & a, const AnyKmer & b) = 0;
  friend bool operator == (const AnyKmer & a, const AnyKmer & b) = 0;
  friend bool operator != (const AnyKmer & a, const AnyKmer & b) = 0;

};


*/


// -------- REC_t class --------

// _bases[0]         is the lowest  significant 2bit
// _bases[_K - 1]    is the highest significant 2bit


class Kmer29
{
  static uint64_t const K_MAX = 29;

  uint64_t _K          :  5;  // 5 bits for K 
  uint64_t _tag        :  1; 
  uint64_t _bases      : 58;  // K_MAX * 2

public:
  typedef Kmer29 kmer_type;
  Kmer29() : _K(0), _tag(0), _bases(0) {}
  explicit Kmer29(const unsigned K) : _K(K), _tag(0), _bases(0) {}

protected:
  uint64_t bases()    const { return _bases; }
  uint64_t and_mask() const { return (1ul << (_K << 1)) - 1; }

public:
  bool     is_valid_kmer()             const { return _K > 0; }
  unsigned K()                         const { return _K; }
  unsigned size()                      const { return _K; }
  unsigned tag()                       const { return _tag; }
  void     tag_set(unsigned const tag = 1)   { _tag = tag; }

  // insert a base at the end, index _K - 1, shifting everything down
  // base assumed to be < 4
  void push_right(const uint64_t base)
  { _bases >>= 2; _bases |= (base << ((_K - 1) << 1)); }
    
  // insert a base at the begining, index 0, shifting everything up
  // base assumed to be < 4
  void push_left(const uint64_t base)
  { _bases <<= 2; _bases &= and_mask(); _bases |= base; }



  // add a base at the end, index _K, incrementing K (if possible)
  // base assumed to be < 4
  void add_right(const uint64_t base)
  {
    ForceAssert(_K < K_MAX);
    _bases |= (base << (_K << 1)); 
    _K++;
  }
    
  // add a base at the beginning, index 0, shifting everything up and 
  // incrementing K (if possible)
  // base assumed to be < 4
  void add_left(const uint64_t base)
  {
    ForceAssert(_K < K_MAX);
    _bases <<= 2; 
    _bases |= base; 
    _K++; 
  }



  // delete the base at the end, index _K - 1, decrementing _K
  // base assumed to be < 4
  void del_right()
  {    
    ForceAssert(_K > 0u);
    _K--;
    _bases &= and_mask();
  }
    
  // delete the base at the begining, index 0, shifting everything down,
  // and decrementing K
  // base assumed to be < 4
  void del_left()
  {
    ForceAssert(_K > 0u);
    _K--;
    _bases >>= 2; 
  }
  

  void set(const unsigned i, const uint64_t base)
  { 
    ForceAssert(i < _K);
    unsigned const shift = (i << 1);
    _bases &= ~(3ul << shift);   // mask
    _bases |= (base << shift);   // imprint
  }

  unsigned operator[](const unsigned i) const 
  {
    ForceAssert(i < _K);
    return 3ul & (_bases >> (i << 1)); 
  }

  uint64_t hash_64bits() const
  {
    const unsigned nbits = (_K << 1);
    uint64_t hash = 0x39743f67bf1a3d0c;
    
    for (unsigned shift = 0; shift < nbits; shift += 8)
      hash ^= hash_mixers[(_bases >> shift) & 255u];
    return hash;
  }

  uint64_t hash_64bits2() const
  {
    const unsigned nbits = (_K << 1);
    uint64_t hash = 0x39745f62be1a8d0c;
    
    for (unsigned shift = 0; shift < nbits; shift += 6) {
      const unsigned shift2 = (_bases >> shift) & 63u;
      if (shift2)
        hash ^= (hash << shift2) | ((~hash) >> (64 - shift2));
    }
    return hash;
  }


  bool match(const Kmer29 & a) const 
  {
    ForceAssert(a._K == _K);
    return a._bases == _bases; 
  }

  friend bool operator<(const Kmer29 & a, const Kmer29 & b) 
  { 
    if (a._bases < b._bases) return true; 
    if (a._bases > b._bases) return false; 
    return (a._tag < b._tag);
  }

  friend bool operator>(const Kmer29 & a, const Kmer29 & b) 
  { 
    if (a._bases > b._bases) return true; 
    if (a._bases < b._bases) return false; 
    return (a._tag > b._tag);
  }

  friend bool operator==(const Kmer29 & a, const Kmer29 & b)
  { 
    return (a._bases == b._bases && a._tag == b._tag);
  }
};






// ------- hollow kmer --------

// _K               is assumed odd
// _kmer[_K / 2]    is neglected for matching purposes
class Kmer29H : public Kmer29
{
public:
  typedef Kmer29H kmer_type;
  explicit Kmer29H(const unsigned K = 0) : Kmer29(K) {}
  Kmer29H(const Kmer29 & kmer) : Kmer29(kmer) {}

protected:
  uint64_t and_mask_flanks() const { return and_mask() & ~(3ul << (K() - 1)); }
  
  uint64_t flanks() const { return bases() & and_mask_flanks(); }
  uint64_t center() const { return 3ul & (bases() >> (K() - 1)); }  // assumes K() odd

public:
  size_t hash_64bits() const 
  {
    const unsigned nbits = (K() << 1);
    const uint64_t bases = flanks();
    uint64_t hash = 0x39743f67bf1a3d0c;
    
    for (unsigned shift = 0; shift < nbits; shift += 8)
      hash ^= hash_mixers[(bases >> shift) & 255u];
    return hash;
  }


  bool match(const Kmer29H & a) const { return a.flanks() == flanks(); }

  friend bool operator==(const Kmer29H & a, const Kmer29H & b) { return a.flanks() == b.flanks(); }
  friend bool operator< (const Kmer29H & a, const Kmer29H & b) { return a.flanks() <  b.flanks(); }
  friend bool operator> (const Kmer29H & a, const Kmer29H & b) { return a.flanks() >  b.flanks(); }
};










// -------- REC_t class --------

template<class UINT_K_t, unsigned NBYTES>
class KmerTemplate
{
  static UINT_K_t const NBYTES_BASES = NBYTES - sizeof(UINT_K_t);
  static UINT_K_t const K_MAX = 4 * NBYTES_BASES;

  UINT_K_t _tag   : 1;
  UINT_K_t _K     : 8 * sizeof(UINT_K_t) - 1;  
  uint8_t  _bases[NBYTES_BASES]; 

  unsigned n_bytes_used() const { return (_K + 3u) >> 2; }

public:
  typedef KmerTemplate kmer_type;
  KmerTemplate() : _tag(0), _K(0) { memset(_bases, 0, NBYTES_BASES); }
  explicit KmerTemplate(const UINT_K_t K = 0) : _tag(0), _K(K) 
  { memset(_bases, 0, NBYTES_BASES); }

  bool     is_valid_kmer() const { return _K > 0; }
  UINT_K_t K()             const { return _K; }
  UINT_K_t size()          const { return _K; }
  unsigned tag()           const { return _tag; }
  void tag_set(unsigned const tag = 1) { _tag = tag; }


  // (0)...(K-1) -> (1)...(K-1)(base)
  //
  // insert a base at the end, index _K - 1, shifting everything down
  // base assumed to be < 4
  void push_right(const uint8_t base)
  {
    const unsigned i1 = n_bytes_used() - 1;

    for (unsigned i = 0; i < i1; i++)
      _bases[i] = (_bases[i] >> 2) | (_bases[i+1] << 6);

    _bases[i1] = (_bases[i1] >> 2) | (base << (((_K - 1u) & 3u) << 1));
  }

    
  // (0)...(K-1) -> (base)(0)...(K-2)
  //
  // insert a base at the begining, index 0, shifting everything up
  // base assumed to be < 4
  void push_left(const uint8_t base)
  {
    const unsigned i1 = n_bytes_used() - 1;

    for (unsigned i = i1; i > 0; i--)
      _bases[i] = (_bases[i] << 2) | (_bases[i-1] >> 6);
    _bases[0] = (_bases[0] << 2) | base;

    if(_K & 3u)
      _bases[i1] &= ~(3u << ((_K & 3u) << 1));  // zero out base to drop
  }




  // (0)...(K-1) -> (0)...(K-1)(base)
  //
  // insert a base at the end, index _K, and increment K
  // base assumed to be < 4
  void add_right(const uint8_t base)
  {
    ForceAssert(_K < KmerTemplate::K_MAX);
    _K++;
    const unsigned i1 = n_bytes_used() - 1;
    _bases[i1] = _bases[i1] | (base << (((_K - 1u) & 3u) << 1));
  }

    
  // (0)...(K-1) -> (base)(0)...(K-1)
  //
  // insert a base at the begining, index 0, shifting everything up
  // and incrementing K
  // base assumed to be < 4
  void add_left(const uint8_t base)
  {
    ForceAssert(_K < KmerTemplate::K_MAX);
    _K++;
    const unsigned i1 = n_bytes_used() - 1;

    for (unsigned i = i1; i > 0; i--)
      _bases[i] = (_bases[i] << 2) | (_bases[i-1] >> 6);
  
    _bases[0] = (_bases[0] << 2) | base;
  }






  // (0)...(K-1) -> (0)...(K-2) 
  //
  // delete the base at the end, index _K - 1, decrementing K
  // base assumed to be < 4
  void del_right()
  {
    ForceAssert(_K > 0u);
    const unsigned i1 = n_bytes_used() - 1;
    _K--;
    _bases[i1] &= ~(3u << ((_K & 3u) << 1)); // zero out last base
  }

    
  // (0)...(K-1) -> (1)...(K-1)
  //
  // insert a base at the begining, index 0, shifting everything up
  // base assumed to be < 4
  void del_left()
  {
    ForceAssert(_K > 0u);
    const unsigned i1 = n_bytes_used() - 1;

    for (unsigned i = 0; i < i1; i++)
      _bases[i] = (_bases[i] >> 2) | (_bases[i+1] << 6);
    
    _K--;
  }


  void set(const UINT_K_t i, const uint8_t base)
  { 
    ForceAssert(i < _K);
    uint8_t * p_i = _bases + (i >> 2);
    unsigned const shift = ((i & 3u) << 1);
    *p_i &= ~(3u << shift);    // mask
    *p_i |= (base << shift);   // imprint
  }

  unsigned operator[](const UINT_K_t i) const 
  {
    ForceAssert(i < _K);
    uint8_t const byte = _bases[i >> 2];
    unsigned const shift = ((i & 3u) << 1);
    return 3ul & (byte >> shift);
  }

  uint64_t hash_64bits() const
  {
    uint64_t hash = 0x39743f67bf1a3d0c;
    const uint8_t * p_hi = _bases + n_bytes_used();
    for (const uint8_t * p = _bases; p != p_hi; p++)
      hash ^= hash_mixers[*p];
    return hash;
  }



  bool match(const KmerTemplate & a) const 
  { 
    ForceAssert(a._K == _K);
    unsigned const nbytes = n_bytes_used();
    for (unsigned i = 0; i != nbytes; i++)
      if (a._bases[i] != _bases[i]) return false;
    return true;
  }

  friend bool operator==(const KmerTemplate & a, const KmerTemplate & b)
  { 
    ForceAssert(a._K == b._K);
    unsigned const nbytes = a.n_bytes_used();
    for (unsigned i = 0; i != nbytes; i++)
      if (a._bases[i] != b._bases[i]) return false;
    return (a._tag == b._tag);
  }

  friend bool operator<(const KmerTemplate & a, const KmerTemplate & b) 
  { 
    ForceAssert(a._K == b._K);
    unsigned const nbytes = a.n_bytes_used();
    for (unsigned i = nbytes; i != 0; i--) {
      if (a._bases[i-1] < b._bases[i-1]) return true;
      if (a._bases[i-1] > b._bases[i-1]) return false;
    }
    return (a._tag < b._tag);
  }

  friend bool operator>(const KmerTemplate & a, const KmerTemplate & b) 
  { 
    ForceAssert(a._K == b._K);
    unsigned const nbytes = a.n_bytes_used();
    for (unsigned i = nbytes; i != 0; i--) {
      if (a._bases[i-1] > b._bases[i-1]) return true;
      if (a._bases[i-1] < b._bases[i-1]) return false;
    }
    return (a._tag > b._tag);
  }
  
 
};



// -------- REC_t class --------


template<class KMER_t>
class KmerKmerQual : public KMER_t
{
  uint8_t _q;

public:
  KmerKmerQual(const unsigned K, const unsigned q = 0) : KMER_t(K), _q(q) {}

  unsigned qual() const { return _q; }
  void qual(const uint8_t q) { _q = q; }

  /*
  friend bool operator<(const KmerKmerQual & a, const KmerKmerQual & b) 
  { 
    if (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b)) return true;
    if (static_cast<const KMER_t &>(b) < static_cast<const KMER_t &>(a)) return false;
    return (a.qual() < b.qual());
  }
  */
};

template<class REC_t>
bool kmer_qual_lt(const REC_t & a, const REC_t & b) { return a.qual() < b.qual(); }

template<class REC_t>
bool kmer_qual_gt(const REC_t & a, const REC_t & b) { return a.qual() > b.qual(); }







// -------- REC_t class --------



template<class KMER_t>
class KmerKmerFreq : public KMER_t
{
  uint32_t _freq;
public:
  KmerKmerFreq(const KMER_t kmer, const unsigned freq = 0) : KMER_t(kmer), _freq(freq) {}
  KmerKmerFreq(const unsigned K, const unsigned freq = 0) : KMER_t(K), _freq(freq) {}
  
  uint32_t freq() const { return _freq; }
  void     set_freq(const unsigned freq) { _freq = freq; }
};


template<class REC_t>
bool kmer_freq_lt(const REC_t & a, const REC_t & b) { return a.freq() < b.freq(); }

template<class REC_t>
bool kmer_freq_gt(const REC_t & a, const REC_t & b) { return a.freq() > b.freq(); }









// ------ a kmer location in a BaseVecVec ------

class BVLoc
{
  uint64_t _ibv : 34;   // up to 2^34 =   16G (at 100b/read => 400+ GB)
  uint64_t _ib  : 28;   // up to 2^28 =  256M
  uint64_t _dir :  2;   // 0: fw  1: rc  2: palindrome  3: invalid
  //uint64_t _fw  :  1;   // 0: fw  1: rc 
  //uint64_t _pal :  1;   // palidrome?

public:
  //BVLoc() : _ibv(0), _ib(0), _fw(0), _pal(0) {} 
  BVLoc() : _ibv(0), _ib(0), _dir(3) {} 

  void set_ibv       (const size_t ibv) { _ibv = ibv; }
  void set_ib        (const size_t ib)  { _ib  = ib; }
  void set_fw        (const bool fw)    { _dir = (fw ? 0 : 1); }
  void set_rc        (const bool rc)    { _dir = (rc ? 1 : 0); }
  void set_palindrome(const bool pal)   { if (pal) _dir = 2; }

  size_t ibv()           const { return _ibv; }
  size_t ib()            const { return _ib; }
  bool   is_fw()         const { return _dir == 0; }
  bool   is_rc()         const { return _dir == 1; }
  bool   is_palindrome() const { return _dir == 2; }
  bool   is_valid_loc()  const { return _dir != 3; }

  String to_str() const
  {
    return (String((is_palindrome() ? "dir= pal" : 
                    is_fw() ? "dir= fw " : 
		    is_rc() ? "dir= rc " : "dir= bad")) +
            " ibv/ib=" + ToString(_ibv) +
            " " + ToString(_ib));
  }
  
  friend
  bool operator < (const BVLoc & a, const BVLoc & b)
  { 
    if (a._ibv < b._ibv) return true;
    if (a._ibv > b._ibv) return false;
    return (a._ib < b._ib);
  }
};




// -------- REC_t class --------


// ------ basically a pair of a kmer and its location ------


template<class KMER_t>  
class KmerBVLoc : public KMER_t, public BVLoc
{
public:
  KmerBVLoc(const KMER_t & kmer) : KMER_t(kmer), BVLoc() {}
  explicit KmerBVLoc(const unsigned K = 0) : KMER_t(K), BVLoc() {}

  friend
  bool operator < (const KmerBVLoc & a, const KmerBVLoc & b)
  { 
    if (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b)) return true;
    if (static_cast<const KMER_t &>(b) < static_cast<const KMER_t &>(a)) return false;
    return (static_cast<const BVLoc &>(a) < static_cast<const BVLoc &>(b));
  }

};



// ------- class that represents a kmer in both FW and RC directions --------
template<class KMER_t>
class KmerFWRC
{
  KMER_t _fw;
  KMER_t _rc;


public:
  KmerFWRC(const KMER_t & kmer) : _fw(kmer), _rc(reverse_complement(kmer)) {}
  explicit KmerFWRC(const unsigned K = 0) : _fw(K), _rc(K) 
  {
    _rc = reverse_complement(_fw);
  }
  
  unsigned K()    const { return _fw.K(); }
  unsigned size() const { return _fw.K(); }
  
  bool is_canonical_fw() const { return _fw < _rc; }
  bool is_canonical_rc() const { return _rc < _fw; }
  bool is_palindrome()   const { return _fw == _rc; }
  
  const KMER_t & fw()        const { return _fw; }
  const KMER_t & rc()        const { return _rc; }
  const KMER_t & canonical() const { return _fw < _rc ? _fw : _rc; }

  KMER_t & fw()        { return _fw; }
  KMER_t & rc()        { return _rc; }
  KMER_t & canonical() { return _fw < _rc ? _fw : _rc; }

  void push_right(const unsigned base)  // no checking if in [0,3]!!!
  {
    _fw.push_right(base);
    _rc.push_left(3ul ^ base); // the complement of base
  }
    
  void push_left(const unsigned base)  // no checking if in [0,3]!!!
  {
    _fw.push_left(base);
    _rc.push_right(3ul ^ base); // the complement of base
  }
};




// -------- class to parse kmers from reads --------
// 
// this looks like an iterator but I chose not to make it an iterator for simplicity.
//
template<class BASE_VEC_t, class KMER_t>
class SubKmers : public KmerFWRC<KMER_t>
{
private:
  const BASE_VEC_t & _bv;
  const unsigned     _nb;
  const unsigned     _K;

  unsigned           _ib;  // index of base after current kmer

public:
  SubKmers(const unsigned K, const BASE_VEC_t & bv, const size_t ib0 = 0) :
    KmerFWRC<KMER_t>(K), _bv(bv), _nb(_bv.size()), _K(K)
  {
    if (_nb < ib0 + _K) // ---- BaseVec too short.  Signal we're done.
      _ib = _nb + 1;  
    else 
      for (_ib = ib0; _ib < ib0 + _K; _ib++) // ---- build 1st kmer.
        push_right(_bv[_ib]);
  }

  bool not_done() const { return (_ib <= _nb && _ib >= _K); }

  void next()
  {
    if (_ib < _nb) push_right(_bv[_ib]);
    _ib++;
  }

  void previous()
  {
    if (_ib > _K) push_left(_bv[_ib - _K - 1]);
    _ib--;
  }

  size_t   index_start()         const { return _ib - _K; }
  size_t   index_stop()          const { return _ib; }
  size_t   index_start_reverse() const { return _nb - _ib; }
  size_t   index_stop_reverse()  const { return _nb - _ib + _K; }
  size_t   n_kmers()             const { return _nb - _K + 1; }
};






// ---- general kmer frequency validator class ----
// a bit obsolete


class ValidatorKmerFreq
{
public:
  const size_t f_min;
  const size_t f_max;
  ValidatorKmerFreq(const size_t f_min = 0, const size_t f_max = 0) : f_min(f_min), f_max(f_max) {}
  template<class REC_t>
  bool operator ()(const REC_t & krec) const { return ((!f_min || krec.freq() >= f_min) && 
						       (!f_max || krec.freq() <= f_max)); }
};


class ValidatorKmerFreqAll : public ValidatorKmerFreq
{
public: 
  ValidatorKmerFreqAll() : ValidatorKmerFreq(0, 0) {}
};




// ---------- Validator class ---------------

class Validator
{
  size_t _min, _max;
public:
  Validator(const size_t min = 0, const size_t max = 0) : _min(min), _max(max) {}
  bool operator()(const size_t x) const 
  { 
    return (((_min == 0 && x > 0) || (_min > 0 && x >= _min)) && 
            (_max == 0 || x <= _max)); 
  }
};


// ------ kmer template functions ------





inline
String hieroglyph(unsigned base)
{
  if (base == 0) return "^";
  if (base == 1) return "(";
  if (base == 2) return "-";
  if (base == 3) return ".";
  ForceAssert(base < 4u);
  return "?";
}


template<class KMER_t>
String hieroglyphs(const KMER_t & kmer)
{
  String h = "";
  unsigned K = kmer.size();
  for (unsigned i = 0; i != K; i++)
    h += hieroglyph(kmer[i]);
  return h;
}




template<class KMER_t>
inline 
KMER_t reverse_complement(const KMER_t & kmerFW) 
{
  const unsigned K = kmerFW.K();
  KMER_t kmerRC = kmerFW;          // copy, because KMER_t might have other stuff 
  for (unsigned i = 0; i != K; i++) 
    kmerRC.push_left(3u ^ kmerFW[i]);   // push bases from the begining
  
  return kmerRC;
}



template<class KMER_t>
inline
KMER_t canonical(const KMER_t & kmerFW) 
{
  const unsigned K = kmerFW.K();
  KMER_t kmerRC = kmerFW;          // copy, because KMER_t might have other stuff 
  for (unsigned i = 0; i != K; i++) 
    kmerRC.push_left(3u ^ kmerFW[i]);   // push bases from the begining
  
  return (kmerFW < kmerRC) ? kmerFW : kmerRC;
}




template <class KMER_t>
inline
bool operator<=(const KMER_t & a, const KMER_t & b) { return !(a > b); }

template <class KMER_t>
inline
bool operator>=(const KMER_t & a, const KMER_t & b) { return !(a < b); }

template <class KMER_t>
inline
bool operator!=(const KMER_t & a, const KMER_t & b) { return !(a == b); }






typedef KmerTemplate< uint8_t,  16>   Kmer60;   //  60 =  64 - 4 * 1
typedef KmerTemplate< uint8_t,  32>   Kmer124;  // 124 = 128 - 4 * 1
typedef KmerTemplate<uint16_t,  64>   Kmer248;  // 248 = 256 - 4 * 2 
typedef KmerTemplate<uint16_t, 128>   Kmer504;  // 504 = 512 - 4 * 2


TRIVIALLY_SERIALIZABLE(Kmer29);
TRIVIALLY_SERIALIZABLE(Kmer60);
TRIVIALLY_SERIALIZABLE(Kmer124);
TRIVIALLY_SERIALIZABLE(Kmer248);
TRIVIALLY_SERIALIZABLE(Kmer504);





#endif
