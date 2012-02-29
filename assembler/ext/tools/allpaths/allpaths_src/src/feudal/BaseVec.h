///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BASEVECTOR_H
#define BASEVECTOR_H

#include "feudal/BinaryStreamTraits.h"
#include "feudal/FieldVec.h"
#include "feudal/Mempool.h"
#include "dna/Bases.h"
#include "system/Assert.h"
#include <algorithm>
#include <iterator>
#include <limits>
#include <utility>
#include <ostream>
#include <vector>

// TODO: it would be really great to get these out of the interface
#include "String.h"
#include "Vec.h"

/// generate an out stream based on the boundaries defined by the iterators
template <class Itr>
inline void PrintBasesIter(std::ostream& out,
                           Itr head, Itr const& tail, int breakCol)
{
    unsigned int cnt = 0;
    for ( ; head != tail; ++head )
    { out << Base::val2Char(*head);
      if ( !(++cnt % breakCol) ) out << '\n'; }
    if ( cnt % breakCol ) out << '\n';
}


class BaseVec : public FieldVec<2, MempoolAllocator<unsigned char> >
{
public:
    class const_rc_iterator
    : public ItrBase,
      public IteratorBase<const_rc_iterator,size_type,difference_type>
    {
    public:
        const_rc_iterator() : mpBaseVec(0), mLast(~0) {}
        const_rc_iterator( BaseVec const* pBaseVec, size_type pos )
        : IteratorBase<const_rc_iterator,size_type,difference_type>(pos),
          mpBaseVec(pBaseVec), mLast(pBaseVec->size()-1) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const
        { return BaseVec::Complement((*mpBaseVec)[mLast-this->pos()]); }

        value_type operator[]( difference_type idx ) const
        { return BaseVec::Complement((*mpBaseVec)[mLast-(this->pos()+idx)]); }

    protected:
        BaseVec const* mpBaseVec;
        size_type mLast;
    };

    typedef const_rc_iterator const_reverse_complement_iterator;
    typedef MempoolAllocator<unsigned char> Alloc;
    typedef FieldVec<2, MempoolAllocator<unsigned char> > Base_;

    //
    // Constructors
    //

    BaseVec() { }
    BaseVec(Alloc const& alloc) : Base_(alloc) { }

    /// SetToSubOf constructor.
    BaseVec( const BaseVec& b, const size_type start_pos, const size_type len )
    { SetToSubOf(b,start_pos,len); }

    /// sized constructor
    explicit BaseVec( size_type sz, size_type extra = 0,
                      Alloc const& alloc = Alloc() ) :
        Base_(sz, value_type(), sz+extra, alloc)
    { }

    /// Initialize from kmer.  Buffer has an odd arrangement:  left-most base
    /// is in lowest two bits.
    BaseVec( unsigned int k, unsigned char const* buf )
    { assignBaseBits(k,buf); }

    /// Copy constructor
    BaseVec(BaseVec const& bv) :
        Base_(bv)
    { }

    template <class Itr>
    BaseVec( Itr first, Itr const& last,
             size_type cap = 0, Alloc const& alloc = Alloc() )
    : Base_(first,last,cap,alloc) {}

    template <class Itr, class Mapper>
    BaseVec( Itr first, Itr const& last, Mapper mapper,
             size_type cap = 0, Alloc const& alloc = Alloc() )
    : Base_(first,last,mapper,cap,alloc) {}

    enum ConstructorFromStringBehavior
    { BASE_REQUIRED, MAP_N_TO_A };

    // Construct from a DNA String
    explicit BaseVec( const String& s,
            const ConstructorFromStringBehavior behavior = BASE_REQUIRED )
    { if ( behavior == BASE_REQUIRED )
        assign(s.begin(),s.end(),CharToBaseMapper());
    else
        assign(s.begin(),s.end(),GenCharToRandomBaseMapper()); }

    //
    // Accessors
    //

    /// returns size (number of bases) as int.  inherently unsafe unless you
    /// have outside knowledge about how big it can possibly get.
    int isize() const {
        AssertLe(size(),static_cast<size_type>(std::numeric_limits<int>::max()));
        return static_cast<int>(size());
    }

    /// MidLength is provided for compatibility with KmerPath's MidLength.
    int MidLength( ) const { return this->isize( ); }

    //
    // Iterators
    //

    const_iterator Begin() const { return cbegin(); }
    const_iterator Begin(size_type pos) const { return cbegin(pos); }
    const_iterator End() const { return cend(); }
    const_reverse_iterator RBegin() const{ return crbegin(); }
    const_reverse_iterator RBegin(size_type pos) const{ return crbegin(pos); }
    const_reverse_iterator REnd() const { return crend(); }
    const_rc_iterator RCBegin() const { return rcbegin(); }
    const_rc_iterator RCBegin(size_type pos) const { return rcbegin(pos); }
    const_rc_iterator RCEnd() const { return rcend(); }
    const_rc_iterator rcbegin() const { return const_rc_iterator(this,0); }
    const_rc_iterator rcbegin(size_type pos) const
    { AssertLe(pos,size()); return const_rc_iterator(this, pos); }
    const_rc_iterator rcend() const { return const_rc_iterator(this, size()); }

    //
    // Mutators
    //

    /// alias for swap()
    void Swap( BaseVec& b ) { swap(b); }

    /// Replace *this by its reverse (NOT reverse complement!).
    BaseVec& Reverse()
    { std::reverse(begin(), end()); return *this; }
    BaseVec& Reverse( const BaseVec& b )
    { *this = b; return Reverse(); }

    /// initialize BaseVec with a DNA String
    void SetFromString( const String& s )
    { assign(s.begin(),s.end(),CharToBaseMapper()); }

    /// Set from a string containing non-ACGTacgt characters.
    /// For any such character a random base from its ambiguity set
    /// will be placed into the basevector.
    void SetFromStringWithNs( const String &s )
    { assign(s.begin(),s.end(),GenCharToRandomBaseMapper()); }

    //
    // Deprecated wrappers of standard methods
    //

    void Set(size_type i, value_type base) { set(i, base); }
    void Setsize(size_type n, size_type extra=0 ) { reserve(n+extra).resize(n); }
    void AppendBase( value_type _base ) { push_back(_base); }
    void SetNoCheck(size_type i, value_type base) { set(i, base); }
    size_type SizeOfDynamicData() const { return physicalSize(size()); }

    /// true if we have a capacity
    bool Initialized() const { return capacity()!=0; }

    /// reduce capacity to 0
    void Reinitialize() { clear().shrink_to_fit(); }

    /// pointless unless you're spying on our internals
    void ZeroOutUnusedBases()
    { size_type sz = size(); resize(capacity()); resize(sz); }

    friend std::ostream& operator<<( std::ostream& os, const BaseVec& bv )
    { std::transform(bv.begin(),bv.end(),std::ostream_iterator<char>(os),BaseToCharMapper());
      return os; }

    // XXX: Used by many modules
    friend void BinaryWrite(int fd, const BaseVec& b);
    friend void BinaryWriteContent(int fd, const BaseVec& b);
    friend void BinaryRead(int fd, BaseVec &b);
    friend void BinaryReadContent(int fd, BaseVec &b);

    //
    // Declarations
    //

    // Utility functions

    /// Cap: in a given basevector, replace any sequence of N > n identical
    /// bases by n of the same base.
    void Cap(unsigned n);

    /// Return the first position of other inside ourselves, or size() if
    /// not found.  Start looking at start, and stop looking when we reach
    /// the minimum of end or size().
    size_type Find(const BaseVec& other, size_type start, size_type end) const;
    size_type Find(const BaseVec& other, size_type start = 0) const
    { return Find(other, start, size()); }

    /// FindAll: Find all start positions of "other" inside "this".
    vec<size_type> FindAll(const BaseVec& other) const;

    /// Return true if two basevectors overlap exactly by r bases.
    /// i.e., the last r bases of this are equal to the first r bases of that.
    bool Overlap(const BaseVec& endswith, size_type r) const;

    /// SetToSubOf(orig_bv, start_pos, len):  Set this to the length len
    /// sub-basevector of orig_bv, starting at position start_pos.  The
    /// case where this == &orig_bv is allowed.  If len is -1, it's
    /// adjusted to mean "the end" of the bvec being copied.
    BaseVec& SetToSubOf(const BaseVec& orig_bv, size_type start_pos,
                            size_type len, size_type extra = 0);

    // Genomic Classifications

    /// Return true if all bases are equal or if empty.
    bool IsHomopolymer() const;

    /// Return the % of highest base and which base it is.
    std::pair<float, unsigned char> HomopolPercent() const;

    /// Return homopolymer count at this base (extends and counts both
    /// left and right from the specified base).
    int Homopol(size_type base_position) const;

    /// Replace *this by its reverse complement.
    BaseVec& ReverseComplement();
    BaseVec& ReverseComplement(const BaseVec &bv)
    { return (*this = bv).ReverseComplement(); }

    // Return type of <Canonicalize()>, which tells you what form of the
    // input basevector was canonical.
    //
    //     rc_form - rc of input was canonical
    //     palindrome - input was palindromic
    //     fw_form - input was canonica
    enum CanonicalForm
    {
        rc_form = -1,
        palindrome = 0,
        fw_form = 1
    };

    // Return the CanonicalForm (without altering the bvec).
    //  - if this is less than its RC, return fw_form
    //  - if RC is less than this, return rc_form
    //  - if this is equal to its RC, return palindrome
    CanonicalForm getCanonicalForm() const;

    /// Turn this basevector into its canonical form.
    /// I.e., reverse-complement this bvec if the CanonicalForm is rc_form.
    CanonicalForm Canonicalize()
    { CanonicalForm result = getCanonicalForm();
      if ( result == rc_form ) ReverseComplement();
      return result; }

    size_type GcBases(size_type start, size_type end) const;
    size_type GcBases(size_type start = 0U) const
    { return GcBases(start,size()); }

    float GcPercent(size_type start, size_type end) const
    {  AssertLe(end,size()); AssertLe(start,end);
       return 100.f*GcBases(start,end) / (end-start); }
    float GcPercent(size_type start = 0U) const
    { return GcPercent(start,size()); }

    /// Make a number out of the k bases at a given offset.
    /// Oddly, the number produced will have the left-most base in the lowest
    /// two bits.  I guess k had better be <= 16.
    unsigned int extractKmer( unsigned int offset, unsigned int k ) const;

    /// Extract internal representation of the bits that represent the bases.
    /// This is meant to be opaque, but it's actually in the same weird format
    /// as described above.  Len must be >= (size()+3)/4.  Unused buffer bits
    /// will be zeroed.
    void extractBaseBits( void* buf, unsigned int len ) const
    { unsigned int physSize = physicalSize(size());
      if ( physSize )
      { AssertLe(physSize,len);
        memcpy(buf,data(),physSize);
        // zero unused bits in final byte
        static unsigned char gMasks[4] = { 0xff, 0x03, 0x0f, 0x3f };
        static_cast<unsigned char*>(buf)[physSize-1] &= gMasks[size()&3]; }
      memset(static_cast<unsigned char*>(buf)+physSize,0,len-physSize); }

    /// Assign an internal representation of the bits that represent the bases
    /// to this this bvec.  Buffer is assumed to be in the same odd order as
    /// above:  each byte has left-most base in the lowest two bits.  If you're
    /// a good citizen, you won't make use of this knowledge, but only assign
    /// base-bit representations that you've obtained via extractBaseBits.
    void assignBaseBits( unsigned int k, void const* buf )
    { resize(k); memcpy(data(),buf,(k+3)/4); }

    /// This byte-oriented method probably needs to be rethought:  it exposes
    /// to much knowledge of internals.
    unsigned int hash( unsigned int byteOffset, unsigned int nBytes ) const;

    // Output

    /// Translate basevector to String of base letters (ACGT)
    String ToString() const;

    /// Print all the bases from this basevector into the stream;
    /// breaks long sequences into 80-character lines.
    void Print(std::ostream& out) const
    { PrintBases(out, 0, size(), false); }

    /// Prints in a fasta format: ">sequence_<id>\n" followed
    /// by the full base sequence stored in the basevector; breaks
    /// long sequences nicely into 80-character lines
    void Print(std::ostream& out, int id) const
    { out << ">sequence_" << id << "\n"; Print(out); }

    void Print(std::ostream& out, String const& string_id) const
    { out << ">" << string_id << "\n"; Print(out); }

    /// Below like above prints, but with specified column breaks.
    void PrintCol(std::ostream& out, int breakCol) const
    { PrintBases(out, 0, size(), false, breakCol); }

    void PrintCol(std::ostream& out, String const& string_id, int breakCol) const
    { out << ">" << string_id << "\n"; PrintCol(out, breakCol); }

    /// Select the "nbases" bases starting at "start".
    /// If rc = False, print them.  Otherwise print the reverse complement
    /// of the bases. Breaks long sequences nicely into 80-character lines.
    void PrintBases(std::ostream& out, size_type start, size_type nbases,
                            Bool rc=False, int breakCol=80 ) const
    { AssertLe(start,size()); AssertLe(nbases,size()-start);
      size_type end = start + nbases;
      if ( rc )
          PrintBasesIter(out,rcbegin(size()-end),rcbegin(size()-start),breakCol);
      else
          PrintBasesIter(out,begin(start),begin(end),breakCol); }

    /// returns true if feudal file is good
    static bool IsGoodFeudalFile(const String& filename, bool verbose=false);

    /// Retrieves the sizes of all the bvec's in a fastb file.
    /// This peers into the fixed-length data in the feudal file, and is therefore
    /// much quicker than instantiating a vecbvec or even iterating through with
    /// a VirtualMasterVec.
    static std::vector<size_type> getSizes( String const& fastbFilename );

    static value_type Complement( value_type base ) { return base ^ 3; }
};

SELF_SERIALIZABLE(BaseVec);

namespace std
{
template<> inline void iter_swap( BaseVec::iterator itr1,
                                  BaseVec::iterator itr2 )
{
    BaseVec::value_type tmp = *itr1;
    itr1.set(*itr2);
    itr2.set(tmp);
}

template<> inline void iter_swap( BaseVec::reverse_iterator itr1,
                                  BaseVec::reverse_iterator itr2 )
{
    BaseVec::value_type tmp = *itr1;
    itr1.set(*itr2);
    itr2.set(tmp);
}
}

inline void swap( BaseVec& bv1, BaseVec& bv2 )
{
    bv1.swap(bv2);
}

// global functions

/// Copy count bases from src at src_start to target at target_start.
/// if rc_from is true, the src will be RC'd before it is copied.
void CopyBases(const BaseVec& src, BaseVec::size_type src_start,
                    BaseVec& target, BaseVec::size_type target_start,
                    BaseVec::size_type count, Bool rc_from=False);

/// copy the RC of seq into rc_seq
inline void StringReverseComplement(const String &seq, String &rc_seq)
{
    rc_seq = seq;
    GeneralizedBase::reverseComplement(rc_seq.begin(), rc_seq.end());
}

inline float GcPercent(const BaseVec& b)
{
    return b.GcPercent();
}

inline float GcPercent( const String& b,
                        BaseVec::size_type start,
                        BaseVec::size_type end )
{ return BaseVec(b).GcPercent(start,end); }

inline float GcPercent( const String& b, BaseVec::size_type start = 0 )
{ return GcPercent(b,start,b.size()); }

inline BaseVec::size_type GcBases( const BaseVec& bv,
                                   BaseVec::size_type start,
                                   BaseVec::size_type end )
{ return bv.GcBases(start, end); }

inline BaseVec::size_type GcBases( const BaseVec& bv,
                                   BaseVec::size_type start = 0U )
{ return GcBases(bv,start,bv.size()); }

inline bool Overlap( const BaseVec& source,
                     const BaseVec& endswith,
                     BaseVec::size_type r)
{ return source.Overlap(endswith, r); }

/// Algorithm: Step through basevector s.  At each location, look for an
/// overlap - i.e., a perfect match with the beginning of t, starting at that
/// location in s, and continuing to the end of s.  If we find an overlap,
/// it must be the largest overlap (because we've been searching from the
/// beginning of s.)
unsigned int LargestOverlap(const BaseVec& s, const BaseVec& t,
                            unsigned int r_max, unsigned int r_min = 0U );

inline unsigned LargestOverlap( BaseVec const& s, BaseVec const& t )
{ return LargestOverlap(s,t,s.size(),0U); }

/// Compute the concatenation of two basevectors
inline BaseVec Cat(const BaseVec& left, const BaseVec& right)
{
    BaseVec join;
    join.reserve(left.size()+right.size());
    join = left;
    join.append(right.begin(),right.end());
    return join;
}

/// Compute the concatenation of three basevectors
inline BaseVec Cat(const BaseVec& b1, const BaseVec& b2, const BaseVec& b3)
{
    BaseVec join;
    join.reserve(b1.size()+b2.size()+b3.size());
    join = b1;
    join.append(b2.begin(),b2.end());
    join.append(b3.begin(),b3.end());
    return join;
}

/// Mirrors Print() for qualvectors so we can use templates.
/// Delegates to basevector::Print, param scores_per_line is not used.
inline void Print(std::ostream &out, const BaseVec &b, const String &name,
    const int /*scores_per_line*/ = 60)
{ b.Print(out, name); }

inline String ToString(const BaseVec& bv)
{ return bv.ToString(); }

#endif
