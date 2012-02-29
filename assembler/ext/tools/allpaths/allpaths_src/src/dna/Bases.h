///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Bases.h
 * \author tsharpe
 * \date Apr 10, 2009
 *
 * \brief Classes representing nucleotides.
 */
#ifndef DNA_BASES_H_
#define DNA_BASES_H_

#include <cstring>
#include <functional>

/// A (possibly ambiguous) base.
class GeneralizedBase
{
public:
    /// The canonical IUB upper-case letter for the ambiguous nucleotide.
    /// The "no base" case's value is '-'.
    char asChar() const
    { return mCharRep; }

    /// Bit mask of nucleotides represented:  A=bit0, C=bit1, G=bit2, T=bit3.
    unsigned char bits() const
    { return mBits; }

    /// Number of choices of nucleotides represented by this one:
    /// 0 for X, 1 for ACGT, 2 for KMRSWY, 3 for BDHV, 4 for N
    size_t getAmbiguityCount() const { return mEndBases-mBases; }

    /// False for ACGT, true otherwise (including the "no base" case).
    bool isAmbiguous() const
    { return mAmbiguous; }

    /// Returns true if the argument has any bases in common with this one.
    bool matches( GeneralizedBase const& base ) const
    { return mBits & base.mBits; }

    /// Returns the complementary generalized base.
    /// Complement of N is N.  Complement of - is -.
    GeneralizedBase const& complement() const
    { return fromChar(mComplement); }

    /// Short-cut for complement().asChar().
    char complementChar() const
    { return mComplement; }

    /// A string giving the bases represented by the ambiguity code, e.g. GeneralizedBase::W.bases() is "AT".
    char const* bases() const
    { return mBases; }

    /// End of above string.
    char const* end() const
    { return mEndBases; }

    /// A random base from within the set of bases represented by the ambiguity code.
    unsigned char random() const;

    /// Singletons, so equality is tested with an address comparison.
    bool operator==( GeneralizedBase const& base ) const
    { return &base == this; }

    /// Singletons, so inequality is tested with an address comparison.
    bool operator!=( GeneralizedBase const& base ) const
    { return &base != this; }

    /// Accepts upper and lower case letters from the set [ACGTRYKMSWBDHVNX.-].
    /// X and - are "no base".  N and . are "any base".
    static GeneralizedBase const& fromChar( char chr )
    { GeneralizedBase const* pBase = gCharToGenBase[static_cast<unsigned char>(chr)];
      if ( !pBase ) { fromCharFatalErr(chr); }
      return *pBase; }

    /// Returns the ambiguity code for chr1 or chr2.
    static char ambiguityCode( char chr1, char chr2 )
    { unsigned char bits = GeneralizedBase::fromChar(chr1).bits() |
                            GeneralizedBase::fromChar(chr2).bits();
      return GeneralizedBase::fromBits(bits).asChar(); }

    /// Accepts values from 0 to 15.
    static GeneralizedBase const& fromBits( unsigned char bits )
    { if ( bits > 15 ) { fromBitsFatalErr(bits); }
      return *gBitsToGenBase[bits]; }

    /// Upper and lower case letters from the set [ACGTRYKMSWBDHVNX.-] return true.
    static bool isGeneralizedBase( char chr )
    { return gCharToGenBase[static_cast<unsigned char>(chr)]; }

    /// Check that all the characters in a sequence are generalized bases.
    /// The Itr type must be an input iterator resolving to a char.
    template <class Itr>
    static bool areGeneralizedBases( Itr begin, Itr const& end )
    { bool result = true;
      for ( ; result && begin != end; ++begin ) result = isGeneralizedBase(*begin);
      return result; }

    static char complementChar( char chr )
    { char result = gRevComp[static_cast<unsigned char>(chr)];
      if ( !result ) fromCharFatalErr(chr);
      return result; }

    /// Reverse-complement a sequence of characters in place, preserving case.
    /// The Itr type must be a bidirectional iterator resolving to a char.
    template <class Itr>
    static void reverseComplement( Itr begin, Itr end )
    { while ( begin != end )
      { char tmp = complementChar(*--end);
        if ( begin == end ) { *begin = tmp; break; }
        *end = complementChar(*begin);
        *begin = tmp;
        ++begin; }
    }

    /// This returns a Base val by choosing an element from the set of bases
    /// represented by chr, which may be ambiguous.  Contrast this with
    /// Base::char2Val, for which the character must be unambiguous.
    static unsigned char char2Val( char chr )
    { return fromChar(chr).random(); }

    /// Returns a 4-bit, generalized base encoding for a character.
    static unsigned char char2Bits( char chr )
    { return fromChar(chr).bits(); }

    /// Inverse of above.
    static char bits2Char( unsigned char bits )
    { return fromBits(bits).asChar(); }

    /// Maps bits (which may be ambiguous) to a random base from the set of bits.
    static unsigned char bits2Val( unsigned char bits )
    { return fromBits(bits).random(); }

    /// Are the bits ambiguous?
    static bool bits2Ambig( unsigned char bits )
    { return fromBits(bits).isAmbiguous(); }

    static GeneralizedBase const R;
    static GeneralizedBase const Y;
    static GeneralizedBase const K;
    static GeneralizedBase const M;
    static GeneralizedBase const S;
    static GeneralizedBase const W;
    static GeneralizedBase const B;
    static GeneralizedBase const V;
    static GeneralizedBase const D;
    static GeneralizedBase const H;
    static GeneralizedBase const N;
    static GeneralizedBase const X; // i.e., no base, or not a base

private:
    GeneralizedBase( GeneralizedBase const& ); // unimplemented -- these are singletons
    GeneralizedBase( char charRep, char complement, bool isAmbiguous, unsigned int bits, char const* bases )
    : mCharRep(charRep), mComplement(complement), mAmbiguous(isAmbiguous), mBits(bits), mBases(bases), mEndBases(bases+strlen(bases))
    {}

    unsigned char chooseRandom() const;

    static void fromBitsFatalErr( unsigned char bits );
    static void fromCharFatalErr( char chr );

    char mCharRep;
    char mComplement;
    bool mAmbiguous;
    unsigned int mBits; // 4-bit representation for figuring out matches
    char const* mBases;
    char const* mEndBases;

    static GeneralizedBase const*const gCharToGenBase[256];
    static GeneralizedBase const*const gBitsToGenBase[16];
    static char gRevComp[256];
    friend class Base;
};

/// An unambiguous base.
class Base : public GeneralizedBase
{
public:
    /// Returns the code for the base:  A=0, C=1, G=2, T=3.
    unsigned char val() const
    { return mVal; }

    /// Returns the complement.
    Base const& complement() const
    { return *static_cast<Base const*>(&GeneralizedBase::complement()); }

    /// Accepts values from 0 to 3.
    static Base const& fromVal( unsigned char val )
    { if ( val & ~3 ) { fromValFatalErr(val); }
      return *gValToBase[val]; }

    /// Accepts upper and lower case letters from the set [ACGT].
    static Base const& fromChar( char chr )
    { Base const* pBase = gCharToBase[static_cast<unsigned char>(chr)];
      if ( !pBase ) { fromCharFatalErr(chr); }
      return *pBase; }

    /// True for upper- and lower-case letters from the set [ACGT].
    static bool isBase( char chr )
    { return gCharToBase[static_cast<unsigned char>(chr)]; }

    /// True for the upper-case letters [ACGT].
    static bool isCanonicalBase( char chr )
    { Base const* pBase = gCharToBase[static_cast<unsigned char>(chr)];
      return pBase && pBase->asChar()==chr; }

    /// Check that all the characters in a sequence are bases.
    /// The Itr type must be an input iterator resolving to a char.
    template <class Itr>
    static bool areBases( Itr begin, Itr const& end )
    { bool result = true;
      for ( ; result && begin != end; ++begin ) result = isBase(*begin);
      return result; }

    /// turns a character into its 2-bit encoded value
    static unsigned char char2Val( char chr )
    { return fromChar(chr).val(); }

    /// turns a 2-bit encoded value into a character
    static char val2Char( unsigned char val )
    { return fromVal(val).asChar(); }

    /// turns a 2-bit encoded value into a 4-bit mask
    static unsigned char val2Bits( unsigned char val )
    { return fromVal(val).bits(); }

    /// reverse complements a byte representing 4 bases, 2-bit encoded
    static unsigned char rcByte( unsigned char val )
    { return gRCByte[val]; }

    /// reverse complements some integer type representing 2-bit encoded bases
    template <class IntType>
    static IntType rc( IntType val )
    { IntType result;
      unsigned char* dest = reinterpret_cast<unsigned char*>(&result);
      unsigned char* end = reinterpret_cast<unsigned char*>(&val);
      unsigned char* src = end + sizeof(val);
      while ( src != end ) *dest++ = rcByte(*--src);
      return result; }

    static Base const A;
    static Base const C;
    static Base const G;
    static Base const T;

    static void fromValFatalErr( unsigned char val );

private:
    Base( char charRep, char complement, bool isAmbiguous, unsigned int bits, char const* bases, unsigned char val )
    : GeneralizedBase(charRep,complement,isAmbiguous,bits,bases), mVal(val)
    {}

    static void fromCharFatalErr( char chr );

    unsigned char mVal;
    static Base const*const gValToBase[4];
    static Base const*const gCharToBase[256];
    static unsigned char const gRCByte[256];
};

typedef unsigned char base_t;

/// BASE_A == Base::A.val()
const unsigned char BASE_A = 0;
/// BASE_C == Base::C.val()
const unsigned char BASE_C = 1;
/// BASE_G == Base::G.val()
const unsigned char BASE_G = 2;
/// BASE_T == Base::T.val()
const unsigned char BASE_T = 3;

inline unsigned char GeneralizedBase::random() const
{ return mEndBases == mBases+1 ? Base::char2Val(*mBases) : chooseRandom(); }

/// Converts the code for a base to the one-character, canonical name of the base.
inline char as_base( unsigned char x )
{ return Base::val2Char(x); }

/// Tests whether a given base (specified as integer code) is a G or a C.
inline bool IsGC( unsigned char x )
{ return x == BASE_C || x == BASE_G; }

/// Convert the one-character name (base letter) of a base to its code.
/// Inverse of as_base().
inline unsigned char as_char( char c )
{ return Base::char2Val(c); }

/// Return the code for the base complementary by Watson-Crick base pairing to the given code.
inline unsigned char GetComplementaryBase( unsigned char base )
{ if ( base > 3 ) Base::fromValFatalErr(base); // only 0,1,2,3 (A,C,T,G) are allowed
  return 3 - base; }

/// Returns the canonical letter for the base complementary to this one, including ambiguous bases.
inline char GetComplementaryBaseChar( char c )
{ return GeneralizedBase::complementChar(c); }

struct BaseToCharMapper
: public std::unary_function<unsigned char,char>
{
    char operator()( unsigned char val ) const
    { return Base::val2Char(val); }
};

struct CharToBaseMapper
: public std::unary_function<char,unsigned char>
{
    unsigned char operator()( char c ) const
    { return Base::char2Val(c); }
};

struct GenCharToRandomBaseMapper
: public std::unary_function<char,unsigned char>
{
    unsigned char operator()( char c ) const
    { return GeneralizedBase::char2Val(c); }
};

struct GenBaseBitsToBaseValMapper
: public std::unary_function<unsigned char,unsigned char>
{
    unsigned char operator()( unsigned char bits ) const
    { return GeneralizedBase::bits2Val(bits); }
};

#endif // DNA_BASES_H
