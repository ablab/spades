///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef TYPES
#define TYPES

#include <climits>
#include <cstddef>
#include <cstdlib>
#include <netinet/in.h>

using namespace std;

// This assumes that all suns run Solaris...
#if __sun == 1
#define __solaris
#endif

#ifdef __solaris
#include <ieeefp.h>
#endif
// ===========================================================================
//
// basic arithmetic type definitions (and a routine to check them)
//
// Presumably on most machines, a "char" is one byte, a "short" is two bytes, and
// and "int" is four bytes.  However, the eight byte arithmetic type might 
// conceivably be called "long" or "long long" or both or neither, depending on the 
// architecture.
//
// ==========================================================================



#ifdef __solaris
inline bool MIN(int x, int y) {
  if (x < y)
    return x;
  else
    return y;
}
/*inline bool MIN(double x, double y) {
  if (x < y)
    return x;
  else
    return y;
    }*/
inline bool MAX(int x, int y) {
  if (x > y)
    return x;
  else
    return y;
}
/*
inline bool MAX(double x, double y) {
  if (x > y)
    return x;
  else
    return y;
}
*/


#endif 

// Typedefs required for Mac OS X, and other systems that don't 
// define a ulong
#ifndef C_ULONG_DEFINED
typedef unsigned long ulong;
#endif

// Mac OS X scandir accepts a callback taking a dirent*, but
// Linux wants a const dirent*.  Typedef the difference away.
#ifdef MACOSX
typedef struct dirent* pScandirCandidate;
#else
typedef const struct dirent* pScandirCandidate;
#endif

// The following should define an eight-byte arithmetic type.
#ifndef LONG_BIT
#error LONG_BIT not defined, cannot define a longlong type
#endif

#if (LONG_BIT == 32)
typedef long long longlong;
typedef unsigned long long ulonglong;
#define LLCONST(num) num ## LL
#define ULLCONST(num) num ## ULL
#define LLFORMAT "%Ld"
#elif (LONG_BIT == 64)
typedef long longlong;
typedef unsigned long ulonglong;
#define LLCONST(num) num ## L
#define ULLCONST(num) num ## UL
#define LLFORMAT "%ld"
#else
#error LONG_BIT has unknown value, cannot define a longlong type
#endif
// To specify compile-time constants over 32 bits (cross-platform), use
// longlong foo = LLCONST(123456789123456789);
// ulonglong ufoo = ULLCONST(12345678901234567890);


// ===========================================================================
//
// The type "bool" is poorly implemented in g++.  Static initialization 
// doesn't seem to work right.  Also I accidentally discovered that on an alpha 
// 21264, sizeof(bool) = 8!  So perhaps it would be better to eliminate the use of 
// bools.
//
// ===========================================================================

#ifndef MYBOOL
#define MYBOOL
     typedef unsigned char Bool;
     const Bool True = 1;
     const Bool False = 0;
#endif

typedef size_t size_type;

typedef unsigned char uchar;

void TestTypes( );

const unsigned int TopBit32 = 0x80000000u; // the top 1 bit of a 32-bit int
const int          Bits31   = 0x7FFFFFFF;  // the low 31 bits of a 32-bit int

const unsigned int infinitely_many = 1000000000;


// ============================================================================
//
// Define exactly one of {Little_Endian, Big_Endian}.  There seems to be a 
// different way to do this on each machine.
//
// ===========================================================================

#include <sys/types.h>
#include <sys/param.h>
#ifdef _BIG_ENDIAN
     #define Big_Endian
#endif
#ifdef _LITTLE_ENDIAN
     #define Little_Endian
#endif
#ifdef __BYTE_ORDER
     #if __BYTE_ORDER == __BIG_ENDIAN
          #define Big_Endian
     #endif
     #if __BYTE_ORDER == __LITTLE_ENDIAN
          #define Little_Endian
     #endif
#endif
#ifdef BYTE_ORDER
     #if BYTE_ORDER == BIG_ENDIAN
          #define Big_Endian
     #endif
     #if BYTE_ORDER == LITTLE_ENDIAN
          #define Little_Endian
     #endif
#endif
#ifndef NDEBUG
     #ifndef Little_Endian
          #ifndef Big_Endian
               #error Endianness of architecture unknown.  See Types.h.
          #endif
     #endif
     #ifdef Little_Endian
          #ifdef Big_Endian
               #error Endianness of architecture undetermined.  See Types.h.
          #endif
     #endif
#endif



///Transform from bigendian if needed, noop #ifdef Big_Endian
inline short int FromBigEndian(short int s) {
#ifdef Little_Endian
  return ntohs(s);
#else
  return s;
#endif
}

inline unsigned short int FromBigEndian(unsigned short int s) {
  return FromBigEndian(short(s));
}



///Transform from bigendian if needed, noop #ifdef Big_Endian
inline int FromBigEndian( int i) {
#ifdef Little_Endian
  return ntohl(i);
#else
  return i;
#endif
}

inline unsigned int FromBigEndian(unsigned int i) {
  return FromBigEndian(int(i));
}



// The below will give a warning in G++ 4.3+ with the -Wconversion flag.
// This is because ntohs is a macro that has another macro in it, __bswap16.
// bswap16 has a bitwise & operator whose result g++ casts as an int.  So,
// there is an implicit type casting from int to unsigned short that occurs
// therein.  
///Transform from bigendian if needed, noop #ifdef Big_Endian
inline longlong FromBigEndian( longlong ll) {
#ifdef Little_Endian
  longlong ret = ntohl(ll & 0xFFFFFFFF);
  ret |= ntohl(ll >> 32 & 0xFFFFFFFF);
  return ret;
#else
  return ll;
#endif
}

inline uint64_t  FromBigEndian(uint64_t ll) {
  return FromBigEndian(longlong(ll));
}

typedef unsigned char unsigned_char_t;
typedef unsigned short unsigned_short_t;
typedef unsigned int unsigned_int_t;

/**
   Macro: FOR_ALL_BUILTIN_TYPES

   Call the given macro for all built-in C++ types.  Useful e.g. for creating
   specialized versions of a function template for working with built-in types.
*/
#define FOR_ALL_BUILTIN_TYPES(M) \
    M(char); \
    M(unsigned_char_t); \
    M(short); \
    M(unsigned_short_t); \
    M(int); \
    M(unsigned_int_t); \
    M(longlong); \
    M(ulonglong)

#endif
