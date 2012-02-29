/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header: StdMethods.h

   Defines macros that take a list of fields of a class, and define standard
   methods of the class (copy constructor, assignment operator, binary read and write
   operators).

   For a class with fields x, y and z, put the line

   > STD_METHODS3(x,y,z);

   into the definition of the class, to create the standard methods.
   Use the appropriate macro for the number of fields in your class.

   *IMPORTANT*: If the class has a parent, use an STD_METHODSPn instead of
   STD_METHODSn macro!  Otherwise, the copy constructor and assignment operator
   will be incorrect.
*/

#ifndef __INCLUDE_StdMethods_h
#define __INCLUDE_StdMethods_h

#include "system/System.h"
#include "Vec.h"

#define CMPLX_BIN_IO(T)                              \
  friend void BinaryWrite( int fd, const vec<T>& s ) \
  {    BinaryWriteComplex( fd, s );    }             \
  friend void BinaryRead( int fd, vec<T>& s )        \
  {    BinaryReadComplex( fd, s );    }              \
  typedef int __ ## T ## _dummy_eat_semicolon_

#define STD_METHODS1(T,f1)                           \
   T( const T& t ): f1(t.f1) { }                     \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)

#define STD_METHODS2(T,f1,f2)                        \
   T( const T& t ): f1(t.f1), f2(t.f2) { }           \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)

#define STD_METHODS3(T,f1,f2,f3)                     \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3) { } \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)

#define STD_METHODS4(T,f1,f2,f3,f4)                  \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4) { }                                  \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
    BinaryWrite( fd, t.f4 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
    BinaryRead( fd, t.f4 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)

#define STD_METHODS5(T,f1,f2,f3,f4,f5)               \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5) { }                        \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
    BinaryWrite( fd, t.f4 );                         \
    BinaryWrite( fd, t.f5 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
    BinaryRead( fd, t.f4 );                          \
    BinaryRead( fd, t.f5 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)

#define STD_METHODS6(T,f1,f2,f3,f4,f5,f6)            \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6) { }              \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
    BinaryWrite( fd, t.f4 );                         \
    BinaryWrite( fd, t.f5 );                         \
    BinaryWrite( fd, t.f6 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
    BinaryRead( fd, t.f4 );                          \
    BinaryRead( fd, t.f5 );                          \
    BinaryRead( fd, t.f6 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)

#define STD_METHODS7(T,f1,f2,f3,f4,f5,f6,f7)         \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6), f7(t.f7) { }    \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     f7 = t.f7;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
    BinaryWrite( fd, t.f4 );                         \
    BinaryWrite( fd, t.f5 );                         \
    BinaryWrite( fd, t.f6 );                         \
    BinaryWrite( fd, t.f7 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
    BinaryRead( fd, t.f4 );                          \
    BinaryRead( fd, t.f5 );                          \
    BinaryRead( fd, t.f6 );                          \
    BinaryRead( fd, t.f7 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)


#define STD_METHODS8(T,f1,f2,f3,f4,f5,f6,f7,f8)      \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6), f7(t.f7),       \
       f8(t.f8) { }                                  \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     f7 = t.f7;                                      \
     f8 = t.f8;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
    BinaryWrite( fd, t.f4 );                         \
    BinaryWrite( fd, t.f5 );                         \
    BinaryWrite( fd, t.f6 );                         \
    BinaryWrite( fd, t.f7 );                         \
    BinaryWrite( fd, t.f8 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
    BinaryRead( fd, t.f4 );                          \
    BinaryRead( fd, t.f5 );                          \
    BinaryRead( fd, t.f6 );                          \
    BinaryRead( fd, t.f7 );                          \
    BinaryRead( fd, t.f8 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)


#define STD_METHODS9(T,f1,f2,f3,f4,f5,f6,f7,f8,f9)   \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6), f7(t.f7),       \
       f8(t.f8), f9(t.f9) { }                        \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     f7 = t.f7;                                      \
     f8 = t.f8;                                      \
     f9 = t.f9;                                      \
     return *this;                                   \
   }                                                 \
  friend void BinaryWrite( int fd, const T& t ) {    \
    BinaryWrite( fd, t.f1 );                         \
    BinaryWrite( fd, t.f2 );                         \
    BinaryWrite( fd, t.f3 );                         \
    BinaryWrite( fd, t.f4 );                         \
    BinaryWrite( fd, t.f5 );                         \
    BinaryWrite( fd, t.f6 );                         \
    BinaryWrite( fd, t.f7 );                         \
    BinaryWrite( fd, t.f8 );                         \
    BinaryWrite( fd, t.f9 );                         \
  }                                                  \
  friend void BinaryRead( int fd, T& t ) {           \
    BinaryRead( fd, t.f1 );                          \
    BinaryRead( fd, t.f2 );                          \
    BinaryRead( fd, t.f3 );                          \
    BinaryRead( fd, t.f4 );                          \
    BinaryRead( fd, t.f5 );                          \
    BinaryRead( fd, t.f6 );                          \
    BinaryRead( fd, t.f7 );                          \
    BinaryRead( fd, t.f8 );                          \
    BinaryRead( fd, t.f9 );                          \
  }                                                  \
  CMPLX_BIN_IO(T)


#define DEF_BIN_IO(T)                                \
   inline void BinaryWrite( int fd, const T& x ) {   \
       WriteBytes( fd, &x, sizeof(x) );              \
   }                                                 \
   inline void BinaryRead( int fd, T& x ) {          \
       ReadBytes( fd, &x, sizeof(x) );               \
   }

DEF_BIN_IO(float)
DEF_BIN_IO(double)
DEF_BIN_IO(char)
DEF_BIN_IO(short)
DEF_BIN_IO(unsigned int)
DEF_BIN_IO(unsigned char)
DEF_BIN_IO(longlong)
DEF_BIN_IO(ulonglong)
// DEF_BIN_IO(int) is defined in System.h

#endif
// #ifndef __INCLUDE_StdMethods_h
