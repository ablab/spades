/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: BinaryIO.h

   Code to help define object serialization.

   For a class T with fields f1, f2, ..., fk
   put DEFINE_BINARY_K(T,f1,...,fk) into the class definition
   to define BinaryWrite() and BinaryRead() for T and vec<T>.

   Besides being more concise, the use of macros ensures that the
   write and read routines are consistent (write/read the same fields
   in the same order) -- just have to make sure once that the macro definitions
   are correct.
*/

#ifndef __INCLUDE_BinaryIO_h
#define __INCLUDE_BinaryIO_h

#include "system/System.h"
#include "Vec.h"

#define DEFINE_BINARY_WRITE_VEC(T)                    \
  friend void BinaryWrite( int fd, const vec<T>& h )  \
  {    BinaryWriteComplex( fd, h );    }              \
  friend void BinaryRead( int fd, vec<T>& h )         \
  {    BinaryReadComplex( fd, h );    }               \
  typedef int __ ## T ## __binaryWriteEatSemicolon__

#define DEFINE_BINARY_IO_1(T,f1)               \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)


#define DEFINE_BINARY_IO_2(T,f1,f2)            \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)

#define DEFINE_BINARY_IO_3(T,f1,f2,f3)         \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)

#define DEFINE_BINARY_IO_4(T,f1,f2,f3,f4)             \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)

#define DEFINE_BINARY_IO_5(T,f1,f2,f3,f4,f5)      \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)

#define DEFINE_BINARY_IO_6(T,f1,f2,f3,f4,f5,f6)   \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)


#define DEFINE_BINARY_IO_7(T,f1,f2,f3,f4,f5,f6,f7)   \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)

#define DEFINE_BINARY_IO_8(T,f1,f2,f3,f4,f5,f6,f7,f8)   \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)


#define DEFINE_BINARY_IO_9(T,f1,f2,f3,f4,f5,f6,f7,f8,f9)   \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)


#define DEFINE_BINARY_IO_10(T,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)   \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
    BinaryWrite( fd, x.f10 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
    BinaryRead( fd, x.f10 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)


#define DEFINE_BINARY_IO_11(T,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11)	\
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
    BinaryWrite( fd, x.f10 );                          \
    BinaryWrite( fd, x.f11 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
    BinaryRead( fd, x.f10 );                           \
    BinaryRead( fd, x.f11 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)



#define DEFINE_BINARY_IO_12(T,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12)	\
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
    BinaryWrite( fd, x.f10 );                          \
    BinaryWrite( fd, x.f11 );                          \
    BinaryWrite( fd, x.f12 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
    BinaryRead( fd, x.f10 );                           \
    BinaryRead( fd, x.f11 );                           \
    BinaryRead( fd, x.f12 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)



#define DEFINE_BINARY_IO_13(T,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13)	\
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
    BinaryWrite( fd, x.f10 );                          \
    BinaryWrite( fd, x.f11 );                          \
    BinaryWrite( fd, x.f12 );                          \
    BinaryWrite( fd, x.f13 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
    BinaryRead( fd, x.f10 );                           \
    BinaryRead( fd, x.f11 );                           \
    BinaryRead( fd, x.f12 );                           \
    BinaryRead( fd, x.f13 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)



#define DEFINE_BINARY_IO_14(T,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14) \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
    BinaryWrite( fd, x.f10 );                          \
    BinaryWrite( fd, x.f11 );                          \
    BinaryWrite( fd, x.f12 );                          \
    BinaryWrite( fd, x.f13 );                          \
    BinaryWrite( fd, x.f14 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
    BinaryRead( fd, x.f10 );                           \
    BinaryRead( fd, x.f11 );                           \
    BinaryRead( fd, x.f12 );                           \
    BinaryRead( fd, x.f13 );                           \
    BinaryRead( fd, x.f14 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)

#define DEFINE_BINARY_IO_15(T,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15) \
  friend void BinaryWrite( int fd, const T& x ) {     \
    BinaryWrite( fd, x.f1 );                          \
    BinaryWrite( fd, x.f2 );                          \
    BinaryWrite( fd, x.f3 );                          \
    BinaryWrite( fd, x.f4 );                          \
    BinaryWrite( fd, x.f5 );                          \
    BinaryWrite( fd, x.f6 );                          \
    BinaryWrite( fd, x.f7 );                          \
    BinaryWrite( fd, x.f8 );                          \
    BinaryWrite( fd, x.f9 );                          \
    BinaryWrite( fd, x.f10 );                          \
    BinaryWrite( fd, x.f11 );                          \
    BinaryWrite( fd, x.f12 );                          \
    BinaryWrite( fd, x.f13 );                          \
    BinaryWrite( fd, x.f14 );                          \
    BinaryWrite( fd, x.f15 );                          \
  }                                                   \
  friend  void BinaryRead( int fd, T& x ) {           \
    BinaryRead( fd, x.f1 );                           \
    BinaryRead( fd, x.f2 );                           \
    BinaryRead( fd, x.f3 );                           \
    BinaryRead( fd, x.f4 );                           \
    BinaryRead( fd, x.f5 );                           \
    BinaryRead( fd, x.f6 );                           \
    BinaryRead( fd, x.f7 );                           \
    BinaryRead( fd, x.f8 );                           \
    BinaryRead( fd, x.f9 );                           \
    BinaryRead( fd, x.f10 );                           \
    BinaryRead( fd, x.f11 );                           \
    BinaryRead( fd, x.f12 );                           \
    BinaryRead( fd, x.f13 );                           \
    BinaryRead( fd, x.f14 );                           \
    BinaryRead( fd, x.f15 );                           \
  }                                                   \
  DEFINE_BINARY_WRITE_VEC(T)



#endif
// #ifndef __INCLUDE_system_BinaryIO_h
