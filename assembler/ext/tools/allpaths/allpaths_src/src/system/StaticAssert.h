///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef STATIC_ASSERT_H
#define STATIC_ASSERT_H

///Compile-time assert helper class
template<bool> struct CompileTimeError;
///Compile-time assert helper class partial specialization
template<> struct CompileTimeError<true> {};
///Compile-time assert macro with no message.
///From Alexandrescu's book Modern C++ Design
#define STATIC_ASSERT(expr) (CompileTimeError<(expr) != 0>())

///Helper class for compile-time assert with message
template<bool> struct CompileTimeChecker;
///Helper class for compile-time assert with message partial specialization.
///This one does not have the generic constructor and will therefore emit
///a compiler error message.
template<> struct CompileTimeChecker<false>  {};
///Helper class for compile-time assert with message partial specialization.
template<> struct CompileTimeChecker<true>  { CompileTimeChecker(...) {} };
///Compile-time assert macro with message (must be valid C++ identifier).
///Adapted from Alexandrescu's book Modern C++ Design
#if __sun == 1
#define STATIC_ASSERT_M(expr,msg) \
{ \
  class ERROR_##msg {}; \
  CompileTimeChecker<(expr) != 0>  foo; \
}
#else //__sun
#define STATIC_ASSERT_M(expr,msg) \
{ \
  class ERROR_##msg {}; \
  CompileTimeChecker<(expr) != 0>  foo((ERROR_##msg())); \
}
#endif
#endif //STATIC_ASSERT_H
