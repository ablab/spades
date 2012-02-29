///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef	SYSTEM_ASSERT_H
#define SYSTEM_ASSERT_H

#include "system/Exit.h"
#include <sstream>

namespace Assert
{

void reportVals( char const* loc, char const* func, char const* vals );
void reportValsAndDie( char const* loc, char const* func, char const* vals )
            __attribute__((__noreturn__));

template <class T, class U>
void reportAndDie( T const& t, U const& u, char const* loc, char const* func )
{ std::ostringstream oss; oss << "arg1 = " << t << " and arg2 = " << u;
  reportValsAndDie(loc,func,oss.str().c_str()); }

template <class T, class U>
inline void eq( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t == u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void ne( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t != u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void gt( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t > u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void ge( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t >= u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void lt( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t < u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void le( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t <= u) ) reportAndDie(t,u,loc,func); }

inline void yes( bool t, char const* loc, char const* func )
{ if ( !t ) reportValsAndDie(loc,func,0); }

inline void no( bool t, char const* loc, char const* func )
{ if ( t ) reportValsAndDie(loc,func,0); }

template <class T, class U>
void report( T const& t, U const& u, char const* loc, char const* func )
{ std::ostringstream oss; oss << "arg1 = " << t << " and arg2 = " << u;
  reportVals(loc,func,oss.str().c_str()); }

template <class T, class U>
inline void eqTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t == u) ) report(t,u,loc,func); }

template <class T, class U>
inline void neTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t != u) ) report(t,u,loc,func); }

template <class T, class U>
inline void gtTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t > u) ) report(t,u,loc,func); }

template <class T, class U>
inline void geTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t >= u) ) report(t,u,loc,func); }

template <class T, class U>
inline void ltTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t < u) ) report(t,u,loc,func); }

template <class T, class U>
inline void leTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t <= u) ) report(t,u,loc,func); }

inline void yesTest( bool t, char const* loc, char const* func )
{ if ( !t ) reportVals(loc,func,0); }

inline void noTest( bool t, char const* loc, char const* func )
{ if ( t ) reportVals(loc,func,0); }

}

#define ASSERT_ARGS_HELP(x) #x
#define ASSERT_ARGS_HELP2(x) ASSERT_ARGS_HELP(x)
#define ASSERT_ARGS(macro,x)\
    #macro "(" #x ") at " __FILE__ ":" ASSERT_ARGS_HELP2(__LINE__),\
    __PRETTY_FUNCTION__
#define ASSERT_ARGS2(macro,x,y)\
    #macro "(" #x "," #y ") at " __FILE__ ":" ASSERT_ARGS_HELP2(__LINE__),\
    __PRETTY_FUNCTION__

// this set of macros halt execution if the condition is false
#define ForceAssertEq(x,y) Assert::eq(x,y,ASSERT_ARGS2(ForceAssertEq,x,y))
#define ForceAssertNe(x,y) Assert::ne(x,y,ASSERT_ARGS2(ForceAssertNe,x,y))
#define ForceAssertGt(x,y) Assert::gt(x,y,ASSERT_ARGS2(ForceAssertGt,x,y))
#define ForceAssertGe(x,y) Assert::ge(x,y,ASSERT_ARGS2(ForceAssertGe,x,y))
#define ForceAssertLt(x,y) Assert::lt(x,y,ASSERT_ARGS2(ForceAssertLt,x,y))
#define ForceAssertLe(x,y) Assert::le(x,y,ASSERT_ARGS2(ForceAssertLe,x,y))
#define ForceAssert(x)     Assert::yes(x,ASSERT_ARGS(ForceAssert,x))
#define ForceAssertNot(x)  Assert::no(x,ASSERT_ARGS(ForceAssertNot,x))

// this set of macros just print a message if the condition is false
#define TestAssertEq(x,y) Assert::eqTest(x,y,ASSERT_ARGS2(AssertEq,x,y))
#define TestAssertNe(x,y) Assert::neTest(x,y,ASSERT_ARGS2(AssertNe,x,y))
#define TestAssertGt(x,y) Assert::gtTest(x,y,ASSERT_ARGS2(AssertGt,x,y))
#define TestAssertGe(x,y) Assert::geTest(x,y,ASSERT_ARGS2(AssertGe,x,y))
#define TestAssertLt(x,y) Assert::ltTest(x,y,ASSERT_ARGS2(AssertLt,x,y))
#define TestAssertLe(x,y) Assert::leTest(x,y,ASSERT_ARGS2(AssertLe,x,y))
#define TestAssert(x)     Assert::yesTest(x,ASSERT_ARGS(Assert,x))
#define TestAssertNot(x)  Assert::noTest(x,ASSERT_ARGS(AssertNot,x))

// this set of macros does nothing if NDEBUG is defined
// they halt execution when the condition is false if NDEBUG isn't defined
#ifdef NDEBUG
#define AssertEq(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertNe(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertGt(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertGe(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertLt(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertLe(x,y) ((void)(sizeof(x)+sizeof(y)))
#define Assert(x)     ((void)sizeof(x))
#define AssertNot(x)  ((void)sizeof(x))
#else
#define AssertEq(x,y) Assert::eq(x,y,ASSERT_ARGS2(AssertEq,x,y))
#define AssertNe(x,y) Assert::ne(x,y,ASSERT_ARGS2(AssertNe,x,y))
#define AssertGt(x,y) Assert::gt(x,y,ASSERT_ARGS2(AssertGt,x,y))
#define AssertGe(x,y) Assert::ge(x,y,ASSERT_ARGS2(AssertGe,x,y))
#define AssertLt(x,y) Assert::lt(x,y,ASSERT_ARGS2(AssertLt,x,y))
#define AssertLe(x,y) Assert::le(x,y,ASSERT_ARGS2(AssertLe,x,y))
#define Assert(x)     Assert::yes(x,ASSERT_ARGS(Assert,x))
#define AssertNot(x)  Assert::no(x,ASSERT_ARGS(AssertNot,x))
#endif

#endif // SYSETEM_ASSERT_H
