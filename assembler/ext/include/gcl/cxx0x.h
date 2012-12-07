// Copyright 2009 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef CXX0X_H
#define CXX0X_H

// See ../testing/cxx0x_test.cc for examples of use.

#if defined(__GNUC__) && __GNUC__ == 4

#ifndef __APPLE__
#define HAS_CXX0X_PRIMITIVE_THREAD_LOCAL
#endif

#if defined(__clang__)
# if __has_feature(cxx_static_assert)
#   define HAS_CXX0X_STATIC_ASSERT
# endif
# if __has_feature(cxx_rvalue_references)
#   define HAS_CXX0X_RVREF
# endif
# if __has_feature(cxx_variadic_templates)
#   define HAS_CXX0X_VARIADIC_TMPL
# endif
# if __has_feature(cxx_default_function_template_args)
#   define HAS_CXX0X_DFLT_FUNC_TMPL_ARGS
# endif
# if __has_feature(cxx_decltype)
#  define HAS_CXX0X_DECL_TYPE
# endif

#define HAS_CXX0X_EXTERN_TEMPLATE
#define HAS_CXX0X_EXPLICIT_AGGR_INIT
#define HAS_CXX0X___FUNC__
#define HAS_CXX0X_C99_PREPROCESSOR
#define HAS_CXX0X_LONG_LONG

# if __has_feature(cxx_defaulted_functions)
#   define HAS_CXX0X_DEFAULTED
# endif
# if __has_feature(cxx_deleted_functions)
#  define HAS_CXX0X_DELETED
# endif

# if __has_feature(cxx_auto_type)
#  define HAS_CXX0X_AUTO_VAR
# endif

# if __has_feature(cxx_nullptr)
#  define HAS_CXX0X_NULLPTR
# endif

# if __has_feature(cxx_constexpr)
#  define HAS_CXX0X_CONSTEXPR_VAR
#  define HAS_CXX0X_CONSTEXPR_FLD
#  define HAS_CXX0X_CONSTEXPR_FUN
#  define HAS_CXX0X_CONSTEXPR_CTOR
# endif

# if __has_feature(cxx_explicit_conversions)
#  define HAS_CXX0X_EXPLICIT_CONV
# endif

# if __has_feature(cxx_strong_enums)
#  define HAS_CXX0X_STRONG_ENUM
# endif

#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__)
#if __GNUC_MINOR__ >= 3
#define HAS_CXX0X_STATIC_ASSERT
#define HAS_CXX0X_RVREF
#define HAS_CXX0X_VARIADIC_TMPL
#define HAS_CXX0X_DFLT_FUNC_TMPL_ARGS
#define HAS_CXX0X_DECL_TYPE
#define HAS_CXX0X_EXTERN_TEMPLATE
#define HAS_CXX0X___FUNC__
#define HAS_CXX0X_C99_PREPROCESSOR
#define HAS_CXX0X_LONG_LONG
#endif
#if __GNUC_MINOR__ >= 4
#define HAS_CXX0X_DEFAULTED
#define HAS_CXX0X_DELETED
#define HAS_CXX0X_AUTO_VAR
#define HAS_CXX0X_NEW_FUNC_SYNTAX
#define HAS_CXX0X_EXPLICIT_AGGR_INIT
#define HAS_CXX0X_STRONG_ENUM
#define HAS_CXX0X_VARIADIC_TMPL_TMPL
#define HAS_CXX0X_INIT_LIST
#define HAS_CXX0X_EXPR_SFINAE
#define HAS_CXX0X_NEW_CHAR_TYPES
#define HAS_CXX0X_EXTENDED_SIZEOF
#define HAS_CXX0X_INLINE_NAMESPACE
#define HAS_CXX0X_ATOMIC_OPS
#define HAS_CXX0X_PROPOGATE_EXCEPT
#endif
#if __GNUC_MINOR__ >= 5
#define HAS_CXX0X_EXPLICIT_CONV
#define HAS_CXX0X_LAMBDA
#define HAS_CXX0X_UNICODE_STRING_LIT
#define HAS_CXX0X_RAW_STRING_LIT
#define HAS_CXX0X_UCN_NAME_LIT
#define HAS_CXX0X_STANDARD_LAYOUT
#define HAS_CXX0X_LOCAL_TMPL_ARGS
#endif
#if __GNUC_MINOR__ >= 6
#define HAS_CXX0X_MOVE_SPECIAL
#define HAS_CXX0X_NULLPTR
#define HAS_CXX0X_FORWARD_ENUM
#define HAS_CXX0X_CONSTEXPR_VAR
#define HAS_CXX0X_CONSTEXPR_FLD
#define HAS_CXX0X_CONSTEXPR_FUN
#define HAS_CXX0X_CONSTEXPR_CTOR
#define HAS_CXX0X_UNRESTRICT_UNION
#define HAS_CXX0X_RANGE_FOR
#define HAS_CXX0X_CORE_NOEXCEPT
#endif
#if __GNUC_MINOR__ >= 7
#define HAS_CXX0X_EXTENDED_FRIEND
#define HAS_CXX0X_EXPLICIT_OVERRIDE
#endif
/* Unimplemented:
// #define HAS_CXX0X_RVREF_THIS
// #define HAS_CXX0X_FIELD_INIT
// #define HAS_CXX0X_TMPL_ALIAS
// #define HAS_CXX0X_USER_DEFINED_LIT
// #define HAS_CXX0X_GARBAGE_COLLECTION
// #define HAS_CXX0X_SEQUENCE_POINTS
// #define HAS_CXX0X_DEPENDENCY_ORDERING
// #define HAS_CXX0X_QUICK_EXIT
// #define HAS_CXX0X_ATOMIC_IN_SIGNAL
// #define HAS_CXX0X_TRIVIAL_THREAD_LOCAL
// #define HAS_CXX0X_NON_TRIVIAL_THREAD_LOCAL
// #define HAS_CXX0X_DYNAMIC_INIT_CONCUR
// #define HAS_CXX0X_EXTENDED_INTEGRAL
// #define HAS_CXX0X_TRIVIAL_PRIVATE
*/
#endif
#endif

namespace std {

// static assert

#ifdef HAS_CXX0X_STATIC_ASSERT
#define CXX0X_STATIC_ASSERT(EXPR, VARNAME) static_assert(EXPR, # VARNAME);
#define CXX0X_CLASS_STATIC_ASSERT(EXPR, VARNAME) static_assert(EXPR, # VARNAME);
#else
template <bool t> class static_asserter;
template<> class static_asserter<true> { };
#define CXX0X_STATIC_ASSERT(EXPR, VARNAME) \
std::static_asserter<EXPR> VARNAME;
#define CXX0X_CLASS_STATIC_ASSERT(EXPR, VARNAME) \
static std::static_asserter<EXPR> VARNAME;
#endif

// defaulted functions

#ifdef HAS_CXX0X_DEFAULTED
#define CXX0X_DEFAULTED_EASY =default;
#define CXX0X_DEFAULTED_HARD( BODY ) =default;
#else
#define CXX0X_DEFAULTED_EASY { }
#define CXX0X_DEFAULTED_HARD( BODY ) BODY
#endif

// deleted functions

#ifdef HAS_CXX0X_DELETED
#define CXX0X_DELETED =delete;
#else
#define CXX0X_DELETED ;
#endif

// aggregate initialization

#ifdef HAS_CXX0X_EXPLICIT_AGGR_INIT
#define CXX0X_AGGR_INIT( DECLS ) DECLS
#define CXX0X_NO_AGGR_INIT( DECLS )
#else
#define CXX0X_AGGR_INIT( DECLS )
#define CXX0X_NO_AGGR_INIT( DECLS ) DECLS
#endif

// private members in trivial classes

#ifdef HAS_CXX0X_TRIVIAL_PRIVATE
#define CXX0X_TRIVIAL_PRIVATE private:
#else
#define CXX0X_TRIVIAL_PRIVATE
#endif

// auto variables

#ifdef HAS_CXX0X_AUTO_VAR
#define CXX0X_AUTO_VAR( VARNAME, EXPR ) auto VARNAME = EXPR
#else
#define CXX0X_AUTO_VAR( VARNAME, EXPR ) __typeof(EXPR) VARNAME = EXPR
#endif

// nullptr

#ifdef HAS_CXX0X_NULLPTR
#define CXX0X_NULLPTR nullptr
#else
const                        // this is a const object
class {
public:
  template<typename T>             // convertible to any type
    operator T*() const            // of null non-member pointer...
    { return 0; }
  template<typename C, typename T> // or any type of null
    operator T C::*() const        // member pointer...
    { return 0; }
private:
  void operator&() const;    // whose address can't be taken
} nullptr = {};              // and whose name is nullptr
#define CXX0X_NULLPTR ::std::nullptr
#endif

// constexpr variables, functions, and constructors

#ifdef HAS_CXX0X_CONSTEXPR_VAR
#define CXX0X_CONSTEXPR_VAR constexpr
#else
#define CXX0X_CONSTEXPR_VAR static const
#endif

#ifdef HAS_CXX0X_CONSTEXPR_FLD
#define CXX0X_CONSTEXPR_FLD constexpr
#else
#define CXX0X_CONSTEXPR_FLD const
#endif

#ifdef HAS_CXX0X_CONSTEXPR_FUN
#define CXX0X_CONSTEXPR_FUN constexpr
#else
#define CXX0X_CONSTEXPR_FUN inline
#endif

#ifdef HAS_CXX0X_CONSTEXPR_CTOR
#define CXX0X_CONSTEXPR_CTOR constexpr
#else
#define CXX0X_CONSTEXPR_CTOR inline
#endif

// explicit conversions

#ifdef HAS_CXX0X_EXPLICIT_CONV
#define CXX0X_EXPLICIT_CONV explicit
#else
#define CXX0X_EXPLICIT_CONV
#endif

// thread local variables

#if defined HAS_CXX0X_NON_TRIVIAL_THREAD_LOCAL
#define CXX0X_NON_TRIVIAL_THREAD_LOCAL thread_local
#define CXX0X_TRIVIAL_THREAD_LOCAL thread_local
#elif defined HAS_CXX0X_TRIVIAL_THREAD_LOCAL
#define CXX0X_TRIVIAL_THREAD_LOCAL thread_local
#elif defined HAS_CXX0X_PRIMITIVE_THREAD_LOCAL
#define HAS_CXX0X_TRIVIAL_THREAD_LOCAL
#define CXX0X_TRIVIAL_THREAD_LOCAL __thread
#endif

// rvalue reference declarations

#ifdef HAS_CXX0X_RVREF
#define CXX0X_RVREF( x ) x
#else
#define CXX0X_RVREF( x )
#endif

#ifdef HAS_CXX0X_MOVE_SPECIAL
#define CXX0X_MOVE_SPECIAL( x ) x
#else
#define CXX0X_MOVE_SPECIAL( x )
#endif

// strong enums and their qualification

#ifdef HAS_CXX0X_STRONG_ENUM
#define CXX0X_ENUM_CLASS enum class
#define CXX0X_ENUM_QUAL( e ) e::
#else
#define CXX0X_ENUM_CLASS enum
#define CXX0X_ENUM_QUAL( e )
#endif

} // namespace std

#endif
