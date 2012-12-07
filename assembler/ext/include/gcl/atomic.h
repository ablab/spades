
#ifndef CXX0X_ATOMIC_H
#define CXX0X_ATOMIC_H


/*
This header implements as much as is currently feasible of
the C++0x atomics,
which is chapter 29 of
http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2010/n3092.pdf,
and the C1x atomics.
*/

#ifdef __cplusplus
#include <cstddef>
#include "cxx0x.h"
namespace std {
#else
#include <stddef.h>
#include <stdbool.h>
#endif


typedef enum memory_order {
memory_order_relaxed, memory_order_consume, memory_order_acquire,
memory_order_release, memory_order_acq_rel, memory_order_seq_cst
} memory_order;


#ifdef __cplusplus
extern "C" {
#endif

extern void atomic_thread_fence( memory_order );
extern void atomic_signal_fence( memory_order );

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus

template <class T> T kill_dependency(T y) { return y; }

#endif


typedef struct atomic_flag atomic_flag;


#ifdef __cplusplus
extern "C" {
#endif

extern bool __atomic_flag_test_and_set_explicit
( volatile atomic_flag*, memory_order );
extern void __atomic_flag_clear_explicit
( volatile atomic_flag*, memory_order );
extern void __atomic_flag_wait_explicit__
( volatile atomic_flag*, memory_order );
extern volatile atomic_flag* __atomic_flag_for_address__
( const volatile void* __z__ )
__attribute__((const));

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus


inline bool atomic_flag_test_and_set
(  atomic_flag* __a__ )
{ return __atomic_flag_test_and_set_explicit( __a__, memory_order_seq_cst ); }

inline bool atomic_flag_test_and_set_explicit
(  atomic_flag* __a__, memory_order __x__ )
{ return __atomic_flag_test_and_set_explicit( __a__, __x__ ); }

inline void atomic_flag_clear
(  atomic_flag* __a__ )
{ return __atomic_flag_clear_explicit( __a__, memory_order_seq_cst ); }

inline void atomic_flag_clear_explicit
(  atomic_flag* __a__, memory_order __x__ )
{ return __atomic_flag_clear_explicit( __a__, __x__ ); }


inline bool atomic_flag_test_and_set
( volatile atomic_flag* __a__ )
{ return __atomic_flag_test_and_set_explicit( __a__, memory_order_seq_cst ); }

inline bool atomic_flag_test_and_set_explicit
( volatile atomic_flag* __a__, memory_order __x__ )
{ return __atomic_flag_test_and_set_explicit( __a__, __x__ ); }

inline void atomic_flag_clear
( volatile atomic_flag* __a__ )
{ return __atomic_flag_clear_explicit( __a__, memory_order_seq_cst ); }

inline void atomic_flag_clear_explicit
( volatile atomic_flag* __a__, memory_order __x__ )
{ return __atomic_flag_clear_explicit( __a__, __x__ ); }


#endif


#ifndef __cplusplus

#define atomic_flag_test_and_set( __a__ ) __atomic_flag_test_and_set_explicit( __a__, memory_order_seq_cst )

#define atomic_flag_test_and_set_explicit( __a__, __x__ ) __atomic_flag_test_and_set_explicit( __a__, __x__ )

#define atomic_flag_clear( __a__ ) __atomic_flag_clear_explicit( __a__, memory_order_seq_cst )

#define atomic_flag_clear_explicit( __a__, __x__ ) __atomic_flag_clear_explicit( __a__, __x__ )

#endif


struct atomic_flag
{
    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_flag() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_flag( bool __v__ ) : __f__( __v__ ) { }
    atomic_flag( const atomic_flag& ) CXX0X_DELETED
    atomic_flag& operator =( const atomic_flag& ) CXX0X_DELETED
    )


    bool test_and_set( memory_order __x__ = memory_order_seq_cst ) 
    { return atomic_flag_test_and_set_explicit( this, __x__ ); }

    void clear( memory_order __x__ = memory_order_seq_cst ) 
    { atomic_flag_clear_explicit( this, __x__ ); }


    bool test_and_set( memory_order __x__ = memory_order_seq_cst ) volatile
    { return atomic_flag_test_and_set_explicit( this, __x__ ); }

    void clear( memory_order __x__ = memory_order_seq_cst ) volatile
    { atomic_flag_clear_explicit( this, __x__ ); }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    bool __f__;
};

#define ATOMIC_FLAG_INIT { false }


#define ATOMIC_CHAR_LOCK_FREE 0
#define ATOMIC_CHAR16_T_LOCK_FREE 0
#define ATOMIC_CHAR32_T_LOCK_FREE 0
#define ATOMIC_WCHAR_T_LOCK_FREE 0
#define ATOMIC_SHORT_LOCK_FREE 0
#define ATOMIC_INT_LOCK_FREE 0
#define ATOMIC_LONG_LOCK_FREE 0
#define ATOMIC_LLONG_LOCK_FREE 0
#define ATOMIC_ADDRESS_LOCK_FREE 0


#define ATOMIC_VAR_INIT( __m__ ) { __m__ }


#define _ATOMIC_LOAD_( __a__, __x__ ) ({ volatile __typeof__((__a__)->__f__)* __p__ = &((__a__)->__f__);    volatile atomic_flag* __g__ = __atomic_flag_for_address__( __p__ );    __atomic_flag_wait_explicit__( __g__, __x__ );    __typeof__((__a__)->__f__) __r__ = *__p__;    atomic_flag_clear_explicit( __g__, __x__ );    __r__; })

#define _ATOMIC_STORE_( __a__, __m__, __x__ ) ({ volatile __typeof__((__a__)->__f__)* __p__ = &((__a__)->__f__);    __typeof__(__m__) __v__ = (__m__);    volatile atomic_flag* __g__ = __atomic_flag_for_address__( __p__ );    __atomic_flag_wait_explicit__( __g__, __x__ );    *__p__ = __v__;    atomic_flag_clear_explicit( __g__, __x__ );    __v__; })

#define _ATOMIC_MODIFY_( __a__, __o__, __m__, __x__ ) ({ volatile __typeof__((__a__)->__f__)* __p__ = &((__a__)->__f__);    __typeof__(__m__) __v__ = (__m__);    volatile atomic_flag* __g__ = __atomic_flag_for_address__( __p__ );    __atomic_flag_wait_explicit__( __g__, __x__ );    __typeof__((__a__)->__f__) __r__ = *__p__;    *__p__ __o__ __v__;    atomic_flag_clear_explicit( __g__, __x__ );    __r__; })

#define _ATOMIC_PTRMOD_( __a__, __o__, __m__, __x__ ) ({ volatile __typeof__((__a__)->__f__)* __p__ = &((__a__)->__f__);    __typeof__(__m__) __v__ = (__m__);    volatile atomic_flag* __g__ = __atomic_flag_for_address__( __p__ );    __atomic_flag_wait_explicit__( __g__, __x__ );    __typeof__((__a__)->__f__) __r__ = *__p__;    *(char*)__p__ __o__ __v__;    atomic_flag_clear_explicit( __g__, __x__ );    __r__; })

#define _ATOMIC_CMPSWP_( __a__, __e__, __m__, __x__ ) ({ volatile __typeof__((__a__)->__f__)* __p__ = &((__a__)->__f__);    __typeof__(__e__) __q__ = (__e__);    __typeof__(__m__) __v__ = (__m__);    bool __r__;    volatile atomic_flag* __g__ = __atomic_flag_for_address__( __p__ );    __atomic_flag_wait_explicit__( __g__, __x__ );    __typeof__((__a__)->__f__) __t__ = *__p__;    if ( __t__ == *__q__ ) { *__p__ = __v__; __r__ = true; }    else { *__q__ = (__typeof__(*__q__))__t__; __r__ = false; }    atomic_flag_clear_explicit( __g__, __x__ );    __r__; })


typedef struct atomic_bool
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_bool() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_bool( bool __v__ )
    : __f__( __v__) { }
    atomic_bool( const atomic_bool& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( bool __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    bool load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (bool)_ATOMIC_LOAD_( this, __x__ ); }

    bool exchange
    ( bool __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (bool)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( bool& __e__, bool __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( bool& __e__, bool __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( bool& __e__, bool __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( bool& __e__, bool __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( bool __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    bool load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (bool)_ATOMIC_LOAD_( this, __x__ ); }

    bool exchange
    ( bool __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (bool)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( bool& __e__, bool __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( bool& __e__, bool __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( bool& __e__, bool __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( bool& __e__, bool __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (bool)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_bool& operator =
    ( const atomic_bool& )  CXX0X_DELETED
    )

    bool operator =( bool __v__ ) 
    { store( __v__ ); return __v__; }

    operator bool() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_bool& operator =
    ( const atomic_bool& ) volatile CXX0X_DELETED
    )

    bool operator =( bool __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator bool() const volatile
    { return load(); }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    bool __f__;
} atomic_bool;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_bool* __a__ )
{ return false; }

inline void atomic_init
(  atomic_bool* __a__, bool __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline bool atomic_load_explicit
(  atomic_bool* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline bool atomic_load(  atomic_bool* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_bool* __a__, bool __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_bool* __a__, bool __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline bool atomic_exchange_explicit
(  atomic_bool* __a__, bool __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline bool atomic_exchange
(  atomic_bool* __a__, bool __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_bool* __a__, bool* __e__, bool __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_bool* __a__, bool* __e__, bool __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_bool* __a__, bool* __e__, bool __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_bool* __a__, bool* __e__, bool __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_bool* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_bool* __a__, bool __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline bool atomic_load_explicit
( volatile atomic_bool* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline bool atomic_load( volatile atomic_bool* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_bool* __a__, bool __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_bool* __a__, bool __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline bool atomic_exchange_explicit
( volatile atomic_bool* __a__, bool __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline bool atomic_exchange
( volatile atomic_bool* __a__, bool __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_bool* __a__, bool* __e__, bool __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_bool* __a__, bool* __e__, bool __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_bool* __a__, bool* __e__, bool __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_bool* __a__, bool* __e__, bool __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


typedef struct atomic_address
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_address() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_address( void* __v__ )
    : __f__( __v__) { }
    atomic_address( const atomic_address& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( void* __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    void* load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (void*)_ATOMIC_LOAD_( this, __x__ ); }

    void* exchange
    ( void* __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (void*)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( void*& __e__, void* __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( void*& __e__, void* __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( void*& __e__, void* __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( void*& __e__, void* __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( void* __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    void* load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (void*)_ATOMIC_LOAD_( this, __x__ ); }

    void* exchange
    ( void* __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (void*)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( void*& __e__, void* __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( void*& __e__, void* __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( void*& __e__, void* __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( void*& __e__, void* __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (void*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_address& operator =
    ( const atomic_address& )  CXX0X_DELETED
    )

    void* operator =( void* __v__ ) 
    { store( __v__ ); return __v__; }

    operator void*() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_address& operator =
    ( const atomic_address& ) volatile CXX0X_DELETED
    )

    void* operator =( void* __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator void*() const volatile
    { return load(); }


    void* fetch_add
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_PTRMOD_( this, +=, __m__, __x__ ); }


    void* fetch_add
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_PTRMOD_( this, +=, __m__, __x__ ); }


    void* fetch_sub
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_PTRMOD_( this, -=, __m__, __x__ ); }


    void* fetch_sub
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_PTRMOD_( this, -=, __m__, __x__ ); }


    void* operator +=( ptrdiff_t __v__ ) 
    { return ((char*)fetch_add( __v__ )) + __v__; }


    void* operator +=( ptrdiff_t __v__ ) volatile
    { return ((char*)fetch_add( __v__ )) + __v__; }


    void* operator -=( ptrdiff_t __v__ ) 
    { return ((char*)fetch_sub( __v__ )) - __v__; }


    void* operator -=( ptrdiff_t __v__ ) volatile
    { return ((char*)fetch_sub( __v__ )) - __v__; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    void* __f__;
} atomic_address;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_address* __a__ )
{ return false; }

inline void atomic_init
(  atomic_address* __a__, void* __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline void* atomic_load_explicit
(  atomic_address* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline void* atomic_load(  atomic_address* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_address* __a__, void* __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_address* __a__, void* __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline void* atomic_exchange_explicit
(  atomic_address* __a__, void* __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline void* atomic_exchange
(  atomic_address* __a__, void* __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_address* __a__, void** __e__, void* __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_address* __a__, void** __e__, void* __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_address* __a__, void** __e__, void* __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_address* __a__, void** __e__, void* __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_address* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_address* __a__, void* __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline void* atomic_load_explicit
( volatile atomic_address* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline void* atomic_load( volatile atomic_address* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_address* __a__, void* __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_address* __a__, void* __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline void* atomic_exchange_explicit
( volatile atomic_address* __a__, void* __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline void* atomic_exchange
( volatile atomic_address* __a__, void* __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_address* __a__, void** __e__, void* __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_address* __a__, void** __e__, void* __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_address* __a__, void** __e__, void* __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_address* __a__, void** __e__, void* __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline void* atomic_fetch_add_explicit
(  atomic_address* __a__, ptrdiff_t __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline void* atomic_fetch_add
(  atomic_address* __a__, ptrdiff_t __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline void* atomic_fetch_add_explicit
( volatile atomic_address* __a__, ptrdiff_t __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline void* atomic_fetch_add
( volatile atomic_address* __a__, ptrdiff_t __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline void* atomic_fetch_sub_explicit
(  atomic_address* __a__, ptrdiff_t __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline void* atomic_fetch_sub
(  atomic_address* __a__, ptrdiff_t __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline void* atomic_fetch_sub_explicit
( volatile atomic_address* __a__, ptrdiff_t __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline void* atomic_fetch_sub
( volatile atomic_address* __a__, ptrdiff_t __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_char
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_char() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_char( char __v__ )
    : __f__( __v__) { }
    atomic_char( const atomic_char& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    char load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (char)_ATOMIC_LOAD_( this, __x__ ); }

    char exchange
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (char)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( char& __e__, char __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( char& __e__, char __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( char& __e__, char __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( char& __e__, char __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    char load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (char)_ATOMIC_LOAD_( this, __x__ ); }

    char exchange
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (char)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( char& __e__, char __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( char& __e__, char __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( char& __e__, char __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( char& __e__, char __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_char& operator =
    ( const atomic_char& )  CXX0X_DELETED
    )

    char operator =( char __v__ ) 
    { store( __v__ ); return __v__; }

    operator char() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_char& operator =
    ( const atomic_char& ) volatile CXX0X_DELETED
    )

    char operator =( char __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator char() const volatile
    { return load(); }


    char fetch_add
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    char fetch_add
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    char fetch_sub
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    char fetch_sub
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    char fetch_and
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    char fetch_and
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    char fetch_or
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    char fetch_or
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    char fetch_xor
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    char fetch_xor
    ( char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    char operator +=( char __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    char operator +=( char __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    char operator -=( char __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    char operator -=( char __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    char operator &=( char __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    char operator &=( char __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    char operator |=( char __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    char operator |=( char __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    char operator ^=( char __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    char operator ^=( char __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    char operator ++( int ) 
    { return fetch_add( 1 ); }
    char operator --( int ) 
    { return fetch_sub( 1 ); }
    char operator ++() 
    { return fetch_add( 1 ) + 1; }
    char operator --() 
    { return fetch_sub( 1 ) - 1; }


    char operator ++( int ) volatile
    { return fetch_add( 1 ); }
    char operator --( int ) volatile
    { return fetch_sub( 1 ); }
    char operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    char operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    char __f__;
} atomic_char;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_char* __a__ )
{ return false; }

inline void atomic_init
(  atomic_char* __a__, char __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline char atomic_load_explicit
(  atomic_char* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline char atomic_load(  atomic_char* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_char* __a__, char __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline char atomic_exchange_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline char atomic_exchange
(  atomic_char* __a__, char __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_char* __a__, char* __e__, char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_char* __a__, char* __e__, char __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_char* __a__, char* __e__, char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_char* __a__, char* __e__, char __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_char* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_char* __a__, char __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline char atomic_load_explicit
( volatile atomic_char* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline char atomic_load( volatile atomic_char* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_char* __a__, char __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline char atomic_exchange_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline char atomic_exchange
( volatile atomic_char* __a__, char __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_char* __a__, char* __e__, char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_char* __a__, char* __e__, char __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_char* __a__, char* __e__, char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_char* __a__, char* __e__, char __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline char atomic_fetch_add_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline char atomic_fetch_add
(  atomic_char* __a__, char __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_add_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline char atomic_fetch_add
( volatile atomic_char* __a__, char __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_sub_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline char atomic_fetch_sub
(  atomic_char* __a__, char __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_sub_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline char atomic_fetch_sub
( volatile atomic_char* __a__, char __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_and_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline char atomic_fetch_and
(  atomic_char* __a__, char __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_and_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline char atomic_fetch_and
( volatile atomic_char* __a__, char __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_or_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline char atomic_fetch_or
(  atomic_char* __a__, char __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_or_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline char atomic_fetch_or
( volatile atomic_char* __a__, char __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_xor_explicit
(  atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline char atomic_fetch_xor
(  atomic_char* __a__, char __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline char atomic_fetch_xor_explicit
( volatile atomic_char* __a__, char __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline char atomic_fetch_xor
( volatile atomic_char* __a__, char __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_schar
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_schar() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_schar( signed char __v__ )
    : __f__( __v__) { }
    atomic_schar( const atomic_schar& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    signed char load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (signed char)_ATOMIC_LOAD_( this, __x__ ); }

    signed char exchange
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (signed char)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( signed char& __e__, signed char __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( signed char& __e__, signed char __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( signed char& __e__, signed char __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( signed char& __e__, signed char __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    signed char load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (signed char)_ATOMIC_LOAD_( this, __x__ ); }

    signed char exchange
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (signed char)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( signed char& __e__, signed char __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( signed char& __e__, signed char __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( signed char& __e__, signed char __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( signed char& __e__, signed char __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (signed char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_schar& operator =
    ( const atomic_schar& )  CXX0X_DELETED
    )

    signed char operator =( signed char __v__ ) 
    { store( __v__ ); return __v__; }

    operator signed char() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_schar& operator =
    ( const atomic_schar& ) volatile CXX0X_DELETED
    )

    signed char operator =( signed char __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator signed char() const volatile
    { return load(); }


    signed char fetch_add
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    signed char fetch_add
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    signed char fetch_sub
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    signed char fetch_sub
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    signed char fetch_and
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    signed char fetch_and
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    signed char fetch_or
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    signed char fetch_or
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    signed char fetch_xor
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    signed char fetch_xor
    ( signed char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    signed char operator +=( signed char __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    signed char operator +=( signed char __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    signed char operator -=( signed char __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    signed char operator -=( signed char __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    signed char operator &=( signed char __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    signed char operator &=( signed char __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    signed char operator |=( signed char __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    signed char operator |=( signed char __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    signed char operator ^=( signed char __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    signed char operator ^=( signed char __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    signed char operator ++( int ) 
    { return fetch_add( 1 ); }
    signed char operator --( int ) 
    { return fetch_sub( 1 ); }
    signed char operator ++() 
    { return fetch_add( 1 ) + 1; }
    signed char operator --() 
    { return fetch_sub( 1 ) - 1; }


    signed char operator ++( int ) volatile
    { return fetch_add( 1 ); }
    signed char operator --( int ) volatile
    { return fetch_sub( 1 ); }
    signed char operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    signed char operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    signed char __f__;
} atomic_schar;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_schar* __a__ )
{ return false; }

inline void atomic_init
(  atomic_schar* __a__, signed char __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline signed char atomic_load_explicit
(  atomic_schar* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline signed char atomic_load(  atomic_schar* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_schar* __a__, signed char __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline signed char atomic_exchange_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline signed char atomic_exchange
(  atomic_schar* __a__, signed char __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_schar* __a__, signed char* __e__, signed char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_schar* __a__, signed char* __e__, signed char __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_schar* __a__, signed char* __e__, signed char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_schar* __a__, signed char* __e__, signed char __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_schar* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_schar* __a__, signed char __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline signed char atomic_load_explicit
( volatile atomic_schar* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline signed char atomic_load( volatile atomic_schar* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_schar* __a__, signed char __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline signed char atomic_exchange_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline signed char atomic_exchange
( volatile atomic_schar* __a__, signed char __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_schar* __a__, signed char* __e__, signed char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_schar* __a__, signed char* __e__, signed char __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_schar* __a__, signed char* __e__, signed char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_schar* __a__, signed char* __e__, signed char __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline signed char atomic_fetch_add_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline signed char atomic_fetch_add
(  atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_add_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline signed char atomic_fetch_add
( volatile atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_sub_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline signed char atomic_fetch_sub
(  atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_sub_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline signed char atomic_fetch_sub
( volatile atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_and_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline signed char atomic_fetch_and
(  atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_and_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline signed char atomic_fetch_and
( volatile atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_or_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline signed char atomic_fetch_or
(  atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_or_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline signed char atomic_fetch_or
( volatile atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_xor_explicit
(  atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline signed char atomic_fetch_xor
(  atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline signed char atomic_fetch_xor_explicit
( volatile atomic_schar* __a__, signed char __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline signed char atomic_fetch_xor
( volatile atomic_schar* __a__, signed char __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_uchar
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_uchar() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_uchar( unsigned char __v__ )
    : __f__( __v__) { }
    atomic_uchar( const atomic_uchar& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned char load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (unsigned char)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned char exchange
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned char)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned char load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (unsigned char)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned char exchange
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned char)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned char& __e__, unsigned char __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned char)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_uchar& operator =
    ( const atomic_uchar& )  CXX0X_DELETED
    )

    unsigned char operator =( unsigned char __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned char() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_uchar& operator =
    ( const atomic_uchar& ) volatile CXX0X_DELETED
    )

    unsigned char operator =( unsigned char __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned char() const volatile
    { return load(); }


    unsigned char fetch_add
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned char fetch_add
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned char fetch_sub
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned char fetch_sub
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned char fetch_and
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned char fetch_and
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned char fetch_or
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned char fetch_or
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned char fetch_xor
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned char fetch_xor
    ( unsigned char __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned char operator +=( unsigned char __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    unsigned char operator +=( unsigned char __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    unsigned char operator -=( unsigned char __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    unsigned char operator -=( unsigned char __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    unsigned char operator &=( unsigned char __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    unsigned char operator &=( unsigned char __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    unsigned char operator |=( unsigned char __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    unsigned char operator |=( unsigned char __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    unsigned char operator ^=( unsigned char __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned char operator ^=( unsigned char __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned char operator ++( int ) 
    { return fetch_add( 1 ); }
    unsigned char operator --( int ) 
    { return fetch_sub( 1 ); }
    unsigned char operator ++() 
    { return fetch_add( 1 ) + 1; }
    unsigned char operator --() 
    { return fetch_sub( 1 ) - 1; }


    unsigned char operator ++( int ) volatile
    { return fetch_add( 1 ); }
    unsigned char operator --( int ) volatile
    { return fetch_sub( 1 ); }
    unsigned char operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    unsigned char operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    unsigned char __f__;
} atomic_uchar;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_uchar* __a__ )
{ return false; }

inline void atomic_init
(  atomic_uchar* __a__, unsigned char __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned char atomic_load_explicit
(  atomic_uchar* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned char atomic_load(  atomic_uchar* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_uchar* __a__, unsigned char __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned char atomic_exchange_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned char atomic_exchange
(  atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_uchar* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned char atomic_load_explicit
( volatile atomic_uchar* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned char atomic_load( volatile atomic_uchar* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned char atomic_exchange_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned char atomic_exchange
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_uchar* __a__, unsigned char* __e__, unsigned char __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline unsigned char atomic_fetch_add_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned char atomic_fetch_add
(  atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_add_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned char atomic_fetch_add
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_sub_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned char atomic_fetch_sub
(  atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_sub_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned char atomic_fetch_sub
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_and_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned char atomic_fetch_and
(  atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_and_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned char atomic_fetch_and
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_or_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned char atomic_fetch_or
(  atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_or_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned char atomic_fetch_or
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_xor_explicit
(  atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned char atomic_fetch_xor
(  atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline unsigned char atomic_fetch_xor_explicit
( volatile atomic_uchar* __a__, unsigned char __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned char atomic_fetch_xor
( volatile atomic_uchar* __a__, unsigned char __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_short
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_short() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_short( short __v__ )
    : __f__( __v__) { }
    atomic_short( const atomic_short& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    short load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (short)_ATOMIC_LOAD_( this, __x__ ); }

    short exchange
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (short)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( short& __e__, short __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( short& __e__, short __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( short& __e__, short __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( short& __e__, short __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    short load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (short)_ATOMIC_LOAD_( this, __x__ ); }

    short exchange
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (short)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( short& __e__, short __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( short& __e__, short __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( short& __e__, short __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( short& __e__, short __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_short& operator =
    ( const atomic_short& )  CXX0X_DELETED
    )

    short operator =( short __v__ ) 
    { store( __v__ ); return __v__; }

    operator short() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_short& operator =
    ( const atomic_short& ) volatile CXX0X_DELETED
    )

    short operator =( short __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator short() const volatile
    { return load(); }


    short fetch_add
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    short fetch_add
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    short fetch_sub
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    short fetch_sub
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    short fetch_and
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    short fetch_and
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    short fetch_or
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    short fetch_or
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    short fetch_xor
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    short fetch_xor
    ( short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    short operator +=( short __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    short operator +=( short __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    short operator -=( short __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    short operator -=( short __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    short operator &=( short __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    short operator &=( short __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    short operator |=( short __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    short operator |=( short __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    short operator ^=( short __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    short operator ^=( short __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    short operator ++( int ) 
    { return fetch_add( 1 ); }
    short operator --( int ) 
    { return fetch_sub( 1 ); }
    short operator ++() 
    { return fetch_add( 1 ) + 1; }
    short operator --() 
    { return fetch_sub( 1 ) - 1; }


    short operator ++( int ) volatile
    { return fetch_add( 1 ); }
    short operator --( int ) volatile
    { return fetch_sub( 1 ); }
    short operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    short operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    short __f__;
} atomic_short;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_short* __a__ )
{ return false; }

inline void atomic_init
(  atomic_short* __a__, short __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline short atomic_load_explicit
(  atomic_short* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline short atomic_load(  atomic_short* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_short* __a__, short __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline short atomic_exchange_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline short atomic_exchange
(  atomic_short* __a__, short __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_short* __a__, short* __e__, short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_short* __a__, short* __e__, short __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_short* __a__, short* __e__, short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_short* __a__, short* __e__, short __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_short* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_short* __a__, short __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline short atomic_load_explicit
( volatile atomic_short* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline short atomic_load( volatile atomic_short* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_short* __a__, short __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline short atomic_exchange_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline short atomic_exchange
( volatile atomic_short* __a__, short __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_short* __a__, short* __e__, short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_short* __a__, short* __e__, short __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_short* __a__, short* __e__, short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_short* __a__, short* __e__, short __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline short atomic_fetch_add_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline short atomic_fetch_add
(  atomic_short* __a__, short __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_add_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline short atomic_fetch_add
( volatile atomic_short* __a__, short __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_sub_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline short atomic_fetch_sub
(  atomic_short* __a__, short __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_sub_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline short atomic_fetch_sub
( volatile atomic_short* __a__, short __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_and_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline short atomic_fetch_and
(  atomic_short* __a__, short __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_and_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline short atomic_fetch_and
( volatile atomic_short* __a__, short __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_or_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline short atomic_fetch_or
(  atomic_short* __a__, short __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_or_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline short atomic_fetch_or
( volatile atomic_short* __a__, short __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_xor_explicit
(  atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline short atomic_fetch_xor
(  atomic_short* __a__, short __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline short atomic_fetch_xor_explicit
( volatile atomic_short* __a__, short __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline short atomic_fetch_xor
( volatile atomic_short* __a__, short __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_ushort
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_ushort() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_ushort( unsigned short __v__ )
    : __f__( __v__) { }
    atomic_ushort( const atomic_ushort& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned short load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (unsigned short)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned short exchange
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned short)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned short load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (unsigned short)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned short exchange
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned short)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned short& __e__, unsigned short __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned short)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_ushort& operator =
    ( const atomic_ushort& )  CXX0X_DELETED
    )

    unsigned short operator =( unsigned short __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned short() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_ushort& operator =
    ( const atomic_ushort& ) volatile CXX0X_DELETED
    )

    unsigned short operator =( unsigned short __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned short() const volatile
    { return load(); }


    unsigned short fetch_add
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned short fetch_add
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned short fetch_sub
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned short fetch_sub
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned short fetch_and
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned short fetch_and
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned short fetch_or
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned short fetch_or
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned short fetch_xor
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned short fetch_xor
    ( unsigned short __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned short operator +=( unsigned short __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    unsigned short operator +=( unsigned short __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    unsigned short operator -=( unsigned short __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    unsigned short operator -=( unsigned short __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    unsigned short operator &=( unsigned short __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    unsigned short operator &=( unsigned short __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    unsigned short operator |=( unsigned short __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    unsigned short operator |=( unsigned short __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    unsigned short operator ^=( unsigned short __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned short operator ^=( unsigned short __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned short operator ++( int ) 
    { return fetch_add( 1 ); }
    unsigned short operator --( int ) 
    { return fetch_sub( 1 ); }
    unsigned short operator ++() 
    { return fetch_add( 1 ) + 1; }
    unsigned short operator --() 
    { return fetch_sub( 1 ) - 1; }


    unsigned short operator ++( int ) volatile
    { return fetch_add( 1 ); }
    unsigned short operator --( int ) volatile
    { return fetch_sub( 1 ); }
    unsigned short operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    unsigned short operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    unsigned short __f__;
} atomic_ushort;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_ushort* __a__ )
{ return false; }

inline void atomic_init
(  atomic_ushort* __a__, unsigned short __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned short atomic_load_explicit
(  atomic_ushort* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned short atomic_load(  atomic_ushort* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_ushort* __a__, unsigned short __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned short atomic_exchange_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned short atomic_exchange
(  atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_ushort* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned short atomic_load_explicit
( volatile atomic_ushort* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned short atomic_load( volatile atomic_ushort* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned short atomic_exchange_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned short atomic_exchange
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_ushort* __a__, unsigned short* __e__, unsigned short __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline unsigned short atomic_fetch_add_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned short atomic_fetch_add
(  atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_add_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned short atomic_fetch_add
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_sub_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned short atomic_fetch_sub
(  atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_sub_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned short atomic_fetch_sub
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_and_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned short atomic_fetch_and
(  atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_and_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned short atomic_fetch_and
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_or_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned short atomic_fetch_or
(  atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_or_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned short atomic_fetch_or
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_xor_explicit
(  atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned short atomic_fetch_xor
(  atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline unsigned short atomic_fetch_xor_explicit
( volatile atomic_ushort* __a__, unsigned short __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned short atomic_fetch_xor
( volatile atomic_ushort* __a__, unsigned short __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_int
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_int() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_int( int __v__ )
    : __f__( __v__) { }
    atomic_int( const atomic_int& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    int load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (int)_ATOMIC_LOAD_( this, __x__ ); }

    int exchange
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (int)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( int& __e__, int __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( int& __e__, int __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( int& __e__, int __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( int& __e__, int __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    int load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (int)_ATOMIC_LOAD_( this, __x__ ); }

    int exchange
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (int)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( int& __e__, int __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( int& __e__, int __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( int& __e__, int __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( int& __e__, int __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_int& operator =
    ( const atomic_int& )  CXX0X_DELETED
    )

    int operator =( int __v__ ) 
    { store( __v__ ); return __v__; }

    operator int() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_int& operator =
    ( const atomic_int& ) volatile CXX0X_DELETED
    )

    int operator =( int __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator int() const volatile
    { return load(); }


    int fetch_add
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    int fetch_add
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    int fetch_sub
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    int fetch_sub
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    int fetch_and
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    int fetch_and
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    int fetch_or
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    int fetch_or
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    int fetch_xor
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    int fetch_xor
    ( int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    int operator +=( int __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    int operator +=( int __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    int operator -=( int __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    int operator -=( int __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    int operator &=( int __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    int operator &=( int __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    int operator |=( int __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    int operator |=( int __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    int operator ^=( int __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    int operator ^=( int __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    int operator ++( int ) 
    { return fetch_add( 1 ); }
    int operator --( int ) 
    { return fetch_sub( 1 ); }
    int operator ++() 
    { return fetch_add( 1 ) + 1; }
    int operator --() 
    { return fetch_sub( 1 ) - 1; }


    int operator ++( int ) volatile
    { return fetch_add( 1 ); }
    int operator --( int ) volatile
    { return fetch_sub( 1 ); }
    int operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    int operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    int __f__;
} atomic_int;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_int* __a__ )
{ return false; }

inline void atomic_init
(  atomic_int* __a__, int __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline int atomic_load_explicit
(  atomic_int* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline int atomic_load(  atomic_int* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_int* __a__, int __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline int atomic_exchange_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline int atomic_exchange
(  atomic_int* __a__, int __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_int* __a__, int* __e__, int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_int* __a__, int* __e__, int __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_int* __a__, int* __e__, int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_int* __a__, int* __e__, int __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_int* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_int* __a__, int __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline int atomic_load_explicit
( volatile atomic_int* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline int atomic_load( volatile atomic_int* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_int* __a__, int __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline int atomic_exchange_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline int atomic_exchange
( volatile atomic_int* __a__, int __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_int* __a__, int* __e__, int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_int* __a__, int* __e__, int __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_int* __a__, int* __e__, int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_int* __a__, int* __e__, int __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline int atomic_fetch_add_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline int atomic_fetch_add
(  atomic_int* __a__, int __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_add_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline int atomic_fetch_add
( volatile atomic_int* __a__, int __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_sub_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline int atomic_fetch_sub
(  atomic_int* __a__, int __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_sub_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline int atomic_fetch_sub
( volatile atomic_int* __a__, int __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_and_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline int atomic_fetch_and
(  atomic_int* __a__, int __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_and_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline int atomic_fetch_and
( volatile atomic_int* __a__, int __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_or_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline int atomic_fetch_or
(  atomic_int* __a__, int __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_or_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline int atomic_fetch_or
( volatile atomic_int* __a__, int __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_xor_explicit
(  atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline int atomic_fetch_xor
(  atomic_int* __a__, int __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline int atomic_fetch_xor_explicit
( volatile atomic_int* __a__, int __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline int atomic_fetch_xor
( volatile atomic_int* __a__, int __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_uint
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_uint() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_uint( unsigned int __v__ )
    : __f__( __v__) { }
    atomic_uint( const atomic_uint& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned int load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (unsigned int)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned int exchange
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned int)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned int load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (unsigned int)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned int exchange
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned int)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned int& __e__, unsigned int __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned int)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_uint& operator =
    ( const atomic_uint& )  CXX0X_DELETED
    )

    unsigned int operator =( unsigned int __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned int() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_uint& operator =
    ( const atomic_uint& ) volatile CXX0X_DELETED
    )

    unsigned int operator =( unsigned int __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned int() const volatile
    { return load(); }


    unsigned int fetch_add
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned int fetch_add
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned int fetch_sub
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned int fetch_sub
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned int fetch_and
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned int fetch_and
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned int fetch_or
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned int fetch_or
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned int fetch_xor
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned int fetch_xor
    ( unsigned int __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned int operator +=( unsigned int __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    unsigned int operator +=( unsigned int __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    unsigned int operator -=( unsigned int __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    unsigned int operator -=( unsigned int __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    unsigned int operator &=( unsigned int __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    unsigned int operator &=( unsigned int __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    unsigned int operator |=( unsigned int __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    unsigned int operator |=( unsigned int __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    unsigned int operator ^=( unsigned int __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned int operator ^=( unsigned int __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned int operator ++( int ) 
    { return fetch_add( 1 ); }
    unsigned int operator --( int ) 
    { return fetch_sub( 1 ); }
    unsigned int operator ++() 
    { return fetch_add( 1 ) + 1; }
    unsigned int operator --() 
    { return fetch_sub( 1 ) - 1; }


    unsigned int operator ++( int ) volatile
    { return fetch_add( 1 ); }
    unsigned int operator --( int ) volatile
    { return fetch_sub( 1 ); }
    unsigned int operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    unsigned int operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    unsigned int __f__;
} atomic_uint;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_uint* __a__ )
{ return false; }

inline void atomic_init
(  atomic_uint* __a__, unsigned int __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned int atomic_load_explicit
(  atomic_uint* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned int atomic_load(  atomic_uint* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_uint* __a__, unsigned int __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned int atomic_exchange_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned int atomic_exchange
(  atomic_uint* __a__, unsigned int __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_uint* __a__, unsigned int* __e__, unsigned int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_uint* __a__, unsigned int* __e__, unsigned int __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_uint* __a__, unsigned int* __e__, unsigned int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_uint* __a__, unsigned int* __e__, unsigned int __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_uint* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_uint* __a__, unsigned int __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned int atomic_load_explicit
( volatile atomic_uint* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned int atomic_load( volatile atomic_uint* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_uint* __a__, unsigned int __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned int atomic_exchange_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned int atomic_exchange
( volatile atomic_uint* __a__, unsigned int __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_uint* __a__, unsigned int* __e__, unsigned int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_uint* __a__, unsigned int* __e__, unsigned int __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_uint* __a__, unsigned int* __e__, unsigned int __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_uint* __a__, unsigned int* __e__, unsigned int __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline unsigned int atomic_fetch_add_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned int atomic_fetch_add
(  atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_add_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned int atomic_fetch_add
( volatile atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_sub_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned int atomic_fetch_sub
(  atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_sub_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned int atomic_fetch_sub
( volatile atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_and_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned int atomic_fetch_and
(  atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_and_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned int atomic_fetch_and
( volatile atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_or_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned int atomic_fetch_or
(  atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_or_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned int atomic_fetch_or
( volatile atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_xor_explicit
(  atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned int atomic_fetch_xor
(  atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline unsigned int atomic_fetch_xor_explicit
( volatile atomic_uint* __a__, unsigned int __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned int atomic_fetch_xor
( volatile atomic_uint* __a__, unsigned int __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_long
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_long() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_long( long __v__ )
    : __f__( __v__) { }
    atomic_long( const atomic_long& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    long load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (long)_ATOMIC_LOAD_( this, __x__ ); }

    long exchange
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long& __e__, long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long& __e__, long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long& __e__, long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long& __e__, long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    long load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (long)_ATOMIC_LOAD_( this, __x__ ); }

    long exchange
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long& __e__, long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long& __e__, long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long& __e__, long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long& __e__, long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_long& operator =
    ( const atomic_long& )  CXX0X_DELETED
    )

    long operator =( long __v__ ) 
    { store( __v__ ); return __v__; }

    operator long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_long& operator =
    ( const atomic_long& ) volatile CXX0X_DELETED
    )

    long operator =( long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator long() const volatile
    { return load(); }


    long fetch_add
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    long fetch_add
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    long fetch_sub
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    long fetch_sub
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    long fetch_and
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    long fetch_and
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    long fetch_or
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    long fetch_or
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    long fetch_xor
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    long fetch_xor
    ( long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    long operator +=( long __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    long operator +=( long __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    long operator -=( long __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    long operator -=( long __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    long operator &=( long __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    long operator &=( long __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    long operator |=( long __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    long operator |=( long __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    long operator ^=( long __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    long operator ^=( long __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    long operator ++( int ) 
    { return fetch_add( 1 ); }
    long operator --( int ) 
    { return fetch_sub( 1 ); }
    long operator ++() 
    { return fetch_add( 1 ) + 1; }
    long operator --() 
    { return fetch_sub( 1 ) - 1; }


    long operator ++( int ) volatile
    { return fetch_add( 1 ); }
    long operator --( int ) volatile
    { return fetch_sub( 1 ); }
    long operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    long operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    long __f__;
} atomic_long;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_long* __a__ )
{ return false; }

inline void atomic_init
(  atomic_long* __a__, long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline long atomic_load_explicit
(  atomic_long* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline long atomic_load(  atomic_long* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_long* __a__, long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline long atomic_exchange_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline long atomic_exchange
(  atomic_long* __a__, long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_long* __a__, long* __e__, long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_long* __a__, long* __e__, long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_long* __a__, long* __e__, long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_long* __a__, long* __e__, long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_long* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_long* __a__, long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline long atomic_load_explicit
( volatile atomic_long* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline long atomic_load( volatile atomic_long* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_long* __a__, long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline long atomic_exchange_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline long atomic_exchange
( volatile atomic_long* __a__, long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_long* __a__, long* __e__, long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_long* __a__, long* __e__, long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_long* __a__, long* __e__, long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_long* __a__, long* __e__, long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline long atomic_fetch_add_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline long atomic_fetch_add
(  atomic_long* __a__, long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_add_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline long atomic_fetch_add
( volatile atomic_long* __a__, long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_sub_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline long atomic_fetch_sub
(  atomic_long* __a__, long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_sub_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline long atomic_fetch_sub
( volatile atomic_long* __a__, long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_and_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline long atomic_fetch_and
(  atomic_long* __a__, long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_and_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline long atomic_fetch_and
( volatile atomic_long* __a__, long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_or_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline long atomic_fetch_or
(  atomic_long* __a__, long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_or_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline long atomic_fetch_or
( volatile atomic_long* __a__, long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_xor_explicit
(  atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline long atomic_fetch_xor
(  atomic_long* __a__, long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline long atomic_fetch_xor_explicit
( volatile atomic_long* __a__, long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline long atomic_fetch_xor
( volatile atomic_long* __a__, long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_ulong
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_ulong() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_ulong( unsigned long __v__ )
    : __f__( __v__) { }
    atomic_ulong( const atomic_ulong& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned long load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (unsigned long)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned long exchange
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned long load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (unsigned long)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned long exchange
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long& __e__, unsigned long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_ulong& operator =
    ( const atomic_ulong& )  CXX0X_DELETED
    )

    unsigned long operator =( unsigned long __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_ulong& operator =
    ( const atomic_ulong& ) volatile CXX0X_DELETED
    )

    unsigned long operator =( unsigned long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned long() const volatile
    { return load(); }


    unsigned long fetch_add
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned long fetch_add
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned long fetch_sub
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned long fetch_sub
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned long fetch_and
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned long fetch_and
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned long fetch_or
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned long fetch_or
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned long fetch_xor
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned long fetch_xor
    ( unsigned long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned long operator +=( unsigned long __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    unsigned long operator +=( unsigned long __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    unsigned long operator -=( unsigned long __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    unsigned long operator -=( unsigned long __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    unsigned long operator &=( unsigned long __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    unsigned long operator &=( unsigned long __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    unsigned long operator |=( unsigned long __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    unsigned long operator |=( unsigned long __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    unsigned long operator ^=( unsigned long __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned long operator ^=( unsigned long __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned long operator ++( int ) 
    { return fetch_add( 1 ); }
    unsigned long operator --( int ) 
    { return fetch_sub( 1 ); }
    unsigned long operator ++() 
    { return fetch_add( 1 ) + 1; }
    unsigned long operator --() 
    { return fetch_sub( 1 ) - 1; }


    unsigned long operator ++( int ) volatile
    { return fetch_add( 1 ); }
    unsigned long operator --( int ) volatile
    { return fetch_sub( 1 ); }
    unsigned long operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    unsigned long operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    unsigned long __f__;
} atomic_ulong;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_ulong* __a__ )
{ return false; }

inline void atomic_init
(  atomic_ulong* __a__, unsigned long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned long atomic_load_explicit
(  atomic_ulong* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned long atomic_load(  atomic_ulong* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_ulong* __a__, unsigned long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned long atomic_exchange_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned long atomic_exchange
(  atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_ulong* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned long atomic_load_explicit
( volatile atomic_ulong* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned long atomic_load( volatile atomic_ulong* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned long atomic_exchange_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned long atomic_exchange
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_ulong* __a__, unsigned long* __e__, unsigned long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline unsigned long atomic_fetch_add_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned long atomic_fetch_add
(  atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_add_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned long atomic_fetch_add
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_sub_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned long atomic_fetch_sub
(  atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_sub_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned long atomic_fetch_sub
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_and_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned long atomic_fetch_and
(  atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_and_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned long atomic_fetch_and
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_or_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned long atomic_fetch_or
(  atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_or_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned long atomic_fetch_or
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_xor_explicit
(  atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned long atomic_fetch_xor
(  atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline unsigned long atomic_fetch_xor_explicit
( volatile atomic_ulong* __a__, unsigned long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned long atomic_fetch_xor
( volatile atomic_ulong* __a__, unsigned long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_llong
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_llong() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_llong( long long __v__ )
    : __f__( __v__) { }
    atomic_llong( const atomic_llong& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    long long load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (long long)_ATOMIC_LOAD_( this, __x__ ); }

    long long exchange
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (long long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long long& __e__, long long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long long& __e__, long long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long long& __e__, long long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long long& __e__, long long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    long long load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (long long)_ATOMIC_LOAD_( this, __x__ ); }

    long long exchange
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (long long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long long& __e__, long long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( long long& __e__, long long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long long& __e__, long long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( long long& __e__, long long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_llong& operator =
    ( const atomic_llong& )  CXX0X_DELETED
    )

    long long operator =( long long __v__ ) 
    { store( __v__ ); return __v__; }

    operator long long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_llong& operator =
    ( const atomic_llong& ) volatile CXX0X_DELETED
    )

    long long operator =( long long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator long long() const volatile
    { return load(); }


    long long fetch_add
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    long long fetch_add
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    long long fetch_sub
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    long long fetch_sub
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    long long fetch_and
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    long long fetch_and
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    long long fetch_or
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    long long fetch_or
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    long long fetch_xor
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    long long fetch_xor
    ( long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    long long operator +=( long long __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    long long operator +=( long long __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    long long operator -=( long long __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    long long operator -=( long long __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    long long operator &=( long long __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    long long operator &=( long long __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    long long operator |=( long long __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    long long operator |=( long long __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    long long operator ^=( long long __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    long long operator ^=( long long __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    long long operator ++( int ) 
    { return fetch_add( 1 ); }
    long long operator --( int ) 
    { return fetch_sub( 1 ); }
    long long operator ++() 
    { return fetch_add( 1 ) + 1; }
    long long operator --() 
    { return fetch_sub( 1 ) - 1; }


    long long operator ++( int ) volatile
    { return fetch_add( 1 ); }
    long long operator --( int ) volatile
    { return fetch_sub( 1 ); }
    long long operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    long long operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    long long __f__;
} atomic_llong;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_llong* __a__ )
{ return false; }

inline void atomic_init
(  atomic_llong* __a__, long long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline long long atomic_load_explicit
(  atomic_llong* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline long long atomic_load(  atomic_llong* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_llong* __a__, long long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline long long atomic_exchange_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline long long atomic_exchange
(  atomic_llong* __a__, long long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_llong* __a__, long long* __e__, long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_llong* __a__, long long* __e__, long long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_llong* __a__, long long* __e__, long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_llong* __a__, long long* __e__, long long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_llong* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_llong* __a__, long long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline long long atomic_load_explicit
( volatile atomic_llong* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline long long atomic_load( volatile atomic_llong* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_llong* __a__, long long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline long long atomic_exchange_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline long long atomic_exchange
( volatile atomic_llong* __a__, long long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_llong* __a__, long long* __e__, long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_llong* __a__, long long* __e__, long long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_llong* __a__, long long* __e__, long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_llong* __a__, long long* __e__, long long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline long long atomic_fetch_add_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline long long atomic_fetch_add
(  atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_add_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline long long atomic_fetch_add
( volatile atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_sub_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline long long atomic_fetch_sub
(  atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_sub_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline long long atomic_fetch_sub
( volatile atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_and_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline long long atomic_fetch_and
(  atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_and_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline long long atomic_fetch_and
( volatile atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_or_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline long long atomic_fetch_or
(  atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_or_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline long long atomic_fetch_or
( volatile atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_xor_explicit
(  atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline long long atomic_fetch_xor
(  atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline long long atomic_fetch_xor_explicit
( volatile atomic_llong* __a__, long long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline long long atomic_fetch_xor
( volatile atomic_llong* __a__, long long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


typedef struct atomic_ullong
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_ullong() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_ullong( unsigned long long __v__ )
    : __f__( __v__) { }
    atomic_ullong( const atomic_ullong& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned long long load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (unsigned long long)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned long long exchange
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned long long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    unsigned long long load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (unsigned long long)_ATOMIC_LOAD_( this, __x__ ); }

    unsigned long long exchange
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned long long)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( unsigned long long& __e__, unsigned long long __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (unsigned long long)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_ullong& operator =
    ( const atomic_ullong& )  CXX0X_DELETED
    )

    unsigned long long operator =( unsigned long long __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned long long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_ullong& operator =
    ( const atomic_ullong& ) volatile CXX0X_DELETED
    )

    unsigned long long operator =( unsigned long long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned long long() const volatile
    { return load(); }


    unsigned long long fetch_add
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned long long fetch_add
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    unsigned long long fetch_sub
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned long long fetch_sub
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    unsigned long long fetch_and
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned long long fetch_and
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    unsigned long long fetch_or
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned long long fetch_or
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    unsigned long long fetch_xor
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned long long fetch_xor
    ( unsigned long long __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    unsigned long long operator +=( unsigned long long __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    unsigned long long operator +=( unsigned long long __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    unsigned long long operator -=( unsigned long long __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    unsigned long long operator -=( unsigned long long __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    unsigned long long operator &=( unsigned long long __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    unsigned long long operator &=( unsigned long long __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    unsigned long long operator |=( unsigned long long __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    unsigned long long operator |=( unsigned long long __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    unsigned long long operator ^=( unsigned long long __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned long long operator ^=( unsigned long long __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    unsigned long long operator ++( int ) 
    { return fetch_add( 1 ); }
    unsigned long long operator --( int ) 
    { return fetch_sub( 1 ); }
    unsigned long long operator ++() 
    { return fetch_add( 1 ) + 1; }
    unsigned long long operator --() 
    { return fetch_sub( 1 ) - 1; }


    unsigned long long operator ++( int ) volatile
    { return fetch_add( 1 ); }
    unsigned long long operator --( int ) volatile
    { return fetch_sub( 1 ); }
    unsigned long long operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    unsigned long long operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    unsigned long long __f__;
} atomic_ullong;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_ullong* __a__ )
{ return false; }

inline void atomic_init
(  atomic_ullong* __a__, unsigned long long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned long long atomic_load_explicit
(  atomic_ullong* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned long long atomic_load(  atomic_ullong* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_ullong* __a__, unsigned long long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned long long atomic_exchange_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned long long atomic_exchange
(  atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_ullong* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline unsigned long long atomic_load_explicit
( volatile atomic_ullong* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline unsigned long long atomic_load( volatile atomic_ullong* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline unsigned long long atomic_exchange_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline unsigned long long atomic_exchange
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_ullong* __a__, unsigned long long* __e__, unsigned long long __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline unsigned long long atomic_fetch_add_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned long long atomic_fetch_add
(  atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_add_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline unsigned long long atomic_fetch_add
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_sub_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned long long atomic_fetch_sub
(  atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_sub_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline unsigned long long atomic_fetch_sub
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_and_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned long long atomic_fetch_and
(  atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_and_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline unsigned long long atomic_fetch_and
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_or_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned long long atomic_fetch_or
(  atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_or_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline unsigned long long atomic_fetch_or
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_xor_explicit
(  atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned long long atomic_fetch_xor
(  atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline unsigned long long atomic_fetch_xor_explicit
( volatile atomic_ullong* __a__, unsigned long long __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline unsigned long long atomic_fetch_xor
( volatile atomic_ullong* __a__, unsigned long long __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


typedef struct atomic_wchar_t
{


    #ifdef __cplusplus
    CXX0X_AGGR_INIT(
    atomic_wchar_t() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic_wchar_t( wchar_t __v__ )
    : __f__( __v__) { }
    atomic_wchar_t( const atomic_wchar_t& ) CXX0X_DELETED
    )


    bool is_lock_free() const 
    { return false; }

    void store
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    wchar_t load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (wchar_t)_ATOMIC_LOAD_( this, __x__ ); }

    wchar_t exchange
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (wchar_t)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    wchar_t load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (wchar_t)_ATOMIC_LOAD_( this, __x__ ); }

    wchar_t exchange
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (wchar_t)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( wchar_t& __e__, wchar_t __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (wchar_t)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_wchar_t& operator =
    ( const atomic_wchar_t& )  CXX0X_DELETED
    )

    wchar_t operator =( wchar_t __v__ ) 
    { store( __v__ ); return __v__; }

    operator wchar_t() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_wchar_t& operator =
    ( const atomic_wchar_t& ) volatile CXX0X_DELETED
    )

    wchar_t operator =( wchar_t __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator wchar_t() const volatile
    { return load(); }


    wchar_t fetch_add
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    wchar_t fetch_add
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, +=, __m__, __x__ ); }


    wchar_t fetch_sub
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    wchar_t fetch_sub
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, -=, __m__, __x__ ); }


    wchar_t fetch_and
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    wchar_t fetch_and
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, &=, __m__, __x__ ); }


    wchar_t fetch_or
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    wchar_t fetch_or
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, |=, __m__, __x__ ); }


    wchar_t fetch_xor
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    wchar_t fetch_xor
    ( wchar_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_MODIFY_( this, ^=, __m__, __x__ ); }


    wchar_t operator +=( wchar_t __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    wchar_t operator +=( wchar_t __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    wchar_t operator -=( wchar_t __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    wchar_t operator -=( wchar_t __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    wchar_t operator &=( wchar_t __v__ ) 
    { return fetch_and( __v__ ) & __v__; }


    wchar_t operator &=( wchar_t __v__ ) volatile
    { return fetch_and( __v__ ) & __v__; }


    wchar_t operator |=( wchar_t __v__ ) 
    { return fetch_or( __v__ ) | __v__; }


    wchar_t operator |=( wchar_t __v__ ) volatile
    { return fetch_or( __v__ ) | __v__; }


    wchar_t operator ^=( wchar_t __v__ ) 
    { return fetch_xor( __v__ ) ^ __v__; }


    wchar_t operator ^=( wchar_t __v__ ) volatile
    { return fetch_xor( __v__ ) ^ __v__; }


    wchar_t operator ++( int ) 
    { return fetch_add( 1 ); }
    wchar_t operator --( int ) 
    { return fetch_sub( 1 ); }
    wchar_t operator ++() 
    { return fetch_add( 1 ) + 1; }
    wchar_t operator --() 
    { return fetch_sub( 1 ) - 1; }


    wchar_t operator ++( int ) volatile
    { return fetch_add( 1 ); }
    wchar_t operator --( int ) volatile
    { return fetch_sub( 1 ); }
    wchar_t operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    wchar_t operator --() volatile
    { return fetch_sub( 1 ) - 1; }


    CXX0X_TRIVIAL_PRIVATE
    #endif
    wchar_t __f__;
} atomic_wchar_t;


#ifdef __cplusplus


inline bool atomic_is_lock_free( const  atomic_wchar_t* __a__ )
{ return false; }

inline void atomic_init
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline wchar_t atomic_load_explicit
(  atomic_wchar_t* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline wchar_t atomic_load(  atomic_wchar_t* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline wchar_t atomic_exchange_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline wchar_t atomic_exchange
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
(  atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
(  atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
(  atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
(  atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


inline bool atomic_is_lock_free( const volatile atomic_wchar_t* __a__ )
{ return false; }

inline void atomic_init
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ __a__->store( __m__, memory_order_relaxed ); }

inline wchar_t atomic_load_explicit
( volatile atomic_wchar_t* __a__, memory_order __x__ )
{ return __a__->load( __x__ );; }

inline wchar_t atomic_load( volatile atomic_wchar_t* __a__ )
{ return __a__->load( memory_order_seq_cst ); }

inline void atomic_store_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ __a__->store( __m__, __x__ ); }

inline void atomic_store
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ __a__->store( __m__, memory_order_seq_cst ); }

inline wchar_t atomic_exchange_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->exchange( __m__, __x__ ); }

inline wchar_t atomic_exchange
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->exchange( __m__, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_weak_explicit
( volatile atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_weak( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_weak
( volatile atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__ )
{ return __a__->compare_exchange_weak( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }

inline bool atomic_compare_exchange_strong_explicit
( volatile atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__,
    memory_order __x__, memory_order __y__ )
{ return __a__->compare_exchange_strong( *__e__, __m__, __x__, __y__ ); }

inline bool atomic_compare_exchange_strong
( volatile atomic_wchar_t* __a__, wchar_t* __e__, wchar_t __m__ )
{ return __a__->compare_exchange_strong( *__e__, __m__,
    memory_order_seq_cst, memory_order_seq_cst ); }


#endif


#ifdef __cplusplus


inline wchar_t atomic_fetch_add_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline wchar_t atomic_fetch_add
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_add_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_add( __m__, __x__ ); }

inline wchar_t atomic_fetch_add
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_add( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_sub_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline wchar_t atomic_fetch_sub
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_sub_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_sub( __m__, __x__ ); }

inline wchar_t atomic_fetch_sub
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_sub( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_and_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline wchar_t atomic_fetch_and
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_and_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_and( __m__, __x__ ); }

inline wchar_t atomic_fetch_and
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_and( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_or_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline wchar_t atomic_fetch_or
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_or_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_or( __m__, __x__ ); }

inline wchar_t atomic_fetch_or
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_or( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_xor_explicit
(  atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline wchar_t atomic_fetch_xor
(  atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


inline wchar_t atomic_fetch_xor_explicit
( volatile atomic_wchar_t* __a__, wchar_t __m__, memory_order __x__ )
{ return __a__->fetch_xor( __m__, __x__ ); }

inline wchar_t atomic_fetch_xor
( volatile atomic_wchar_t* __a__, wchar_t __m__ )
{ return __a__->fetch_xor( __m__, memory_order_seq_cst ); }


#endif


#endif


typedef atomic_schar atomic_int_least8_t;
typedef atomic_uchar atomic_uint_least8_t;
typedef atomic_short atomic_int_least16_t;
typedef atomic_ushort atomic_uint_least16_t;
typedef atomic_int atomic_int_least32_t;
typedef atomic_uint atomic_uint_least32_t;
typedef atomic_llong atomic_int_least64_t;
typedef atomic_ullong atomic_uint_least64_t;

typedef atomic_schar atomic_int_fast8_t;
typedef atomic_uchar atomic_uint_fast8_t;
typedef atomic_short atomic_int_fast16_t;
typedef atomic_ushort atomic_uint_fast16_t;
typedef atomic_int atomic_int_fast32_t;
typedef atomic_uint atomic_uint_fast32_t;
typedef atomic_llong atomic_int_fast64_t;
typedef atomic_ullong atomic_uint_fast64_t;

typedef atomic_long atomic_intptr_t;
typedef atomic_ulong atomic_uintptr_t;

typedef atomic_long atomic_ssize_t;
typedef atomic_ulong atomic_size_t;

typedef atomic_long atomic_ptrdiff_t;

typedef atomic_llong atomic_intmax_t;
typedef atomic_ullong atomic_uintmax_t;

#ifndef __cplusplus

typedef atomic_int_least16_t atomic_char16_t;
typedef atomic_int_least32_t atomic_char32_t;
typedef atomic_int_least32_t atomic_wchar_t;

#endif


#ifdef __cplusplus


template< typename T >
struct atomic
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( T* __v__ )
    CXX0X_AGGR_INIT( : __f__( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    bool is_lock_free() const 
    { return false; }

    void store
    ( T __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    T load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (T)_ATOMIC_LOAD_( this, __x__ ); }

    T exchange
    ( T __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (T)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T& __e__, T __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T& __e__, T __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T& __e__, T __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T& __e__, T __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( T __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    T load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (T)_ATOMIC_LOAD_( this, __x__ ); }

    T exchange
    ( T __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (T)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T& __e__, T __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T& __e__, T __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T& __e__, T __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T& __e__, T __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (T)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic& operator =
    ( const atomic& )  CXX0X_DELETED
    )

    T operator =( T __v__ ) 
    { store( __v__ ); return __v__; }

    operator T() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic& operator =
    ( const atomic& ) volatile CXX0X_DELETED
    )

    T operator =( T __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator T() const volatile
    { return load(); }


    CXX0X_TRIVIAL_PRIVATE
    T __f__;
};


template< typename T >
struct atomic< T* >
: atomic_address
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( T* __v__ )
    CXX0X_AGGR_INIT( : atomic_address( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    bool is_lock_free() const 
    { return false; }

    void store
    ( T* __m__, memory_order __x__ = memory_order_seq_cst ) 
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    T* load
    ( memory_order __x__ = memory_order_seq_cst ) const 
    { return (T*)_ATOMIC_LOAD_( this, __x__ ); }

    T* exchange
    ( T* __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return (T*)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T*& __e__, T* __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T*& __e__, T* __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T*& __e__, T* __m__,
      memory_order __x__, memory_order __y__ ) 
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T*& __e__, T* __m__,
      memory_order __x__ = memory_order_seq_cst ) 
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    bool is_lock_free() const volatile
    { return false; }

    void store
    ( T* __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { (void)_ATOMIC_STORE_( this, __m__, __x__ ); }

    T* load
    ( memory_order __x__ = memory_order_seq_cst ) const volatile
    { return (T*)_ATOMIC_LOAD_( this, __x__ ); }

    T* exchange
    ( T* __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return (T*)_ATOMIC_MODIFY_( this, =, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T*& __e__, T* __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_weak
    ( T*& __e__, T* __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T*& __e__, T* __m__,
      memory_order __x__, memory_order __y__ ) volatile
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }

    bool compare_exchange_strong
    ( T*& __e__, T* __m__,
      memory_order __x__ = memory_order_seq_cst ) volatile
    { return (T*)_ATOMIC_CMPSWP_( this, &__e__, __m__, __x__ ); }


    CXX0X_AGGR_INIT(
    atomic_address& operator =
    ( const atomic_address& )  CXX0X_DELETED
    )

    T* operator =( T* __v__ ) 
    { store( __v__ ); return __v__; }

    operator T*() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_address& operator =
    ( const atomic_address& ) volatile CXX0X_DELETED
    )

    T* operator =( T* __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator T*() const volatile
    { return load(); }


    T* fetch_add
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_PTRMOD_( this, +=, __m__ * sizeof (T), __x__ ); }

    T* fetch_sub
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) 
    { return _ATOMIC_PTRMOD_( this, -=, __m__ * sizeof (T), __x__ ); }


    T* fetch_add
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_PTRMOD_( this, +=, __m__ * sizeof (T), __x__ ); }

    T* fetch_sub
    ( ptrdiff_t __m__, memory_order __x__ = memory_order_seq_cst ) volatile
    { return _ATOMIC_PTRMOD_( this, -=, __m__ * sizeof (T), __x__ ); }


    T* operator +=( ptrdiff_t __v__ ) 
    { return fetch_add( __v__ ) + __v__; }


    T* operator +=( ptrdiff_t __v__ ) volatile
    { return fetch_add( __v__ ) + __v__; }


    T* operator -=( ptrdiff_t __v__ ) 
    { return fetch_sub( __v__ ) - __v__; }


    T* operator -=( ptrdiff_t __v__ ) volatile
    { return fetch_sub( __v__ ) - __v__; }


    T* operator ++( int ) 
    { return fetch_add( 1 ); }
    T* operator --( int ) 
    { return fetch_sub( 1 ); }
    T* operator ++() 
    { return fetch_add( 1 ) + 1; }
    T* operator --() 
    { return fetch_sub( 1 ) - 1; }


    T* operator ++( int ) volatile
    { return fetch_add( 1 ); }
    T* operator --( int ) volatile
    { return fetch_sub( 1 ); }
    T* operator ++() volatile
    { return fetch_add( 1 ) + 1; }
    T* operator --() volatile
    { return fetch_sub( 1 ) - 1; }


};


template<>
struct atomic< bool >
: atomic_bool
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( bool __v__ )
    CXX0X_AGGR_INIT( : atomic_bool( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_bool& operator =
    ( const atomic_bool& )  CXX0X_DELETED
    )

    bool operator =( bool __v__ ) 
    { store( __v__ ); return __v__; }

    operator bool() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_bool& operator =
    ( const atomic_bool& ) volatile CXX0X_DELETED
    )

    bool operator =( bool __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator bool() const volatile
    { return load(); }


};


template<>
struct atomic< void* >
: atomic_address
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( void* __v__ )
    CXX0X_AGGR_INIT( : atomic_address( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_address& operator =
    ( const atomic_address& )  CXX0X_DELETED
    )

    void* operator =( void* __v__ ) 
    { store( __v__ ); return __v__; }

    operator void*() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_address& operator =
    ( const atomic_address& ) volatile CXX0X_DELETED
    )

    void* operator =( void* __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator void*() const volatile
    { return load(); }


};


template<>
struct atomic< char >
: atomic_char
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( char __v__ )
    CXX0X_AGGR_INIT( : atomic_char( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_char& operator =
    ( const atomic_char& )  CXX0X_DELETED
    )

    char operator =( char __v__ ) 
    { store( __v__ ); return __v__; }

    operator char() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_char& operator =
    ( const atomic_char& ) volatile CXX0X_DELETED
    )

    char operator =( char __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator char() const volatile
    { return load(); }


};


template<>
struct atomic< signed char >
: atomic_schar
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( signed char __v__ )
    CXX0X_AGGR_INIT( : atomic_schar( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_schar& operator =
    ( const atomic_schar& )  CXX0X_DELETED
    )

    signed char operator =( signed char __v__ ) 
    { store( __v__ ); return __v__; }

    operator signed char() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_schar& operator =
    ( const atomic_schar& ) volatile CXX0X_DELETED
    )

    signed char operator =( signed char __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator signed char() const volatile
    { return load(); }


};


template<>
struct atomic< unsigned char >
: atomic_uchar
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( unsigned char __v__ )
    CXX0X_AGGR_INIT( : atomic_uchar( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_uchar& operator =
    ( const atomic_uchar& )  CXX0X_DELETED
    )

    unsigned char operator =( unsigned char __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned char() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_uchar& operator =
    ( const atomic_uchar& ) volatile CXX0X_DELETED
    )

    unsigned char operator =( unsigned char __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned char() const volatile
    { return load(); }


};


template<>
struct atomic< short >
: atomic_short
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( short __v__ )
    CXX0X_AGGR_INIT( : atomic_short( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_short& operator =
    ( const atomic_short& )  CXX0X_DELETED
    )

    short operator =( short __v__ ) 
    { store( __v__ ); return __v__; }

    operator short() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_short& operator =
    ( const atomic_short& ) volatile CXX0X_DELETED
    )

    short operator =( short __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator short() const volatile
    { return load(); }


};


template<>
struct atomic< unsigned short >
: atomic_ushort
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( unsigned short __v__ )
    CXX0X_AGGR_INIT( : atomic_ushort( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_ushort& operator =
    ( const atomic_ushort& )  CXX0X_DELETED
    )

    unsigned short operator =( unsigned short __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned short() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_ushort& operator =
    ( const atomic_ushort& ) volatile CXX0X_DELETED
    )

    unsigned short operator =( unsigned short __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned short() const volatile
    { return load(); }


};


template<>
struct atomic< int >
: atomic_int
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( int __v__ )
    CXX0X_AGGR_INIT( : atomic_int( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_int& operator =
    ( const atomic_int& )  CXX0X_DELETED
    )

    int operator =( int __v__ ) 
    { store( __v__ ); return __v__; }

    operator int() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_int& operator =
    ( const atomic_int& ) volatile CXX0X_DELETED
    )

    int operator =( int __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator int() const volatile
    { return load(); }


};


template<>
struct atomic< unsigned int >
: atomic_uint
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( unsigned int __v__ )
    CXX0X_AGGR_INIT( : atomic_uint( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_uint& operator =
    ( const atomic_uint& )  CXX0X_DELETED
    )

    unsigned int operator =( unsigned int __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned int() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_uint& operator =
    ( const atomic_uint& ) volatile CXX0X_DELETED
    )

    unsigned int operator =( unsigned int __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned int() const volatile
    { return load(); }


};


template<>
struct atomic< long >
: atomic_long
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( long __v__ )
    CXX0X_AGGR_INIT( : atomic_long( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_long& operator =
    ( const atomic_long& )  CXX0X_DELETED
    )

    long operator =( long __v__ ) 
    { store( __v__ ); return __v__; }

    operator long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_long& operator =
    ( const atomic_long& ) volatile CXX0X_DELETED
    )

    long operator =( long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator long() const volatile
    { return load(); }


};


template<>
struct atomic< unsigned long >
: atomic_ulong
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( unsigned long __v__ )
    CXX0X_AGGR_INIT( : atomic_ulong( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_ulong& operator =
    ( const atomic_ulong& )  CXX0X_DELETED
    )

    unsigned long operator =( unsigned long __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_ulong& operator =
    ( const atomic_ulong& ) volatile CXX0X_DELETED
    )

    unsigned long operator =( unsigned long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned long() const volatile
    { return load(); }


};


template<>
struct atomic< long long >
: atomic_llong
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( long long __v__ )
    CXX0X_AGGR_INIT( : atomic_llong( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_llong& operator =
    ( const atomic_llong& )  CXX0X_DELETED
    )

    long long operator =( long long __v__ ) 
    { store( __v__ ); return __v__; }

    operator long long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_llong& operator =
    ( const atomic_llong& ) volatile CXX0X_DELETED
    )

    long long operator =( long long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator long long() const volatile
    { return load(); }


};


template<>
struct atomic< unsigned long long >
: atomic_ullong
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( unsigned long long __v__ )
    CXX0X_AGGR_INIT( : atomic_ullong( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_ullong& operator =
    ( const atomic_ullong& )  CXX0X_DELETED
    )

    unsigned long long operator =( unsigned long long __v__ ) 
    { store( __v__ ); return __v__; }

    operator unsigned long long() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_ullong& operator =
    ( const atomic_ullong& ) volatile CXX0X_DELETED
    )

    unsigned long long operator =( unsigned long long __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator unsigned long long() const volatile
    { return load(); }


};


template<>
struct atomic< wchar_t >
: atomic_wchar_t
{


    atomic() CXX0X_DEFAULTED_EASY
    CXX0X_CONSTEXPR_CTOR atomic( wchar_t __v__ )
    CXX0X_AGGR_INIT( : atomic_wchar_t( __v__ ) )
    { CXX0X_NO_AGGR_INIT( __f__ = __v__; ) }
    atomic( const atomic& ) CXX0X_DELETED


    CXX0X_AGGR_INIT(
    atomic_wchar_t& operator =
    ( const atomic_wchar_t& )  CXX0X_DELETED
    )

    wchar_t operator =( wchar_t __v__ ) 
    { store( __v__ ); return __v__; }

    operator wchar_t() const 
    { return load(); }


    CXX0X_AGGR_INIT(
    atomic_wchar_t& operator =
    ( const atomic_wchar_t& ) volatile CXX0X_DELETED
    )

    wchar_t operator =( wchar_t __v__ ) volatile
    { store( __v__ ); return __v__; }

    operator wchar_t() const volatile
    { return load(); }


};


#endif


#ifndef __cplusplus

#define atomic_is_lock_free( __a__ ) false

#define atomic_load( __a__ ) _ATOMIC_LOAD_( __a__, memory_order_seq_cst )

#define atomic_load_explicit( __a__, __x__ ) _ATOMIC_LOAD_( __a__, __x__ )

#define atomic_store( __a__, __m__ ) _ATOMIC_STORE_( __a__, __m__, memory_order_seq_cst )

#define atomic_store_explicit( __a__, __m__, __x__ ) _ATOMIC_STORE_( __a__, __m__, __x__ )

#define atomic_exchange( __a__, __m__ ) _ATOMIC_MODIFY_( __a__, =, __m__, memory_order_seq_cst )

#define atomic_exchange_explicit( __a__, __m__, __x__ ) _ATOMIC_MODIFY_( __a__, =, __m__, __x__ )

#define atomic_compare_exchange_strong( __a__, __e__, __m__ ) _ATOMIC_CMPSWP_( __a__, __e__, __m__, memory_order_seq_cst )

#define atomic_compare_exchange_strong_explicit( __a__, __e__, __m__, __x__, __y__ ) _ATOMIC_CMPSWP_( __a__, __e__, __m__, __x__ )

#define atomic_compare_exchange_weak( __a__, __e__, __m__ ) _ATOMIC_CMPSWP_( __a__, __e__, __m__, memory_order_seq_cst )

#define atomic_compare_exchange_weak_explicit( __a__, __e__, __m__, __x__, __y__ ) _ATOMIC_CMPSWP_( __a__, __e__, __m__, __x__ )


#define atomic_fetch_add( __a__, __m__ ) _ATOMIC_MODIFY_( __a__, +=, __m__, memory_order_seq_cst )

#define atomic_fetch_add_explicit( __a__, __m__, __x__ ) _ATOMIC_MODIFY_( __a__, +=, __m__, __x__ )


#define atomic_fetch_sub( __a__, __m__ ) _ATOMIC_MODIFY_( __a__, -=, __m__, memory_order_seq_cst )

#define atomic_fetch_sub_explicit( __a__, __m__, __x__ ) _ATOMIC_MODIFY_( __a__, -=, __m__, __x__ )


#define atomic_fetch_and( __a__, __m__ ) _ATOMIC_MODIFY_( __a__, &=, __m__, memory_order_seq_cst )

#define atomic_fetch_and_explicit( __a__, __m__, __x__ ) _ATOMIC_MODIFY_( __a__, &=, __m__, __x__ )


#define atomic_fetch_or( __a__, __m__ ) _ATOMIC_MODIFY_( __a__, |=, __m__, memory_order_seq_cst )

#define atomic_fetch_or_explicit( __a__, __m__, __x__ ) _ATOMIC_MODIFY_( __a__, |=, __m__, __x__ )


#define atomic_fetch_xor( __a__, __m__ ) _ATOMIC_MODIFY_( __a__, ^=, __m__, memory_order_seq_cst )

#define atomic_fetch_xor_explicit( __a__, __m__, __x__ ) _ATOMIC_MODIFY_( __a__, ^=, __m__, __x__ )


#endif


#ifdef __cplusplus
} // namespace std
#endif


#endif // CXX0X_ATOMIC_H

