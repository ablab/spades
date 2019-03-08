#pragma once

#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <llvm/ADT/IntrusiveRefCntPtr.h>

// Work around MSVC not having alignof
#if defined(_MSC_VER) && _MSC_VER < 1900
#define CEREAL_ALIGNOF __alignof
#else // not MSVC 2013 or older
#define CEREAL_ALIGNOF alignof
#endif // end MSVC check

namespace cereal
{
  //! Saving llvm::IntrusiveRefCntPtr for non polymorphic types
  template <class Archive, class T> inline
  typename std::enable_if<!std::is_polymorphic<T>::value, void>::type
  CEREAL_SAVE_FUNCTION_NAME( Archive & ar, llvm::IntrusiveRefCntPtr<T> const & ptr )
  {
    ar( CEREAL_NVP_("ptr_wrapper", memory_detail::make_ptr_wrapper( ptr )) );
  }

  //! Loading llvm::IntrusiveRefCntPtr, case when no user load and construct for non polymorphic types
  template <class Archive, class T> inline
  typename std::enable_if<!std::is_polymorphic<T>::value, void>::type
  CEREAL_LOAD_FUNCTION_NAME( Archive & ar, llvm::IntrusiveRefCntPtr<T> & ptr )
  {
    ar( CEREAL_NVP_("ptr_wrapper", memory_detail::make_ptr_wrapper( ptr )) );
  }

  // ######################################################################
  // Pointer wrapper implementations follow below

  //! Saving llvm::IntrusiveRefCntPtr (wrapper implementation)
  /*! @internal */
  template <class Archive, class T> inline
  void CEREAL_SAVE_FUNCTION_NAME( Archive & ar, memory_detail::PtrWrapper<llvm::IntrusiveRefCntPtr<T> const &> const & wrapper )
  {
    auto & ptr = wrapper.ptr;

    uint32_t id = ar.registerSharedPointer( ptr.get() );
    ar( CEREAL_NVP_("id", id) );

    if( id & detail::msb_32bit )
    {
      ar( CEREAL_NVP_("data", *ptr) );
    }
  }

  //! Loading llvm::IntrusiveRefCntPtr, case when user load and construct (wrapper implementation)
  /*! @internal */
  template <class Archive, class T> inline
  typename std::enable_if<traits::has_load_and_construct<T, Archive>::value, void>::type
  CEREAL_LOAD_FUNCTION_NAME( Archive & ar, memory_detail::PtrWrapper<llvm::IntrusiveRefCntPtr<T> &> & wrapper )
  {
    auto & ptr = wrapper.ptr;

    uint32_t id;

    ar( CEREAL_NVP_("id", id) );

    if( id & detail::msb_32bit )
    {
      // Storage type for the pointer - since we can't default construct this type,
      // we'll allocate it using std::aligned_storage and use a custom deleter
      using ST = typename std::aligned_storage<sizeof(T), CEREAL_ALIGNOF(T)>::type;

      // Valid flag - set to true once construction finishes
      //  This prevents us from calling the destructor on
      //  uninitialized data.
      auto valid = std::make_shared<bool>( false );

      // Allocate our storage, which we will treat as
      //  uninitialized until initialized with placement new
      ptr.reset( reinterpret_cast<T *>( new ST() ),
          [=]( T * t )
          {
            if( *valid )
              t->~T();

            delete reinterpret_cast<ST *>( t );
          } );

      // Register the pointer
      ar.registerSharedPointer( id, ptr );

      // Perform the actual loading and allocation
      memory_detail::loadAndConstructSharedPtr( ar, ptr.get(), typename ::cereal::traits::has_shared_from_this<T>::type() );

      // Mark pointer as valid (initialized)
      *valid = true;
    }
    else
      ptr = std::static_pointer_cast<T>(ar.getSharedPointer(id));
  }

  //! Loading llvm::IntrusiveRefCntPtr, case when no user load and construct (wrapper implementation)
  /*! @internal */
  template <class Archive, class T> inline
  typename std::enable_if<!traits::has_load_and_construct<T, Archive>::value, void>::type
  CEREAL_LOAD_FUNCTION_NAME( Archive & ar, memory_detail::PtrWrapper<llvm::IntrusiveRefCntPtr<T> &> & wrapper )
  {
    auto & ptr = wrapper.ptr;

    uint32_t id;

    ar( CEREAL_NVP_("id", id) );

    if( id & detail::msb_32bit )
    {
      ptr.reset( detail::Construct<T, Archive>::load_andor_construct() );
      ar.registerSharedPointer( id, ptr );
      ar( CEREAL_NVP_("data", *ptr) );
    }
    else
      ptr = std::static_pointer_cast<T>(ar.getSharedPointer(id));
  }

} // namespace cereal

// TODO Add polymorphic intrusive ptr support

#undef CEREAL_ALIGNOF
