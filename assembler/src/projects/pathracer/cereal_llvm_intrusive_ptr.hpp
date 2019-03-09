#pragma once

#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <llvm/ADT/IntrusiveRefCntPtr.h>

//FIXME --- remove it
// #include "utils/logger/logger.hpp"

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
  CEREAL_LOAD_FUNCTION_NAME( Archive & ar, memory_detail::PtrWrapper<llvm::IntrusiveRefCntPtr<T> &> & wrapper );
  // It makes no sense for intrusive_ptr. We have to have valid pointer BEFORE the actual deserialization (since the deserialized
  // object could contain pointers to itself). But it's impossible (or, at least, hardly possible) to consturuct intrusive_ptr referring to
  // a raw memory block

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
      ptr = detail::Construct<T, Archive>::load_andor_construct();
      auto sptr = std::make_shared<llvm::IntrusiveRefCntPtr<T>>(ptr);
      // INFO("Cached pointer: " << sptr.get() << "id " << id);
      ar.registerSharedPointer( id, sptr );
      ar( CEREAL_NVP_("data", *ptr) );
    }
    else
    {
      // INFO("Loading existing pointer: " << ar.getSharedPointer(id).get() << " id " << id);
      llvm::IntrusiveRefCntPtr<T> *pptr = static_cast<llvm::IntrusiveRefCntPtr<T>*>(ar.getSharedPointer(id).get());
      ptr = pptr ? *pptr : nullptr;
    }
  }

} // namespace cereal

// TODO Add polymorphic intrusive ptr support

#undef CEREAL_ALIGNOF
