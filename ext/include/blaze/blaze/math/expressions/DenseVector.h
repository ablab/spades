//=================================================================================================
/*!
//  \file blaze/math/expressions/DenseVector.h
//  \brief Header file for the DenseVector base class
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_EXPRESSIONS_DENSEVECTOR_H_
#define _BLAZE_MATH_EXPRESSIONS_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Vector.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/system/Inline.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_vector Dense Vectors
// \ingroup vector
*/
/*!\defgroup dense_vector_expression Expressions
// \ingroup dense_vector
*/
/*!\brief Base class for N-dimensional dense vectors.
// \ingroup dense_vector
//
// The DenseVector class is a base class for all arbitrarily sized (N-dimensional) dense
// vectors. It provides an abstraction from the actual type of the dense vector, but enables
// a conversion back to this type via the Vector base class.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DenseVector
   : public Vector<VT,TF>
{
 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   DenseVector() = default;
   DenseVector( const DenseVector& ) = default;
   DenseVector( DenseVector&& ) = default;
   ~DenseVector() = default;
   DenseVector& operator=( const DenseVector& ) = default;
   DenseVector& operator=( DenseVector&& ) = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseVector global functions */
//@{
template< typename VT, bool TF >
typename VT::ElementType* data( DenseVector<VT,TF>& dv ) noexcept;

template< typename VT, bool TF >
const typename VT::ElementType* data( const DenseVector<VT,TF>& dv ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for vectors without mutable data access.
// \ingroup dense_vector
//
// \param dv The given dense vector.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense vector without mutable data access,
// which is represented by \c nullptr.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto data_backend( DenseVector<VT,TF>& dv ) noexcept
   -> DisableIf_t< HasMutableDataAccess_v<VT>, typename VT::ElementType* >
{
   MAYBE_UNUSED( dv );

   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for vectors with mutable data access.
// \ingroup dense_vector
//
// \param dv The given dense vector.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense vector with mutable data access.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto data_backend( DenseVector<VT,TF>& dv ) noexcept
   -> EnableIf_t< HasMutableDataAccess_v<VT>, typename VT::ElementType* >
{
   return (*dv).data();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the dense vector elements.
// \ingroup dense_vector
//
// \param dv The given dense vector.
// \return Pointer to the internal element storage.
//
// This function provides a unified interface to access the given dense vector's internal
// element storage. In contrast to the \c data() member function, which is only available
// in case the vector has some internal storage, this function can be used on all kinds of
// dense vectors. In case the given dense vector does not provide low-level data access,
// the function returns \c nullptr.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ElementType* data( DenseVector<VT,TF>& dv ) noexcept
{
   return data_backend( *dv );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for vectors without constant data access.
// \ingroup dense_vector
//
// \param dv The given dense vector.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense vector without constant data access,
// which is represented by \c nullptr.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto data_backend( const DenseVector<VT,TF>& dv ) noexcept
   -> DisableIf_t< HasConstDataAccess_v<VT>, const typename VT::ElementType* >
{
   MAYBE_UNUSED( dv );

   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for vectors with constant data access.
// \ingroup dense_vector
//
// \param dv The given dense vector.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense vector with constant data access.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto data_backend( const DenseVector<VT,TF>& dv ) noexcept
   -> EnableIf_t< HasConstDataAccess_v<VT>, const typename VT::ElementType* >
{
   return (*dv).data();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the dense vector elements.
// \ingroup dense_vector
//
// \param dv The given dense vector.
// \return Pointer to the internal element storage.
//
// This function provides a unified interface to access the given dense vector's internal
// element storage. In contrast to the \c data() member function, which is only available
// in case the vector has some internal storage, this function can be used on all kinds of
// dense vectors. In case the given dense vector does not provide low-level data access,
// the function returns \c nullptr.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE const typename VT::ElementType* data( const DenseVector<VT,TF>& dv ) noexcept
{
   return data_backend( *dv );
}
//*************************************************************************************************

} // namespace blaze

#endif
