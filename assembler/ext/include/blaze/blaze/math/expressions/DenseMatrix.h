//=================================================================================================
/*!
//  \file blaze/math/expressions/DenseMatrix.h
//  \brief Header file for the DenseMatrix base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DENSEMATRIX_H_
#define _BLAZE_MATH_EXPRESSIONS_DENSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/system/Inline.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_matrix Dense Matrices
// \ingroup matrix
*/
/*!\defgroup dense_matrix_expression Expressions
// \ingroup dense_matrix
*/
/*!\brief Base class for dense matrices.
// \ingroup dense_matrix
//
// The DenseMatrix class is a base class for all dense matrix classes. It provides an
// abstraction from the actual type of the dense matrix, but enables a conversion back
// to this type via the Matrix base class.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
class DenseMatrix
   : public Matrix<MT,SO>
{
 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   DenseMatrix() = default;
   DenseMatrix( const DenseMatrix& ) = default;
   DenseMatrix( DenseMatrix&& ) = default;
   ~DenseMatrix() = default;
   DenseMatrix& operator=( const DenseMatrix& ) = default;
   DenseMatrix& operator=( DenseMatrix&& ) = default;
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
/*!\name DenseMatrix global functions */
//@{
template< typename MT, bool SO >
typename MT::ElementType* data( DenseMatrix<MT,SO>& dm ) noexcept;

template< typename MT, bool SO >
const typename MT::ElementType* data( const DenseMatrix<MT,SO>& dm ) noexcept;

template< typename MT, bool SO >
size_t spacing( const DenseMatrix<MT,SO>& dm ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for matrices without mutable data access.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense matrix without mutable data access,
// which is represented by \c nullptr.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE auto data_backend( DenseMatrix<MT,SO>& dm ) noexcept
   -> DisableIf_t< HasMutableDataAccess_v<MT>, typename MT::ElementType* >
{
   MAYBE_UNUSED( dm );

   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for matrices with mutable data access.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense matrix with mutable data access.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE auto data_backend( DenseMatrix<MT,SO>& dm ) noexcept
   -> EnableIf_t< HasMutableDataAccess_v<MT>, typename MT::ElementType* >
{
   return (*dm).data();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the dense matrix elements.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return Pointer to the internal element storage.
//
// This function provides a unified interface to access the given dense matrix's internal
// element storage. In contrast to the \c data() member function, which is only available
// in case the matrix has some internal storage, this function can be used on all kinds of
// dense matrices. In case the given dense matrix does not provide low-level data access,
// the function returns \c nullptr.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE typename MT::ElementType* data( DenseMatrix<MT,SO>& dm ) noexcept
{
   return data_backend( *dm );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for matrices without constant data access.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense matrix without constant data access,
// which is represented by \c nullptr.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE auto data_backend( const DenseMatrix<MT,SO>& dm ) noexcept
   -> DisableIf_t< HasConstDataAccess_v<MT>, const typename MT::ElementType* >
{
   MAYBE_UNUSED( dm );

   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c data() function for matrices with constant data access.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return Pointer to the internal element storage.
//
// This function returns the internal storage of a dense matrix with constant data access.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE auto data_backend( const DenseMatrix<MT,SO>& dm ) noexcept
   -> EnableIf_t< HasConstDataAccess_v<MT>, const typename MT::ElementType* >
{
   return (*dm).data();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the dense matrix elements.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return Pointer to the internal element storage.
//
// This function provides a unified interface to access the given dense matrix's internal
// element storage. In contrast to the \c data() member function, which is only available
// in case the matrix has some internal storage, this function can be used on all kinds of
// dense matrices. In case the given dense matrix does not provide low-level data access,
// the function returns \c nullptr.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE const typename MT::ElementType* data( const DenseMatrix<MT,SO>& dm ) noexcept
{
   return data_backend( *dm );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows/columns.
// \ingroup dense_matrix
//
// \param dm The given matrix.
// \return The spacing between the beginning of two rows/columns.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
BLAZE_ALWAYS_INLINE size_t spacing( const DenseMatrix<MT,SO>& dm ) noexcept
{
   return (*dm).spacing();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for row-major dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given row-major dense
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetLower_backend( DenseMatrix<MT,false>& dm )
   -> DisableIf_t< IsUniform_v<MT> || IsUpper_v<MT> >
{
   using blaze::reset;

   const size_t m( (*dm).rows()    );
   const size_t n( (*dm).columns() );

   for( size_t i=1UL; i<m; ++i ) {
      const size_t jend( min( i, n ) );
      for( size_t j=0UL; j<jend; ++j ) {
         reset( (*dm)(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for column-major dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given column-major dense
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetLower_backend( DenseMatrix<MT,true>& dm )
   -> DisableIf_t< IsUniform_v<MT> || IsUpper_v<MT> >
{
   using blaze::reset;

   const size_t m   ( (*dm).rows()    );
   const size_t n   ( (*dm).columns() );
   const size_t jend( min( m, n ) );

   for( size_t j=0UL; j<jend; ++j ) {
      for( size_t i=j+1UL; i<m; ++i ) {
         reset( (*dm)(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for uniform dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given uniform dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline auto resetLower_backend( DenseMatrix<MT,SO>& dm )
   -> EnableIf_t< IsUniform_v<MT> && !IsUpper_v<MT> >
{
   using blaze::reset;

   reset( *dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for upper dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given upper dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline auto resetLower_backend( DenseMatrix<MT,SO>& dm )
   -> EnableIf_t< IsUpper_v<MT> >
{
   MAYBE_UNUSED( dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the lower part of the given dense matrix.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline void resetLower( DenseMatrix<MT,SO>& dm )
{
   resetLower_backend( *dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for row-major dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given row-major dense
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetUpper_backend( DenseMatrix<MT,false>& dm )
   -> DisableIf_t< IsUniform_v<MT> || IsLower_v<MT> >
{
   using blaze::reset;

   const size_t m   ( (*dm).rows()    );
   const size_t n   ( (*dm).columns() );
   const size_t iend( min( m, n ) );

   for( size_t i=0UL; i<iend; ++i ) {
      for( size_t j=i+1UL; j<n; ++j ) {
         reset( (*dm)(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for column-major dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given column-major dense
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetUpper_backend( DenseMatrix<MT,true>& dm )
   -> DisableIf_t< IsUniform_v<MT> || IsLower_v<MT> >
{
   using blaze::reset;

   const size_t m( (*dm).rows()    );
   const size_t n( (*dm).columns() );

   for( size_t j=1UL; j<n; ++j ) {
      const size_t iend( min( j, m ) );
      for( size_t i=0UL; i<iend; ++i ) {
         reset( (*dm)(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for uniform dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given uniform dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline auto resetUpper_backend( DenseMatrix<MT,SO>& dm )
   -> EnableIf_t< IsUniform_v<MT> || !IsLower_v<MT> >
{
   using blaze::reset;

   reset( *dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for lower dense matrices.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given lower dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline auto resetUpper_backend( DenseMatrix<MT,SO>& dm )
   -> EnableIf_t< IsLower_v<MT> >
{
   MAYBE_UNUSED( dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the upper part of the given dense matrix.
// \ingroup dense_matrix
//
// \param matrix The given dense matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline void resetUpper( DenseMatrix<MT,SO>& dm )
{
   resetUpper_backend( *dm );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
