//=================================================================================================
/*!
//  \file blaze/math/expressions/SparseVector.h
//  \brief Header file for the SparseVector base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SPARSEVECTOR_H_
#define _BLAZE_MATH_EXPRESSIONS_SPARSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Vector.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup sparse_vector Sparse Vectors
// \ingroup vector
*/
/*!\defgroup sparse_vector_expression Expressions
// \ingroup sparse_vector
*/
/*!\brief Base class for sparse vectors.
// \ingroup sparse_vector
//
// The SparseVector class is a base class for all arbitrarily sized (N-dimensional) sparse
// vectors. It provides an abstraction from the actual type of the sparse vector, but enables
// a conversion back to this type via the Vector base class.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
class SparseVector
   : public Vector<VT,TF>
{
 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   SparseVector() = default;
   SparseVector( const SparseVector& ) = default;
   SparseVector( SparseVector&& ) = default;
   ~SparseVector() = default;
   SparseVector& operator=( const SparseVector& ) = default;
   SparseVector& operator=( SparseVector&& ) = default;
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
/*!\name SparseVector global functions */
//@{
template< typename VT, bool TF >
typename VT::Iterator find( SparseVector<VT,TF>& sv, size_t index );

template< typename VT, bool TF >
typename VT::ConstIterator find( const SparseVector<VT,TF>& sv, size_t index );

template< typename VT, bool TF >
typename VT::Iterator lowerBound( SparseVector<VT,TF>& sv, size_t index );

template< typename VT, bool TF >
typename VT::ConstIterator lowerBound( const SparseVector<VT,TF>& sv, size_t index );

template< typename VT, bool TF >
typename VT::Iterator upperBound( SparseVector<VT,TF>& sv, size_t index );

template< typename VT, bool TF >
typename VT::ConstIterator upperBound( const SparseVector<VT,TF>& sv, size_t index );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific sparse vector element.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse vector (the end() iterator) is returned. Note that
// the returned sparse vector iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::Iterator
   find( SparseVector<VT,TF>& sv, size_t index )
{
   return (*sv).find( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific sparse vector element.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse vector (the end() iterator) is returned. Note that
// the returned sparse vector iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator
   find( const SparseVector<VT,TF>& sv, size_t index )
{
   return (*sv).find( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse vector iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::Iterator
   lowerBound( SparseVector<VT,TF>& sv, size_t index )
{
   return (*sv).lowerBound( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse vector iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator
   lowerBound( const SparseVector<VT,TF>& sv, size_t index )
{
   return (*sv).lowerBound( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse vector iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::Iterator
   upperBound( SparseVector<VT,TF>& sv, size_t index )
{
   return (*sv).upperBound( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse vector iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator
   upperBound( const SparseVector<VT,TF>& sv, size_t index )
{
   return (*sv).upperBound( index );
}
//*************************************************************************************************

} // namespace blaze

#endif
