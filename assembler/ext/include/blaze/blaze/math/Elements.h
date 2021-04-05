//=================================================================================================
/*!
//  \file blaze/math/Elements.h
//  \brief Header file for the complete Elements implementation
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

#ifndef _BLAZE_MATH_ELEMENTS_H_
#define _BLAZE_MATH_ELEMENTS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Elements.h>
#include <blaze/math/Exception.h>
#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/views/Elements.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION FOR DENSE ELEMENT SELECTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for dense element selections.
// \ingroup random
//
// This specialization of the Rand class randomizes dense element selections.
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Rand< Elements<VT,TF,true,CEAs...> >
{
 public:
   //**********************************************************************************************
   /*!\brief Randomization of a dense element selection.
   //
   // \param elements The element selection to be randomized.
   // \return void
   */
   template< typename ET >  // Type of the element selection
   inline void randomize( ET&& elements ) const
   {
      using blaze::randomize;

      using ElementsType = RemoveReference_t<ET>;

      BLAZE_CONSTRAINT_MUST_BE_ELEMENTS_TYPE( ElementsType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ElementsType );

      for( size_t i=0UL; i<elements.size(); ++i ) {
         randomize( elements[i] );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a dense element selection.
   //
   // \param elements The element selection to be randomized.
   // \param min The smallest possible value for an element.
   // \param max The largest possible value for an element.
   // \return void
   */
   template< typename ET     // Type of the element selection
           , typename Arg >  // Min/max argument type
   inline void randomize( ET&& elements, const Arg& min, const Arg& max ) const
   {
      using blaze::randomize;

      using ElementsType = RemoveReference_t<ET>;

      BLAZE_CONSTRAINT_MUST_BE_ELEMENTS_TYPE( ElementsType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ElementsType );

      for( size_t i=0UL; i<elements.size(); ++i ) {
         randomize( elements[i], min, max );
      }
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION FOR SPARSE ELEMENT SELECTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for sparse element selections.
// \ingroup random
//
// This specialization of the Rand class randomizes sparse element selections.
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Rand< Elements<VT,TF,false,CEAs...> >
{
 public:
   //**********************************************************************************************
   /*!\brief Randomization of a sparse element selection.
   //
   // \param elements The element selection to be randomized.
   // \return void
   */
   template< typename ET >  // Type of the element selection
   inline void randomize( ET&& elements ) const
   {
      using ElementsType = RemoveReference_t<ET>;
      using ElementType  = ElementType_t<ElementsType>;

      BLAZE_CONSTRAINT_MUST_BE_ELEMENTS_TYPE( ElementsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ElementsType );

      const size_t size( elements.size() );

      if( size == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

      elements.reset();
      elements.reserve( nonzeros );

      while( elements.nonZeros() < nonzeros ) {
         elements[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse element selection.
   //
   // \param elements The element selection to be randomized.
   // \param nonzeros The number of non-zero elements of the random element selection.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename ET >  // Type of the element selection
   inline void randomize( ET&& elements, size_t nonzeros ) const
   {
      using ElementsType = RemoveReference_t<ET>;
      using ElementType  = ElementType_t<ElementsType>;

      BLAZE_CONSTRAINT_MUST_BE_ELEMENTS_TYPE( ElementsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ElementsType );

      const size_t size( elements.size() );

      if( nonzeros > size ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( size == 0UL ) return;

      elements.reset();
      elements.reserve( nonzeros );

      while( elements.nonZeros() < nonzeros ) {
         elements[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse element selection.
   //
   // \param elements The element selection to be randomized.
   // \param min The smallest possible value for an element.
   // \param max The largest possible value for an element.
   // \return void
   */
   template< typename ET     // Type of the element selection
           , typename Arg >  // Min/max argument type
   inline void randomize( ET&& elements, const Arg& min, const Arg& max ) const
   {
      using ElementsType = RemoveReference_t<ET>;
      using ElementType  = ElementType_t<ElementsType>;

      BLAZE_CONSTRAINT_MUST_BE_ELEMENTS_TYPE( ElementsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ElementsType );

      const size_t size( elements.size() );

      if( size == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

      elements.reset();
      elements.reserve( nonzeros );

      while( elements.nonZeros() < nonzeros ) {
         elements[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>( min, max );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse element selection.
   //
   // \param elements The element selection to be randomized.
   // \param nonzeros The number of non-zero elements of the random element selection.
   // \param min The smallest possible value for an element.
   // \param max The largest possible value for an element.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename ET     // Type of the element selection
           , typename Arg >  // Min/max argument type
   inline void randomize( ET&& elements, size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      using ElementsType = RemoveReference_t<ET>;
      using ElementType  = ElementType_t<ElementsType>;

      BLAZE_CONSTRAINT_MUST_BE_ELEMENTS_TYPE( ElementsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ElementsType );

      const size_t size( elements.size() );

      if( nonzeros > size ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( size == 0UL ) return;

      elements.reset();
      elements.reserve( nonzeros );

      while( elements.nonZeros() < nonzeros ) {
         elements[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>( min, max );
      }
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
