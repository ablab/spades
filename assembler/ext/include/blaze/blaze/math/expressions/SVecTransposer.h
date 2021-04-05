//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecTransposer.h
//  \brief Header file for the sparse vector transposer
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECTRANSPOSER_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECTRANSPOSER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SVECTRANSPOSER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the transposition of a sparse vector.
// \ingroup sparse_vector_expression
//
// The SVecTransposer class is a wrapper object for the temporary transposition of a sparse vector.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
class SVecTransposer
   : public SparseVector< SVecTransposer<VT,TF>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This           = SVecTransposer<VT,TF>;  //!< Type of this SVecTransposer instance.
   using BaseType       = SparseVector<This,TF>;  //!< Base type of this SVecTransposer instance.
   using ResultType     = TransposeType_t<VT>;    //!< Result type for expression template evaluations.
   using TransposeType  = ResultType_t<VT>;       //!< Transpose type for expression template evaluations.
   using ElementType    = ElementType_t<VT>;      //!< Resulting element type.
   using ReturnType     = ReturnType_t<VT>;       //!< Return type for expression template evaluations.
   using CompositeType  = const This&;            //!< Data type for composite expression templates.
   using Reference      = Reference_t<VT>;        //!< Reference to a non-constant matrix value.
   using ConstReference = ConstReference_t<VT>;   //!< Reference to a constant matrix value.
   using Iterator       = Iterator_t<VT>;         //!< Iterator over non-constant elements.
   using ConstIterator  = ConstIterator_t<VT>;    //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = VT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SVecTransposer class.
   //
   // \param sv The sparse vector operand.
   */
   explicit inline SVecTransposer( VT& sv ) noexcept
      : sv_( sv )  // The sparse vector operand
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline ConstReference operator[]( size_t index ) const {
      BLAZE_USER_ASSERT( index < sv_.size(), "Invalid vector access index" );
      return const_cast<const VT&>( sv_ )[index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ConstReference at( size_t index ) const {
      if( index >= sv_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the sparse vector.
   //
   // \return Iterator to the first non-zero element of the sparse vector.
   */
   inline Iterator begin() {
      return sv_.begin();
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the sparse vector.
   //
   // \return Iterator to the first non-zero element of the sparse vector.
   */
   inline ConstIterator begin() const {
      return sv_.cbegin();
   }
   //**********************************************************************************************

   //**Cbegin function*****************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the sparse vector.
   //
   // \return Iterator to the first non-zero element of the sparse vector.
   */
   inline ConstIterator cbegin() const {
      return sv_.cbegin();
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the sparse vector.
   //
   // \return Iterator just past the last non-zero element of the sparse vector.
   */
   inline Iterator end() {
      return sv_.end();
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the sparse vector.
   //
   // \return Iterator just past the last non-zero element of the sparse vector.
   */
   inline ConstIterator end() const {
      return sv_.cend();
   }
   //**********************************************************************************************

   //**Cend function*******************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the sparse vector.
   //
   // \return Iterator just past the last non-zero element of the sparse vector.
   */
   inline ConstIterator cend() const {
      return sv_.cend();
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return sv_.size();
   }
   //**********************************************************************************************

   //**Capacity function***************************************************************************
   /*!\brief Returns the maximum capacity of the vector.
   //
   // \return The capacity of the vector.
   */
   inline size_t capacity() const noexcept {
      return sv_.capacity();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the vector.
   //
   // \return The number of non-zero elements in the vector.
   */
   inline size_t nonZeros() const {
      return sv_.nonZeros();
   }
   //**********************************************************************************************

   //**Reset function******************************************************************************
   /*!\brief Resets the vector elements.
   //
   // \return void
   */
   inline void reset() {
      return sv_.reset();
   }
   //**********************************************************************************************

   //**Insert function*****************************************************************************
   /*!\brief Inserting an element into the sparse vector.
   //
   // \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
   // \param value The value of the element to be inserted.
   // \return Iterator to the inserted element.
   // \exception std::invalid_argument Invalid sparse vector access index.
   //
   // This function inserts a new element into the sparse vector. However, duplicate elements are
   // not allowed. In case the sparse matrix already contains an element with row index \a index,
   // a \a std::invalid_argument exception is thrown.
   */
   inline Iterator insert( size_t index, const ElementType& value ) {
      return sv_.insert( index, value );
   }
   //**********************************************************************************************

   //**Find function*******************************************************************************
   /*!\brief Inserting an element into the sparse vector.
   //
   // \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
   // \return Iterator to the element in case the index is found, end() iterator otherwise.
   //
   // This function can be used to check whether a specific element is contained in the sparse
   // vector. It specifically searches for the element with index \a index. In case the element
   // is found, the function returns an iterator to the element. Otherwise an iterator just past
   // the last non-zero element of the sparse vector (the end() iterator) is returned. Note that
   // the returned sparse vector iterator is subject to invalidation due to inserting operations
   // via the subscript operator or the insert() function!
   */
   inline Iterator find( size_t index ) {
      return sv_.find( index );
   }
   //**********************************************************************************************

   //**Reserve function****************************************************************************
   /*!\brief Setting the minimum capacity of the sparse vector.
   //
   // \param nonzeros The new minimum capacity of the sparse vector.
   // \return void
   //
   // This function increases the capacity of the sparse vector to at least \a nonzeros elements.
   // The current values of the vector elements are preserved.
   */
   inline void reserve( size_t nonzeros ) {
      sv_.reserve( nonzeros );
   }
   //**********************************************************************************************

   //**Append function*****************************************************************************
   /*!\brief Appending an element to the sparse vector.
   //
   // \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
   // \param value The value of the element to be appended.
   // \param check \a true if the new value should be checked for default values, \a false if not.
   // \return void
   //
   // This function provides a very efficient way to fill a sparse vector with elements. It
   // appends a new element to the end of the sparse vector without any additional memory
   // allocation. Therefore it is strictly necessary to keep the following preconditions in
   // mind:
   //
   //  - the index of the new element must be strictly larger than the largest index of non-zero
   //    elements in the sparse vector
   //  - the current number of non-zero elements must be smaller than the capacity of the vector
   //
   // Ignoring these preconditions might result in undefined behavior! The optional \a check
   // parameter specifies whether the new value should be tested for a default value. If the new
   // value is a default value (for instance 0 in case of an integral element type) the value is
   // not appended. Per default the values are not tested.
   //
   // \note Although append() does not allocate new memory, it still invalidates all iterators
   // returned by the end() functions!
   */
   inline void append( size_t index, const ElementType& value, bool check=false ) {
      sv_.append( index, value, check );
   }
   //**********************************************************************************************

   //**IsIntact function***************************************************************************
   /*!\brief Returns whether the invariants of the vector are intact.
   //
   // \return \a true in case the vector's invariants are intact, \a false otherwise.
   */
   inline bool isIntact() const noexcept
   {
      using blaze::isIntact;
      return isIntact( sv_ );
   }
   //**********************************************************************************************

   //**CanAlias function***************************************************************************
   /*!\brief Returns whether the vector can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this vector, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool canAlias( const Other* alias ) const noexcept
   {
      return sv_.canAlias( alias );
   }
   //**********************************************************************************************

   //**IsAliased function**************************************************************************
   /*!\brief Returns whether the vector is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this vector, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool isAliased( const Other* alias ) const noexcept
   {
      return sv_.isAliased( alias );
   }
   //**********************************************************************************************

   //**CanSMPAssign function***********************************************************************
   /*!\brief Returns whether the vector can be used in SMP assignments.
   //
   // \return \a true in case the vector can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept
   {
      return sv_.canSMPAssign();
   }
   //**********************************************************************************************

   //**Transpose assignment of vectors*************************************************************
   /*!\brief Implementation of the transpose assignment of a vector.
   //
   // \param rhs The right-hand side vector to be assigned.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side vector
   inline void assign( const Vector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      sv_.assign( trans( *rhs ) );
   }
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*!\brief Calculating a new vector capacity.
   //
   // \return The new compressed vector capacity.
   //
   // This function calculates a new vector capacity based on the current capacity of the sparse
   // vector. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
   */
   inline size_t extendCapacity() const noexcept
   {
      using blaze::max;
      using blaze::min;

      size_t nonzeros( 2UL*sv_.capacity()+1UL );
      nonzeros = max( nonzeros, 7UL );
      nonzeros = min( nonzeros, sv_.size() );

      BLAZE_INTERNAL_ASSERT( nonzeros > sv_.capacity(), "Invalid capacity value" );

      return nonzeros;
   }
   //**********************************************************************************************

   //**Member variables****************************************************************************
   VT& sv_;  //!< The sparse vector operand.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, !TF );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the sparse vector contained in a SVecTransposer.
// \ingroup sparse_vector_expression
//
// \param v The sparse vector to be resetted.
// \return void
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline void reset( SVecTransposer<VT,TF>& v )
{
   v.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given SVecTransposer are intact.
// \ingroup sparse_vector_expression
//
// \param v The sparse vector to be tested.
// \return \a true in caes the given vector's invariants are intact, \a false otherwise.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline bool isIntact( const SVecTransposer<VT,TF>& v ) noexcept
{
   return v.isIntact();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct Size< SVecTransposer<VT,TF>, 0UL >
   : public Size<VT,0UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAXSIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct MaxSize< SVecTransposer<VT,TF>, 0UL >
   : public MaxSize<VT,0UL>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
