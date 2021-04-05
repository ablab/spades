//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecTransposer.h
//  \brief Header file for the dense vector transposer
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECTRANSPOSER_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECTRANSPOSER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECTRANSPOSER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the transposition of a dense vector.
// \ingroup dense_vector_expression
//
// The DVecTransposer class is a wrapper object for the temporary transposition of a dense vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DVecTransposer
   : public DenseVector< DVecTransposer<VT,TF>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This           = DVecTransposer<VT,TF>;     //!< Type of this DVecTransposer instance.
   using BaseType       = DenseVector<This,TF>;      //!< Base type of this DVecTransposer instance.
   using ResultType     = TransposeType_t<VT>;       //!< Result type for expression template evaluations.
   using TransposeType  = ResultType_t<VT>;          //!< Transpose type for expression template evaluations.
   using ElementType    = ElementType_t<VT>;         //!< Type of the vector elements.
   using SIMDType       = SIMDTrait_t<ElementType>;  //!< SIMD type of the vector elements.
   using ReturnType     = ReturnType_t<VT>;          //!< Return type for expression template evaluations.
   using CompositeType  = const This&;               //!< Data type for composite expression templates.
   using Reference      = Reference_t<VT>;           //!< Reference to a non-constant vector value.
   using ConstReference = ConstReference_t<VT>;      //!< Reference to a constant vector value.
   using Pointer        = Pointer_t<VT>;             //!< Pointer to a non-constant vector value.
   using ConstPointer   = ConstPointer_t<VT>;        //!< Pointer to a constant vector value.
   using Iterator       = Iterator_t<VT>;            //!< Iterator over non-constant elements.
   using ConstIterator  = ConstIterator_t<VT>;       //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the vector is involved
       in can be optimized via SIMD operations. In case the dense vector operand is vectorizable,
       the \a simdEnabled compilation flag is set to \a true, otherwise it is set to \a false. */
   static constexpr bool simdEnabled = VT::simdEnabled;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = VT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecTransposer class.
   //
   // \param dv The dense vector operand.
   */
   explicit inline DVecTransposer( VT& dv ) noexcept
      : dv_( dv )  // The dense vector operand
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline Reference operator[]( size_t index ) {
      BLAZE_USER_ASSERT( index < dv_.size(), "Invalid vector access index" );
      return dv_[index];
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline ConstReference operator[]( size_t index ) const {
      BLAZE_USER_ASSERT( index < dv_.size(), "Invalid vector access index" );
      return const_cast<const VT&>( dv_ )[index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline Reference at( size_t index ) {
      if( index >= dv_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
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
      if( index >= dv_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the vector elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline Pointer data() noexcept {
      return dv_.data();
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the vector elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline ConstPointer data() const noexcept {
      return dv_.data();
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first element of the dense vector.
   //
   // \return Iterator to the first element of the dense vector.
   */
   inline Iterator begin() {
      return dv_.begin();
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first element of the dense vector.
   //
   // \return Iterator to the first element of the dense vector.
   */
   inline ConstIterator begin() const {
      return dv_.cbegin();
   }
   //**********************************************************************************************

   //**Cbegin function*****************************************************************************
   /*!\brief Returns an iterator to the first element of the dense vector.
   //
   // \return Iterator to the first element of the dense vector.
   */
   inline ConstIterator cbegin() const {
      return dv_.cbegin();
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last element of the dense vector.
   //
   // \return Iterator just past the last element of the dense vector.
   */
   inline Iterator end() {
      return dv_.end();
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last element of the dense vector.
   //
   // \return Iterator just past the last element of the dense vector.
   */
   inline ConstIterator end() const {
      return dv_.cend();
   }
   //**********************************************************************************************

   //**Cend function*******************************************************************************
   /*!\brief Returns an iterator just past the last element of the dense vector.
   //
   // \return Iterator just past the last element of the dense vector.
   */
   inline ConstIterator cend() const {
      return dv_.cend();
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return dv_.size();
   }
   //**********************************************************************************************

   //**Reset function******************************************************************************
   /*!\brief Resets the vector elements.
   //
   // \return void
   */
   inline void reset() {
      return dv_.reset();
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
      return isIntact( dv_ );
   }
   //**********************************************************************************************

   //**CanAliased function*************************************************************************
   /*!\brief Returns whether the vector can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this vector, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool canAlias( const Other* alias ) const noexcept
   {
      return dv_.canAlias( alias );
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
      return dv_.isAliased( alias );
   }
   //**********************************************************************************************

   //**IsAligned function**************************************************************************
   /*!\brief Returns whether the vector is properly aligned in memory.
   //
   // \return \a true in case the vector is aligned, \a false if not.
   */
   inline bool isAligned() const noexcept
   {
      return dv_.isAligned();
   }
   //**********************************************************************************************

   //**CanSMPAssign function***********************************************************************
   /*!\brief Returns whether the vector can be used in SMP assignments.
   //
   // \return \a true in case the vector can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept
   {
      return dv_.canSMPAssign();
   }
   //**********************************************************************************************

   //**Load function*******************************************************************************
   /*!\brief Load of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \return The loaded SIMD element.
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE SIMDType load( size_t index ) const noexcept
   {
      return dv_.load( index );
   }
   //**********************************************************************************************

   //**Loada function******************************************************************************
   /*!\brief Aligned load of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \return The loaded SIMD element.
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept
   {
      return dv_.loada( index );
   }
   //**********************************************************************************************

   //**Loadu function******************************************************************************
   /*!\brief Unaligned load of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \return The loaded SIMD element.
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept
   {
      return dv_.loadu( index );
   }
   //**********************************************************************************************

   //**Store function******************************************************************************
   /*!\brief Store of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void store( size_t index, const SIMDType& value ) noexcept
   {
      dv_.store( index, value );
   }
   //**********************************************************************************************

   //**Storea function******************************************************************************
   /*!\brief Aligned store of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept
   {
      dv_.storea( index, value );
   }
   //**********************************************************************************************

   //**Storeu function*****************************************************************************
   /*!\brief Unaligned store of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept
   {
      dv_.storeu( index, value );
   }
   //**********************************************************************************************

   //**Stream function*****************************************************************************
   /*!\brief Aligned, non-temporal store of a SIMD element of the vector.
   //
   // \param index Access index. The index must be smaller than the number of vector elements.
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept
   {
      dv_.stream( index, value );
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

      dv_.assign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose addition assignment of vectors****************************************************
   /*!\brief Implementation of the transpose addition assignment of a vector.
   //
   // \param rhs The right-hand side vector to be added.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side vector
   inline void addAssign( const Vector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      dv_.addAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose subtraction assignment of vectors*************************************************
   /*!\brief Implementation of the transpose subtraction assignment of a vector.
   //
   // \param rhs The right-hand side vector to be subtracted.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side vector
   inline void subAssign( const Vector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      dv_.subAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose multiplication assignment of vectors**********************************************
   /*!\brief Implementation of the transpose multiplication assignment of a vector.
   //
   // \param rhs The right-hand side vector to be multiplied.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side vector
   inline void multAssign( const Vector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      dv_.multAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose division assignment of dense vectors**********************************************
   /*!\brief Implementation of the transpose division assignment of a dense vector.
   //
   // \param rhs The right-hand side dense vector divisor.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side dense vector
   inline void divAssign( const DenseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      dv_.divAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   VT& dv_;  //!< The dense vector operand.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
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
/*!\brief Resetting the dense vector contained in a DVecTransposer.
// \ingroup dense_vector_expression
//
// \param v The dense vector to be resetted.
// \return void
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void reset( DVecTransposer<VT,TF>& v )
{
   v.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given DVecTransposer are intact.
// \ingroup dense_vector_expression
//
// \param v The dense vector to be tested.
// \return \a true in caes the given vector's invariants are intact, \a false otherwise.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline bool isIntact( const DVecTransposer<VT,TF>& v ) noexcept
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
struct Size< DVecTransposer<VT,TF>, 0UL >
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
struct MaxSize< DVecTransposer<VT,TF>, 0UL >
   : public MaxSize<VT,0UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct HasConstDataAccess< DVecTransposer<VT,TF> >
   : public HasConstDataAccess<VT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct HasMutableDataAccess< DVecTransposer<VT,TF> >
   : public HasMutableDataAccess<VT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct IsAligned< DVecTransposer<VT,TF> >
   : public IsAligned<VT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGUOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct IsContiguous< DVecTransposer<VT,TF> >
   : public IsContiguous<VT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct IsPadded< DVecTransposer<VT,TF> >
   : public IsPadded<VT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
