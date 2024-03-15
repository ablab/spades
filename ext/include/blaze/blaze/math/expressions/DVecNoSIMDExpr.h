//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecNoSIMDExpr.h
//  \brief Header file for the dense vector no-SIMD expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECNOSIMDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECNOSIMDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/Modification.h>
#include <blaze/math/expressions/VecNoSIMDExpr.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/GetMemberType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECNOSIMDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the non-SIMD evaluation of dense vectors.
// \ingroup dense_vector_expression
//
// The DVecNoSIMDExpr class represents the compile time expression for the non-SIMD evaluation
// of a dense vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DVecNoSIMDExpr
   : public VecNoSIMDExpr< DenseVector< DVecNoSIMDExpr<VT,TF>, TF > >
   , private Modification<VT>
{
 private:
   //**Type definitions****************************************************************************
   //! Definition of the GetConstIterator type trait.
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetConstIterator, ConstIterator, INVALID_TYPE );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the no-SIMD expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the no-SIMD expression. In case the given dense vector
       expression of type \a VT requires an intermediate evaluation, \a useAssign will be set
       to 1 and the no-SIMD expression will be evaluated via the \a assign function family.
       Otherwise \a useAssign will be set to 0 and the expression will be evaluated via the
       subscript operator. */
   static constexpr bool useAssign = RequiresEvaluation_v<VT>;

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the target vector is SMP assignable and the dense vector operand requires an
       intermediate evaluation, the variable is set to 1 and the expression specific evaluation
       strategy is selected. Otherwise the variable is set to 0 and the default strategy is
       chosen. */
   template< typename VT2 >
   static constexpr bool UseSMPAssign_v = ( VT2::smpAssignable && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecNoSIMDExpr instance.
   using This = DVecNoSIMDExpr<VT,TF>;

   //! Base type of this DVecNoSIMDExpr instance.
   using BaseType = VecNoSIMDExpr< DenseVector<This,TF> >;

   using ResultType    = ResultType_t<VT>;     //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<VT>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<VT>;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DVecNoSIMDExpr& >;

   //! Iterator over the elements of the dense vector.
   using ConstIterator = GetConstIterator_t<VT>;

   //! Composite data type of the dense vector expression.
   using Operand = If_t< IsExpression_v<VT>, const VT, const VT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = VT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecNoSIMDExpr class.
   //
   // \param dv The dense vector operand of the no-SIMD expression.
   */
   explicit inline DVecNoSIMDExpr( const VT& dv ) noexcept
      : dv_( dv )  // Dense vector of the no-SIMD expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      return dv_[index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      return dv_.at( index );
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the vector elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline const ElementType* data() const noexcept {
      return dv_.data();
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( dv_.begin() );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( dv_.end() );
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

   //**Operand access******************************************************************************
   /*!\brief Returns the dense vector operand.
   //
   // \return The dense vector operand.
   */
   inline Operand operand() const noexcept {
      return dv_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return dv_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return dv_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return dv_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return dv_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand dv_;  //!< Dense vector of the no-SIMD expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector no-SIMD
   // expression to a dense vector. Due to the explicit application of the SFINAE principle, this
   // function can only be selected by the compiler in case the operand requires an intermediate
   // evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector no-SIMD expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side no-SIMD expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector no-SIMD
   // expression to a sparse vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto assign( SparseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector
   // no-SIMD expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto addAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector no-SIMD expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto subAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector no-SIMD expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto multAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      multAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector
   // no-SIMD expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto divAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      divAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector
   // no-SIMD expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector no-SIMD expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side no-SIMD expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector
   // no-SIMD expression to a sparse vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector no-SIMD expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAddAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a
   // dense vector no-SIMD expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpSubAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // dense vector no-SIMD expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpMultAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP division assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a dense vector no-SIMD expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side no-SIMD expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a dense
   // vector no-SIMD expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT2,TF>& lhs, const DVecNoSIMDExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpDivAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to sparse vectors***************************************************
   // No special implementation for the SMP division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Disables the SIMD evaluation of the given dense vector expression \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The SIMD-disabled dense vector.
//
// The \a nosimd function disables the SIMD evaluation of the given dense vector expression
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a nosimd function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = nosimd( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) nosimd( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecNoSIMDExpr<VT,TF>;
   return ReturnType( *dv );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct HasConstDataAccess< DVecNoSIMDExpr<VT,TF> >
   : public HasConstDataAccess<VT>
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
struct IsAligned< DVecNoSIMDExpr<VT,TF> >
   : public IsAligned<VT>
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
struct IsPadded< DVecNoSIMDExpr<VT,TF> >
   : public IsPadded<VT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
