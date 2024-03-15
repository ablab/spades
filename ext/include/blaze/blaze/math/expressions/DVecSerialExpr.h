//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecSerialExpr.h
//  \brief Header file for the dense vector serial evaluation expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSERIALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSERIALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecSerialExpr.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSERIALEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the forced serial evaluation of dense vectors.
// \ingroup dense_vector_expression
//
// The DVecSerialExpr class represents the compile time expression for the forced serial
// evaluation of a dense vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DVecSerialExpr
   : public VecSerialExpr< DenseVector< DVecSerialExpr<VT,TF>, TF > >
   , private Computation
{
 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecSerialExpr instance.
   using This = DVecSerialExpr<VT,TF>;

   //! Base type of this DVecSerialExpr instance.
   using BaseType = VecSerialExpr< DenseVector<This,TF> >;

   using ResultType    = ResultType_t<VT>;     //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<VT>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<VT>;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

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
   /*!\brief Constructor for the DVecSerialExpr class.
   //
   // \param dv The dense vector operand of the serial evaluation expression.
   */
   explicit inline DVecSerialExpr( const VT& dv ) noexcept
      : dv_( dv )  // Dense vector of the serial evaluation expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < dv_.size(), "Invalid vector access index" );
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
      if( index >= dv_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
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

   //**Conversion operator*************************************************************************
   /*!\brief Conversion to the type of the dense vector operand.
   //
   // \return The dense vector operand.
   */
   inline operator Operand() const noexcept {
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
   Operand dv_;  //!< Dense vector of the serial evaluation expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector serial evaluation expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector serial
   // evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector serial evaluation expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector serial
   // evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector serial evaluation expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector
   // serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector serial evaluation expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector
   // serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void addAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector serial evaluation expression to a dense
             vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector serial evaluation expression to a sparse
             vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void subAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector serial evaluation expression to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      multAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector serial evaluation expression to a sparse
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void multAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      multAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense vector serial evaluation expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector
   // serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      divAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense vector serial evaluation expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector
   // serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void divAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      divAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector serial evaluation expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector
   // serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void smpAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector serial evaluation expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector
   // serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void smpAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector serial evaluation expression to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void smpAddAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector serial evaluation expression to a sparse
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void smpAddAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector serial evaluation expression to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // vector serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void smpSubAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector serial evaluation expression to a sparse
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // vector serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void smpSubAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector serial evaluation expression to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a dense
   // vector serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void smpMultAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      multAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector serial evaluation expression to a
   //        sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a dense
   // vector serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void smpMultAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      multAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a dense vector serial evaluation expression to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side serial evaluation expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a dense vector
   // serial evaluation expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void smpDivAssign( DenseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      divAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to sparse vectors***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a dense vector serial evaluation expression to a sparse
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side serial evaluation expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a dense vector
   // serial evaluation expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void smpDivAssign( SparseVector<VT2,TF>& lhs, const DVecSerialExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      divAssign( *lhs, rhs.dv_ );
   }
   /*! \endcond */
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
/*!\brief Forces the serial evaluation of the given dense vector expression \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The evaluated dense vector.
//
// The \a serial function forces the serial evaluation of the given dense vector expression
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a serial function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = serial( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) serial( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecSerialExpr<VT,TF>;
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
struct IsAligned< DVecSerialExpr<VT,TF> >
   : public IsAligned<VT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
