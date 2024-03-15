//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecRepeatExpr.h
//  \brief Header file for the dense vector repeat expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECREPEATEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECREPEATEXPR_H_


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
#include <blaze/math/expressions/RepeatExprData.h>
#include <blaze/math/expressions/Transformation.h>
#include <blaze/math/expressions/VecRepeatExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/subvector/Dense.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECREPEATEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense vector repeat() function.
// \ingroup dense_vector_expression
//
// The DVecRepeatExpr class represents the compile time expression for repeating a dense vector
// via the repeat() function.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CRAs >  // Compile time repeater arguments
class DVecRepeatExpr
   : public VecRepeatExpr< DenseVector< DVecRepeatExpr<VT,TF,CRAs...>, TF >, CRAs... >
   , private If_t< IsComputation_v<VT>, Computation, Transformation >
   , private RepeatExprData<1UL,CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RepeatExprData<1UL,CRAs...>;  //!< The type of the RepeatExprData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecRepeatExpr instance.
   using This = DVecRepeatExpr<VT,TF,CRAs...>;

   //! Base type of this DVecRepeatExpr instance.
   using BaseType = VecRepeatExpr< DenseVector<This,TF>, CRAs... >;

   using ResultType    = RepeatTrait_t<VT,CRAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;            //!< Resulting element type.
   using ReturnType    = ReturnType_t<VT>;             //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite data type of the dense vector expression.
   using Operand = If_t< IsExpression_v<VT>, const VT, const VT& >;
   //**********************************************************************************************

 public:
   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecRepeatExpr class.
   //
   // \param dv The dense vector operand of the repeater expression.
   // \param args The number of repetitions.
   */
   template< typename... RRAs >  // Runtime repeater arguments
   explicit inline DVecRepeatExpr( const VT& dv, RRAs... args ) noexcept
      : DataType( args... )  // Base class initialization
      , dv_     ( dv )       // Dense vector of the repeater expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return dv_[index%dv_.size()];
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
      if( index >= size() ) {
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
      return dv_.size() * this->template repetitions<0UL>();
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
   using DataType::repetitions;
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return IsExpression_v<VT> && dv_.canAlias( alias );
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
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand dv_;  //!< Dense vector of the repeater expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector repeater expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side repeater expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector repeater
   // expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT2,TF>& lhs, const DVecRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CompositeType_t<VT> x( serial( rhs.dv_ ) );  // Evaluation of the dense vector operand

      const size_t reps( rhs.template repetitions<0UL>() );
      const size_t size( x.size() );

      for( size_t rep=0UL; rep<reps; ++rep ) {
         subvector( *lhs, rep*size, size, unchecked ) = serial( x );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector repeater expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side repeater expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector repeater
   // expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT2,TF>& lhs, const DVecRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CompositeType_t<VT> x( serial( rhs.dv_ ) );  // Evaluation of the dense vector operand

      const size_t reps( rhs.template repetitions<0UL>() );
      const size_t size( x.size() );

      (*lhs).reset();
      (*lhs).reserve( reps*size );

      for( size_t rep=0UL; rep<reps; ++rep ) {
         for( size_t i=0UL; i<size; ++i ) {
            (*lhs).append( rep*size+i, x[i], true );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector repeater expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side repeater expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector
   // repeater expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT2,TF>& lhs, const DVecRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CompositeType_t<VT> x( serial( rhs.dv_ ) );  // Evaluation of the dense vector operand

      const size_t reps( rhs.template repetitions<0UL>() );
      const size_t size( x.size() );

      for( size_t rep=0UL; rep<reps; ++rep ) {
         subvector( *lhs, rep*size, size, unchecked ) += serial( x );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector repeater expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side repeater expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector repeater expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT2,TF>& lhs, const DVecRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CompositeType_t<VT> x( serial( rhs.dv_ ) );  // Evaluation of the dense vector operand

      const size_t reps( rhs.template repetitions<0UL>() );
      const size_t size( x.size() );

      for( size_t rep=0UL; rep<reps; ++rep ) {
         subvector( *lhs, rep*size, size, unchecked ) -= serial( x );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector repeater expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side repeater expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector repeater expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT2,TF>& lhs, const DVecRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CompositeType_t<VT> x( serial( rhs.dv_ ) );  // Evaluation of the dense vector operand

      const size_t reps( rhs.template repetitions<0UL>() );
      const size_t size( x.size() );

      for( size_t rep=0UL; rep<reps; ++rep ) {
         subvector( *lhs, rep*size, size, unchecked ) *= serial( x );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense vector repeater expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side repeater expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector
   // repeater expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT2,TF>& lhs, const DVecRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CompositeType_t<VT> x( serial( rhs.dv_ ) );  // Evaluation of the dense vector operand

      const size_t reps( rhs.template repetitions<0UL>() );
      const size_t size( x.size() );

      for( size_t rep=0UL; rep<reps; ++rep ) {
         subvector( *lhs, rep*size, size, unchecked ) /= serial( x );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
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
/*!\brief Repeats the given dense vector.
// \ingroup dense_vector
//
// \param dv The dense vector to be repeated.
// \param repetitions The number of repetitions.
// \return The repeated dense vector.
//
// This function returns an expression representing the repeated dense vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<int,columnVector> a1{ 1, 0, -2 };
   blaze::DynamicVector<int,rowVector> b1{ 0, -1, 7 };

   blaze::DynamicVector<int,columnVector> a2;
   blaze::DynamicVector<int,rowVector> b2;

   // ... Resizing and initialization

   a2 = repeat( a1, 3UL );  // Results in ( 1, 0, -2, 1, 0, -2, 1, 0, -2 )
   b2 = repeat( b1, 3UL );  // Results in ( 0, -1, 7, 0, -1, 7, 0, -1, 7 )
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) repeat( const DenseVector<VT,TF>& dv, size_t repetitions )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecRepeatExpr<VT,TF>;
   return ReturnType( *dv, repetitions );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Repeats the given dense vector.
// \ingroup dense_vector
//
// \param dv The dense vector to be repeated.
// \return The repeated dense vector.
//
// This function returns an expression representing the repeated dense vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<int,columnVector> a1{ 1, 0, -2 };
   blaze::DynamicVector<int,rowVector> b1{ 0, -1, 7 };

   blaze::DynamicVector<int,columnVector> a2;
   blaze::DynamicVector<int,rowVector> b2;

   // ... Resizing and initialization

   a2 = repeat<3UL>( a1 );  // Results in ( 1, 0, -2, 1, 0, -2, 1, 0, -2 )
   b2 = repeat<3UL>( b1 );  // Results in ( 0, -1, 7, 0, -1, 7, 0, -1, 7 )
   \endcode
*/
template< size_t R0    // Compile time repetitions
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) repeat( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecRepeatExpr<VT,TF,R0>;
   return ReturnType( *dv );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Repeating the given dense vector.
// \ingroup dense_vector
//
// \param dv The dense vector to be repeated.
// \param repetitions The number of repetitions.
// \return The repeated dense vector.
//
// This auxiliary overload of the \c repeat() function accepts both a compile time and a runtime
// expansion. The runtime argument is discarded in favor of the compile time argument.
*/
template< size_t R0    // Compile time repetitions
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) repeat( const DenseVector<VT,TF>& dv, size_t repetitions )
{
   MAYBE_UNUSED( repetitions );

   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecRepeatExpr<VT,TF,R0>;
   return ReturnType( *dv );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
