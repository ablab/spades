//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecExpandExpr.h
//  \brief Header file for the sparse vector expansion expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECEXPANDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECEXPANDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/ExpandExprData.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/Transformation.h>
#include <blaze/math/expressions/VecExpandExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/ExpandTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/GetMemberType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SVECEXPANDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse vector expansion.
// \ingroup sparse_vector_expression
//
// The SVecExpandExpr class represents the compile time expression for expansions of
// sparse vectors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CEAs >  // Compile time expansion arguments
class SVecExpandExpr
   : public VecExpandExpr< SparseMatrix< SVecExpandExpr<VT,TF,CEAs...>, !TF >, CEAs... >
   , private If_t< IsComputation_v<VT>, Computation, Transformation >
   , private ExpandExprData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<VT>;  //!< Result type of the sparse vector expression.

   using DataType = ExpandExprData<CEAs...>;  //!< The type of the ExpandExprData base class.

   //! Definition of the GetConstIterator type trait.
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetConstIterator, ConstIterator, INVALID_TYPE );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the expansion expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the expansion expression. In case the given sparse
       vector expression of type \a VT is a computation or requires an intermediate evaluation,
       \a useAssign will be set to 1 and the expansion expression will be evaluated via the
       \a assign function family. Otherwise \a useAssign will be set to 0 and the expression
       will be evaluated via the subscript operator. */
   static constexpr bool useAssign = IsComputation_v<VT> || RequiresEvaluation_v<VT>;

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the target matrix is SMP assignable and the sparse vector operand requires an
       intermediate evaluation, the variable is set to 1 and the expression specific evaluation
       strategy is selected. Otherwise the variable is set to 0 and the default strategy is
       chosen. */
   template< typename MT >
    static constexpr bool UseSMPAssign_v = ( MT::smpAssignable && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SVecExpandExpr instance.
   using This = SVecExpandExpr<VT,TF,CEAs...>;

   //! Base type of this SVecExpandExpr instance.
   using BaseType = VecExpandExpr< SparseMatrix<This,!TF>, CEAs... >;

   using ResultType    = ExpandTrait_t<VT,CEAs...>;    //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;            //!< Resulting element type.
   using ReturnType    = ReturnType_t<VT>;             //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const SVecExpandExpr& >;

   //! Iterator over the elements of the sparse matrix.
   using ConstIterator = GetConstIterator_t<VT>;

   //! Composite data type of the sparse matrix expression.
   using Operand = If_t< IsExpression_v<VT>, const VT, const VT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = VT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SVecExpandExpr class.
   //
   // \param sv The sparse vector operand of the expansion expression.
   // \param args The runtime expansion expression arguments.
   */
   template< typename... REAs >  // Runtime expansion arguments
   explicit inline SVecExpandExpr( const VT& sv, REAs... args ) noexcept
      : DataType( args... )  // Base class initialization
      , sv_     ( sv )       // Sparse vector of the expansion expression
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      if( TF ) {
         BLAZE_INTERNAL_ASSERT( i < expansion(), "Invalid row access index"    );
         BLAZE_INTERNAL_ASSERT( j < sv_.size() , "Invalid column access index" );
         return sv_[j];
      }
      else {
         BLAZE_INTERNAL_ASSERT( i < sv_.size() , "Invalid row access index"    );
         BLAZE_INTERNAL_ASSERT( j < expansion(), "Invalid column access index" );
         return sv_[i];
      }
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid matrix access index.
   */
   inline ReturnType at( size_t i, size_t j ) const {
      if( i >= ( TF ? expansion() : sv_.size() ) ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= ( TF ? sv_.size() : expansion() ) ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      MAYBE_UNUSED( i );
      return sv_.begin();
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      MAYBE_UNUSED( i );
      return sv_.end();
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return ( TF ? expansion() : sv_.size() );
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return ( TF ? sv_.size() : expansion() );
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return sv_.nonZeros() * expansion();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row/column.
   //
   // \param i The index of the row/column.
   // \return The number of non-zero elements of row/column \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      MAYBE_UNUSED( i );
      return sv_.nonZeros();
   }
   //**********************************************************************************************

   //**Find function*******************************************************************************
   /*!\brief Searches for a specific matrix element.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the element in case the index is found, end() iterator otherwise.
   */
   inline ConstIterator find( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );
      return sv_.find( TF ? j : i );
   }
   //**********************************************************************************************

   //**LowerBound function*************************************************************************
   /*!\brief Returns an iterator to the first index not less then the given index.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the first index not less then the given index, end() iterator otherwise.
   */
   inline ConstIterator lowerBound( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );
      return sv_.lowerBound( TF ? j : i );
   }
   //**********************************************************************************************

   //**UpperBound function*************************************************************************
   /*!\brief Returns an iterator to the first index greater then the given index.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the first index greater then the given index, end() iterator otherwise.
   */
   inline ConstIterator upperBound( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );
      return sv_.upperBound( TF ? j : i );
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the sparse vector operand.
   //
   // \return The sparse vector operand.
   */
   inline Operand operand() const noexcept {
      return sv_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::expansion;
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return sv_.isAliased( alias );
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
      return sv_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return sv_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand sv_;  //!< Sparse vector of the expansion expression.
   //**********************************************************************************************

   //**Assignment to matrices**********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector expansion
   // expression to a matrix. Due to the explicit application of the SFINAE principle, this
   // function can only be selected by the compiler in case the operand requires an intermediate
   // evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( *rhs.sv_ ) );

      assign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to matrices*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse vector
   // expansion expression to a matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto addAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( *rhs.sv_ ) );

      addAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to matrices**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse vector
   // expansion expression to a matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto subAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( *rhs.sv_ ) );

      subAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to matrices********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // vector expansion expression to a matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand requires
   // an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto schurAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( *rhs.sv_ ) );

      schurAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // vector expansion expression to a matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand requires
   // an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto multAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( *rhs.sv_ ) );

      multAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to matrices******************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expansion to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse vector
   // expansion expansion to a matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the expression specific parallel
   // evaluation strategy is selected.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline auto smpAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( *rhs.sv_ );

      smpAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to matrices*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expansion to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // vector expansion expansion to a matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline auto smpAddAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( *rhs.sv_ );

      smpAddAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expansion to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // vector expansion expansion to a matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline auto smpSubAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( *rhs.sv_ );

      smpSubAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expansion for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // sparse vector expansion expansion to a matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline auto smpSchurAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( *rhs.sv_ );

      smpSchurAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a sparse vector expansion expression to a matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side expansion expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // sparse matrix expansion expression to a matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline auto smpMultAssign( Matrix<MT,SO>& lhs, const SVecExpandExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      using blaze::expand;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( *rhs.sv_ );

      smpMultAssign( *lhs, expand<CEAs...>( tmp, rhs.expansion() ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
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
/*!\brief Expansion of the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be expanded.
// \param expansion The expansion.
// \return The expansion of the vector.
//
// This function returns an expression representing the expansion of the given sparse vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedVector<int,columnVector> a{ 1, 0, -2, 0 }
   blaze::CompressedVector<int,rowVector> b{ 0, -1, 7, 0 }

   blaze::CompressedMatrix<double,columnMajor> A;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization

   // Expansion of the column vector 'a' to 4x3 column-major matrix
   //
   //    (  1  1  1 )
   //    (  0  0  0 )
   //    ( -2 -2 -2 )
   //    (  0  0  0 )
   //
   A = expand( a, 3UL );

   // Expansion of the row vector 'b' to a 3x4 row-major matrix
   //
   //    ( 0, -1, 7, 0 )
   //    ( 0, -1, 7, 0 )
   //    ( 0, -1, 7, 0 )
   //
   B = expand( b, 3UL );
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) expand( const SparseVector<VT,TF>& sv, size_t expansion )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const SVecExpandExpr<VT,TF>;
   return ReturnType( *sv, expansion );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be expanded.
// \return The expansion of the vector.
//
// This function returns an expression representing the expansion of the given sparse vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedVector<int,columnVector> a{ 1, 0, -2, 0 }
   blaze::CompressedVector<int,rowVector> b{ 0, -1, 7, 0 }

   blaze::CompressedMatrix<double,columnMajor> A;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization

   // Expansion of the column vector 'a' to 4x3 column-major matrix
   //
   //    (  1  1  1 )
   //    (  0  0  0 )
   //    ( -2 -2 -2 )
   //    (  0  0  0 )
   //
   A = expand<3UL>( a );

   // Expansion of the row vector 'b' to a 3x4 row-major matrix
   //
   //    ( 0, -1, 7, 0 )
   //    ( 0, -1, 7, 0 )
   //    ( 0, -1, 7, 0 )
   //
   B = expand<3UL>( b );
   \endcode
*/
template< size_t E     // Compile time expansion argument
        , typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) expand( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const SVecExpandExpr<VT,TF,E>;
   return ReturnType( *sv );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expansion of the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be expanded.
// \param expansion The expansion.
// \return The expansion of the vector.
//
// This auxiliary overload of the \c expand() function accepts both a compile time and a runtime
// expansion. The runtime argument is discarded in favor of the compile time argument.
*/
template< size_t E     // Compile time expansion argument
        , typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) expand( const SparseVector<VT,TF>& sv, size_t expansion )
{
   MAYBE_UNUSED( expansion );

   BLAZE_FUNCTION_TRACE;

   using ReturnType = const SVecExpandExpr<VT,TF,E>;
   return ReturnType( *sv );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
