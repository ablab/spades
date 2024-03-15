//=================================================================================================
/*!
//  \file blaze/math/views/columns/Dense.h
//  \brief Columns specialization for dense matrices
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

#ifndef _BLAZE_MATH_VIEWS_COLUMNS_DENSE_H_
#define _BLAZE_MATH_VIEWS_COLUMNS_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/Columns.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnsTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/columns/BaseTemplate.h>
#include <blaze/math/views/columns/ColumnsData.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Columns for column selections on column-major dense matrices.
// \ingroup columns
//
// This specialization of Columns adapts the class template to the requirements of column-major
// dense matrices.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
class Columns<MT,true,true,SF,CCAs...>
   : public View< DenseMatrix< Columns<MT,true,true,SF,CCAs...>, true > >
   , private ColumnsData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnsData<CCAs...>;                 //!< The type of the ColumnsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT1, typename MT2 >
   static constexpr bool EnforceEvaluation_v =
      ( IsRestricted_v<MT1> && RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Columns instance.
   using This = Columns<MT,true,true,SF,CCAs...>;

   //! Base type of this Columns instance.
   using BaseType = View< DenseMatrix<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Columns instance.
   using ResultType    = ColumnsTrait_t<MT,N>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const Columns&;               //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant column value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant column value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<MT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, Iterator_t<MT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Columns( MT& matrix, RCAs... args );

   Columns( const Columns& ) = default;
   Columns( Columns&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Columns() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t j ) noexcept;
   inline ConstPointer   data  ( size_t j ) const noexcept;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Columns& operator=( const ElementType& rhs );
   inline Columns& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Columns& operator=( const Columns& rhs );

   template< typename MT2, bool SO2 >
   inline Columns& operator=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t rows() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Columns& transpose();
   inline Columns& ctranspose();

   template< typename Other > inline Columns& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<MT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<MT2> > );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CCAs2 >
   inline bool canAlias( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CCAs2 >
   inline bool isAliased( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the columns.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CCAs2 > friend class Columns;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMNS_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for column selections on column-major dense matrices.
//
// \param matrix The matrix containing the columns.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Columns<MT,true,true,SF,CCAs...>::Columns( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the columns
{
   if( isChecked( args... ) ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( matrix_.columns() <= idx(j) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::Reference
   Columns<MT,true,true,SF,CCAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstReference
   Columns<MT,true,true,SF,CCAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::Reference
   Columns<MT,true,true,SF,CCAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstReference
   Columns<MT,true,true,SF,CCAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column selection. Note
// that you can NOT assume that all matrix elements lie adjacent to each other! The underlying
// matrix may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::Pointer
   Columns<MT,true,true,SF,CCAs...>::data() noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column selection. Note
// that you can NOT assume that all matrix elements lie adjacent to each other! The underlying
// matrix may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstPointer
   Columns<MT,true,true,SF,CCAs...>::data() const noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of column \a i.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::Pointer
   Columns<MT,true,true,SF,CCAs...>::data( size_t j ) noexcept
{
   return matrix_.data( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of column \a i.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstPointer
   Columns<MT,true,true,SF,CCAs...>::data( size_t j ) const noexcept
{
   return matrix_.data( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::Iterator
   Columns<MT,true,true,SF,CCAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.begin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstIterator
   Columns<MT,true,true,SF,CCAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstIterator
   Columns<MT,true,true,SF,CCAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::Iterator
   Columns<MT,true,true,SF,CCAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.end( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstIterator
   Columns<MT,true,true,SF,CCAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cend( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,true,SF,CCAs...>::ConstIterator
   Columns<MT,true,true,SF,CCAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cend( idx(j) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Homogenous assignment to all column selection elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column selection.
//
// This function homogeneously assigns the given value to all column selection elements. Note
// that in case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal
// elements of the underlying matrix are modified.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::operator=( const ElementType& rhs )
{
   for( size_t j=0UL; j<columns(); ++j ) {
      column( matrix_, idx(j), unchecked ) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column selection elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column selection.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the column
// selection by means of an initializer list. The column selection elements are assigned the
// values from the given initializer list. Missing values are initialized as default. Note that
// if the size of the top-level initializer list does not match the number of rows of the column
// selection or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column selection" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerMatrix<ElementType> tmp( list, columns() );
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      size_t j( 0UL );
      for( const auto& element : rowList ) {
         matrix_(i,idx(j)) = element;
         ++j;
      }
      for( ; j<columns(); ++j ) {
         matrix_(i,idx(j)) = ElementType();
      }
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Columns.
//
// \param rhs Dense column selection to be copied.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense column selection is initialized as a copy of the given dense column selection. In
// case the current sizes of the two column selections don't match, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper
// triangular, or symmetric matrix and the assignment would violate its lower, upper, or
// symmetry property, respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::operator=( const Columns& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && compareIndices( *this, rhs ) ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense column selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::operator=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<MT2>, const MT2& >;
   Right right( *rhs );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( right, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<MT2> ) {
      reset();
   }

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,true,true,SF,CCAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAddAssign( column( matrix_, idx(j), unchecked ), column( *rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( (*rhs).canAlias( this ) ) {
      const AddType tmp( *this + (*rhs) );
      smpAssign( left, tmp );
   }
   else {
      smpAddAssign( left, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,true,true,SF,CCAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( tmp, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,true,true,SF,CCAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !trySubAssign( column( matrix_, idx(j), unchecked ), column( *rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( (*rhs).canAlias( this ) ) {
      const SubType tmp( *this - (*rhs ) );
      smpAssign( left, tmp );
   }
   else {
      smpSubAssign( left, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,true,true,SF,CCAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( tmp, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,true,true,SF,CCAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryMultAssign( column( matrix_, idx(j), unchecked ), column( *rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( (*rhs).canAlias( this ) ) {
      const SchurType tmp( *this % (*rhs) );
      if( IsSparseMatrix_v<SchurType> )
         reset();
      smpAssign( left, tmp );
   }
   else {
      smpSchurAssign( left, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,true,true,SF,CCAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( tmp, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      reset();
   }

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline MT& Columns<MT,true,true,SF,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline const MT& Columns<MT,true,true,SF,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the column selection.
//
// \return The number of rows of the column selection.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,true,SF,CCAs...>::rows() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total number
// of elements of a column.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,true,SF,CCAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column selection.
//
// \return The capacity of the dense column selection.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,true,SF,CCAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
//
// This function returns the current capacity of the specified column.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,true,SF,CCAs...>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense column selection.
//
// \return The number of non-zero elements in the dense column selection.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,true,SF,CCAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns(); ++j ) {
      nonzeros += matrix_.nonZeros( idx(j) );
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
//
// This function returns the current number of non-zero elements in the specified column.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,true,SF,CCAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.nonZeros( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,true,SF,CCAs...>::reset()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      matrix_.reset( idx(j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function resets the values in the specified column to their default value.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,true,SF,CCAs...>::reset( size_t j )
{
   matrix_.reset( idx(j) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense column selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the dense column selection.
//
// This function scales the column selection by applying the given scalar value \a scalar to each
// element of the column selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to
// scale a column selection on a lower or upper unitriangular matrix. The attempt to scale such
// a column selection results in a compile time error!
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the scalar value
inline Columns<MT,true,true,SF,CCAs...>&
   Columns<MT,true,true,SF,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index ( idx(j) );
      const size_t ibegin( IsLower<MT>::value ? ( IsStrictlyLower_v<MT> ? index+1UL : index ) : 0UL );
      const size_t iend  ( IsUpper<MT>::value ? ( IsStrictlyUpper_v<MT> ? index : index+1UL ) : rows() );

      for( size_t i=ibegin; i<iend; ++i ) {
         matrix_(i,index) *= scalar;
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address can alias with the dense column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,true,true,SF,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can alias with the given dense column selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address can alias with the dense column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , bool SF              // Symmetry flag
        , typename... CCAs >   // Compile time column arguments
template< typename MT2         // Data type of the foreign dense column selection
        , bool SO2             // Storage order of the foreign dense column selection
        , bool SF2             // Symmetry flag of the foreign dense column selection
        , typename... CCAs2 >  // Compile time column arguments of the foreign dense column selection
inline bool
   Columns<MT,true,true,SF,CCAs...>::canAlias( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense column selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,true,true,SF,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is aliased with the given dense column selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense column selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , bool SF              // Symmetry flag
        , typename... CCAs >   // Compile time column arguments
template< typename MT2         // Data type of the foreign dense column selection
        , bool SO2             // Storage order of the foreign dense column selection
        , bool SF2             // Symmetry flag of the foreign dense column selection
        , typename... CCAs2 >  // Compile time column arguments of the foreign dense column selection
inline bool
   Columns<MT,true,true,SF,CCAs...>::isAliased( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is properly aligned in memory.
//
// \return \a true in case the dense column selection is aligned, \a false if not.
//
// This function returns whether the dense column selection is guaranteed to be properly aligned
// in memory, i.e. whether the beginning and the end of the dense column selection are guaranteed
// to conform to the alignment restrictions of the element type \a Type.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,true,true,SF,CCAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can be used in SMP assignments.
//
// \return \a true in case the dense column selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense column selection can be used in SMP assignments. In
// contrast to the \a smpAssignable member enumeration, which is based solely on compile time
// information, this function additionally provides runtime information (as for instance the
// current number of rows and/or columns of the dense column selection).
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,true,true,SF,CCAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense column selection. The
// row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Columns<MT,true,true,SF,CCAs...>::SIMDType
   Columns<MT,true,true,SF,CCAs...>::load( size_t i, size_t j ) const noexcept
{
   return matrix_.load( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense column selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Columns<MT,true,true,SF,CCAs...>::SIMDType
   Columns<MT,true,true,SF,CCAs...>::loada( size_t i, size_t j ) const noexcept
{
   return matrix_.loada( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense column
// selection. The row index must be smaller than the number of rows and the column index must be
// smaller than the number of columns. Additionally, the row index must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Columns<MT,true,true,SF,CCAs...>::SIMDType
   Columns<MT,true,true,SF,CCAs...>::loadu( size_t i, size_t j ) const noexcept
{
   return matrix_.loadu( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense column selection. The
// row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Columns<MT,true,true,SF,CCAs...>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.store( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense column selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Columns<MT,true,true,SF,CCAs...>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.storea( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense column
// selection. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the row index must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Columns<MT,true,true,SF,CCAs...>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.storeu( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the dense
// column selection. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a multiple of
// the number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Columns<MT,true,true,SF,CCAs...>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.stream( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::assign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(i    ,index) = (*rhs)(i    ,j);
         matrix_(i+1UL,index) = (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(ipos,index) = (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::assign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   if( useStreaming &&
       rows()*columns() > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t j=0UL; j<columns(); ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT2> right( (*rhs).begin(j) );

         for( ; i<ipos; i+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<rows(); ++i ) {
            *left = *right;
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<columns(); ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT2> right( (*rhs).begin(j) );

         for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<ipos; i+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<rows(); ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,true,true,SF,CCAs...>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) = (*rhs)(i    ,j);
            matrix_(i+1UL,index) = (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() ) {
            matrix_(ipos,index) = (*rhs)(ipos,j);
         }
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) = (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(element->index(),index) = element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(i,idx(element->index())) = element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );
      if( IsDiagonal_v<MT2> ) {
         matrix_(j,index) += (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) += (*rhs)(i    ,j);
            matrix_(i+1UL,index) += (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(ipos,index) += (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT2> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT2> )
                           ?( IsStrictlyUpper_v<MT2> ? j : j+1UL )
                           :( rows() ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( prevMultiple( iend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT2> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<iend; ++i ) {
         *left += *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,true,true,SF,CCAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) += (*rhs)(i    ,j);
            matrix_(i+1UL,index) += (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() )
            matrix_(ipos,index) += (*rhs)(ipos,j);
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) += (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(element->index(),index) += element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(i,idx(element->index())) += element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );

      if( IsDiagonal_v<MT2> ) {
         matrix_(j,index) -= (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) -= (*rhs)(i    ,j);
            matrix_(i+1UL,index) -= (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(ipos,index) -= (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT2> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT2> )
                           ?( IsStrictlyUpper_v<MT2> ? j : j+1UL )
                           :( rows() ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( prevMultiple( iend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT2> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<iend; ++i ) {
         *left -= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,true,true,SF,CCAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) -= (*rhs)(i    ,j);
            matrix_(i+1UL,index) -= (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() )
            matrix_(ipos,index) -= (*rhs)(ipos,j);
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) -= (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(element->index(),index) -= element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(i,idx(element->index())) -= element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(i    ,index) *= (*rhs)(i    ,j);
         matrix_(i+1UL,index) *= (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(ipos,index) *= (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a column-major dense
//        matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Columns<MT,true,true,SF,CCAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ipos( prevMultiple( rows(), SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

      size_t i( 0UL );
      Iterator left( begin(j) );
      ConstIterator_t<MT2> right( (*rhs).begin(j) );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<rows(); ++i ) {
         *left *= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,true,true,SF,CCAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) *= (*rhs)(i    ,j);
            matrix_(i+1UL,index) *= (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() )
            matrix_(ipos,index) *= (*rhs)(ipos,j);
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) *= (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(i,index) );
         matrix_(i,index) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(i,index) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,true,SF,CCAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(i,idx(j)) );
         matrix_(i,idx(j)) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(i,idx(j)) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL ROW-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Columns for column selections on general row-major dense matrices.
// \ingroup columns
//
// This specialization of Columns adapts the class template to the requirements of general
// row-major dense matrices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
class Columns<MT,false,true,false,CCAs...>
   : public View< DenseMatrix< Columns<MT,false,true,false,CCAs...>, true > >
   , private ColumnsData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnsData<CCAs...>;                 //!< The type of the ColumnsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT1, typename MT2 >
   static constexpr bool EnforceEvaluation_v =
      ( IsRestricted_v<MT1> && RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Columns instance.
   using This = Columns<MT,false,true,false,CCAs...>;

   //! Base type of this Columns instance.
   using BaseType = View< DenseMatrix<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Columns instance.
   using ResultType    = ColumnsTrait_t<MT,N>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const Columns&;               //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant column value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant column value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**ColumnsIterator class definition************************************************************
   /*!\brief Iterator over the elements of a selected dense column.
   */
   template< typename MatrixType      // Type of the dense matrix
           , typename IteratorType >  // Type of the dense matrix iterator
   class ColumnsIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! The iterator category.
      using IteratorCategory = typename std::iterator_traits<IteratorType>::iterator_category;

      //! Type of the underlying elements.
      using ValueType = typename std::iterator_traits<IteratorType>::value_type;

      //! Pointer return type.
      using PointerType = typename std::iterator_traits<IteratorType>::pointer;

      //! Reference return type.
      using ReferenceType = typename std::iterator_traits<IteratorType>::reference;

      //! Difference between two iterators.
      using DifferenceType = typename std::iterator_traits<IteratorType>::difference_type;

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the ColumnsIterator class.
      */
      inline ColumnsIterator() noexcept
         : matrix_( nullptr )  // The dense matrix containing the column
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   (     )      // Iterator to the current dense element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ColumnsIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnsIterator( MatrixType& matrix, size_t row, size_t column ) noexcept
         : matrix_( &matrix )  // The dense matrix containing the selected column
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   (         )  // Iterator to the current dense element
      {
         if( row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ColumnsIterator instances.
      //
      // \param it The column iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline ColumnsIterator( const ColumnsIterator<MatrixType2,IteratorType2>& it ) noexcept
         : matrix_( it.matrix_ )  // The dense matrix containing the seleted column
         , row_   ( it.row_    )  // The current row index
         , column_( it.column_ )  // The current column index
         , pos_   ( it.pos_    )  // Iterator to the current dense element
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ColumnsIterator& operator+=( size_t inc ) noexcept {
         using blaze::reset;
         row_ += inc;
         if( row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline ColumnsIterator& operator-=( size_t dec ) noexcept {
         using blaze::reset;
         row_ -= dec;
         if( row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ColumnsIterator& operator++() noexcept {
         using blaze::reset;
         ++row_;
         if( row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ColumnsIterator operator++( int ) noexcept {
         const ColumnsIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ColumnsIterator& operator--() noexcept {
         using blaze::reset;
         --row_;
         if( row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ColumnsIterator operator--( int ) noexcept {
         const ColumnsIterator tmp( *this );
         --(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Subscript operator***********************************************************************
      /*!\brief Direct access to the dense row elements.
      //
      // \param index Access index.
      // \return Reference to the accessed value.
      */
      inline ReferenceType operator[]( size_t index ) const {
         BLAZE_USER_ASSERT( row_+index < matrix_->rows(), "Invalid access index detected" );
         const IteratorType pos( matrix_->begin( row_+index ) + column_ );
         return *pos;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense column element at the current iterator position.
      //
      // \return Reference to the current value.
      */
      inline ReferenceType operator*() const {
         return *pos_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense column element at the current iterator position.
      //
      // \return Pointer to the dense column element at the current iterator position.
      */
      inline PointerType operator->() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ == rhs.row_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ < rhs.row_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ > rhs.row_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<=( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ <= rhs.row_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>=( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ >= rhs.row_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two column iterators.
      //
      // \param rhs The right-hand side column iterator.
      // \return The number of elements between the two column iterators.
      */
      inline DifferenceType operator-( const ColumnsIterator& rhs ) const noexcept {
         return row_ - rhs.row_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ColumnsIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ColumnsIterator operator+( const ColumnsIterator& it, size_t inc ) noexcept {
         return ColumnsIterator( *it.matrix_, it.row_+inc, it.column_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ColumnsIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ColumnsIterator operator+( size_t inc, const ColumnsIterator& it ) noexcept {
         return ColumnsIterator( *it.matrix_, it.row_+inc, it.column_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ColumnsIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ColumnsIterator operator-( const ColumnsIterator& it, size_t dec ) noexcept {
         return ColumnsIterator( *it.matrix_, it.row_-dec, it.column_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The dense matrix containing the selected column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current dense element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class ColumnsIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = ColumnsIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, ColumnsIterator< MT, Iterator_t<MT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Columns( MT& matrix, RCAs... args );

   Columns( const Columns& ) = default;
   Columns( Columns&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Columns() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t j ) noexcept;
   inline ConstPointer   data  ( size_t j ) const noexcept;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Columns& operator=( const ElementType& rhs );
   inline Columns& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Columns& operator=( const Columns& rhs );

   template< typename MT2, bool SO2 >
   inline Columns& operator=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t rows() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Columns& transpose();
   inline Columns& ctranspose();

   template< typename Other > inline Columns& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CCAs2 >
   inline bool canAlias( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CCAs2 >
   inline bool isAliased( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the columns.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CCAs2 > friend class Columns;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMNS_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE       ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for column selections on general row-major dense matrices.
//
// \param matrix The matrix containing the columns.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Columns<MT,false,true,false,CCAs...>::Columns( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the columns
{
   if( isChecked( args... ) ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( matrix_.columns() <= idx(j) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::Reference
   Columns<MT,false,true,false,CCAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstReference
   Columns<MT,false,true,false,CCAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::Reference
   Columns<MT,false,true,false,CCAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstReference
   Columns<MT,false,true,false,CCAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column selection. Note
// that you can NOT assume that all matrix elements lie adjacent to each other! The underlying
// matrix may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::Pointer
   Columns<MT,false,true,false,CCAs...>::data() noexcept
{
   return matrix_.data() + idx(0UL);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column selection. Note
// that you can NOT assume that all matrix elements lie adjacent to each other! The underlying
// matrix may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstPointer
   Columns<MT,false,true,false,CCAs...>::data() const noexcept
{
   return matrix_.data() + idx(0UL);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::Pointer
   Columns<MT,false,true,false,CCAs...>::data( size_t j ) noexcept
{
   return matrix_.data() + idx(j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstPointer
   Columns<MT,false,true,false,CCAs...>::data( size_t j ) const noexcept
{
   return matrix_.data() + idx(j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::Iterator
   Columns<MT,false,true,false,CCAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return Iterator( matrix_, 0UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstIterator
   Columns<MT,false,true,false,CCAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return ConstIterator( matrix_, 0UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstIterator
   Columns<MT,false,true,false,CCAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return ConstIterator( matrix_, 0UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::Iterator
   Columns<MT,false,true,false,CCAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return Iterator( matrix_, rows(), idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstIterator
   Columns<MT,false,true,false,CCAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return ConstIterator( matrix_, rows(), idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,false,CCAs...>::ConstIterator
   Columns<MT,false,true,false,CCAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return ConstIterator( matrix_, rows(), idx(j) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Homogenous assignment to all column selection elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column selection.
//
// This function homogeneously assigns the given value to all column selection elements. Note
// that in case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal
// elements of the underlying matrix are modified.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::operator=( const ElementType& rhs )
{
   for( size_t j=0UL; j<columns(); ++j ) {
      column( matrix_, idx(j), unchecked ) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column selection elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid initializer list dimension.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the column
// selection by means of an initializer list. The column selection elements are assigned the
// values from the given initializer list. Missing values are initialized as default. Note that
// if the size of the top-level initializer list does not match the number of rows of the column
// selection or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column selection" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerMatrix<ElementType> tmp( list, columns() );
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      size_t j( 0UL );
      for( const auto& element : rowList ) {
         matrix_(i,idx(j)) = element;
         ++j;
      }
      for( ; j<columns(); ++j ) {
         matrix_(i,idx(j)) = ElementType();
      }
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Columns.
//
// \param rhs Dense column selection to be copied.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense column selection is initialized as a copy of the given dense column selection. In
// case the current sizes of the two column selections don't match, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper
// triangular, or symmetric matrix and the assignment would violate its lower, upper, or
// symmetry property, respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::operator=( const Columns& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && compareIndices( *this, rhs ) ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense column selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::operator=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<MT2>, const MT2& >;
   Right right( *rhs );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( right, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<MT2> ) {
      reset();
   }

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,false,true,false,CCAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAddAssign( column( matrix_, idx(j), unchecked ), column( *rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( (*rhs).canAlias( this ) ) {
      const AddType tmp( *this + (*rhs) );
      smpAssign( left, tmp );
   }
   else {
      smpAddAssign( left, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,false,true,false,CCAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( tmp, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,false,true,false,CCAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !trySubAssign( column( matrix_, idx(j), unchecked ), column( *rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( (*rhs).canAlias( this ) ) {
      const SubType tmp( *this - (*rhs ) );
      smpAssign( left, tmp );
   }
   else {
      smpSubAssign( left, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the column selection.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,false,true,false,CCAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( tmp, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,false,true,false,CCAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryMultAssign( column( matrix_, idx(j), unchecked ), column( *rhs, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( (*rhs).canAlias( this ) ) {
      const SchurType tmp( *this % (*rhs) );
      if( IsSparseMatrix_v<SchurType> )
         reset();
      smpAssign( left, tmp );
   }
   else {
      smpSchurAssign( left, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Columns<MT,false,true,false,CCAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Columns& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( column( matrix_, idx(j), unchecked ), column( tmp, j, unchecked ), 0UL ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      reset();
   }

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline MT& Columns<MT,false,true,false,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline const MT& Columns<MT,false,true,false,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the column selection.
//
// \return The number of rows of the column selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,false,CCAs...>::rows() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total number
// of elements of a column.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,false,CCAs...>::spacing() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column selection.
//
// \return The capacity of the dense column selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,false,CCAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
//
// This function returns the current capacity of the specified column.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,false,CCAs...>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense column selection.
//
// \return The number of non-zero elements in the dense column selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,false,CCAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns(); ++j ) {
      nonzeros += nonZeros( j );
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
//
// This function returns the current number of non-zero elements in the specified column.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,false,CCAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   size_t nonzeros( 0UL );

   const size_t index( idx(j) );
   for( size_t i=0UL; i<rows(); ++i ) {
      if( !isDefault( matrix_( i, index ) ) )
         ++nonzeros;
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,true,false,CCAs...>::reset()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      reset( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function resets the values in the specified column to their default value.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,true,false,CCAs...>::reset( size_t j )
{
   using blaze::reset;

   const size_t index( idx(j) );
   for( size_t i=0UL; i<rows(); ++i ) {
      reset( matrix_( i, index ) );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense column selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the dense column selection.
//
// This function scales the column selection by applying the given scalar value \a scalar to each
// element of the column selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to
// scale a column selection on a lower or upper unitriangular matrix. The attempt to scale such
// a column selection results in a compile time error!
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the scalar value
inline Columns<MT,false,true,false,CCAs...>&
   Columns<MT,false,true,false,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index ( idx(j) );
      const size_t ibegin( IsLower<MT>::value ? ( IsStrictlyLower_v<MT> ? index+1UL : index ) : 0UL );
      const size_t iend  ( IsUpper<MT>::value ? ( IsStrictlyUpper_v<MT> ? index : index+1UL ) : rows() );

      for( size_t i=ibegin; i<iend; ++i ) {
         matrix_(i,index) *= scalar;
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address can alias with the dense column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,true,false,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can alias with the given dense column
//        selection \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address can alias with the dense column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CCAs >   // Compile time column arguments
template< typename MT2         // Data type of the foreign dense column selection
        , bool SO2             // Storage order of the foreign dense column selection
        , bool SF2             // Symmetry flag of the foreign dense column selection
        , typename... CCAs2 >  // Compile time column arguments of the foreign dense column selection
inline bool
   Columns<MT,false,true,false,CCAs...>::canAlias( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense column selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,true,false,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is aliased with the given dense column
//        selection \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense column selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CCAs >   // Compile time column arguments
template< typename MT2         // Data type of the foreign dense column selection
        , bool SO2             // Storage order of the foreign dense column selection
        , bool SF2             // Symmetry flag of the foreign dense column selection
        , typename... CCAs2 >  // Compile time column arguments of the foreign dense column selection
inline bool
   Columns<MT,false,true,false,CCAs...>::isAliased( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is properly aligned in memory.
//
// \return \a true in case the dense column selection is aligned, \a false if not.
//
// This function returns whether the dense column selection is guaranteed to be properly aligned
// in memory, i.e. whether the beginning and the end of the dense column selection are guaranteed
// to conform to the alignment restrictions of the element type \a Type.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,false,true,false,CCAs...>::isAligned() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can be used in SMP assignments.
//
// \return \a true in case the dense column selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense column selection can be used in SMP assignments. In
// contrast to the \a smpAssignable member enumeration, which is based solely on compile time
// information, this function additionally provides runtime information (as for instance the
// current number of rows and/or columns of the dense column selection).
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,false,true,false,CCAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(i    ,index) = (*rhs)(i    ,j);
         matrix_(i+1UL,index) = (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(ipos,index) = (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) = (*rhs)(i    ,j);
            matrix_(i+1UL,index) = (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() ) {
            matrix_(ipos,index) = (*rhs)(ipos,j);
         }
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) = (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(element->index(),index) = element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(i,idx(element->index())) = element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );
      if( IsDiagonal_v<MT2> ) {
         matrix_(j,index) += (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) += (*rhs)(i    ,j);
            matrix_(i+1UL,index) += (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(ipos,index) += (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) += (*rhs)(i    ,j);
            matrix_(i+1UL,index) += (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() )
            matrix_(ipos,index) += (*rhs)(ipos,j);
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) += (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(element->index(),index) += element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(i,idx(element->index())) += element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );

      if( IsDiagonal_v<MT2> ) {
         matrix_(j,index) -= (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) -= (*rhs)(i    ,j);
            matrix_(i+1UL,index) -= (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(ipos,index) -= (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) -= (*rhs)(i    ,j);
            matrix_(i+1UL,index) -= (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() )
            matrix_(ipos,index) -= (*rhs)(ipos,j);
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) -= (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(element->index(),index) -= element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(i,idx(element->index())) -= element->value();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(i    ,index) *= (*rhs)(i    ,j);
         matrix_(i+1UL,index) *= (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(ipos,index) *= (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Columns<MT,false,true,false,CCAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t ipos( prevMultiple( (*rhs).rows(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).rows(), "Invalid end calculation" );

      for( size_t j=0UL; j<columns(); ++j ) {
         const size_t index( idx(j) );
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(i    ,index) *= (*rhs)(i    ,j);
            matrix_(i+1UL,index) *= (*rhs)(i+1UL,j);
         }
         if( ipos < (*rhs).rows() )
            matrix_(ipos,index) *= (*rhs)(ipos,j);
      }
   }
   else
   {
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t ii=0UL; ii<rows(); ii+=block ) {
            const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
            for( size_t j=jj; j<jend; ++j ) {
               const size_t index( idx(j) );
               for( size_t i=ii; i<iend; ++i ) {
                  matrix_(i,index) *= (*rhs)(i,j);
               }
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(i,index) );
         matrix_(i,index) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(i,index) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,true,false,CCAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(i,idx(j)) );
         matrix_(i,idx(j)) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(i,idx(j)) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC ROW-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Columns for column selections on symmetric row-major dense matrices.
// \ingroup columns
//
// This specialization of Columns adapts the class template to the requirements of symmetric
// row-major dense matrices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
class Columns<MT,false,true,true,CCAs...>
   : public View< DenseMatrix< Columns<MT,false,true,true,CCAs...>, true > >
   , private ColumnsData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnsData<CCAs...>;                 //!< The type of the ColumnsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Columns instance.
   using This = Columns<MT,false,true,true,CCAs...>;

   //! Base type of this Columns instance.
   using BaseType = View< DenseMatrix<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Columns instance.
   using ResultType    = ColumnsTrait_t<MT,N>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const Columns&;               //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant column value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant column value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<MT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, Iterator_t<MT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Columns( MT& matrix, RCAs... args );

   Columns( const Columns& ) = default;
   Columns( Columns&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Columns() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t j ) noexcept;
   inline ConstPointer   data  ( size_t j ) const noexcept;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Columns& operator=( const ElementType& rhs );

   Columns& operator=( const Columns& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t rows() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CCAs2 >
   inline bool canAlias( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CCAs2 >
   inline bool isAliased( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the columns.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CCAs2 > friend class Columns;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMNS_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for column selections on symmetric row-major dense matrices.
//
// \param matrix The matrix containing the columns.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Columns<MT,false,true,true,CCAs...>::Columns( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the columns
{
   if( isChecked( args... ) ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( matrix_.columns() <= idx(j) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::Reference
   Columns<MT,false,true,true,CCAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(idx(j),i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstReference
   Columns<MT,false,true,true,CCAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(idx(j),i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::Reference
   Columns<MT,false,true,true,CCAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstReference
   Columns<MT,false,true,true,CCAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column selection. Note
// that you can NOT assume that all matrix elements lie adjacent to each other! The underlying
// matrix may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::Pointer
   Columns<MT,false,true,true,CCAs...>::data() noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column selection. Note
// that you can NOT assume that all matrix elements lie adjacent to each other! The underlying
// matrix may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstPointer
   Columns<MT,false,true,true,CCAs...>::data() const noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of column \a i.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::Pointer
   Columns<MT,false,true,true,CCAs...>::data( size_t j ) noexcept
{
   return matrix_.data( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of column \a i.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstPointer
   Columns<MT,false,true,true,CCAs...>::data( size_t j ) const noexcept
{
   return matrix_.data( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::Iterator
   Columns<MT,false,true,true,CCAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.begin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstIterator
   Columns<MT,false,true,true,CCAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns an iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstIterator
   Columns<MT,false,true,true,CCAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::Iterator
   Columns<MT,false,true,true,CCAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.end( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstIterator
   Columns<MT,false,true,true,CCAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cend( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,true,true,CCAs...>::ConstIterator
   Columns<MT,false,true,true,CCAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return matrix_.cend( idx(j) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Homogenous assignment to all column selection elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column selection.
//
// This function homogeneously assigns the given value to all column selection elements. Note
// that in case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal
// elements of the underlying matrix are modified.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,true,true,CCAs...>&
   Columns<MT,false,true,true,CCAs...>::operator=( const ElementType& rhs )
{
   for( size_t j=0UL; j<columns(); ++j ) {
      row( matrix_, idx(j), unchecked ) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline MT& Columns<MT,false,true,true,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline const MT& Columns<MT,false,true,true,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the column selection.
//
// \return The number of rows of the column selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,true,CCAs...>::rows() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total number
// of elements of a column.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,true,CCAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column selection.
//
// \return The capacity of the dense column selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,true,CCAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
//
// This function returns the current capacity of the specified column.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,true,CCAs...>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense column selection.
//
// \return The number of non-zero elements in the dense column selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,true,CCAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns(); ++j ) {
      nonzeros += matrix_.nonZeros( idx(j) );
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
//
// This function returns the current number of non-zero elements in the specified column.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,true,true,CCAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.nonZeros( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,true,true,CCAs...>::reset()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      matrix_.reset( idx(j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function resets the values in the specified column to their default value.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,true,true,CCAs...>::reset( size_t j )
{
   matrix_.reset( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address can alias with the dense column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,true,true,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can alias with the given dense column selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address can alias with the dense column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CCAs >   // Compile time column arguments
template< typename MT2         // Data type of the foreign dense column selection
        , bool SO2             // Storage order of the foreign dense column selection
        , bool SF2             // Symmetry flag of the foreign dense column selection
        , typename... CCAs2 >  // Compile time column arguments of the foreign dense column selection
inline bool
   Columns<MT,false,true,true,CCAs...>::canAlias( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense column selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,true,true,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is aliased with the given dense column selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense column selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CCAs >   // Compile time column arguments
template< typename MT2         // Data type of the foreign dense column selection
        , bool SO2             // Storage order of the foreign dense column selection
        , bool SF2             // Symmetry flag of the foreign dense column selection
        , typename... CCAs2 >  // Compile time column arguments of the foreign dense column selection
inline bool
   Columns<MT,false,true,true,CCAs...>::isAliased( const Columns<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection is properly aligned in memory.
//
// \return \a true in case the dense column selection is aligned, \a false if not.
//
// This function returns whether the dense column selection is guaranteed to be properly aligned
// in memory, i.e. whether the beginning and the end of the dense column selection are guaranteed
// to conform to the alignment restrictions of the element type \a Type.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,false,true,true,CCAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column selection can be used in SMP assignments.
//
// \return \a true in case the dense column selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense column selection can be used in SMP assignments. In
// contrast to the \a smpAssignable member enumeration, which is based solely on compile time
// information, this function additionally provides runtime information (as for instance the
// current number of rows and/or columns of the dense column selection).
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,false,true,true,CCAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense column selection. The
// row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Columns<MT,false,true,true,CCAs...>::SIMDType
   Columns<MT,false,true,true,CCAs...>::load( size_t i, size_t j ) const noexcept
{
   return matrix_.load( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense column selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Columns<MT,false,true,true,CCAs...>::SIMDType
   Columns<MT,false,true,true,CCAs...>::loada( size_t i, size_t j ) const noexcept
{
   return matrix_.loada( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense column selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense column
// selection. The row index must be smaller than the number of rows and the column index must be
// smaller than the number of columns. Additionally, the row index must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Columns<MT,false,true,true,CCAs...>::SIMDType
   Columns<MT,false,true,true,CCAs...>::loadu( size_t i, size_t j ) const noexcept
{
   return matrix_.loadu( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
