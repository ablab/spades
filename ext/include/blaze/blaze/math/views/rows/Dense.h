//=================================================================================================
/*!
//  \file blaze/math/views/rows/Dense.h
//  \brief Rows specialization for dense matrices
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

#ifndef _BLAZE_MATH_VIEWS_ROWS_DENSE_H_
#define _BLAZE_MATH_VIEWS_ROWS_DENSE_H_


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
#include <blaze/math/constraints/Rows.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/RowsTrait.h>
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
#include <blaze/math/views/rows/BaseTemplate.h>
#include <blaze/math/views/rows/RowsData.h>
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
//  CLASS TEMPLATE SPECIALIZATION FOR ROW-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Rows for row selections on row-major dense matrices.
// \ingroup rows
//
// This specialization of Rows adapts the class template to the requirements of row-major
// dense matrices.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
class Rows<MT,true,true,SF,CRAs...>
   : public View< DenseMatrix< Rows<MT,true,true,SF,CRAs...>, false > >
   , private RowsData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowsData<CRAs...>;                    //!< The type of the RowsData base class.
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
   //! Type of this Rows instance.
   using This = Rows<MT,true,true,SF,CRAs...>;

   //! Base type of this Rows instance.
   using BaseType = View< DenseMatrix<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Rows instance.
   using ResultType    = RowsTrait_t<MT,N>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const Rows&;                  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant row value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant row value.
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
   template< typename... RRAs >
   explicit inline Rows( MT& matrix, RRAs... args );

   Rows( const Rows& ) = default;
   Rows( Rows&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Rows() = default;
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
   inline Pointer        data  ( size_t i ) noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Rows& operator=( const ElementType& rhs );
   inline Rows& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Rows& operator=( const Rows& rhs );

   template< typename MT2, bool SO2 >
   inline Rows& operator=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::rows;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t columns() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Rows& transpose();
   inline Rows& ctranspose();

   template< typename Other > inline Rows& scale( const Other& scalar );
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

   template< typename MT2, bool SO2, bool SF2, typename... CRAs2 >
   inline bool canAlias( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CRAs2 >
   inline bool isAliased( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

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
   inline auto assign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>&  rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the rows.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CRAs2 > friend class Rows;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROWS_TYPE        ( MT );
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
/*!\brief Constructor for row selections on row-major dense matrices.
//
// \param matrix The matrix containing the rows.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Rows<MT,true,true,SF,CRAs...>::Rows( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the rows
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( matrix_.rows() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief 2D-access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::Reference
   Rows<MT,true,true,SF,CRAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstReference
   Rows<MT,true,true,SF,CRAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::Reference
   Rows<MT,true,true,SF,CRAs...>::at( size_t i, size_t j )
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
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstReference
   Rows<MT,true,true,SF,CRAs...>::at( size_t i, size_t j ) const
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
/*!\brief Low-level data access to the selected row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row selection. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The underlying matrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::Pointer
   Rows<MT,true,true,SF,CRAs...>::data() noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row selection. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The underlying matrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstPointer
   Rows<MT,true,true,SF,CRAs...>::data() const noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::Pointer
   Rows<MT,true,true,SF,CRAs...>::data( size_t i ) noexcept
{
   return matrix_.data( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstPointer
   Rows<MT,true,true,SF,CRAs...>::data( size_t i ) const noexcept
{
   return matrix_.data( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::Iterator
   Rows<MT,true,true,SF,CRAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.begin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstIterator
   Rows<MT,true,true,SF,CRAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstIterator
   Rows<MT,true,true,SF,CRAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::Iterator
   Rows<MT,true,true,SF,CRAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.end( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstIterator
   Rows<MT,true,true,SF,CRAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cend( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,true,SF,CRAs...>::ConstIterator
   Rows<MT,true,true,SF,CRAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cend( idx(i) );
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
/*!\brief Homogenous assignment to all row selection elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row selection.
//
// This function homogeneously assigns the given value to all row selection elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::operator=( const ElementType& rhs )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      row( matrix_, idx(i), unchecked ) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all row selection elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row selection.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the row
// selection by means of an initializer list. The row selection elements are assigned the values
// from the given initializer list. Missing values are initialized as default. Note that if the
// size of the top-level initializer list does not match the number of rows of the row selection
// or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row selection" );
   }

   if( IsRestricted_v<MT> ) {
      size_t i( 0UL );
      for( const auto& rowList : list ) {
         const InitializerVector<ElementType> tmp( rowList, columns() );
         if( !tryAssign( row( matrix_, idx(i), unchecked ), tmp, 0UL ) ){
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
         ++i;
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), left.begin(i) ), left.end(i), ElementType() );
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Rows.
//
// \param rhs Dense row selection to be copied.
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense row selection is initialized as a copy of the given dense row selection. In case
// the current sizes of the two row selections don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::operator=( const Rows& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( rhs, i, unchecked ), 0UL ) ) {
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
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense row selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::operator=( const Matrix<MT2,SO2>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( right, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be added to the row selection.
// \return Reference to the dense row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,true,true,SF,CRAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAddAssign( row( matrix_, idx(i), unchecked ), row( *rhs, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be added to the row selection.
// \return Reference to the dense row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,true,true,SF,CRAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( tmp, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be subtracted from the row selection.
// \return Reference to the dense row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,true,true,SF,CRAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !trySubAssign( row( matrix_, idx(i), unchecked ), row( *rhs, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be subtracted from the row selection.
// \return Reference to the dense row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,true,true,SF,CRAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( tmp, i, unchecked ), 0UL ) ) {
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
// \return Reference to the dense row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,true,true,SF,CRAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryMultAssign( row( matrix_, idx(i), unchecked ), row( *rhs, i, unchecked ), 0UL ) ) {
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
// \return Reference to the dense row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,true,true,SF,CRAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( tmp, i, unchecked ), 0UL ) ) {
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
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline MT& Rows<MT,true,true,SF,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline const MT& Rows<MT,true,true,SF,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the row selection.
//
// \return The number of columns of the row selection.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,true,SF,CRAs...>::columns() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two rows.
//
// \return The spacing between the beginning of two rows.
//
// This function returns the spacing between the beginning of two rows, i.e. the total number of
// elements of a row.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,true,SF,CRAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row selection.
//
// \return The capacity of the dense row selection.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,true,SF,CRAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,true,SF,CRAs...>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense row selection.
//
// \return The number of non-zero elements in the dense row selection.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,true,SF,CRAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i ) {
      nonzeros += matrix_.nonZeros( idx(i) );
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,true,SF,CRAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.nonZeros( idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,true,SF,CRAs...>::reset()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      matrix_.reset( idx(i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row to the default initial values.
//
// \param i The index of the row.
// \return void
//
// This function resets the values in the specified row to their default value.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,true,SF,CRAs...>::reset( size_t i )
{
   matrix_.reset( idx(i) );
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
/*!\brief In-place transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i ), idx(i), 0UL ) ) {
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
/*!\brief In-place conjugate transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i ), idx(i), 0UL ) ) {
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
/*!\brief Scaling of the dense row selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the dense row selection.
//
// This function scales the row selection by applying the given scalar value \a scalar to each
// element of the row selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to
// scale a row selection on a lower or upper unitriangular matrix. The attempt to scale such a
// row selection results in a compile time error!
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the scalar value
inline Rows<MT,true,true,SF,CRAs...>&
   Rows<MT,true,true,SF,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index ( idx(i) );
      const size_t jbegin( IsUpper<MT>::value ? ( IsStrictlyUpper_v<MT> ? index+1UL : index ) : 0UL );
      const size_t jend  ( IsLower<MT>::value ? ( IsStrictlyLower_v<MT> ? index : index+1UL ) : columns() );

      for( size_t j=jbegin; j<jend; ++j ) {
         matrix_(index,j) *= scalar;
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
/*!\brief Returns whether the dense row selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address can alias with the dense row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,true,true,SF,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection can alias with the given dense row selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address can alias with the dense row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , bool SF              // Symmetry flag
        , typename... CRAs >   // Compile time row arguments
template< typename MT2         // Data type of the foreign dense row selection
        , bool SO2             // Storage order of the foreign dense row selection
        , bool SF2             // Symmetry flag of the foreign dense row selection
        , typename... CRAs2 >  // Compile time row arguments of the foreign dense row selection
inline bool
   Rows<MT,true,true,SF,CRAs...>::canAlias( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense row selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,true,true,SF,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is aliased with the given dense row selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense row selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , bool SF              // Symmetry flag
        , typename... CRAs >   // Compile time row arguments
template< typename MT2         // Data type of the foreign dense row selection
        , bool SO2             // Storage order of the foreign dense row selection
        , bool SF2             // Symmetry flag of the foreign dense row selection
        , typename... CRAs2 >  // Compile time row arguments of the foreign dense row selection
inline bool
   Rows<MT,true,true,SF,CRAs...>::isAliased( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is properly aligned in memory.
//
// \return \a true in case the dense row selection is aligned, \a false if not.
//
// This function returns whether the dense row selection is guaranteed to be properly aligned in
// memory, i.e. whether the beginning and the end of the dense row selection are guaranteed to
// conform to the alignment restrictions of the element type \a Type.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,true,true,SF,CRAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection can be used in SMP assignments.
//
// \return \a true in case the dense row selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row selection can be used in SMP assignments. In
// contrast to the \a smpAssignable member enumeration, which is based solely on compile time
// information, this function additionally provides runtime information (as for instance the
// current number of rows and/or columns of the dense row selection).
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,true,true,SF,CRAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense row selection. The row
// index must be smaller than the number of rows and the column index must be smaller than the
// number of columns. Additionally, the column index must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Rows<MT,true,true,SF,CRAs...>::SIMDType
   Rows<MT,true,true,SF,CRAs...>::load( size_t i, size_t j ) const noexcept
{
   return matrix_.load( idx(i), j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense row selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Rows<MT,true,true,SF,CRAs...>::SIMDType
   Rows<MT,true,true,SF,CRAs...>::loada( size_t i, size_t j ) const noexcept
{
   return matrix_.loada( idx(i), j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense row selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Rows<MT,true,true,SF,CRAs...>::SIMDType
   Rows<MT,true,true,SF,CRAs...>::loadu( size_t i, size_t j ) const noexcept
{
   return matrix_.loadu( idx(i), j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense row selection. The row
// index must be smaller than the number of rows and the column index must be smaller than the
// number of columns. Additionally, the column index must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Rows<MT,true,true,SF,CRAs...>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.store( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense row selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Rows<MT,true,true,SF,CRAs...>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.storea( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense row selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Rows<MT,true,true,SF,CRAs...>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.storeu( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the dense
// row selection. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the column index must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Rows<MT,true,true,SF,CRAs...>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   matrix_.stream( idx(i), j, value );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::assign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(index,j    ) = (*rhs)(i,j    );
         matrix_(index,j+1UL) = (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(index,jpos) = (*rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a row-major dense matrix.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::assign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   if( useStreaming &&
       rows()*columns() > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<rows(); ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT2> right( (*rhs).begin(i) );

         for( ; j<jpos; j+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<columns(); ++j ) {
            *left = *right;
         }
      }
   }
   else
   {
      for( size_t i=0UL; i<rows(); ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT2> right( (*rhs).begin(i) );

         for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<jpos; j+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<columns(); ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,true,true,SF,CRAs...>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) = (*rhs)(i,j    );
            matrix_(index,j+1UL) = (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() ) {
            matrix_(index,jpos) = (*rhs)(i,jpos);
         }
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) = (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(index,element->index()) = element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(idx(element->index()),j) = element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );
      if( IsDiagonal_v<MT2> ) {
         matrix_(index,i) += (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) += (*rhs)(i,j    );
            matrix_(index,j+1UL) += (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(index,jpos) += (*rhs)(i,jpos);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a row-major dense matrix.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT2> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT2> )
                           ?( IsStrictlyLower_v<MT2> ? i : i+1UL )
                           :( columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( prevMultiple( jend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT2> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jend; ++j ) {
         *left += *right; ++left; ++right;
      }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,true,true,SF,CRAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) += (*rhs)(i,j    );
            matrix_(index,j+1UL) += (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() ) {
            matrix_(index,jpos) += (*rhs)(i,jpos);
         }
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) += (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(index,element->index()) += element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(idx(element->index()),j) += element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );

      if( IsDiagonal_v<MT2> ) {
         matrix_(index,i) -= (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) -= (*rhs)(i,j    );
            matrix_(index,j+1UL) -= (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(index,jpos) -= (*rhs)(i,jpos);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a row-major dense matrix.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT2> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT2> )
                           ?( IsStrictlyLower_v<MT2> ? i : i+1UL )
                           :( columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( prevMultiple( jend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT2> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jend; ++j ) {
         *left -= *right; ++left; ++right;
      }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,true,true,SF,CRAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) -= (*rhs)(i,j    );
            matrix_(index,j+1UL) -= (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() ) {
            matrix_(index,jpos) -= (*rhs)(i,jpos);
         }
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) -= (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(index,element->index()) -= element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(idx(element->index()),j) -= element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(index,j    ) *= (*rhs)(i,j    );
         matrix_(index,j+1UL) *= (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(index,jpos) *= (*rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a row-major dense matrix.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline auto Rows<MT,true,true,SF,CRAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jpos( prevMultiple( columns(), SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

      size_t j( 0UL );
      Iterator left( begin(i) );
      ConstIterator_t<MT2> right( (*rhs).begin(i) );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<columns(); ++j ) {
         *left *= *right; ++left; ++right;
      }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,true,true,SF,CRAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) *= (*rhs)(i,j    );
            matrix_(index,j+1UL) *= (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() ) {
            matrix_(index,jpos) *= (*rhs)(i,jpos);
         }
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) *= (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(index,j) );
         matrix_(index,j) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(index,j) );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,true,SF,CRAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(idx(i),j) );
         matrix_(idx(i),j) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(idx(i),j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL COLUMN-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Rows for row selections on general column-major dense matrices.
// \ingroup rows
//
// This specialization of Rows adapts the class template to the requirements of general
// column-major dense matrices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
class Rows<MT,false,true,false,CRAs...>
   : public View< DenseMatrix< Rows<MT,false,true,false,CRAs...>, false > >
   , private RowsData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowsData<CRAs...>;                    //!< The type of the RowsData base class.
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
   //! Type of this Rows instance.
   using This = Rows<MT,false,true,false,CRAs...>;

   //! Base type of this Rows instance.
   using BaseType = View< DenseMatrix<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Rows instance.
   using ResultType    = RowsTrait_t<MT,N>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const Rows&;                  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant row value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant row value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**RowsIterator class definition***************************************************************
   /*!\brief Iterator over the elements of a selected dense row.
   */
   template< typename MatrixType      // Type of the dense matrix
           , typename IteratorType >  // Type of the dense matrix iterator
   class RowsIterator
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
      /*!\brief Default constructor of the RowsIterator class.
      */
      inline RowsIterator() noexcept
         : matrix_( nullptr )  // The dense matrix containing the row
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   (     )      // Iterator to the current dense element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the RowsIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      */
      inline RowsIterator( MatrixType& matrix, size_t row, size_t column ) noexcept
         : matrix_( &matrix )  // The dense matrix containing the selected row
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   (         )  // Iterator to the current dense element
      {
         if( column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different RowsIterator instances.
      //
      // \param it The row iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline RowsIterator( const RowsIterator<MatrixType2,IteratorType2>& it ) noexcept
         : matrix_( it.matrix_ )  // The dense matrix containing the seleted row
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
      inline RowsIterator& operator+=( size_t inc ) noexcept {
         using blaze::reset;
         column_ += inc;
         if( column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
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
      inline RowsIterator& operator-=( size_t dec ) noexcept {
         using blaze::reset;
         column_ -= dec;
         if( column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline RowsIterator& operator++() noexcept {
         using blaze::reset;
         ++column_;
         if( column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const RowsIterator operator++( int ) noexcept {
         const RowsIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline RowsIterator& operator--() noexcept {
         using blaze::reset;
         --column_;
         if( column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const RowsIterator operator--( int ) noexcept {
         const RowsIterator tmp( *this );
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
         BLAZE_USER_ASSERT( column_+index < matrix_->columns(), "Invalid access index detected" );
         const IteratorType pos( matrix_->begin( column_+index ) + row_ );
         return *pos;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense row element at the current iterator position.
      //
      // \return Reference to the current value.
      */
      inline ReferenceType operator*() const {
         return *pos_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense row element at the current iterator position.
      //
      // \return Pointer to the dense row element at the current iterator position.
      */
      inline PointerType operator->() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ == rhs.column_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ < rhs.column_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ > rhs.column_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<=( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ <= rhs.column_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>=( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ >= rhs.column_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two row iterators.
      //
      // \param rhs The right-hand side row iterator.
      // \return The number of elements between the two row iterators.
      */
      inline DifferenceType operator-( const RowsIterator& rhs ) const noexcept {
         return column_ - rhs.column_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a RowsIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const RowsIterator operator+( const RowsIterator& it, size_t inc ) noexcept {
         return RowsIterator( *it.matrix_, it.row_, it.column_+inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a RowsIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const RowsIterator operator+( size_t inc, const RowsIterator& it ) noexcept {
         return RowsIterator( *it.matrix_, it.row_, it.column_+inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a RowsIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const RowsIterator operator-( const RowsIterator& it, size_t dec ) noexcept {
         return RowsIterator( *it.matrix_, it.row_, it.column_-dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The dense matrix containing the selected row.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current dense element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class RowsIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = RowsIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, RowsIterator< MT, Iterator_t<MT> > >;
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
   template< typename... RRAs >
   explicit inline Rows( MT& matrix, RRAs... args );

   Rows( const Rows& ) = default;
   Rows( Rows&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Rows() = default;
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
   inline Pointer        data  ( size_t i ) noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Rows& operator=( const ElementType& rhs );
   inline Rows& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Rows& operator=( const Rows& rhs );

   template< typename MT2, bool SO2 >
   inline Rows& operator=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::rows;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t columns() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Rows& transpose();
   inline Rows& ctranspose();

   template< typename Other > inline Rows& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CRAs2 >
   inline bool canAlias( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CRAs2 >
   inline bool isAliased( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,true>&   rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>&  rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the rows.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CRAs2 > friend class Rows;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROWS_TYPE            ( MT );
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
/*!\brief Constructor for row selections on general column-major dense matrices.
//
// \param matrix The matrix containing the rows.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Rows<MT,false,true,false,CRAs...>::Rows( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the rows
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( matrix_.rows() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::Reference
   Rows<MT,false,true,false,CRAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstReference
   Rows<MT,false,true,false,CRAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::Reference
   Rows<MT,false,true,false,CRAs...>::at( size_t i, size_t j )
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
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstReference
   Rows<MT,false,true,false,CRAs...>::at( size_t i, size_t j ) const
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
/*!\brief Low-level data access to the selected row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row selection. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The underlying matrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::Pointer
   Rows<MT,false,true,false,CRAs...>::data() noexcept
{
   return matrix_.data() + idx(0UL);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row selection. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The underlying matrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstPointer
   Rows<MT,false,true,false,CRAs...>::data() const noexcept
{
   return matrix_.data() + idx(0UL);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::Pointer
   Rows<MT,false,true,false,CRAs...>::data( size_t i ) noexcept
{
   return matrix_.data() + idx(i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstPointer
   Rows<MT,false,true,false,CRAs...>::data( size_t i ) const noexcept
{
   return matrix_.data() + idx(i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::Iterator
   Rows<MT,false,true,false,CRAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return Iterator( matrix_, idx(i), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstIterator
   Rows<MT,false,true,false,CRAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return ConstIterator( matrix_, idx(i), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstIterator
   Rows<MT,false,true,false,CRAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return ConstIterator( matrix_, idx(i), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::Iterator
   Rows<MT,false,true,false,CRAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return Iterator( matrix_, idx(i), columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstIterator
   Rows<MT,false,true,false,CRAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return ConstIterator( matrix_, idx(i), columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,false,CRAs...>::ConstIterator
   Rows<MT,false,true,false,CRAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return ConstIterator( matrix_, idx(i), columns() );
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
/*!\brief Homogenous assignment to all row selection elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row selection.
//
// This function homogeneously assigns the given value to all row selection elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::operator=( const ElementType& rhs )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      row( matrix_, idx(i), unchecked ) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all row selection elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid initializer list dimension.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the row
// selection by means of an initializer list. The row selection elements are assigned the values
// from the given initializer list. Missing values are initialized as default. Note that if the
// size of the top-level initializer list does not match the number of rows of the row selection
// or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row selection" );
   }

   if( IsRestricted_v<MT> ) {
      size_t i( 0UL );
      for( const auto& rowList : list ) {
         const InitializerVector<ElementType> tmp( rowList, columns() );
         if( !tryAssign( row( matrix_, idx(i), unchecked ), tmp, 0UL ) ){
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
         ++i;
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), left.begin(i) ), left.end(i), ElementType() );
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Rows.
//
// \param rhs Dense row selection to be copied.
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense row selection is initialized as a copy of the given dense row selection. In case
// the current sizes of the two row selections don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::operator=( const Rows& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( rhs, i, unchecked ), 0UL ) ) {
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
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense row selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::operator=( const Matrix<MT2,SO2>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( right, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be added to the row selection.
// \return Reference to the dense row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,false,true,false,CRAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAddAssign( row( matrix_, idx(i), unchecked ), row( *rhs, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be added to the row selection.
// \return Reference to the dense row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,false,true,false,CRAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( tmp, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be subtracted from the row selection.
// \return Reference to the dense row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,false,true,false,CRAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !trySubAssign( row( matrix_, idx(i), unchecked ), row( *rhs, i, unchecked ), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be subtracted from the row selection.
// \return Reference to the dense row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,false,true,false,CRAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( tmp, i, unchecked ), 0UL ) ) {
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
// \return Reference to the dense row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,false,true,false,CRAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryMultAssign( row( matrix_, idx(i), unchecked ), row( *rhs, i, unchecked ), 0UL ) ) {
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
// \return Reference to the dense row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO2 >          // Storage order of the right-hand side matrix
inline auto Rows<MT,false,true,false,CRAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Rows& >
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( row( matrix_, idx(i), unchecked ), row( tmp, i, unchecked ), 0UL ) ) {
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
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline MT& Rows<MT,false,true,false,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline const MT& Rows<MT,false,true,false,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the row selection.
//
// \return The number of columns of the row selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,false,CRAs...>::columns() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two rows.
//
// \return The spacing between the beginning of two rows.
//
// This function returns the spacing between the beginning of two rows, i.e. the total number of
// elements of a row.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,false,CRAs...>::spacing() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row selection.
//
// \return The capacity of the dense row selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,false,CRAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,false,CRAs...>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense row selection.
//
// \return The number of non-zero elements in the dense row selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,false,CRAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i ) {
      nonzeros += nonZeros( i );
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,false,CRAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   size_t nonzeros( 0UL );

   const size_t index( idx(i) );
   for( size_t j=0UL; j<columns(); ++j ) {
      if( !isDefault( matrix_( index, j ) ) )
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
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,true,false,CRAs...>::reset()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      reset( i );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row to the default initial values.
//
// \param i The index of the row.
// \return void
//
// This function resets the values in the specified row to their default value.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,true,false,CRAs...>::reset( size_t i )
{
   using blaze::reset;

   const size_t index( idx(i) );
   for( size_t j=0UL; j<columns(); ++j ) {
      reset( matrix_( index, j ) );
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
/*!\brief In-place transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i ), idx(i), 0UL ) ) {
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
/*!\brief In-place conjugate transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i ), idx(i), 0UL ) ) {
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
/*!\brief Scaling of the dense row selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the dense row selection.
//
// This function scales the row selection by applying the given scalar value \a scalar to each
// element of the row selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to
// scale a row selection on a lower or upper unitriangular matrix. The attempt to scale such a
// row selection results in a compile time error!
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the scalar value
inline Rows<MT,false,true,false,CRAs...>&
   Rows<MT,false,true,false,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index ( idx(i) );
      const size_t jbegin( IsUpper<MT>::value ? ( IsStrictlyUpper_v<MT> ? index+1UL : index ) : 0UL );
      const size_t jend  ( IsLower<MT>::value ? ( IsStrictlyLower_v<MT> ? index : index+1UL ) : columns() );

      for( size_t j=jbegin; j<jend; ++j ) {
         matrix_(index,j) *= scalar;
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
/*!\brief Returns whether the dense row selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address can alias with the dense row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,true,false,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection can alias with the given dense row selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address can alias with the dense row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CRAs >   // Compile time row arguments
template< typename MT2         // Data type of the foreign dense row selection
        , bool SO2             // Storage order of the foreign dense row selection
        , bool SF2             // Symmetry flag of the foreign dense row selection
        , typename... CRAs2 >  // Compile time row arguments of the foreign dense row selection
inline bool
   Rows<MT,false,true,false,CRAs...>::canAlias( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense row selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,true,false,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is aliased with the given dense row selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense row selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CRAs >   // Compile time row arguments
template< typename MT2         // Data type of the foreign dense row selection
        , bool SO2             // Storage order of the foreign dense row selection
        , bool SF2             // Symmetry flag of the foreign dense row selection
        , typename... CRAs2 >  // Compile time row arguments of the foreign dense row selection
inline bool
   Rows<MT,false,true,false,CRAs...>::isAliased( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is properly aligned in memory.
//
// \return \a true in case the dense row selection is aligned, \a false if not.
//
// This function returns whether the dense row selection is guaranteed to be properly aligned in
// memory, i.e. whether the beginning and the end of the dense row selection are guaranteed to
// conform to the alignment restrictions of the element type \a Type.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,false,true,false,CRAs...>::isAligned() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection can be used in SMP assignments.
//
// \return \a true in case the dense row selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row selection can be used in SMP assignments. In
// contrast to the \a smpAssignable member enumeration, which is based solely on compile time
// information, this function additionally provides runtime information (as for instance the
// current number of rows and/or columns of the dense row selection).
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,false,true,false,CRAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() > SMP_DMATASSIGN_THRESHOLD );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(index,j    ) = (*rhs)(i,j    );
         matrix_(index,j+1UL) = (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(index,jpos) = (*rhs)(i,jpos);
      }
   }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) = (*rhs)(i,j    );
            matrix_(index,j+1UL) = (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() )
            matrix_(index,jpos) = (*rhs)(i,jpos);
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) = (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(index,element->index()) = element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(idx(element->index()),j) = element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );
      if( IsDiagonal_v<MT2> ) {
         matrix_(index,i) += (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) += (*rhs)(i,j    );
            matrix_(index,j+1UL) += (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(index,jpos) += (*rhs)(i,jpos);
         }
      }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) += (*rhs)(i,j    );
            matrix_(index,j+1UL) += (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() )
            matrix_(index,jpos) += (*rhs)(i,jpos);
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) += (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(index,element->index()) += element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(idx(element->index()),j) += element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );

      if( IsDiagonal_v<MT2> ) {
         matrix_(index,i) -= (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) -= (*rhs)(i,j    );
            matrix_(index,j+1UL) -= (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(index,jpos) -= (*rhs)(i,jpos);
         }
      }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) -= (*rhs)(i,j    );
            matrix_(index,j+1UL) -= (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() )
            matrix_(index,jpos) -= (*rhs)(i,jpos);
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) -= (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(index,element->index()) -= element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(idx(element->index()),j) -= element->value();
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(index,j    ) *= (*rhs)(i,j    );
         matrix_(index,j+1UL) *= (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(index,jpos) *= (*rhs)(i,jpos);
      }
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side dense matrix
inline void Rows<MT,false,true,false,CRAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   if( rows() < block && columns() < block )
   {
      const size_t jpos( prevMultiple( (*rhs).columns(), 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).columns(), "Invalid end calculation" );

      for( size_t i=0UL; i<rows(); ++i ) {
         const size_t index( idx(i) );
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(index,j    ) *= (*rhs)(i,j    );
            matrix_(index,j+1UL) *= (*rhs)(i,j+1UL);
         }
         if( jpos < (*rhs).columns() )
            matrix_(index,jpos) *= (*rhs)(i,jpos);
      }
   }
   else
   {
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t jj=0UL; jj<columns(); jj+=block ) {
            const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               const size_t index( idx(i) );
               for( size_t j=jj; j<jend; ++j ) {
                  matrix_(index,j) *= (*rhs)(i,j);
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(index,j) );
         matrix_(index,j) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(index,j) );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,true,false,CRAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(idx(i),j) );
         matrix_(idx(i),j) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(idx(i),j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC COLUMN-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Rows for row selections on symmetric column-major dense matrices.
// \ingroup rows
//
// This specialization of Rows adapts the class template to the requirements of symmetric
// column-major dense matrices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
class Rows<MT,false,true,true,CRAs...>
   : public View< DenseMatrix< Rows<MT,false,true,true,CRAs...>, false > >
   , private RowsData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowsData<CRAs...>;                    //!< The type of the RowsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Rows instance.
   using This = Rows<MT,false,true,true,CRAs...>;

   //! Base type of this Rows instance.
   using BaseType = View< DenseMatrix<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Rows instance.
   using ResultType    = RowsTrait_t<MT,N>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const Rows&;                  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant row value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant row value.
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
   template< typename... RRAs >
   explicit inline Rows( MT& matrix, RRAs... args );

   Rows( const Rows& ) = default;
   Rows( Rows&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Rows() = default;
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
   inline Pointer        data  ( size_t i ) noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Rows& operator=( const ElementType& rhs );

   Rows& operator=( const Rows& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::rows;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t columns() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CRAs2 >
   inline bool canAlias( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, typename... CRAs2 >
   inline bool isAliased( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

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
   Operand matrix_;  //!< The matrix containing the rows.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CRAs2 > friend class Rows;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROWS_TYPE           ( MT );
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
/*!\brief Constructor for row selections on symmetric column-major dense matrices.
//
// \param matrix The matrix containing the rows.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Rows<MT,false,true,true,CRAs...>::Rows( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the rows
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( matrix_.rows() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::Reference
   Rows<MT,false,true,true,CRAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(j,idx(i));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstReference
   Rows<MT,false,true,true,CRAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(j,idx(i));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::Reference
   Rows<MT,false,true,true,CRAs...>::at( size_t i, size_t j )
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
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstReference
   Rows<MT,false,true,true,CRAs...>::at( size_t i, size_t j ) const
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
/*!\brief Low-level data access to the selected row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row selection. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The underlying matrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::Pointer
   Rows<MT,false,true,true,CRAs...>::data() noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the selected row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row selection. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The underlying matrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstPointer
   Rows<MT,false,true,true,CRAs...>::data() const noexcept
{
   return matrix_.data( idx(0UL) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::Pointer
   Rows<MT,false,true,true,CRAs...>::data( size_t i ) noexcept
{
   return matrix_.data( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstPointer
   Rows<MT,false,true,true,CRAs...>::data( size_t i ) const noexcept
{
   return matrix_.data( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::Iterator
   Rows<MT,false,true,true,CRAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.begin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstIterator
   Rows<MT,false,true,true,CRAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns an iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstIterator
   Rows<MT,false,true,true,CRAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::Iterator
   Rows<MT,false,true,true,CRAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.end( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstIterator
   Rows<MT,false,true,true,CRAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cend( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,true,true,CRAs...>::ConstIterator
   Rows<MT,false,true,true,CRAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return matrix_.cend( idx(i) );
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
/*!\brief Homogenous assignment to all row selection elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row selection.
//
// This function homogeneously assigns the given value to all row selection elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,true,true,CRAs...>&
   Rows<MT,false,true,true,CRAs...>::operator=( const ElementType& rhs )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      column( matrix_, idx(i), unchecked ) = rhs;
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
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline MT& Rows<MT,false,true,true,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline const MT& Rows<MT,false,true,true,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the row selection.
//
// \return The number of columns of the row selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,true,CRAs...>::columns() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two rows.
//
// \return The spacing between the beginning of two rows.
//
// This function returns the spacing between the beginning of two rows, i.e. the total number of
// elements of a row.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,true,CRAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row selection.
//
// \return The capacity of the dense row selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,true,CRAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,true,CRAs...>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense row selection.
//
// \return The number of non-zero elements in the dense row selection.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,true,CRAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i ) {
      nonzeros += matrix_.nonZeros( idx(i) );
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,true,true,CRAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.nonZeros( idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,true,true,CRAs...>::reset()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      matrix_.reset( idx(i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row to the default initial values.
//
// \param i The index of the row.
// \return void
//
// This function resets the values in the specified row to their default value.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,true,true,CRAs...>::reset( size_t i )
{
   matrix_.reset( idx(i) );
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
/*!\brief Returns whether the dense row selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address can alias with the dense row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,true,true,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection can alias with the given dense row selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address can alias with the dense row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CRAs >   // Compile time row arguments
template< typename MT2         // Data type of the foreign dense row selection
        , bool SO2             // Storage order of the foreign dense row selection
        , bool SF2             // Symmetry flag of the foreign dense row selection
        , typename... CRAs2 >  // Compile time row arguments of the foreign dense row selection
inline bool
   Rows<MT,false,true,true,CRAs...>::canAlias( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense row selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,true,true,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is aliased with the given dense row selection
//        \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row selection, \a false if not.
//
// This function returns whether the given address is aliased with the dense row selection. In
// contrast to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , typename... CRAs >   // Compile time row arguments
template< typename MT2         // Data type of the foreign dense row selection
        , bool SO2             // Storage order of the foreign dense row selection
        , bool SF2             // Symmetry flag of the foreign dense row selection
        , typename... CRAs2 >  // Compile time row arguments of the foreign dense row selection
inline bool
   Rows<MT,false,true,true,CRAs...>::isAliased( const Rows<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection is properly aligned in memory.
//
// \return \a true in case the dense row selection is aligned, \a false if not.
//
// This function returns whether the dense row selection is guaranteed to be properly aligned in
// memory, i.e. whether the beginning and the end of the dense row selection are guaranteed to
// conform to the alignment restrictions of the element type \a Type.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,false,true,true,CRAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row selection can be used in SMP assignments.
//
// \return \a true in case the dense row selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row selection can be used in SMP assignments. In
// contrast to the \a smpAssignable member enumeration, which is based solely on compile time
// information, this function additionally provides runtime information (as for instance the
// current number of rows and/or columns of the dense row selection).
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,false,true,true,CRAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense row selection. The row
// index must be smaller than the number of rows and the column index must be smaller than the
// number of columns. Additionally, the column index must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Rows<MT,false,true,true,CRAs...>::SIMDType
   Rows<MT,false,true,true,CRAs...>::load( size_t i, size_t j ) const noexcept
{
   return matrix_.load( j, idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense row selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Rows<MT,false,true,true,CRAs...>::SIMDType
   Rows<MT,false,true,true,CRAs...>::loada( size_t i, size_t j ) const noexcept
{
   return matrix_.loada( j, idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense row selection.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense row selection.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the dense matrix
        , typename... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Rows<MT,false,true,true,CRAs...>::SIMDType
   Rows<MT,false,true,true,CRAs...>::loadu( size_t i, size_t j ) const noexcept
{
   return matrix_.loadu( j, idx(i) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
