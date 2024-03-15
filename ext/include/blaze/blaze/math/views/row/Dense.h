//=================================================================================================
/*!
//  \file blaze/math/views/row/Dense.h
//  \brief Row specialization for dense matrices
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

#ifndef _BLAZE_MATH_VIEWS_ROW_DENSE_H_
#define _BLAZE_MATH_VIEWS_ROW_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/row/BaseTemplate.h>
#include <blaze/math/views/row/RowData.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/EnableIf.h>
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
/*!\brief Specialization of Row for rows on row-major dense matrices.
// \ingroup row
//
// This specialization of Row adapts the class template to the requirements of row-major
// dense matrices.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
class Row<MT,true,true,SF,CRAs...>
   : public View< DenseVector< Row<MT,true,true,SF,CRAs...>, true > >
   , private RowData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowData<CRAs...>;                     //!< The type of the RowData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Row instance.
   using This = Row<MT,true,true,SF,CRAs...>;

   //! Base type of this Row instance.
   using BaseType = View< DenseVector<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Row instance.
   using ResultType    = RowTrait_t<MT,CRAs...>;       //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Row&;                   //!< Data type for composite expression templates.

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
   explicit inline Row( MT& matrix, RRAs... args );

   Row( const Row& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Row() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Row& operator=( const ElementType& rhs );
   inline Row& operator=( initializer_list<ElementType> list );
   inline Row& operator=( const Row& rhs );

   template< typename VT > inline Row& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator*=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator/=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator%=( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t  size() const noexcept;
   inline size_t  spacing() const noexcept;
   inline size_t  capacity() const noexcept;
   inline size_t  nonZeros() const;
   inline void    reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Row& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDAdd_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDSub_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDMult_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDDiv_v< ElementType, ElementType_t<VT> > );
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

   template< typename MT2, bool SO2, bool SF2, size_t... CRAs2 >
   inline bool canAlias( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CRAs2 >
   inline bool isAliased( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT >
   inline auto assign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CRAs2 > friend class Row;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
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
/*!\brief Constructor for rows on row-major dense matrices.
//
// \param matrix The matrix containing the row.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , size_t... CRAs >    // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Row<MT,true,true,SF,CRAs...>::Row( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the row
{
   if( isChecked( args... ) ) {
      if( matrix_.rows() <= row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row() < matrix_.rows(), "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::Reference
   Row<MT,true,true,SF,CRAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstReference
   Row<MT,true,true,SF,CRAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::Reference
   Row<MT,true,true,SF,CRAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstReference
   Row<MT,true,true,SF,CRAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::Pointer
   Row<MT,true,true,SF,CRAs...>::data() noexcept
{
   return matrix_.data( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstPointer
   Row<MT,true,true,SF,CRAs...>::data() const noexcept
{
   return matrix_.data( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::Iterator
   Row<MT,true,true,SF,CRAs...>::begin()
{
   return matrix_.begin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstIterator
   Row<MT,true,true,SF,CRAs...>::begin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstIterator
   Row<MT,true,true,SF,CRAs...>::cbegin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::Iterator
   Row<MT,true,true,SF,CRAs...>::end()
{
   return matrix_.end( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstIterator
   Row<MT,true,true,SF,CRAs...>::end() const
{
   return matrix_.cend( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,true,SF,CRAs...>::ConstIterator
   Row<MT,true,true,SF,CRAs...>::cend() const
{
   return matrix_.cend( row() );
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
/*!\brief Homogenous assignment to all row elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row.
//
// This function homogeneously assigns the given value to all elements of the row. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, row(), j, rhs ) )
         left(row(),j) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all row elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// row by means of an initializer list. The row elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the row, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate
// an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,true> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );

   std::fill( std::copy( list.begin(), list.end(), left.begin() ), left.end(), ElementType() );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Row.
//
// \param rhs Dense row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator=( const Row& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsExpression_v<MT> && rhs.canAlias( this ) ) {
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
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector_v<VT> )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator*=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator/=( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpDivAssign( left, tmp );
   }
   else {
      smpDivAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the assigned row.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::operator%=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, right );

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
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline MT& Row<MT,true,true,SF,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline const MT& Row<MT,true,true,SF,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the row.
//
// \return The size of the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,true,SF,CRAs...>::size() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the row.
//
// \return The minimum capacity of the row.
//
// This function returns the minimum capacity of the row, which corresponds to the current size
// plus padding.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,true,SF,CRAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row.
//
// \return The maximum capacity of the dense row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,true,SF,CRAs...>::capacity() const noexcept
{
   return matrix_.capacity( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,true,SF,CRAs...>::nonZeros() const
{
   return matrix_.nonZeros( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,true,true,SF,CRAs...>::reset()
{
   matrix_.reset( row() );
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
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
//
// This function scales the row by applying the given scalar value \a scalar to each element
// of the row. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a row
// on a lower or upper unitriangular matrix. The attempt to scale such a row results in a
// compile time error!
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the scalar value
inline Row<MT,true,true,SF,CRAs...>&
   Row<MT,true,true,SF,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsStrictlyUpper_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsStrictlyLower_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_(row(),j) *= scalar;
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
/*!\brief Returns whether the dense row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,true,true,SF,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can alias with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , bool SF            // Symmetry flag
        , size_t... CRAs >   // Compile time row arguments
template< typename MT2       // Data type of the foreign dense row
        , bool SO2           // Storage order of the foreign dense row
        , bool SF2           // Symmetry flag of the foreign dense row
        , size_t... CRAs2 >  // Compile time row arguments of the foreign dense row
inline bool
   Row<MT,true,true,SF,CRAs...>::canAlias( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row() == alias->row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,true,true,SF,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , bool SF            // Symmetry flag
        , size_t... CRAs >   // Compile time row arguments
template< typename MT2       // Data type of the foreign dense row
        , bool SO2           // Storage order of the foreign dense row
        , bool SF2           // Symmetry flag of the foreign dense row
        , size_t... CRAs2 >  // Compile time row arguments of the foreign dense row
inline bool
   Row<MT,true,true,SF,CRAs...>::isAliased( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row() == alias->row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is properly aligned in memory.
//
// \return \a true in case the dense row is aligned, \a false if not.
//
// This function returns whether the dense row is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense row are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool Row<MT,true,true,SF,CRAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can be used in SMP assignments.
//
// \return \a true in case the dense row can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense row).
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool Row<MT,true,true,SF,CRAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Row<MT,true,true,SF,CRAs...>::SIMDType
   Row<MT,true,true,SF,CRAs...>::load( size_t index ) const noexcept
{
   return matrix_.load( row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Row<MT,true,true,SF,CRAs...>::SIMDType
   Row<MT,true,true,SF,CRAs...>::loada( size_t index ) const noexcept
{
   return matrix_.loada( row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Row<MT,true,true,SF,CRAs...>::SIMDType
   Row<MT,true,true,SF,CRAs...>::loadu( size_t index ) const noexcept
{
   return matrix_.loadu( row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store a specific SIMD element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,true,true,SF,CRAs...>::store( size_t index, const SIMDType& value ) noexcept
{
   matrix_.store( row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense row. The
// index must be smaller than the number of matrix columns. This function must \b NOT be
// called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,true,true,SF,CRAs...>::storea( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storea( row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unligned store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store a specific SIMD element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,true,true,SF,CRAs...>::storeu( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storeu( row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific SIMD element of the dense
// row. The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,true,true,SF,CRAs...>::stream( size_t index, const SIMDType& value ) noexcept
{
   matrix_.stream( row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::assign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) = (*rhs)[j    ];
      matrix_(row(),j+1UL) = (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) = (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::assign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t columns( size() );

   const size_t jpos( remainder ? prevMultiple( columns, SIMDSIZE ) : columns );
   BLAZE_INTERNAL_ASSERT( jpos <= columns, "Invalid end calculation" );

   size_t j( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   if( useStreaming && columns > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(*rhs).isAliased( this ) )
   {
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && j<columns; ++j ) {
         *left = *right; ++left; ++right;
      }
   }
   else
   {
      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && j<columns; ++j ) {
         *left = *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,true,SF,CRAs...>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(row(),element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) += (*rhs)[j    ];
      matrix_(row(),j+1UL) += (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) += (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t columns( size() );

   const size_t jpos( remainder ? prevMultiple( columns, SIMDSIZE ) : columns );
   BLAZE_INTERNAL_ASSERT( jpos <= columns, "Invalid end calculation" );

   size_t j( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; j<jpos; j+=SIMDSIZE ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && j<columns; ++j ) {
      *left += *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,true,SF,CRAs...>::addAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(row(),element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) -= (*rhs)[j    ];
      matrix_(row(),j+1UL) -= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) -= (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t columns( size() );

   const size_t jpos( remainder ? prevMultiple( columns, SIMDSIZE ) : columns );
   BLAZE_INTERNAL_ASSERT( jpos <= columns, "Invalid end calculation" );

   size_t j( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; j<jpos; j+=SIMDSIZE ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && j<columns; ++j ) {
      *left -= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,true,SF,CRAs...>::subAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(row(),element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::multAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) *= (*rhs)[j    ];
      matrix_(row(),j+1UL) *= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) *= (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::multAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t columns( size() );

   const size_t jpos( remainder ? prevMultiple( columns, SIMDSIZE ) : columns );
   BLAZE_INTERNAL_ASSERT( jpos <= columns, "Invalid end calculation" );

   size_t j( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; j<jpos; j+=SIMDSIZE ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && j<columns; ++j ) {
      *left *= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,true,SF,CRAs...>::multAssign( const SparseVector<VT,true>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t j( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; j<index; ++j )
         reset( matrix_(row(),j) );
      matrix_(row(),j) *= element->value();
      ++j;
   }

   for( ; j<size(); ++j ) {
      reset( matrix_(row(),j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::divAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) /= (*rhs)[j    ];
      matrix_(row(),j+1UL) /= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) /= (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,true,true,SF,CRAs...>::divAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t columns( size() );

   const size_t jpos( prevMultiple( columns, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns, "Invalid end calculation" );

   size_t j( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; j<jpos; j+=SIMDSIZE ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; j<columns; ++j ) {
      *left /= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Row for general column-major dense matrices.
// \ingroup row
//
// This specialization of Row adapts the class template to the requirements of general
// column-major dense matrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
class Row<MT,false,true,false,CRAs...>
   : public View< DenseVector< Row<MT,false,true,false,CRAs...>, true > >
   , private RowData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowData<CRAs...>;                     //!< The type of the RowData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Row instance.
   using This = Row<MT,false,true,false,CRAs...>;

   //! Base type of this Row instance.
   using BaseType = View< DenseVector<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Row instance.
   using ResultType    = RowTrait_t<MT,CRAs...>;       //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ElementType_t<MT>;            //!< Return type for expression template evaluations
   using CompositeType = const Row&;                   //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant row value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant row value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**RowIterator class definition****************************************************************
   /*!\brief Iterator over the elements of the dense row.
   */
   template< typename MatrixType      // Type of the dense matrix
           , typename IteratorType >  // Type of the dense matrix iterator
   class RowIterator
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
      /*!\brief Default constructor of the RowIterator class.
      */
      inline RowIterator() noexcept
         : matrix_( nullptr )  // The dense matrix containing the row
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   (     )      // Iterator to the current dense element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the RowIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      */
      inline RowIterator( MatrixType& matrix, size_t row, size_t column ) noexcept
         : matrix_( &matrix )  // The dense matrix containing the row
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   (         )  // Iterator to the current dense element
      {
         if( column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different RowIterator instances.
      //
      // \param it The row iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline RowIterator( const RowIterator<MatrixType2,IteratorType2>& it ) noexcept
         : matrix_( it.matrix_ )  // The dense matrix containing the row
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
      inline RowIterator& operator+=( size_t inc ) noexcept {
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
      inline RowIterator& operator-=( size_t dec ) noexcept {
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
      inline RowIterator& operator++() noexcept {
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
      inline const RowIterator operator++( int ) noexcept {
         const RowIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline RowIterator& operator--() noexcept {
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
      inline const RowIterator operator--( int ) noexcept {
         const RowIterator tmp( *this );
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
      /*!\brief Equality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ == rhs.column_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ < rhs.column_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ > rhs.column_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<=( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ <= rhs.column_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>=( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ >= rhs.column_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two row iterators.
      //
      // \param rhs The right-hand side row iterator.
      // \return The number of elements between the two row iterators.
      */
      inline DifferenceType operator-( const RowIterator& rhs ) const noexcept {
         return column_ - rhs.column_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a RowIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const RowIterator operator+( const RowIterator& it, size_t inc ) noexcept {
         return RowIterator( *it.matrix_, it.row_, it.column_+inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a RowIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const RowIterator operator+( size_t inc, const RowIterator& it ) noexcept {
         return RowIterator( *it.matrix_, it.row_, it.column_+inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a RowIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const RowIterator operator-( const RowIterator& it, size_t dec ) noexcept {
         return RowIterator( *it.matrix_, it.row_, it.column_-dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The dense matrix containing the row.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current dense element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class RowIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = RowIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, RowIterator< MT, Iterator_t<MT> > >;
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
   explicit inline Row( MT& matrix, RRAs... args );

   Row( const Row& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Row() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Row& operator=( const ElementType& rhs );
   inline Row& operator=( initializer_list<ElementType> list );
   inline Row& operator=( const Row& rhs );

   template< typename VT > inline Row& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator*=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator/=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator%=( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t  size() const noexcept;
   inline size_t  spacing() const noexcept;
   inline size_t  capacity() const noexcept;
   inline size_t  nonZeros() const;
   inline void    reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Row& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CRAs2 >
   inline bool canAlias( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CRAs2 >
   inline bool isAliased( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   template< typename VT > inline void assign    ( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void assign    ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void addAssign ( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void addAssign ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void subAssign ( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void subAssign ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void multAssign( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void multAssign( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void divAssign ( const DenseVector <VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CRAs2 > friend class Row;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE       ( MT );
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
/*!\brief Constructor for rows on column-major dense matrices.
//
// \param matrix The matrix containing the row.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CRAs >    // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Row<MT,false,true,false,CRAs...>::Row( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the row
{
   if( isChecked( args... ) ) {
      if( matrix_.rows() <= row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row() < matrix_.rows(), "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::Reference
   Row<MT,false,true,false,CRAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstReference
   Row<MT,false,true,false,CRAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::Reference
   Row<MT,false,true,false,CRAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstReference
   Row<MT,false,true,false,CRAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::Pointer
   Row<MT,false,true,false,CRAs...>::data() noexcept
{
   return matrix_.data() + row();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstPointer
   Row<MT,false,true,false,CRAs...>::data() const noexcept
{
   return matrix_.data() + row();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::Iterator
   Row<MT,false,true,false,CRAs...>::begin()
{
   return Iterator( matrix_, row(), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstIterator
   Row<MT,false,true,false,CRAs...>::begin() const
{
   return ConstIterator( matrix_, row(), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstIterator
   Row<MT,false,true,false,CRAs...>::cbegin() const
{
   return ConstIterator( matrix_, row(), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::Iterator
   Row<MT,false,true,false,CRAs...>::end()
{
   return Iterator( matrix_, row(), size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstIterator
   Row<MT,false,true,false,CRAs...>::end() const
{
   return ConstIterator( matrix_, row(), size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,false,CRAs...>::ConstIterator
   Row<MT,false,true,false,CRAs...>::cend() const
{
   return ConstIterator( matrix_, row(), size() );
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
/*!\brief Homogenous assignment to all row elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row.
//
// This function homogeneously assigns the given value to all elements of the row. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, row(), j, rhs ) )
         left(row(),j) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all row elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// row by means of an initializer list. The row elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the row, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate
// an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,true> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );

   std::fill( std::copy( list.begin(), list.end(), left.begin() ), left.end(), ElementType() );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Row.
//
// \param rhs Dense row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator=( const Row& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsExpression_v<MT> && rhs.canAlias( this ) ) {
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
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector_v<VT> )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator*=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator/=( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpDivAssign( left, tmp );
   }
   else {
      smpDivAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the assigned row.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::operator%=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, right );

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
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline MT& Row<MT,false,true,false,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline const MT& Row<MT,false,true,false,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the row.
//
// \return The size of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,false,CRAs...>::size() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the row.
//
// \return The minimum capacity of the row.
//
// This function returns the minimum capacity of the row, which corresponds to the current size
// plus padding.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,false,CRAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row.
//
// \return The maximum capacity of the dense row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,false,CRAs...>::capacity() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,false,CRAs...>::nonZeros() const
{
   const size_t columns( size() );
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns; ++j )
      if( !isDefault( matrix_(row(),j) ) )
         ++nonzeros;

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
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,true,false,CRAs...>::reset()
{
   using blaze::clear;

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j )
      clear( matrix_(row(),j) );
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
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
//
// This function scales the row by applying the given scalar value \a scalar to each element
// of the row. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a row
// on a lower or upper unitriangular matrix. The attempt to scale such a row results in a
// compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the scalar value
inline Row<MT,false,true,false,CRAs...>&
   Row<MT,false,true,false,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsStrictlyUpper_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsStrictlyLower_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_(row(),j) *= scalar;
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
/*!\brief Returns whether the dense row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,true,false,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can alias with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CRAs >   // Compile time row arguments
template< typename MT2       // Data type of the foreign dense row
        , bool SO2           // Storage order of the foreign dense row
        , bool SF2           // Symmetry flag of the foreign dense row
        , size_t... CRAs2 >  // Compile time row arguments of the foreign dense row
inline bool
   Row<MT,false,true,false,CRAs...>::canAlias( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row() == alias->row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,true,false,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CRAs >   // Compile time row arguments
template< typename MT2       // Data type of the foreign dense row
        , bool SO2           // Storage order of the foreign dense row
        , bool SF2           // Symmetry flag of the foreign dense row
        , size_t... CRAs2 >  // Compile time row arguments of the foreign dense row
inline bool
   Row<MT,false,true,false,CRAs...>::isAliased( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row() == alias->row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is properly aligned in memory.
//
// \return \a true in case the dense row is aligned, \a false if not.
//
// This function returns whether the dense row is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense row are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline bool Row<MT,false,true,false,CRAs...>::isAligned() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can be used in SMP assignments.
//
// \return \a true in case the dense row can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense row).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline bool Row<MT,false,true,false,CRAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,true,false,CRAs...>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) = (*rhs)[j    ];
      matrix_(row(),j+1UL) = (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) = (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,false,CRAs...>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(row(),element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,true,false,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) += (*rhs)[j    ];
      matrix_(row(),j+1UL) += (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) += (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,false,CRAs...>::addAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(row(),element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,true,false,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) -= (*rhs)[j    ];
      matrix_(row(),j+1UL) -= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) -= (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,false,CRAs...>::subAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(row(),element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,true,false,CRAs...>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) *= (*rhs)[j    ];
      matrix_(row(),j+1UL) *= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) *= (*rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,false,CRAs...>::multAssign( const SparseVector<VT,true>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t j( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; j<index; ++j )
         reset( matrix_(row(),j) );
      matrix_(row(),j) *= element->value();
      ++j;
   }

   for( ; j<size(); ++j ) {
      reset( matrix_(row(),j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,true,false,CRAs...>::divAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row(),j    ) /= (*rhs)[j    ];
      matrix_(row(),j+1UL) /= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(row(),jpos) /= (*rhs)[jpos];
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
/*!\brief Specialization of Row for symmetric column-major dense matrices.
// \ingroup row
//
// This specialization of Row adapts the class template to the requirements of symmetric
// column-major dense matrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
class Row<MT,false,true,true,CRAs...>
   : public View< DenseVector< Row<MT,false,true,true,CRAs...>, true > >
   , private RowData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowData<CRAs...>;                     //!< The type of the RowData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Row instance.
   using This = Row<MT,false,true,true,CRAs...>;

   //! Base type of this Row instance.
   using BaseType = View< DenseVector<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Row instance.
   using ResultType    = RowTrait_t<MT,CRAs...>;       //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the row elements.
   using ReturnType    = ElementType_t<MT>;            //!< Return type for expression template evaluations
   using CompositeType = const Row&;                   //!< Data type for composite expression templates.

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
   explicit inline Row( MT& matrix, RRAs... args );

   Row( const Row& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Row() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Row& operator=( const ElementType& rhs );
   inline Row& operator=( initializer_list<ElementType> list );
   inline Row& operator=( const Row& rhs );

   template< typename VT > inline Row& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator*=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator/=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator%=( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t  size() const noexcept;
   inline size_t  spacing() const noexcept;
   inline size_t  capacity() const noexcept;
   inline size_t  nonZeros() const;
   inline void    reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Row& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDAdd_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDSub_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDMult_v< ElementType, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDDiv_v< ElementType, ElementType_t<VT> > );
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

   template< typename MT2, bool SO2, bool SF2, size_t... CRAs2 >
   inline bool canAlias( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CRAs2 >
   inline bool isAliased( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT >
   inline auto assign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,true>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,true>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CRAs2 > friend class Row;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
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
/*!\brief Constructor for rows on column-major symmetric dense matrices.
//
// \param matrix The matrix containing the row.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CRAs >    // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Row<MT,false,true,true,CRAs...>::Row( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the row
{
   if( isChecked( args... ) ) {
      if( matrix_.rows() <= row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row() < matrix_.rows(), "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::Reference
   Row<MT,false,true,true,CRAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(index,row());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstReference
   Row<MT,false,true,true,CRAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(index,row());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::Reference
   Row<MT,false,true,true,CRAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstReference
   Row<MT,false,true,true,CRAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::Pointer
   Row<MT,false,true,true,CRAs...>::data() noexcept
{
   return matrix_.data( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstPointer
   Row<MT,false,true,true,CRAs...>::data() const noexcept
{
   return matrix_.data( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::Iterator
   Row<MT,false,true,true,CRAs...>::begin()
{
   return matrix_.begin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstIterator
   Row<MT,false,true,true,CRAs...>::begin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstIterator
   Row<MT,false,true,true,CRAs...>::cbegin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::Iterator
   Row<MT,false,true,true,CRAs...>::end()
{
   return matrix_.end( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstIterator
   Row<MT,false,true,true,CRAs...>::end() const
{
   return matrix_.cend( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,true,true,CRAs...>::ConstIterator
   Row<MT,false,true,true,CRAs...>::cend() const
{
   return matrix_.cend( row() );
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
/*!\brief Homogenous assignment to all row elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, row(), rhs ) )
         left(i,row()) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all row elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// row by means of an initializer list. The row elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the row, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate
// an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,true> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );

   std::fill( std::copy( list.begin(), list.end(), left.begin() ), left.end(), ElementType() );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Row.
//
// \param rhs Dense row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator=( const Row& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsExpression_v<MT> && rhs.canAlias( this ) ) {
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
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector_v<VT> )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator*=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator/=( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      smpDivAssign( left, tmp );
   }
   else {
      smpDivAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the assigned row.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::operator%=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, right );

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
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline MT& Row<MT,false,true,true,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline const MT& Row<MT,false,true,true,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the row.
//
// \return The size of the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,true,CRAs...>::size() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the row.
//
// \return The minimum capacity of the row.
//
// This function returns the minimum capacity of the row, which corresponds to the current size
// plus padding.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,true,CRAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row.
//
// \return The maximum capacity of the dense row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,true,CRAs...>::capacity() const noexcept
{
   return matrix_.capacity( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,true,true,CRAs...>::nonZeros() const
{
   return matrix_.nonZeros( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,true,true,CRAs...>::reset()
{
   matrix_.reset( row() );
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
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
//
// This function scales the row by applying the given scalar value \a scalar to each element
// of the row. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a row
// on a lower or upper unitriangular matrix. The attempt to scale such a row results in a
// compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the scalar value
inline Row<MT,false,true,true,CRAs...>&
   Row<MT,false,true,true,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsStrictlyLower_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsStrictlyUpper_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_(i,row()) *= scalar;
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
/*!\brief Returns whether the dense row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,true,true,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can alias with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CRAs >   // Compile time row arguments
template< typename MT2       // Data type of the foreign dense row
        , bool SO2           // Storage order of the foreign dense row
        , bool SF2           // Symmetry flag of the foreign dense row
        , size_t... CRAs2 >  // Compile time row arguments of the foreign dense row
inline bool
   Row<MT,false,true,true,CRAs...>::canAlias( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row() == alias->row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,true,true,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CRAs >   // Compile time row arguments
template< typename MT2       // Data type of the foreign dense row
        , bool SO2           // Storage order of the foreign dense row
        , bool SF2           // Symmetry flag of the foreign dense row
        , size_t... CRAs2 >  // Compile time row arguments of the foreign dense row
inline bool
   Row<MT,false,true,true,CRAs...>::isAliased( const Row<MT2,SO2,true,SF2,CRAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row() == alias->row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is properly aligned in memory.
//
// \return \a true in case the dense row is aligned, \a false if not.
//
// This function returns whether the dense row is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense row are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline bool Row<MT,false,true,true,CRAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can be used in SMP assignments.
//
// \return \a true in case the dense row can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense row).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
inline bool Row<MT,false,true,true,CRAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Row<MT,false,true,true,CRAs...>::SIMDType
   Row<MT,false,true,true,CRAs...>::load( size_t index ) const noexcept
{
   return matrix_.load( index, row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Row<MT,false,true,true,CRAs...>::SIMDType
   Row<MT,false,true,true,CRAs...>::loada( size_t index ) const noexcept
{
   return matrix_.loada( index, row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE typename Row<MT,false,true,true,CRAs...>::SIMDType
   Row<MT,false,true,true,CRAs...>::loadu( size_t index ) const noexcept
{
   return matrix_.loadu( index, row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store a specific SIMD element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,false,true,true,CRAs...>::store( size_t index, const SIMDType& value ) noexcept
{
   matrix_.store( index, row(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense row. The
// index must be smaller than the number of matrix columns. This function must \b NOT be
// called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,false,true,true,CRAs...>::storea( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storea( index, row(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unligned store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store a specific SIMD element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,false,true,true,CRAs...>::storeu( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storeu( index, row(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific SIMD element of the dense
// row. The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE void
   Row<MT,false,true,true,CRAs...>::stream( size_t index, const SIMDType& value ) noexcept
{
   matrix_.stream( index, row(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::assign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row()) = (*rhs)[i    ];
      matrix_(i+1UL,row()) = (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,row()) = (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::assign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t rows( size() );

   const size_t ipos( remainder ? prevMultiple( rows, SIMDSIZE ) : rows );
   BLAZE_INTERNAL_ASSERT( ipos <= rows, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   if( useStreaming && rows > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(*rhs).isAliased( this ) )
   {
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<rows; ++i ) {
         *left = *right; ++left; ++right;
      }
   }
   else
   {
      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<rows; ++i ) {
         *left = *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,true,CRAs...>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),row()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row()) += (*rhs)[i    ];
      matrix_(i+1UL,row()) += (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,row()) += (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t rows( size() );

   const size_t ipos( remainder ? prevMultiple( rows, SIMDSIZE ) : rows );
   BLAZE_INTERNAL_ASSERT( ipos <= rows, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && i<rows; ++i ) {
      *left += *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,true,CRAs...>::addAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),row()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row()) -= (*rhs)[i    ];
      matrix_(i+1UL,row()) -= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,row()) -= (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t rows( size() );

   const size_t ipos( remainder ? prevMultiple( rows, SIMDSIZE ) : rows );
   BLAZE_INTERNAL_ASSERT( ipos <= rows, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && i<rows; ++i ) {
      *left -= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,true,CRAs...>::subAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),row()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::multAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row()) *= (*rhs)[i    ];
      matrix_(i+1UL,row()) *= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,row()) *= (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::multAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<MT> || !IsPadded_v<VT> );

   const size_t rows( size() );

   const size_t ipos( remainder ? prevMultiple( rows, SIMDSIZE ) : rows );
   BLAZE_INTERNAL_ASSERT( ipos <= rows, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && i<rows; ++i ) {
      *left *= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,true,true,CRAs...>::multAssign( const SparseVector<VT,true>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( matrix_(i,row()) );
      matrix_(i,row()) *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( matrix_(i,row()) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::divAssign( const DenseVector<VT,true>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row()) /= (*rhs)[i    ];
      matrix_(i+1UL,row()) /= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,row()) /= (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Row<MT,false,true,true,CRAs...>::divAssign( const DenseVector<VT,true>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t rows( size() );

   const size_t ipos( prevMultiple( rows, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<rows; ++i ) {
      *left /= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
