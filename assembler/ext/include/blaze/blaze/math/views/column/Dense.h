//=================================================================================================
/*!
//  \file blaze/math/views/column/Dense.h
//  \brief Column specialization for dense matrices
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_DENSE_H_
#define _BLAZE_MATH_VIEWS_COLUMN_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
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
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
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
#include <blaze/math/views/column/BaseTemplate.h>
#include <blaze/math/views/column/ColumnData.h>
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
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Column for columns on column-major dense matrices.
// \ingroup column
//
// This specialization of Column adapts the class template to the requirements of column-major
// dense matrices.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
class Column<MT,true,true,SF,CCAs...>
   : public View< DenseVector< Column<MT,true,true,SF,CCAs...>, false > >
   , private ColumnData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnData<CCAs...>;                  //!< The type of the ColumnData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Column instance.
   using This = Column<MT,true,true,SF,CCAs...>;

   //! Base type of this Column instance.
   using BaseType = View< DenseVector<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Column instance.
   using ResultType    = ColumnTrait_t<MT,CCAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Column&;                //!< Data type for composite expression templates.

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
   explicit inline Column( MT& matrix, RCAs... args );

   Column( const Column& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Column() = default;
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
   inline Column& operator=( const ElementType& rhs );
   inline Column& operator=( initializer_list<ElementType> list );
   inline Column& operator=( const Column& rhs );

   template< typename VT > inline Column& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator*=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator/=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator%=( const Vector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::column;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Column& scale( const Other& scalar );
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

   template< typename MT2, bool SO2, bool SF2, size_t... CCAs2 >
   inline bool canAlias( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CCAs2 >
   inline bool isAliased( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

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
   inline auto assign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CCAs2 > friend class Column;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
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
/*!\brief Constructor for columns on column-major dense matrices.
//
// \param matrix The matrix containing the column.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , bool SF             // Symmetry flag
        , size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Column<MT,true,true,SF,CCAs...>::Column( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the column
{
   if( isChecked( args... ) ) {
      if( matrix_.columns() <= column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( column() < matrix_.columns(), "Invalid column access index" );
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
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::Reference
   Column<MT,true,true,SF,CCAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstReference
   Column<MT,true,true,SF,CCAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::Reference
   Column<MT,true,true,SF,CCAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstReference
   Column<MT,true,true,SF,CCAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column. Note that in case
// of a row-major matrix you can NOT assume that the column elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::Pointer
   Column<MT,true,true,SF,CCAs...>::data() noexcept
{
   return matrix_.data( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column. Note that in case
// of a row-major matrix you can NOT assume that the column elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstPointer
   Column<MT,true,true,SF,CCAs...>::data() const noexcept
{
   return matrix_.data( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::Iterator
   Column<MT,true,true,SF,CCAs...>::begin()
{
   return matrix_.begin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstIterator
   Column<MT,true,true,SF,CCAs...>::begin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstIterator
   Column<MT,true,true,SF,CCAs...>::cbegin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::Iterator
   Column<MT,true,true,SF,CCAs...>::end()
{
   return matrix_.end( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstIterator
   Column<MT,true,true,SF,CCAs...>::end() const
{
   return matrix_.cend( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,true,SF,CCAs...>::ConstIterator
   Column<MT,true,true,SF,CCAs...>::cend() const
{
   return matrix_.cend( column() );
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
/*!\brief Homogenous assignment to all column elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column.
//
// This function homogeneously assigns the given value to all elements of the column. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, column(), rhs ) )
         left(i,column()) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// column by means of an initializer list. The column elements are assigned the values from the
// given initializer list. Missing values are reset to their default state. Note that in case
// the size of the initializer list exceeds the size of the column, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,false> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
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
/*!\brief Copy assignment operator for Column.
//
// \param rhs Dense column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator=( const Column& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be added to the dense column.
// \return Reference to the assigned column.
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator+=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the dense column.
// \return Reference to the assigned column.
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator-=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator*=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator/=( const DenseVector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::operator%=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
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
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline MT& Column<MT,true,true,SF,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline const MT& Column<MT,true,true,SF,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,true,SF,CCAs...>::size() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the column.
//
// \return The minimum capacity of the column.
//
// This function returns the minimum capacity of the column, which corresponds to the current
// size plus padding.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,true,SF,CCAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column.
//
// \return The maximum capacity of the dense column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,true,SF,CCAs...>::capacity() const noexcept
{
   return matrix_.capacity( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,true,SF,CCAs...>::nonZeros() const
{
   return matrix_.nonZeros( column() );
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
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,true,true,SF,CCAs...>::reset()
{
   matrix_.reset( column() );
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
/*!\brief Scaling of the column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the dense column.
//
// This function scales the column by applying the given scalar value \a scalar to each element
// of the column. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a column
// on a lower or upper unitriangular matrix. The attempt to scale such a column results in a
// compile time error!
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the scalar value
inline Column<MT,true,true,SF,CCAs...>&
   Column<MT,true,true,SF,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsStrictlyLower_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsStrictlyUpper_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_(i,column()) *= scalar;
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
/*!\brief Returns whether the dense column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In
// contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,true,true,SF,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column can alias with the given dense column \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In
// contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , bool SF            // Symmetry flag
        , size_t... CCAs >   // Compile time column arguments
template< typename MT2       // Data type of the foreign dense column
        , bool SO2           // Storage order of the foreign dense column
        , bool SF2           // Symmetry flag of the foreign dense column
        , size_t... CCAs2 >  // Compile time column arguments of the foreign dense column
inline bool
   Column<MT,true,true,SF,CCAs...>::canAlias( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( column() == alias->column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In
// contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,true,true,SF,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given dense column \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In
// contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , bool SF            // Symmetry flag
        , size_t... CCAs >   // Compile time column arguments
template< typename MT2       // Data type of the foreign dense column
        , bool SO2           // Storage order of the foreign dense column
        , bool SF2           // Symmetry flag of the foreign dense column
        , size_t... CCAs2 >  // Compile time column arguments of the foreign dense column
inline bool
   Column<MT,true,true,SF,CCAs...>::isAliased( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( column() == alias->column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is properly aligned in memory.
//
// \return \a true in case the dense column is aligned, \a false if not.
//
// This function returns whether the dense column is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense column are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool Column<MT,true,true,SF,CCAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column can be used in SMP assignments.
//
// \return \a true in case the dense column can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense column can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size of
// the dense column).
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool Column<MT,true,true,SF,CCAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense column. The index must
// be smaller than the number of matrix rows. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Column<MT,true,true,SF,CCAs...>::SIMDType
   Column<MT,true,true,SF,CCAs...>::load( size_t index ) const noexcept
{
   return matrix_.load( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Column<MT,true,true,SF,CCAs...>::SIMDType
   Column<MT,true,true,SF,CCAs...>::loada( size_t index ) const noexcept
{
   return matrix_.loada( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Column<MT,true,true,SF,CCAs...>::SIMDType
   Column<MT,true,true,SF,CCAs...>::loadu( size_t index ) const noexcept
{
   return matrix_.loadu( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store a specific SIMD element of the dense column. The index must
// be smaller than the number of matrix rows. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,true,true,SF,CCAs...>::store( size_t index, const SIMDType& value ) noexcept
{
   matrix_.store( index, column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,true,true,SF,CCAs...>::storea( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storea( index, column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,true,true,SF,CCAs...>::storeu( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storeu( index, column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific SIMD element of the dense
// column. The index must be smaller than the number of matrix rows. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,true,true,SF,CCAs...>::stream( size_t index, const SIMDType& value ) noexcept
{
   matrix_.stream( index, column(), value );
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::assign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) = (*rhs)[i    ];
      matrix_(i+1UL,column()) = (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) = (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::assign( const DenseVector<VT,false>& rhs )
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
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,true,SF,CCAs...>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),column()) = element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) += (*rhs)[i    ];
      matrix_(i+1UL,column()) += (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) += (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
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
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,true,SF,CCAs...>::addAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),column()) += element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) -= (*rhs)[i    ];
      matrix_(i+1UL,column()) -= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) -= (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
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
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,true,SF,CCAs...>::subAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),column()) -= element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::multAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) *= (*rhs)[i    ];
      matrix_(i+1UL,column()) *= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) *= (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::multAssign( const DenseVector<VT,false>& rhs )
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
      *left *= *right; ++left, ++right;
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,true,SF,CCAs...>::multAssign( const SparseVector<VT,false>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( matrix_(i,column()) );
      matrix_(i,column()) *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( matrix_(i,column()) );
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::divAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) /= (*rhs)[i    ];
      matrix_(i+1UL,column()) /= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) /= (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,true,true,SF,CCAs...>::divAssign( const DenseVector<VT,false>& rhs )
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








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL ROW-MAJOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Column for general row-major dense matrices.
// \ingroup column
//
// This specialization of Column adapts the class template to the requirements of general
// row-major dense matrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
class Column<MT,false,true,false,CCAs...>
   : public View< DenseVector< Column<MT,false,true,false,CCAs...>, false > >
   , private ColumnData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnData<CCAs...>;                  //!< The type of the ColumnData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Column instance.
   using This = Column<MT,false,true,false,CCAs...>;

   //! Base type of this Column instance.
   using BaseType = View< DenseVector<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Column instance.
   using ResultType    = ColumnTrait_t<MT,CCAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Column&;                //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant column value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant column value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**ColumnIterator class definition*************************************************************
   /*!\brief Iterator over the elements of the dense column.
   */
   template< typename MatrixType      // Type of the dense matrix
           , typename IteratorType >  // Type of the dense matrix iterator
   class ColumnIterator
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
      /*!\brief Default constructor of the ColumnIterator class.
      */
      inline ColumnIterator() noexcept
         : matrix_( nullptr )  // The dense matrix containing the column
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   (     )      // Iterator to the current dense element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column ) noexcept
         : matrix_( &matrix )  // The dense matrix containing the column
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   (         )  // Iterator to the current dense element
      {
         if( row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ColumnIterator instances.
      //
      // \param it The column iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline ColumnIterator( const ColumnIterator<MatrixType2,IteratorType2>& it ) noexcept
         : matrix_( it.matrix_ )  // The dense matrix containing the column
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
      inline ColumnIterator& operator+=( size_t inc ) noexcept {
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
      inline ColumnIterator& operator-=( size_t dec ) noexcept {
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
      inline ColumnIterator& operator++() noexcept {
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
      inline const ColumnIterator operator++( int ) noexcept {
         const ColumnIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ColumnIterator& operator--() noexcept {
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
      inline const ColumnIterator operator--( int ) noexcept {
         const ColumnIterator tmp( *this );
         --(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Subscript operator***********************************************************************
      /*!\brief Direct access to the dense column elements.
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
      /*!\brief Direct access to the dense vector element at the current iterator position.
      //
      // \return The current value of the dense element.
      */
      inline ReferenceType operator*() const {
         return *pos_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense vector element at the current iterator position.
      //
      // \return Reference to the dense vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ == rhs.row_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ < rhs.row_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ > rhs.row_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<=( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ <= rhs.row_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>=( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ >= rhs.row_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two column iterators.
      //
      // \param rhs The right-hand side column iterator.
      // \return The number of elements between the two column iterators.
      */
      inline DifferenceType operator-( const ColumnIterator& rhs ) const noexcept {
         return row_ - rhs.row_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ColumnIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ColumnIterator operator+( const ColumnIterator& it, size_t inc ) noexcept {
         return ColumnIterator( *it.matrix_, it.row_+inc, it.column_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ColumnIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ColumnIterator operator+( size_t inc, const ColumnIterator& it ) noexcept {
         return ColumnIterator( *it.matrix_, it.row_+inc, it.column_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ColumnIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ColumnIterator operator-( const ColumnIterator& it, size_t dec ) noexcept {
         return ColumnIterator( *it.matrix_, it.row_-dec, it.column_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The dense matrix containing the column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current dense element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class ColumnIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = ColumnIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, ColumnIterator< MT, Iterator_t<MT> > >;
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
   explicit inline Column( MT& matrix, RCAs... args );

   Column( const Column& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Column() = default;
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
   inline Column& operator=( const ElementType& rhs );
   inline Column& operator=( initializer_list<ElementType> list );
   inline Column& operator=( const Column& rhs );

   template< typename VT > inline Column& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator*=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator/=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator%=( const Vector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::column;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Column& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias ( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CCAs2 >
   inline bool canAlias ( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CCAs2 >
   inline bool isAliased( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   template< typename VT > inline void assign    ( const DenseVector <VT,false>& rhs );
   template< typename VT > inline void assign    ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline void addAssign ( const DenseVector <VT,false>& rhs );
   template< typename VT > inline void addAssign ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline void subAssign ( const DenseVector <VT,false>& rhs );
   template< typename VT > inline void subAssign ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline void multAssign( const DenseVector <VT,false>& rhs );
   template< typename VT > inline void multAssign( const SparseVector<VT,false>& rhs );
   template< typename VT > inline void divAssign ( const DenseVector <VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CCAs2 > friend class Column;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( MT );
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
/*!\brief Constructor for columns on row-major dense matrices.
//
// \param matrix The matrix containing the column.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Column<MT,false,true,false,CCAs...>::Column( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the column
{
   if( isChecked( args... ) ) {
      if( matrix_.columns() <= column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( column() < matrix_.columns(), "Invalid column access index" );
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
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::Reference
   Column<MT,false,true,false,CCAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstReference
   Column<MT,false,true,false,CCAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::Reference
   Column<MT,false,true,false,CCAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstReference
   Column<MT,false,true,false,CCAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column. Note that in case
// of a row-major matrix you can NOT assume that the column elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::Pointer
   Column<MT,false,true,false,CCAs...>::data() noexcept
{
   return matrix_.data() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column. Note that in case
// of a row-major matrix you can NOT assume that the column elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstPointer
   Column<MT,false,true,false,CCAs...>::data() const noexcept
{
   return matrix_.data() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::Iterator
   Column<MT,false,true,false,CCAs...>::begin()
{
   return Iterator( matrix_, 0UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstIterator
   Column<MT,false,true,false,CCAs...>::begin() const
{
   return ConstIterator( matrix_, 0UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstIterator
   Column<MT,false,true,false,CCAs...>::cbegin() const
{
   return ConstIterator( matrix_, 0UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::Iterator
   Column<MT,false,true,false,CCAs...>::end()
{
   return Iterator( matrix_, size(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstIterator
   Column<MT,false,true,false,CCAs...>::end() const
{
   return ConstIterator( matrix_, size(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,false,CCAs...>::ConstIterator
   Column<MT,false,true,false,CCAs...>::cend() const
{
   return ConstIterator( matrix_, size(), column() );
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
/*!\brief Homogenous assignment to all column elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column.
//
// This function homogeneously assigns the given value to all elements of the column. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, column(), rhs ) )
         left(i,column()) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// column by means of an initializer list. The column elements are assigned the values from the
// given initializer list. Missing values are reset to their default state. Note that in case
// the size of the initializer list exceeds the size of the column, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,false> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
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
/*!\brief Copy assignment operator for Column.
//
// \param rhs Dense column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator=( const Column& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be added to the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator+=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator-=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator*=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator/=( const DenseVector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::operator%=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
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
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline MT& Column<MT,false,true,false,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline const MT& Column<MT,false,true,false,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,false,CCAs...>::size() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the column.
//
// \return The minimum capacity of the column.
//
// This function returns the minimum capacity of the column, which corresponds to the current
// size plus padding.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,false,CCAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column.
//
// \return The maximum capacity of the dense column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,false,CCAs...>::capacity() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,false,CCAs...>::nonZeros() const
{
   const size_t rows( size() );
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows; ++i )
      if( !isDefault( matrix_(i,column()) ) )
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
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,true,false,CCAs...>::reset()
{
   using blaze::clear;

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i )
      clear( matrix_(i,column()) );
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
/*!\brief Scaling of the column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the dense column.
//
// This function scales the column by applying the given scalar value \a scalar to each element
// of the column. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a column
// on a lower or upper unitriangular matrix. The attempt to scale such a column results in a
// compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the scalar value
inline Column<MT,false,true,false,CCAs...>&
   Column<MT,false,true,false,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsStrictlyLower_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsStrictlyUpper_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_(i,column()) *= scalar;
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
/*!\brief Returns whether the dense column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,true,false,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column can alias with the given dense column \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CCAs >   // Compile time column arguments
template< typename MT2       // Data type of the foreign dense column
        , bool SO2           // Storage order of the foreign dense column
        , bool SF2           // Symmetry flag of the foreign dense column
        , size_t... CCAs2 >  // Compile time column arguments of the foreign dense column
inline bool
   Column<MT,false,true,false,CCAs...>::canAlias( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( column() == alias->column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,true,false,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given dense column \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CCAs >   // Compile time column arguments
template< typename MT2       // Data type of the foreign dense column
        , bool SO2           // Storage order of the foreign dense column
        , bool SF2           // Symmetry flag of the foreign dense column
        , size_t... CCAs2 >  // Compile time column arguments of the foreign dense column
inline bool
   Column<MT,false,true,false,CCAs...>::isAliased( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( column() == alias->column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is properly aligned in memory.
//
// \return \a true in case the dense column is aligned, \a false if not.
//
// This function returns whether the dense column is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense column are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline bool Column<MT,false,true,false,CCAs...>::isAligned() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column can be used in SMP assignments.
//
// \return \a true in case the dense column can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense column can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size of
// the dense column).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline bool Column<MT,false,true,false,CCAs...>::canSMPAssign() const noexcept
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,true,false,CCAs...>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) = (*rhs)[i    ];
      matrix_(i+1UL,column()) = (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) = (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,false,CCAs...>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),column()) = element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,true,false,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) += (*rhs)[i    ];
      matrix_(i+1UL,column()) += (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) += (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,false,CCAs...>::addAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),column()) += element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,true,false,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) -= (*rhs)[i    ];
      matrix_(i+1UL,column()) -= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) -= (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,false,CCAs...>::subAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(element->index(),column()) -= element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,true,false,CCAs...>::multAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) *= (*rhs)[i    ];
      matrix_(i+1UL,column()) *= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) *= (*rhs)[ipos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,false,CCAs...>::multAssign( const SparseVector<VT,false>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( matrix_(i,column()) );
      matrix_(i,column()) *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( matrix_(i,column()) );
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,true,false,CCAs...>::divAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,column()) /= (*rhs)[i    ];
      matrix_(i+1UL,column()) /= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      matrix_(ipos,column()) /= (*rhs)[ipos];
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
/*!\brief Specialization of Column for symmetric row-major dense matrices.
// \ingroup column
//
// This specialization of Column adapts the class template to the requirements of symmetric
// row-major dense matrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
class Column<MT,false,true,true,CCAs...>
   : public View< DenseVector< Column<MT,false,true,true,CCAs...>, false > >
   , private ColumnData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnData<CCAs...>;                  //!< The type of the ColumnData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Column instance.
   using This = Column<MT,false,true,true,CCAs...>;

   //! Base type of this Column instance.
   using BaseType = View< DenseVector<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Column instance.
   using ResultType    = ColumnTrait_t<MT,CCAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< SIMD type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Column&;                //!< Data type for composite expression templates.

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
   explicit inline Column( MT& matrix, RCAs... args );

   Column( const Column& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Column() = default;
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
   inline Column& operator=( const ElementType& rhs );
   inline Column& operator=( initializer_list<ElementType> list );
   inline Column& operator=( const Column& rhs );

   template< typename VT > inline Column& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator*=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator/=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator%=( const Vector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::column;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Column& scale( const Other& scalar );
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

   template< typename MT2, bool SO2, bool SF2, size_t... CCAs2 >
   inline bool canAlias( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, bool SO2, bool SF2, size_t... CCAs2 >
   inline bool isAliased( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept;

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
   inline auto assign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,false>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,false>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CCAs2 > friend class Column;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE( MT );
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
/*!\brief Constructor for columns on row-major symmetric dense matrices.
//
// \param matrix The matrix containing the column.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Column<MT,false,true,true,CCAs...>::Column( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the column
{
   if( isChecked( args... ) ) {
      if( matrix_.columns() <= column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( column() < matrix_.columns(), "Invalid column access index" );
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
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::Reference
   Column<MT,false,true,true,CCAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(column(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstReference
   Column<MT,false,true,true,CCAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(column(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::Reference
   Column<MT,false,true,true,CCAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstReference
   Column<MT,false,true,true,CCAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column. Note that in case
// of a row-major matrix you can NOT assume that the column elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::Pointer
   Column<MT,false,true,true,CCAs...>::data() noexcept
{
   return matrix_.data( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the column elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense column. Note that in case
// of a row-major matrix you can NOT assume that the column elements lie adjacent to each other!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstPointer
   Column<MT,false,true,true,CCAs...>::data() const noexcept
{
   return matrix_.data( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::Iterator
   Column<MT,false,true,true,CCAs...>::begin()
{
   return matrix_.begin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstIterator
   Column<MT,false,true,true,CCAs...>::begin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstIterator
   Column<MT,false,true,true,CCAs...>::cbegin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::Iterator
   Column<MT,false,true,true,CCAs...>::end()
{
   return matrix_.end( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstIterator
   Column<MT,false,true,true,CCAs...>::end() const
{
   return matrix_.cend( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,true,true,CCAs...>::ConstIterator
   Column<MT,false,true,true,CCAs...>::cend() const
{
   return matrix_.cend( column() );
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
/*!\brief Homogenous assignment to all column elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, column(), j, rhs ) )
         left(column(),j) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// column by means of an initializer list. The column elements are assigned the values from the
// given initializer list. Missing values are reset to their default state. Note that in case
// the size of the initializer list exceeds the size of the column, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,false> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
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
/*!\brief Copy assignment operator for Column.
//
// \param rhs Dense column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator=( const Column& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be added to the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator+=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator-=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, 0UL, column() ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator*=( const Vector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator/=( const DenseVector<VT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, 0UL, column() ) ) {
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
// \return Reference to the assigned column.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::operator%=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
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
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline MT& Column<MT,false,true,true,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline const MT& Column<MT,false,true,true,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,true,CCAs...>::size() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the column.
//
// \return The minimum capacity of the column.
//
// This function returns the minimum capacity of the column, which corresponds to the current
// size plus padding.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,true,CCAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column.
//
// \return The maximum capacity of the dense column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,true,CCAs...>::capacity() const noexcept
{
   return matrix_.capacity( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,true,true,CCAs...>::nonZeros() const
{
   return matrix_.nonZeros( column() );
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
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,true,true,CCAs...>::reset()
{
   matrix_.reset( column() );
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
/*!\brief Scaling of the column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the dense column.
//
// This function scales the column by applying the given scalar value \a scalar to each element
// of the column. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a column
// on a lower or upper unitriangular matrix. The attempt to scale such a column results in a
// compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the scalar value
inline Column<MT,false,true,true,CCAs...>&
   Column<MT,false,true,true,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsStrictlyUpper_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsStrictlyLower_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_(column(),j) *= scalar;
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
/*!\brief Returns whether the dense column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In
// contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,true,true,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column can alias with the given dense column \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In
// contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CCAs >   // Compile time column arguments
template< typename MT2       // Data type of the foreign dense column
        , bool SO2           // Storage order of the foreign dense column
        , bool SF2           // Symmetry flag of the foreign dense column
        , size_t... CCAs2 >  // Compile time column arguments of the foreign dense column
inline bool
   Column<MT,false,true,true,CCAs...>::canAlias( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( column() == alias->column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In
// contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,true,true,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given dense column \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In
// contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CCAs >   // Compile time column arguments
template< typename MT2       // Data type of the foreign dense column
        , bool SO2           // Storage order of the foreign dense column
        , bool SF2           // Symmetry flag of the foreign dense column
        , size_t... CCAs2 >  // Compile time column arguments of the foreign dense column
inline bool
   Column<MT,false,true,true,CCAs...>::isAliased( const Column<MT2,SO2,true,SF2,CCAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( column() == alias->column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is properly aligned in memory.
//
// \return \a true in case the dense column is aligned, \a false if not.
//
// This function returns whether the dense column is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense column are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline bool Column<MT,false,true,true,CCAs...>::isAligned() const noexcept
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column can be used in SMP assignments.
//
// \return \a true in case the dense column can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense column can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size of
// the dense column).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
inline bool Column<MT,false,true,true,CCAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense column. The index must
// be smaller than the number of matrix rows. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Column<MT,false,true,true,CCAs...>::SIMDType
   Column<MT,false,true,true,CCAs...>::load( size_t index ) const noexcept
{
   return matrix_.load( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Column<MT,false,true,true,CCAs...>::SIMDType
   Column<MT,false,true,true,CCAs...>::loada( size_t index ) const noexcept
{
   return matrix_.loada( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE typename Column<MT,false,true,true,CCAs...>::SIMDType
   Column<MT,false,true,true,CCAs...>::loadu( size_t index ) const noexcept
{
   return matrix_.loadu( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store a specific SIMD element of the dense column. The index must
// be smaller than the number of matrix rows. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,false,true,true,CCAs...>::store( size_t index, const SIMDType& value ) noexcept
{
   matrix_.store( column(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,false,true,true,CCAs...>::storea( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storea( column(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store a specific SIMD element of the dense column. The
// index must be smaller than the number of matrix rows. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,false,true,true,CCAs...>::storeu( size_t index, const SIMDType& value ) noexcept
{
   matrix_.storeu( column(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific SIMD element of the dense
// column. The index must be smaller than the number of matrix rows. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE void
   Column<MT,false,true,true,CCAs...>::stream( size_t index, const SIMDType& value ) noexcept
{
   matrix_.stream( column(), index, value );
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::assign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(column(),j    ) = (*rhs)[j    ];
      matrix_(column(),j+1UL) = (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(column(),jpos) = (*rhs)[jpos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::assign( const DenseVector<VT,false>& rhs )
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,true,CCAs...>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(column(),element->index()) = element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(column(),j    ) += (*rhs)[j    ];
      matrix_(column(),j+1UL) += (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(column(),jpos) += (*rhs)[jpos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,true,CCAs...>::addAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(column(),element->index()) += element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(column(),j    ) -= (*rhs)[j    ];
      matrix_(column(),j+1UL) -= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(column(),jpos) -= (*rhs)[jpos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,true,CCAs...>::subAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      matrix_(column(),element->index()) -= element->value();
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::multAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(column(),j    ) *= (*rhs)[j    ];
      matrix_(column(),j+1UL) *= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(column(),jpos) *= (*rhs)[jpos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::multAssign( const DenseVector<VT,false>& rhs )
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,true,true,CCAs...>::multAssign( const SparseVector<VT,false>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t j( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; j<index; ++j )
         reset( matrix_(column(),j) );
      matrix_(column(),j) *= element->value();
      ++j;
   }

   for( ; j<size(); ++j ) {
      reset( matrix_(column(),j) );
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::divAssign( const DenseVector<VT,false>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t jpos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(column(),j    ) /= (*rhs)[j    ];
      matrix_(column(),j+1UL) /= (*rhs)[j+1UL];
   }
   if( jpos < (*rhs).size() )
      matrix_(column(),jpos) /= (*rhs)[jpos];
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
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline auto Column<MT,false,true,true,CCAs...>::divAssign( const DenseVector<VT,false>& rhs )
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

} // namespace blaze

#endif
