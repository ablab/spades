//=================================================================================================
/*!
//  \file blaze/math/views/band/Dense.h
//  \brief Band specialization for dense matrices
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

#ifndef _BLAZE_MATH_VIEWS_BAND_DENSE_H_
#define _BLAZE_MATH_VIEWS_BAND_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/band/BandData.h>
#include <blaze/math/views/band/BaseTemplate.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DENSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Band for dense matrices.
// \ingroup band
//
// This specialization of Band adapts the class template to the requirements of dense matrices.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
class Band<MT,TF,true,false,CBAs...>
   : public View< DenseVector< Band<MT,TF,true,false,CBAs...>, TF > >
   , private BandData<CBAs...>
{
 private:
   //**Type definitions****************************************************************************
   using RT       = ResultType_t<MT>;                     //!< Result type of the dense matrix expression.
   using DataType = BandData<CBAs...>;                    //!< The type of the BandData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Band instance.
   using This = Band<MT,TF,true,false,CBAs...>;

   //! Base type of this Band instance.
   using BaseType = View< DenseVector<This,TF> >;

   using ViewedType    = MT;                           //!< The type viewed by this Band instance.
   using ResultType    = BandTrait_t<RT,CBAs...>;      //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the band elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations

   //! Data type for composite expression templates.
   using CompositeType = If_t< RequiresEvaluation_v<MT>, const ResultType, const Band& >;

   //! Reference to a constant band value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant band value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant band value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant band value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**BandIterator class definition***************************************************************
   /*!\brief Iterator over the elements of the dense band.
   */
   template< typename MatrixType      // Type of the dense matrix
           , typename IteratorType >  // Type of the dense matrix iterator
   class BandIterator
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
      /*!\brief Default constructor of the BandIterator class.
      */
      inline BandIterator() noexcept
         : matrix_( nullptr )  // The dense matrix containing the band
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   (     )      // Iterator to the current dense element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the BandIterator class.
      //
      // \param matrix The matrix containing the band.
      // \param rowIndex The initial row index.
      // \param columnIndex The initial column index.
      */
      inline BandIterator( MatrixType& matrix, size_t rowIndex, size_t columnIndex ) noexcept
         : matrix_( &matrix     )  // The dense matrix containing the band
         , row_   ( rowIndex    )  // The current row index
         , column_( columnIndex )  // The current column index
         , pos_   (             )  // Iterator to the current dense element
      {
         if( IsRowMajorMatrix_v<MatrixType> && row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else if( IsColumnMajorMatrix_v<MatrixType> && column_ != matrix_->columns() )
            pos_ = matrix_->begin( column_ ) + row_;
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different BandIterator instances.
      //
      // \param it The band iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline BandIterator( const BandIterator<MatrixType2,IteratorType2>& it ) noexcept
         : matrix_( it.matrix_ )  // The dense matrix containing the band
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
      inline BandIterator& operator+=( size_t inc ) noexcept {
         using blaze::reset;

         row_    += inc;
         column_ += inc;

         if( IsRowMajorMatrix_v<MatrixType> && row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else if( IsColumnMajorMatrix_v<MatrixType> && column_ != matrix_->columns() )
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
      inline BandIterator& operator-=( size_t dec ) noexcept {
         using blaze::reset;

         row_    -= dec;
         column_ -= dec;

         if( IsRowMajorMatrix_v<MatrixType> && row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else if( IsColumnMajorMatrix_v<MatrixType> && column_ != matrix_->columns() )
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
      inline BandIterator& operator++() noexcept {
         using blaze::reset;

         ++row_;
         ++column_;

         if( IsRowMajorMatrix_v<MatrixType> && row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else if( IsColumnMajorMatrix_v<MatrixType> && column_ != matrix_->columns() )
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
      inline const BandIterator operator++( int ) noexcept {
         const BandIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BandIterator& operator--() noexcept {
         using blaze::reset;

         --row_;
         --column_;

         if( IsRowMajorMatrix_v<MatrixType> && row_ != matrix_->rows() )
            pos_ = matrix_->begin( row_ ) + column_;
         else if( IsColumnMajorMatrix_v<MatrixType> && column_ != matrix_->columns() )
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
      inline const BandIterator operator--( int ) noexcept {
         const BandIterator tmp( *this );
         --(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Subscript operator***********************************************************************
      /*!\brief Direct access to the dense band elements.
      //
      // \param index Access index.
      // \return Reference to the accessed value.
      */
      inline ReferenceType operator[]( size_t index ) const {
         BLAZE_USER_ASSERT( row_   +index < matrix_->rows()   , "Invalid access index detected" );
         BLAZE_USER_ASSERT( column_+index < matrix_->columns(), "Invalid access index detected" );
         const IteratorType pos( IsRowMajorMatrix_v<MatrixType>
                               ? matrix_->begin( row_+index ) + column_ + index
                               : matrix_->begin( column_+index ) + row_ + index );
         return *pos;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense band element at the current iterator position.
      //
      // \return Reference to the current value.
      */
      inline ReferenceType operator*() const {
         return *pos_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense band element at the current iterator position.
      //
      // \return Pointer to the dense band element at the current iterator position.
      */
      inline PointerType operator->() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two BandIterator objects.
      //
      // \param rhs The right-hand side band iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const BandIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ == rhs.row_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two BandIterator objects.
      //
      // \param rhs The right-hand side band iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const BandIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two BandIterator objects.
      //
      // \param rhs The right-hand side band iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<( const BandIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ < rhs.row_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two BandIterator objects.
      //
      // \param rhs The right-hand side band iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>( const BandIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ > rhs.row_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two BandIterator objects.
      //
      // \param rhs The right-hand side band iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator<=( const BandIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ <= rhs.row_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two BandIterator objects.
      //
      // \param rhs The right-hand side band iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator>=( const BandIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ >= rhs.row_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two band iterators.
      //
      // \param rhs The right-hand side band iterator.
      // \return The number of elements between the two band iterators.
      */
      inline DifferenceType operator-( const BandIterator& rhs ) const noexcept {
         return row_ - rhs.row_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a BandIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const BandIterator operator+( const BandIterator& it, size_t inc ) noexcept {
         return BandIterator( *it.matrix_, it.row_+inc, it.column_+inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a BandIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const BandIterator operator+( size_t inc, const BandIterator& it ) noexcept {
         return BandIterator( *it.matrix_, it.row_+inc, it.column_+inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a BandIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const BandIterator operator-( const BandIterator& it, size_t dec ) noexcept {
         return BandIterator( *it.matrix_, it.row_-dec, it.column_-dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The dense matrix containing the band.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current dense element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class BandIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = BandIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, BandIterator< MT, Iterator_t<MT> > >;
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
   template< typename... RBAs >
   explicit inline Band( MT& matrix, RBAs... args );

   Band( const Band& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Band() = default;
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
   inline Band& operator=( const ElementType& rhs );
   inline Band& operator=( initializer_list<ElementType> list );
   inline Band& operator=( const Band& rhs );

   template< typename VT > inline Band& operator= ( const Vector<VT,TF>& rhs );
   template< typename VT > inline Band& operator+=( const Vector<VT,TF>& rhs );
   template< typename VT > inline Band& operator-=( const Vector<VT,TF>& rhs );
   template< typename VT > inline Band& operator*=( const Vector<VT,TF>& rhs );
   template< typename VT > inline Band& operator/=( const DenseVector<VT,TF>&  rhs );
   template< typename VT > inline Band& operator%=( const Vector<VT,TF>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::band;
   using DataType::row;
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
   template< typename Other > inline Band& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename MT2, ptrdiff_t... CBAs2 >
   inline bool canAlias( const Band<MT2,TF,true,false,CBAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, ptrdiff_t... CBAs2 >
   inline bool isAliased( const Band<MT2,TF,true,false,CBAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   template< typename VT > inline void assign    ( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void assign    ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void addAssign ( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void addAssign ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void subAssign ( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void subAssign ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void multAssign( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void multAssign( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void divAssign ( const DenseVector <VT,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the band.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool TF2, bool DF2, bool MF2, ptrdiff_t... CBAs2 > friend class Band;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE  ( MT );
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
/*!\brief Constructor for bands on dense matrices.
//
// \param matrix The matrix containing the band.
// \param args The runtime band arguments.
// \exception std::invalid_argument Invalid band access index.
//
// By default, the provided band arguments are checked at runtime. In case the band is not
// properly specified (i.e. if the specified index does not correspond to a valid band in the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename... RBAs >   // Runtime band arguments
inline Band<MT,TF,true,false,CBAs...>::Band( MT& matrix, RBAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the band
{
   if( isChecked( args... ) ) {
      if( ( band() > 0L && column() >= matrix.columns() ) ||
          ( band() < 0L && row() >= matrix.rows() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( band() <= 0L || column() < matrix.columns(), "Invalid band access index" );
      BLAZE_USER_ASSERT( band() >= 0L || row() < matrix.rows(), "Invalid band access index" );
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
/*!\brief Subscript operator for the direct access to the band elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::Reference
   Band<MT,TF,true,false,CBAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid band access index" );
   return matrix_(row()+index,column()+index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the band elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstReference
   Band<MT,TF,true,false,CBAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid band access index" );
   return const_cast<const MT&>( matrix_ )(row()+index,column()+index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the band elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid band access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::Reference
   Band<MT,TF,true,false,CBAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid band access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the band elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid band access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstReference
   Band<MT,TF,true,false,CBAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid band access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the band elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense band. Note that you
// can NOT assume that the band elements lie adjacent to each other!
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::Pointer
   Band<MT,TF,true,false,CBAs...>::data() noexcept
{
   if( IsRowMajorMatrix_v<MT> )
      return matrix_.data() + row() * matrix_.spacing() + column();
   else
      return matrix_.data() + row() + column() * matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense band. Note that you
// can NOT assume that the band elements lie adjacent to each other!
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstPointer
   Band<MT,TF,true,false,CBAs...>::data() const noexcept
{
   if( IsRowMajorMatrix_v<MT> )
      return matrix_.data() + row() * matrix_.spacing() + column();
   else
      return matrix_.data() + row() + column() * matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the band.
//
// \return Iterator to the first element of the band.
//
// This function returns an iterator to the first element of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::Iterator
   Band<MT,TF,true,false,CBAs...>::begin()
{
   return Iterator( matrix_, row(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the band.
//
// \return Iterator to the first element of the band.
//
// This function returns an iterator to the first element of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstIterator
   Band<MT,TF,true,false,CBAs...>::begin() const
{
   return ConstIterator( matrix_, row(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the band.
//
// \return Iterator to the first element of the band.
//
// This function returns an iterator to the first element of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstIterator
   Band<MT,TF,true,false,CBAs...>::cbegin() const
{
   return ConstIterator( matrix_, row(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the band.
//
// \return Iterator just past the last element of the band.
//
// This function returns an iterator just past the last element of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::Iterator
   Band<MT,TF,true,false,CBAs...>::end()
{
   const size_t n( size() );
   return Iterator( matrix_, row()+n, column()+n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the band.
//
// \return Iterator just past the last element of the band.
//
// This function returns an iterator just past the last element of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstIterator
   Band<MT,TF,true,false,CBAs...>::end() const
{
   const size_t n( size() );
   return ConstIterator( matrix_, row()+n, column()+n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the band.
//
// \return Iterator just past the last element of the band.
//
// This function returns an iterator just past the last element of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,true,false,CBAs...>::ConstIterator
   Band<MT,TF,true,false,CBAs...>::cend() const
{
   const size_t n( size() );
   return ConstIterator( matrix_, row()+n, column()+n );
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
/*!\brief Homogenous assignment to all band elements.
//
// \param rhs Scalar value to be assigned to all band elements.
// \return Reference to the assigned band.
//
// This function homogeneously assigns the given value to all elements of the band. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( matrix_ ) );

   if( ( IsLower_v<MT> && column() > 0UL ) ||
       ( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> ) && row() == 0UL ) ||
       ( IsUpper_v<MT> && row() > 0UL ) ||
       ( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> ) && column() == 0UL ) )
      return *this;

   const size_t n( size() );
   for( size_t i=0UL; i<n; ++i ) {
      if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, row()+i, column()+i, rhs ) )
         left(row()+i,column()+i) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all band elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to band.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the dense
// band by means of an initializer list. The band elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the band, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate an
// invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to band" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerVector<ElementType,false> tmp( list, size() );
      if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
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
/*!\brief Copy assignment operator for Band.
//
// \param rhs Dense band to be copied.
// \return Reference to the assigned band.
// \exception std::invalid_argument Band sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two bands don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator=( const Band& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Band sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, band(), row(), column() ) ) {
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
// \return Reference to the assigned band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, band(), row(), column() ) ) {
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
// \param rhs The right-hand side vector to be added to the dense band.
// \return Reference to the assigned band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator+=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAddAssign( matrix_, right, band(), row(), column() ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the dense band.
// \return Reference to the assigned band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator-=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !trySubAssign( matrix_, right, band(), row(), column() ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the dense band.
// \return Reference to the assigned band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator*=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryMultAssign( matrix_, right, band(), row(), column() ) ) {
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator/=( const DenseVector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryDivAssign( matrix_, right, band(), row(), column() ) ) {
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::operator%=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType right( *this % (*rhs) );

   if( !tryAssign( matrix_, right, band(), row(), column() ) ) {
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
/*!\brief Returns the matrix containing the band.
//
// \return The matrix containing the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline MT& Band<MT,TF,true,false,CBAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the band.
//
// \return The matrix containing the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline const MT& Band<MT,TF,true,false,CBAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the band.
//
// \return The size of the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,true,false,CBAs...>::size() const noexcept
{
   return min( matrix_.rows() - row(), matrix_.columns() - column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the band.
//
// \return The minimum capacity of the band.
//
// This function returns the minimum capacity of the band, which corresponds to the current size
// plus padding.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,true,false,CBAs...>::spacing() const noexcept
{
   return size();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense band.
//
// \return The maximum capacity of the dense band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,true,false,CBAs...>::capacity() const noexcept
{
   return size();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the band.
//
// \return The number of non-zero elements in the band.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the band.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,true,false,CBAs...>::nonZeros() const
{
   const size_t n( size() );
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<n; ++i )
      if( !isDefault( matrix_(row()+i,column()+i) ) )
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void Band<MT,TF,true,false,CBAs...>::reset()
{
   using blaze::clear;

   if( ( IsLower_v<MT> && column() > 0UL ) ||
       ( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> ) && row() == 0UL ) ||
       ( IsUpper_v<MT> && row() > 0UL ) ||
       ( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> ) && column() == 0UL ) )
      return;

   const size_t n( size() );
   for( size_t i=0UL; i<n; ++i )
      clear( matrix_(row()+i,column()+i) );
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
/*!\brief Scaling of the band by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the band scaling.
// \return Reference to the dense band.
//
// This function scales the band by applying the given scalar value \a scalar to each element
// of the band. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a band
// on a lower or upper unitriangular matrix. The attempt to scale such a band results in a
// compile time error!
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Other >     // Data type of the scalar value
inline Band<MT,TF,true,false,CBAs...>&
   Band<MT,TF,true,false,CBAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( ( IsLower_v<MT>         && column() >  0UL ) ||
       ( IsStrictlyLower_v<MT> && row()    == 0UL ) ||
       ( IsUpper_v<MT>         && row()    >  0UL ) ||
       ( IsStrictlyUpper_v<MT> && column() == 0UL ) )
      return *this;

   const size_t n( size() );
   for( size_t i=0UL; i<n; ++i ) {
      matrix_(row()+i,column()+i) *= scalar;
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
/*!\brief Returns whether the dense band can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense band, \a false if not.
//
// This function returns whether the given address can alias with the dense band. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Other >     // Data type of the foreign expression
inline bool Band<MT,TF,true,false,CBAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense band can alias with the given dense band \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense band, \a false if not.
//
// This function returns whether the given address can alias with the dense band. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT           // Type of the dense matrix
        , bool TF               // Transpose flag
        , ptrdiff_t... CBAs >   // Compile time band arguments
template< typename MT2          // Matrix type of the foreign dense band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the foreign dense band
inline bool
   Band<MT,TF,true,false,CBAs...>::canAlias( const Band<MT2,TF,true,false,CBAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( band() == alias->band() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense band is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense band, \a false if not.
//
// This function returns whether the given address is aliased with the dense band. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Other >     // Data type of the foreign expression
inline bool Band<MT,TF,true,false,CBAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense band is aliased with the given dense band \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense band, \a false if not.
//
// This function returns whether the given address is aliased with the dense band. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT           // Type of the dense matrix
        , bool TF               // Transpose flag
        , ptrdiff_t... CBAs >   // Compile time band arguments
template< typename MT2          // Matrix type of the foreign dense band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the foreign dense band
inline bool
   Band<MT,TF,true,false,CBAs...>::isAliased( const Band<MT2,TF,true,false,CBAs2...>* alias ) const noexcept
{
   return matrix_.isAliased( &alias->matrix_ ) && ( band() == alias->band() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense band is properly aligned in memory.
//
// \return \a true in case the dense band is aligned, \a false if not.
//
// This function returns whether the dense band is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense band are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline bool Band<MT,TF,true,false,CBAs...>::isAligned() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense band can be used in SMP assignments.
//
// \return \a true in case the dense band can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense band can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense band).
*/
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline bool Band<MT,TF,true,false,CBAs...>::canSMPAssign() const noexcept
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline void Band<MT,TF,true,false,CBAs...>::assign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(row()+i    ,column()+i    ) = (*rhs)[i    ];
      matrix_(row()+i+1UL,column()+i+1UL) = (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() ) {
      matrix_(row()+ipos,column()+ipos) = (*rhs)[ipos];
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side sparse vector
inline void Band<MT,TF,true,false,CBAs...>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      matrix_(row()+index,column()+index) = element->value();
   }
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline void Band<MT,TF,true,false,CBAs...>::addAssign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(row()+i    ,column()+i    ) += (*rhs)[i    ];
      matrix_(row()+i+1UL,column()+i+1UL) += (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() ) {
      matrix_(row()+ipos,column()+ipos) += (*rhs)[ipos];
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side sparse vector
inline void Band<MT,TF,true,false,CBAs...>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      matrix_(row()+index,column()+index) += element->value();
   }
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline void Band<MT,TF,true,false,CBAs...>::subAssign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(row()+i    ,column()+i    ) -= (*rhs)[i    ];
      matrix_(row()+i+1UL,column()+i+1UL) -= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() ) {
      matrix_(row()+ipos,column()+ipos) -= (*rhs)[ipos];
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side sparse vector
inline void Band<MT,TF,true,false,CBAs...>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      matrix_(row()+index,column()+index) -= element->value();
   }
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline void Band<MT,TF,true,false,CBAs...>::multAssign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(row()+i    ,column()+i    ) *= (*rhs)[i    ];
      matrix_(row()+i+1UL,column()+i+1UL) *= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() ) {
      matrix_(row()+ipos,column()+ipos) *= (*rhs)[ipos];
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side sparse vector
inline void Band<MT,TF,true,false,CBAs...>::multAssign( const SparseVector<VT,TF>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( matrix_(row()+i,column()+i) );
      matrix_(row()+index,column()+index) *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( matrix_(row()+i,column()+i) );
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
template< typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline void Band<MT,TF,true,false,CBAs...>::divAssign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( (*rhs).size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= (*rhs).size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(row()+i    ,column()+i    ) /= (*rhs)[i    ];
      matrix_(row()+i+1UL,column()+i+1UL) /= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() ) {
      matrix_(row()+ipos,column()+ipos) /= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DENSE MATRIX MULTIPLICATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Band for dense matrix multiplications.
// \ingroup band
//
// This specialization of Band adapts the class template to the requirements of dense matrix
// multiplications.
*/
template< typename MT          // Type of the dense matrix multiplication
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
class Band<MT,TF,true,true,CBAs...>
   : public View< DenseVector< Band<MT,TF,true,true,CBAs...>, TF > >
   , private BandData<CBAs...>
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   //! The type of the BandData base class.
   using DataType = BandData<CBAs...>;

   //! The type of the left-hand side matrix operand.
   using LeftOperand = RemoveReference_t< LeftOperand_t<MT> >;

   //! The type of the right-hand side matrix operand.
   using RightOperand = RemoveReference_t< RightOperand_t<MT> >;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Band instance.
   using This = Band<MT,TF,true,true,CBAs...>;

   //! Base type of this Band instance.
   using BaseType = View< DenseVector<This,TF> >;

   //! The type viewed by this Band instance.
   using ViewedType = MT;

   //! Result type for expression template evaluations.
   using ResultType = BandTrait_t<ResultType_t<MT>,CBAs...>;

   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< RequiresEvaluation_v<LeftOperand> ||
                               RequiresEvaluation_v<RightOperand>
                             , const ResultType, const Band& >;

   //! Type for the assignment of the left-hand side matrix operand.
   using LT = If_t< IsSparseMatrix_v<LeftOperand> && IsColumnMajorMatrix_v<LeftOperand>
                  , OppositeType_t<LeftOperand>
                  , CompositeType_t<LeftOperand> >;

   //! Type for the assignment of the right-hand side matrix operand.
   using RT = If_t< IsSparseMatrix_v<RightOperand> && IsRowMajorMatrix_v<RightOperand>
                  , OppositeType_t<RightOperand>
                  , CompositeType_t<RightOperand> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for bands on dense matrix multiplications.
   //
   // \param mmm The matrix multiplication containing the band.
   // \param args The runtime band arguments.
   // \exception std::invalid_argument Invalid band access index.
   */
   template< typename... RBAs >  // Runtime band arguments
   explicit inline Band( const MT& mmm, RBAs... args )
      : DataType( args... )  // Base class initialization
      , matrix_ ( mmm )      // The matrix multiplication containing the band
   {
      if( isChecked( args... ) ) {
         if( ( band() > 0L && column() >= mmm.columns() ) ||
             ( band() < 0L && row() >= mmm.rows() ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
         }
      }
      else {
         BLAZE_USER_ASSERT( band() <= 0L || column() < mmm.columns(), "Invalid band access index" );
         BLAZE_USER_ASSERT( band() >= 0L || row() < mmm.rows(), "Invalid band access index" );
      }
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return matrix_(row()+index,column()+index);
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
      return min( matrix_.rows() - row(), matrix_.columns() - column() );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::band;
   using DataType::row;
   using DataType::column;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the matrix multiplication expression containing the band.
   //
   // \return The matrix multiplication expression containing the band.
   */
   inline const MT& operand() const noexcept {
      return matrix_;
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
      return matrix_.isAliased( &unview( *alias ) );
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
      return matrix_.isAliased( &unview( *alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   constexpr bool isAligned() const noexcept {
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   MT matrix_;  //!< The matrix multiplication containing the band.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a band view on a dense matrix multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a band view on a dense
   // matrix multiplication to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      using blaze::row;
      using blaze::column;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LT A( serial( rhs.operand().leftOperand()  ) );
      RT B( serial( rhs.operand().rightOperand() ) );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] = row( A, rhs.row()+i, unchecked ) * column( B, rhs.column()+i, unchecked );
      }
   }
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a band view on a dense matrix multiplication to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side band view to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a band view on a dense
   // matrix multiplication to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const Band& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a band view on a dense matrix multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a band view on a
   // dense matrix multiplication to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      using blaze::row;
      using blaze::column;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LT A( serial( rhs.operand().leftOperand()  ) );
      RT B( serial( rhs.operand().rightOperand() ) );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] += row( A, rhs.row()+i, unchecked ) * column( B, rhs.column()+i, unchecked );
      }
   }
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a band view on a dense matrix multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a band view
   // on a dense matrix multiplication to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      using blaze::row;
      using blaze::column;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LT A( serial( rhs.operand().leftOperand()  ) );
      RT B( serial( rhs.operand().rightOperand() ) );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] -= row( A, rhs.row()+i, unchecked ) * column( B, rhs.column()+i, unchecked );
      }
   }
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a band view on a dense matrix multiplication to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a band view
   // on a dense matrix multiplication to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      using blaze::row;
      using blaze::column;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LT A( serial( rhs.operand().leftOperand()  ) );
      RT B( serial( rhs.operand().rightOperand() ) );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] *= row( A, rhs.row()+i, unchecked ) * column( B, rhs.column()+i, unchecked );
      }
   }
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*!\brief Division assignment of a band view on a dense matrix multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a band view on a
   // dense matrix multiplication to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      using blaze::row;
      using blaze::column;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LT A( serial( rhs.operand().leftOperand()  ) );
      RT B( serial( rhs.operand().rightOperand() ) );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] /= row( A, rhs.row()+i, unchecked ) * column( B, rhs.column()+i, unchecked );
      }
   }
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATMATMULTEXPR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DENSE MATRIX REPEATER EXPRESSIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Band for dense matrix repeater expressions.
// \ingroup band
//
// This specialization of Band adapts the class template to the requirements of dense matrix
// repeater expressions.
*/
template< typename MT          // Type of the dense matrix
        , bool SO              // Storage order
        , size_t... CRAs       // Compile time repeater arguments
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
class Band< DMatRepeatExpr<MT,SO,CRAs...>, TF, true, false, CBAs... >
   : public View< DenseVector< Band< DMatRepeatExpr<MT,SO,CRAs...>, TF, true, false, CBAs... >, TF > >
   , private BandData<CBAs...>
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense matrix repeater expression.
   using RE = DMatRepeatExpr<MT,SO,CRAs...>;

   //! The type of the BandData base class.
   using DataType = BandData<CBAs...>;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Band instance.
   using This = Band<RE,TF,true,false,CBAs...>;

   //! Base type of this Band instance.
   using BaseType = View< DenseVector<This,TF> >;

   //! The type viewed by this Band instance.
   using ViewedType = RE;

   //! Result type for expression template evaluations.
   using ResultType = BandTrait_t<ResultType_t<RE>,CBAs...>;

   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for bands on dense matrix repeater expressions.
   //
   // \param re The matrix repeater expression containing the band.
   // \param args The runtime band arguments.
   // \exception std::invalid_argument Invalid band access index.
   */
   template< typename... RBAs >  // Runtime band arguments
   explicit inline Band( const RE& re, RBAs... args )
      : DataType( args... )  // Base class initialization
      , matrix_ ( re )       // The matrix repeater expression containing the band
   {
      if( isChecked( args... ) ) {
         if( ( band() > 0L && column() >= re.columns() ) ||
             ( band() < 0L && row() >= re.rows() ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
         }
      }
      else {
         BLAZE_USER_ASSERT( band() <= 0L || column() < re.columns(), "Invalid band access index" );
         BLAZE_USER_ASSERT( band() >= 0L || row() < re.rows(), "Invalid band access index" );
      }
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return matrix_(row()+index,column()+index);
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
      return min( matrix_.rows() - row(), matrix_.columns() - column() );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::band;
   using DataType::row;
   using DataType::column;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the matrix repeater expression containing the band.
   //
   // \return The matrix repeater expression containing the band.
   */
   inline const MT& operand() const noexcept {
      return matrix_;
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
      return matrix_.isAliased( &unview( *alias ) );
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
      return matrix_.isAliased( &unview( *alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   constexpr bool isAligned() const noexcept {
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   RE matrix_;  //!< The matrix repeater expression containing the band.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
