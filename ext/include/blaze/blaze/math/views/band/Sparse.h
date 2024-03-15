//=================================================================================================
/*!
//  \file blaze/math/views/band/Sparse.h
//  \brief Band specialization for sparse matrices
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

#ifndef _BLAZE_MATH_VIEWS_BAND_SPARSE_H_
#define _BLAZE_MATH_VIEWS_BAND_SPARSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/views/band/BandData.h>
#include <blaze/math/views/band/BaseTemplate.h>
#include <blaze/math/views/Check.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Band for sparse matrices.
// \ingroup band
//
// This specialization of Band adapts the class template to the requirements of sparse matrices.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
class Band<MT,TF,false,false,CBAs...>
   : public View< SparseVector< Band<MT,TF,false,false,CBAs...>, TF > >
   , private BandData<CBAs...>
{
 private:
   //**Type definitions****************************************************************************
   using RT       = ResultType_t<MT>;                     //!< Result type of the sparse matrix expression.
   using DataType = BandData<CBAs...>;                    //!< The type of the BandData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Band instance.
   using This = Band<MT,TF,false,false,CBAs...>;

   //! Base type of this Band instance.
   using BaseType = View< SparseVector<This,TF> >;

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
   //**********************************************************************************************

   //**BandElement class definition****************************************************************
   /*!\brief Access proxy for a specific element of the sparse band.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class BandElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the BandElement class.
      //
      // \param pos Iterator to the current position within the sparse band.
      // \param index The index of the element within the band.
      */
      inline BandElement( IteratorType pos, size_t index )
         : pos_  ( pos   )  // Iterator to the current position within the sparse band
         , index_( index )  // Index of the element
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse band element.
      //
      // \param value The new value of the sparse band element.
      // \return Reference to the sparse band element.
      */
      template< typename T > inline BandElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse band element.
      //
      // \param value The right-hand side value for the addition.
      // \return Reference to the sparse band element.
      */
      template< typename T > inline BandElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse band element.
      //
      // \param value The right-hand side value for the subtraction.
      // \return Reference to the sparse band element.
      */
      template< typename T > inline BandElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse band element.
      //
      // \param value The right-hand side value for the multiplication.
      // \return Reference to the sparse band element.
      */
      template< typename T > inline BandElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse band element.
      //
      // \param value The right-hand side value for the division.
      // \return Reference to the sparse band element.
      */
      template< typename T > inline BandElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const BandElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse band element.
      //
      // \return The current value of the sparse band element.
      */
      inline decltype(auto) value() const {
         return pos_->value();
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return index_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse band.
      size_t index_;      //!< Index of the element.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**BandIterator class definition***************************************************************
   /*!\brief Iterator over the elements of the sparse band.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class BandIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;             //!< The iterator category.
      using ValueType        = BandElement<MatrixType,IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                             //!< Pointer return type.
      using ReferenceType    = ValueType;                             //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                             //!< Difference between two iterators.

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
      inline BandIterator()
         : matrix_( nullptr )  // The sparse matrix containing the band
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the BandIterator class.
      //
      // \param matrix The matrix containing the band.
      // \param rowIndex The initial row index.
      // \param columnIndex The initial column index.
      */
      inline BandIterator( MatrixType& matrix, size_t rowIndex, size_t columnIndex )
         : matrix_( &matrix     )  // The sparse matrix containing the band
         , row_   ( rowIndex    )  // The current row index
         , column_( columnIndex )  // The current column index
         , pos_   ()               // Iterator to the current sparse element
      {
         for( ; row_ < matrix_->rows() && column_ < matrix_->columns(); ++row_, ++column_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( IsRowMajorMatrix_v<MatrixType> ? row_ : column_ ) )
               break;
         }
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the BandIterator class.
      //
      // \param matrix The matrix containing the band.
      // \param rowIndex The initial row index.
      // \param columnIndex The initial column index.
      // \param pos Initial position of the iterator
      */
      inline BandIterator( MatrixType& matrix, size_t rowIndex, size_t columnIndex, IteratorType pos )
         : matrix_( &matrix     )  // The sparse matrix containing the band
         , row_   ( rowIndex    )  // The current row index
         , column_( columnIndex )  // The current column index
         , pos_   ( pos         )  // Iterator to the current sparse element
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row_, column_ ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different BandIterator instances.
      //
      // \param it The band iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline BandIterator( const BandIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the band
         , row_   ( it.row_    )  // The current row index
         , column_( it.column_ )  // The current column index
         , pos_   ( it.pos_    )  // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BandIterator& operator++() {
         ++row_;
         ++column_;

         for( ; row_ < matrix_->rows() && column_ < matrix_->columns(); ++row_, ++column_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( IsRowMajorMatrix_v<MatrixType> ? row_ : column_ ) )
               break;
         }

         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const BandIterator operator++( int ) {
         const BandIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, min( row_, column_ ) );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, min( row_, column_ ) );
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

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two band iterators.
      //
      // \param rhs The right-hand side band iterator.
      // \return The number of elements between the two band iterators.
      */
      inline DifferenceType operator-( const BandIterator& rhs ) const {
         size_t counter( 0UL );
         size_t row( rhs.row_ );
         size_t column( rhs.column_ );

         for( ; row<row_; ++row, ++column ) {
            const auto end( matrix_->end( IsRowMajorMatrix_v<MatrixType> ? row : column ) );
            if( matrix_->find( row, column ) != end )
               ++counter;
         }

         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The sparse matrix containing the band.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class BandIterator;
      template< typename MT2, bool DF2, bool TF2, bool MF2, ptrdiff_t... CBAs2 > friend class Band;
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
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

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
                           inline Band& operator= ( initializer_list<ElementType> list );
                           inline Band& operator= ( const Band& rhs );
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
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   reserve( size_t n );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set   ( size_t index, const ElementType& value );
   inline Iterator insert( size_t index, const ElementType& value );
   inline void     append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t index );
   inline Iterator erase( Iterator pos );
   inline Iterator erase( Iterator first, Iterator last );

   template< typename Pred, typename = DisableIf_t< IsIntegral_v<Pred> > >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
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
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT >    inline void assign   ( const DenseVector <VT,TF>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,TF>& rhs );
   template< typename VT >    inline void addAssign( const Vector<VT,TF>& rhs );
   template< typename VT >    inline void subAssign( const Vector<VT,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the band.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE  ( MT );
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
/*!\brief Constructor for bands on sparse matrices.
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename... RBAs >   // Runtime band arguments
inline Band<MT,TF,false,false,CBAs...>::Band( MT& matrix, RBAs... args )
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Reference
   Band<MT,TF,false,false,CBAs...>::operator[]( size_t index )
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstReference
   Band<MT,TF,false,false,CBAs...>::operator[]( size_t index ) const
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Reference
   Band<MT,TF,false,false,CBAs...>::at( size_t index )
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstReference
   Band<MT,TF,false,false,CBAs...>::at( size_t index ) const
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
/*!\brief Returns an iterator to the first element of the band.
//
// \return Iterator to the first element of the band.
//
// This function returns an iterator to the first element of the band.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::begin()
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::begin() const
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::cbegin() const
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::end()
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::end() const
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::cend() const
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
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to band" );
   }

   const InitializerVector<ElementType,false> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Band.
//
// \param rhs Sparse band to be copied.
// \return Reference to the assigned band.
// \exception std::invalid_argument Band sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two bands don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator=( const Band& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && band() == rhs.band() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      assign( left, tmp );
   }
   else {
      assign( left, rhs );
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const CompositeType_t<VT> tmp( *rhs );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the sparse band.
// \return Reference to the sparse band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator+=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the sparse band.
// \return Reference to the sparse band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator-=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

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
// \param rhs The right-hand side vector to be multiplied with the sparse band.
// \return Reference to the sparse band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator*=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( MultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

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
// \return Reference to the sparse band.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator/=( const DenseVector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( DivType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( DivType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

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
// \return Reference to the sparse band.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::operator%=( const Vector<VT,TF>& rhs )
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

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, band(), row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline MT& Band<MT,TF,false,false,CBAs...>::operand() noexcept
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline const MT& Band<MT,TF,false,false,CBAs...>::operand() const noexcept
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,false,false,CBAs...>::size() const noexcept
{
   return min( matrix_.rows() - row(), matrix_.columns() - column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse band.
//
// \return The maximum capacity of the sparse band.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,false,false,CBAs...>::capacity() const noexcept
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline size_t Band<MT,TF,false,false,CBAs...>::nonZeros() const
{
   return end() - begin();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void Band<MT,TF,false,false,CBAs...>::reset()
{
   using blaze::clear;

   if( ( IsLower_v<MT> && column() > 0UL ) ||
       ( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> ) && row() == 0UL ) ||
       ( IsUpper_v<MT> && row() > 0UL ) ||
       ( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> ) && column() == 0UL ) )
      return;

   const size_t n( size() );
   for( size_t i=0UL; i<n; ++i )
      matrix_.erase( row()+i, column()+i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse band.
//
// \param n The new minimum capacity of the sparse band.
// \return void
//
// This function increases the capacity of the sparse band to at least \a n elements. The
// current values of the band elements are preserved.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
void Band<MT,TF,false,false,CBAs...>::reserve( size_t n )
{
   MAYBE_UNUSED( n );

   return;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse band.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse band. In case the sparse band already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::set( size_t index, const ElementType& value )
{
   const size_t rowIndex   ( row()    + index );
   const size_t columnIndex( column() + index );
   return Iterator( matrix_, rowIndex, columnIndex, matrix_.set( rowIndex, columnIndex, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse band.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse band access index.
//
// This function inserts a new element into the sparse band. However, duplicate elements
// are not allowed. In case the sparse band already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::insert( size_t index, const ElementType& value )
{
   const size_t rowIndex   ( row()    + index );
   const size_t columnIndex( column() + index );
   return Iterator( matrix_, rowIndex, columnIndex, matrix_.insert( rowIndex, columnIndex, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse band.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse band with elements. It appends
// a new element to the end of the sparse band without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse band
//  - the current number of non-zero elements must be smaller than the capacity of the band
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void Band<MT,TF,false,false,CBAs...>::append( size_t index, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( row()+index, column()+index, value );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse band.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse band.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void Band<MT,TF,false,false,CBAs...>::erase( size_t index )
{
   matrix_.erase( row()+index, column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse band.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse band.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::erase( Iterator pos )
{
   const size_t rowIndex   ( pos.row_    );
   const size_t columnIndex( pos.column_ );

   if( rowIndex == matrix_.rows() || columnIndex == matrix_.columns() )
      return pos;

   matrix_.erase( ( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ), pos.pos_ );
   return Iterator( matrix_, rowIndex+1UL, columnIndex+1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse band.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse band.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::erase( Iterator first, Iterator last )
{
   for( ; first!=last; ++first ) {
      const size_t index( IsRowMajorMatrix<MT>::value ? first.row_ : first.column_ );
      matrix_.erase( index, first.pos_ );
   }
   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse band.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse band. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto band2 = band( A, 2UL );
   band2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Pred        // Type of the unary predicate
        , typename >           // Type restriction on the unary predicate
inline void Band<MT,TF,false,false,CBAs...>::erase( Pred predicate )
{
   for( Iterator element=begin(); element!=end(); ++element ) {
      if( predicate( element->value() ) ) {
         const size_t index( IsRowMajorMatrix<MT>::value ? element.row_ : element.column_ );
         matrix_.erase( index, element.pos_ );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse band.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse band. The
// elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto band2 = band( A, 2UL );
   band2.erase( band2.begin(), band2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Pred >      // Type of the unary predicate
inline void Band<MT,TF,false,false,CBAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   for( ; first!=last; ++first ) {
      if( predicate( first->value() ) ) {
         const size_t index( IsRowMajorMatrix<MT>::value ? first.row_ : first.column_ );
         matrix_.erase( index, first.pos_ );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific band element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// band. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse band (the end() iterator) is returned. Note that
// the returned sparse band iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::find( size_t index )
{
   const size_t rowIndex   ( row()+index    );
   const size_t columnIndex( column()+index );
   const Iterator_t<MT> pos( matrix_.find( rowIndex, columnIndex ) );

   if( pos != matrix_.end( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ) )
      return Iterator( matrix_, rowIndex, columnIndex, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific band element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// band. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse band (the end() iterator) is returned. Note that
// the returned sparse band iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::find( size_t index ) const
{
   const size_t rowIndex   ( row()+index    );
   const size_t columnIndex( column()+index );
   const ConstIterator_t<MT> pos( matrix_.find( rowIndex, columnIndex ) );

   if( pos != matrix_.end( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ) )
      return ConstIterator( matrix_, rowIndex, columnIndex, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse band iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::lowerBound( size_t index )
{
   for( size_t i=index; i<size(); ++i )
   {
      const size_t rowIndex   ( row()+i    );
      const size_t columnIndex( column()+i );
      const Iterator_t<MT> pos( matrix_.find( rowIndex, columnIndex ) );

      if( pos != matrix_.end( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ) )
         return Iterator( matrix_, rowIndex, columnIndex, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse band iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::lowerBound( size_t index ) const
{
   for( size_t i=index; i<size(); ++i )
   {
      const size_t rowIndex   ( row()+i    );
      const size_t columnIndex( column()+i );
      const ConstIterator_t<MT> pos( matrix_.find( rowIndex, columnIndex ) );

      if( pos != matrix_.end( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ) )
         return ConstIterator( matrix_, rowIndex, columnIndex, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse band iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::Iterator
   Band<MT,TF,false,false,CBAs...>::upperBound( size_t index )
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const size_t rowIndex   ( row()+i    );
      const size_t columnIndex( column()+i );
      const Iterator_t<MT> pos( matrix_.find( rowIndex, columnIndex ) );

      if( pos != matrix_.end( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ) )
         return Iterator( matrix_, rowIndex, columnIndex, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse band iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline typename Band<MT,TF,false,false,CBAs...>::ConstIterator
   Band<MT,TF,false,false,CBAs...>::upperBound( size_t index ) const
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const size_t rowIndex   ( row()+i    );
      const size_t columnIndex( column()+i );
      const ConstIterator_t<MT> pos( matrix_.find( rowIndex, columnIndex ) );

      if( pos != matrix_.end( IsRowMajorMatrix_v<MT> ? rowIndex : columnIndex ) )
         return ConstIterator( matrix_, rowIndex, columnIndex, pos );
   }

   return end();
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Other >     // Data type of the scalar value
inline Band<MT,TF,false,false,CBAs...>&
   Band<MT,TF,false,false,CBAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( ( IsLower_v<MT>         && column() >  0UL ) ||
       ( IsStrictlyLower_v<MT> && row()    == 0UL ) ||
       ( IsUpper_v<MT>         && row()    >  0UL ) ||
       ( IsStrictlyUpper_v<MT> && column() == 0UL ) )
      return *this;

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;

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
/*!\brief Returns whether the sparse band can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse band, \a false if not.
//
// This function returns whether the given address can alias with the sparse band. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Other >     // Data type of the foreign expression
inline bool Band<MT,TF,false,false,CBAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse band is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse band, \a false if not.
//
// This function returns whether the given address is aliased with the sparse band. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename Other >     // Data type of the foreign expression
inline bool Band<MT,TF,false,false,CBAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side dense vector
inline void Band<MT,TF,false,false,CBAs...>::assign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      matrix_(row()+i,column()+i) = (*rhs)[i];
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
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side sparse vector
inline void Band<MT,TF,false,false,CBAs...>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      for( ; i<element->index(); ++i )
         matrix_.erase( row()+i, column()+i );
      matrix_(row()+i,column()+i) = element->value();
      ++i;
   }
   for( ; i<size(); ++i ) {
      matrix_.erase( row()+i, column()+i );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a vector.
//
// \param rhs The right-hand side vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline void Band<MT,TF,false,false,CBAs...>::addAssign( const Vector<VT,TF>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a vector.
//
// \param rhs The right-hand side vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
template< typename VT >        // Type of the right-hand side vector
inline void Band<MT,TF,false,false,CBAs...>::subAssign( const Vector<VT,TF>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SPARSE MATRIX MULTIPLICATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Band for sparse matrix multiplications.
// \ingroup band
//
// This specialization of Band adapts the class template to the requirements of sparse matrix
// multiplications.
*/
template< typename MT          // Type of the sparse matrix multiplication
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
class Band<MT,TF,false,true,CBAs...>
   : public View< SparseVector< Band<MT,TF,false,true,CBAs...>, TF > >
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
   using This = Band<MT,TF,false,true,CBAs...>;

   //! Base type of this Band instance.
   using BaseType = View< SparseVector<This,TF> >;

   //! The type viewed by this Band instance.
   using ViewedType = MT;

   //! Result type for expression template evaluations.
   using ResultType = BandTrait_t<ResultType_t<MT>,CBAs...>;

   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

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
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for bands on sparse matrix multiplications.
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

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the vector.
   //
   // \return The number of non-zero elements in the vector.
   */
   inline size_t nonZeros() const {
      return 0UL;
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

 private:
   //**Member variables****************************************************************************
   MT matrix_;  //!< The matrix multiplication containing the band.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a band view on a sparse matrix multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a band view on a sparse
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
   /*!\brief Assignment of a band view on a sparse matrix multiplication to a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side band view to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a band view on a sparse
   // matrix multiplication to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const Band& rhs )
   {
      using blaze::row;
      using blaze::column;

      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LT A( serial( rhs.operand().leftOperand()  ) );
      RT B( serial( rhs.operand().rightOperand() ) );

      const size_t n( rhs.size() );
      ElementType_t<VT> tmp{};
      size_t nonzeros( 0UL );

      for( size_t i=0UL; i<n; ++i ) {
         tmp = row( A, rhs.row()+i, unchecked ) * column( B, rhs.column()+i, unchecked );
         if( !isDefault<strict>( tmp ) ) {
            if( (*lhs).capacity() <= nonzeros ) {
               (*lhs).reserve( min( max( 2UL*(*lhs).capacity(), 7UL ), (*lhs).size() ) );
            }
            (*lhs).append( i, tmp, false );
            ++nonzeros;
         }
      }
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a band view on a dense matrix multiplication to a dense vector.
   // \ingroup sparse_vector
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
   // \ingroup sparse_vector
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
   // \ingroup sparse_vector
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

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATMATMULTEXPR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SPARSE MATRIX REPEATER EXPRESSIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Band for sparse matrix repeater expressions.
// \ingroup band
//
// This specialization of Band adapts the class template to the requirements of sparse matrix
// repeater expressions.
*/
template< typename MT          // Type of the sparse matrix
        , bool SO              // Storage order
        , size_t... CRAs       // Compile time repeater arguments
        , bool TF              // Transpose flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
class Band< SMatRepeatExpr<MT,SO,CRAs...>, TF, false, false, CBAs... >
   : public View< SparseVector< Band< SMatRepeatExpr<MT,SO,CRAs...>, TF, false, false, CBAs... >, TF > >
   , private BandData<CBAs...>
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the sparse matrix repeater expression.
   using RE = SMatRepeatExpr<MT,SO,CRAs...>;

   //! The type of the BandData base class.
   using DataType = BandData<CBAs...>;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Band instance.
   using This = Band<RE,TF,false,false,CBAs...>;

   //! Base type of this Band instance.
   using BaseType = View< SparseVector<This,TF> >;

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
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for bands on sparse matrix repeater expressions.
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

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the vector.
   //
   // \return The number of non-zero elements in the vector.
   */
   inline size_t nonZeros() const {
      return 0UL;
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

 private:
   //**Member variables****************************************************************************
   RE matrix_;  //!< The matrix repeater expression containing the band.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a band view on a sparse matrix repeater expression to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a band view on a sparse
   // matrix repeater expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] = rhs[i];
      }
   }
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a band view on a sparse matrix repeater expression to a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side band view to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a band view on a sparse
   // matrix repeater expression to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const Band& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const size_t n( rhs.size() );
      ElementType_t<VT> tmp{};
      size_t nonzeros( 0UL );

      for( size_t i=0UL; i<n; ++i ) {
         tmp = rhs[i];
         if( !isDefault<strict>( tmp ) ) {
            if( (*lhs).capacity() <= nonzeros ) {
               (*lhs).reserve( min( max( 2UL*(*lhs).capacity(), 7UL ), (*lhs).size() ) );
            }
            (*lhs).append( i, tmp, false );
            ++nonzeros;
         }
      }
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a band view on a dense matrix repeater expression to a dense
   //        vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a band view on a
   // dense matrix repeater expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] += rhs[i];
      }
   }
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a band view on a dense matrix repeater expression to a
   //        dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a band view
   // on a dense matrix repeater expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] -= rhs[i];
      }
   }
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a band view on a dense matrix repeater expression to a
   //        dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side band view to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a band view
   // on a dense matrix repeater expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT,TF>& lhs, const Band& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const size_t n( rhs.size() );
      for( size_t i=0UL; i<n; ++i ) {
         (*lhs)[i] *= rhs[i];
      }
   }
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
