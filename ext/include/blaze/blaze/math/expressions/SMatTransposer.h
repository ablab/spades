//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatTransposer.h
//  \brief Header file for the sparse matrix transposer
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATTRANSPOSER_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATTRANSPOSER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATTRANSPOSER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the transposition of a sparse matrix.
// \ingroup sparse_matrix_expression
//
// The SMatTransposer class is a wrapper object for the temporary transposition of a sparse matrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
class SMatTransposer
   : public SparseMatrix< SMatTransposer<MT,SO>, SO >
{
 public:
   //**Type definitions****************************************************************************
   using This           = SMatTransposer<MT,SO>;  //!< Type of this SMatTransposer instance.
   using BaseType       = SparseMatrix<This,SO>;  //!< Base type of this SMatTransposer instance.
   using ResultType     = TransposeType_t<MT>;    //!< Result type for expression template evaluations.
   using OppositeType   = MT;                     //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType  = ResultType_t<MT>;       //!< Transpose type for expression template evaluations.
   using ElementType    = ElementType_t<MT>;      //!< Resulting element type.
   using ReturnType     = ReturnType_t<MT>;       //!< Return type for expression template evaluations.
   using CompositeType  = const This&;            //!< Data type for composite expression templates.
   using Reference      = Reference_t<MT>;        //!< Reference to a non-constant matrix value.
   using ConstReference = ConstReference_t<MT>;   //!< Reference to a constant matrix value.
   using Iterator       = Iterator_t<MT>;         //!< Iterator over non-constant elements.
   using ConstIterator  = ConstIterator_t<MT>;    //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatTransposer class.
   //
   // \param sm The sparse matrix operand.
   */
   explicit inline SMatTransposer( MT& sm ) noexcept
      : sm_( sm )  // The sparse matrix operand
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline ConstReference operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < sm_.columns(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < sm_.rows()   , "Invalid column access index" );
      return const_cast<const MT&>( sm_ )(j,i);
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
   inline ConstReference at( size_t i, size_t j ) const {
      if( i >= sm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= sm_.rows() ) {
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
   //
   // This function returns a row/column iterator to the first non-zero element of row/column \a i.
   // In case the storage order is set to \a rowMajor the function returns an iterator to the first
   // non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
   // returns an iterator to the first non-zero element of column \a i.
   */
   inline Iterator begin( size_t i ) {
      return sm_.begin(i);
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   //
   // This function returns a row/column iterator to the first non-zero element of row/column \a i.
   // In case the storage order is set to \a rowMajor the function returns an iterator to the first
   // non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
   // returns an iterator to the first non-zero element of column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return sm_.cbegin(i);
   }
   //**********************************************************************************************

   //**Cbegin function*****************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   //
   // This function returns a row/column iterator to the first non-zero element of row/column \a i.
   // In case the storage order is set to \a rowMajor the function returns an iterator to the first
   // non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
   // returns an iterator to the first non-zero element of column \a i.
   */
   inline ConstIterator cbegin( size_t i ) const {
      return sm_.cbegin(i);
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   //
   // This function returns an row/column iterator just past the last non-zero element of row/column
   // \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
   // past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
   // the function returns an iterator just past the last non-zero element of column \a i.
   */
   inline Iterator end( size_t i ) {
      return sm_.end(i);
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   //
   // This function returns an row/column iterator just past the last non-zero element of row/column
   // \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
   // past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
   // the function returns an iterator just past the last non-zero element of column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return sm_.cend(i);
   }
   //**********************************************************************************************

   //**Cend function*******************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   //
   // This function returns an row/column iterator just past the last non-zero element of row/column
   // \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
   // past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
   // the function returns an iterator just past the last non-zero element of column \a i.
   */
   inline ConstIterator cend( size_t i ) const {
      return sm_.cend(i);
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return sm_.columns();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return sm_.rows();
   }
   //**********************************************************************************************

   //**Capacity function***************************************************************************
   /*!\brief Returns the maximum capacity of the matrix.
   //
   // \return The capacity of the matrix.
   */
   inline size_t capacity() const noexcept {
      return sm_.capacity();
   }
   //**********************************************************************************************

   //**Capacity function***************************************************************************
   /*!\brief Returns the current capacity of the specified row/column.
   //
   // \param i The index of the row/column.
   // \return The current capacity of row/column \a i.
   */
   inline size_t capacity( size_t i ) const noexcept {
      return sm_.capacity( i );
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the matrix.
   //
   // \return The number of non-zero elements in the matrix.
   */
   inline size_t nonZeros() const {
      return sm_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row/column.
   //
   // \param i The index of the row/column.
   // \return The number of non-zero elements of row/column \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return sm_.nonZeros( i );
   }
   //**********************************************************************************************

   //**Reset function******************************************************************************
   /*!\brief Resets the matrix elements.
   //
   // \return void
   */
   inline void reset() {
      return sm_.reset();
   }
   //**********************************************************************************************

   //**Insert function*****************************************************************************
   /*!\brief Inserting an element into the sparse matrix.
   //
   // \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
   // \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
   // \param value The value of the element to be inserted.
   // \return Iterator to the newly inserted element.
   // \exception std::invalid_argument Invalid sparse matrix access index.
   //
   // This function insert a new element into the sparse matrix. However, duplicate elements are
   // not allowed. In case the sparse matrix already contains an element with row index \a i and
   // column index \a j, a \a std::invalid_argument exception is thrown.
   */
   inline Iterator insert( size_t i, size_t j, const ElementType& value ) {
      return sm_.insert( j, i, value );
   }
   //**********************************************************************************************

   //**Reserve function****************************************************************************
   /*!\brief Setting the minimum capacity of the sparse matrix.
   //
   // \param nonzeros The new minimum capacity of the sparse matrix.
   // \return void
   //
   // This function increases the capacity of the sparse matrix to at least \a nonzeros elements.
   // The current values of the matrix elements and the individual capacities of the matrix rows
   // are preserved.
   */
   inline void reserve( size_t nonzeros ) {
      sm_.reserve( nonzeros );
   }
   //**********************************************************************************************

   //**Reserve function****************************************************************************
   /*!\brief Setting the minimum capacity of a specific row/column of the sparse matrix.
   //
   // \param i The row/column index of the new element \f$[0..M-1]\f$ or \f$[0..N-1]\f$.
   // \param nonzeros The new minimum capacity of the specified row.
   // \return void
   //
   // This function increases the capacity of row/column \a i of the sparse matrix to at least
   // \a nonzeros elements. The current values of the sparse matrix and all other individual
   // row/column capacities are preserved. In case the storage order is set to \a rowMajor, the
   // function reserves capacity for row \a i and the index has to be in the range \f$[0..M-1]\f$.
   // In case the storage order is set to \a columnMajor, the function reserves capacity for column
   // \a i and the index has to be in the range \f$[0..N-1]\f$.
   */
   inline void reserve( size_t i, size_t nonzeros ) {
      sm_.reserve( i, nonzeros );
   }
   //**********************************************************************************************

   //**Append function*****************************************************************************
   /*!\brief Appending an element to the specified row/column of the sparse matrix.
   //
   // \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
   // \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
   // \param value The value of the element to be appended.
   // \param check \a true if the new value should be checked for default values, \a false if not.
   // \return void
   //
   // This function provides a very efficient way to fill a sparse matrix with elements. It
   // appends a new element to the end of the specified row/column without any additional
   // memory allocation. Therefore it is strictly necessary to keep the following preconditions
   // in mind:
   //
   //  - the index of the new element must be strictly larger than the largest index of
   //    non-zero elements in the specified row/column of the sparse matrix
   //  - the current number of non-zero elements in row/column \a i must be smaller than
   //    the capacity of row/column \a i.
   //
   // Ignoring these preconditions might result in undefined behavior! The optional \a check
   // parameter specifies whether the new value should be tested for a default value. If the new
   // value is a default value (for instance 0 in case of an integral element type) the value is
   // not appended. Per default the values are not tested.
   //
   // \note Although append() does not allocate new memory, it still invalidates all iterators
   // returned by the end() functions!
   */
   inline void append( size_t i, size_t j, const ElementType& value, bool check=false ) {
      sm_.append( j, i, value, check );
   }
   //**********************************************************************************************

   //**Finalize function***************************************************************************
   /*!\brief Finalizing the element insertion of a row/column.
   //
   // \param i The index of the row/column to be finalized \f$[0..M-1]\f$.
   // \return void
   //
   // This function is part of the low-level interface to efficiently fill the matrix with elements.
   // After completion of row/column \a i via the append() function, this function can be called to
   // finalize row/column \a i and prepare the next row/column for insertion process via append().
   //
   // \note Although finalize() does not allocate new memory, it still invalidates all iterators
   // returned by the end() functions!
   */
   inline void finalize( size_t i ) {
      sm_.finalize( i );
   }
   //**********************************************************************************************

   //**IsIntact function***************************************************************************
   /*!\brief Returns whether the invariants of the matrix are intact.
   //
   // \return \a true in case the matrix's invariants are intact, \a false otherwise.
   */
   inline bool isIntact() const noexcept {
      using blaze::isIntact;
      return isIntact( sm_ );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the matrix can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this matrix, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool canAlias( const Other* alias ) const noexcept
   {
      return sm_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the matrix is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this matrix, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool isAliased( const Other* alias ) const noexcept
   {
      return sm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**CanSMPAssign function***********************************************************************
   /*!\brief Returns whether the matrix can be used in SMP assignments.
   //
   // \return \a true in case the matrix can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept
   {
      return sm_.canSMPAssign();
   }
   //**********************************************************************************************

   //**Transpose assignment of matrices************************************************************
   /*!\brief Implementation of the transpose assignment of a matrix.
   //
   // \param rhs The right-hand side matrix to be assigned.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename MT2  // Type of the right-hand side matrix
           , bool SO2 >    // Storage order of the right-hand side matrix
   inline void assign( const Matrix<MT2,SO2>& rhs )
   {
      sm_.assign( trans( *rhs ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   MT& sm_;  //!< The sparse matrix operand.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE  ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the sparse matrix contained in a SMatTransposer.
// \ingroup sparse_matrix_expression
//
// \param m The sparse matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline void reset( SMatTransposer<MT,SO>& m )
{
   m.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given SMatTransposer are intact.
// \ingroup sparse_matrix_expression
//
// \param m The sparse matrix to be tested.
// \return \a true in caes the given matrix's invariants are intact, \a false otherwise.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline bool isIntact( const SMatTransposer<MT,SO>& m ) noexcept
{
   return m.isIntact();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct Size< SMatTransposer<MT,SO>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO >
struct Size< SMatTransposer<MT,SO>, 1UL >
   : public Size<MT,1UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAXSIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct MaxSize< SMatTransposer<MT,SO>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO >
struct MaxSize< SMatTransposer<MT,SO>, 1UL >
   : public MaxSize<MT,1UL>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
