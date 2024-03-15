//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatTransposer.h
//  \brief Header file for the dense matrix transposer
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATTRANSPOSER_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATTRANSPOSER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATTRANSPOSER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the transposition of a dense matrix.
// \ingroup dense_matrix_expression
//
// The DMatTransposer class is a wrapper object for the temporary transposition of a dense matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
class DMatTransposer
   : public DenseMatrix< DMatTransposer<MT,SO>, SO >
{
 public:
   //**Type definitions****************************************************************************
   using This           = DMatTransposer<MT,SO>;     //!< Type of this DMatTransposer instance.
   using BaseType       = DenseMatrix<This,SO>;      //!< Base type of this DMatTransposer instance.
   using ResultType     = TransposeType_t<MT>;       //!< Result type for expression template evaluations.
   using OppositeType   = OppositeType_t<MT>;        //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType  = ResultType_t<MT>;          //!< Transpose type for expression template evaluations.
   using ElementType    = ElementType_t<MT>;         //!< Type of the matrix elements.
   using SIMDType       = SIMDTrait_t<ElementType>;  //!< SIMD type of the matrix elements.
   using ReturnType     = ReturnType_t<MT>;          //!< Return type for expression template evaluations.
   using CompositeType  = const This&;               //!< Data type for composite expression templates.
   using Reference      = Reference_t<MT>;           //!< Reference to a non-constant matrix value.
   using ConstReference = ConstReference_t<MT>;      //!< Reference to a constant matrix value.
   using Pointer        = Pointer_t<MT>;             //!< Pointer to a non-constant matrix value.
   using ConstPointer   = ConstPointer_t<MT>;        //!< Pointer to a constant matrix value.
   using Iterator       = Iterator_t<MT>;            //!< Iterator over non-constant elements.
   using ConstIterator  = ConstIterator_t<MT>;       //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the matrix is involved
       in can be optimized via SIMD operations. In case the dense matrix operand is vectorizable,
       the \a simdEnabled compilation flag is set to \a true, otherwise it is set to \a false. */
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatTransposer class.
   //
   // \param dm The dense matrix operand.
   */
   explicit inline DMatTransposer( MT& dm ) noexcept
      : dm_( dm )  // The dense matrix operand
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline Reference operator()( size_t i, size_t j ) {
      BLAZE_INTERNAL_ASSERT( i < dm_.columns(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.rows()   , "Invalid column access index" );
      return dm_(j,i);
   }
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline ConstReference operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < dm_.columns(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.rows()   , "Invalid column access index" );
      return const_cast<const MT&>( dm_ )(j,i);
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
   inline Reference at( size_t i, size_t j ) {
      if( i >= dm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= dm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
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
      if( i >= dm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= dm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the matrix elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline Pointer data() noexcept {
      return dm_.data();
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the matrix elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline ConstPointer data() const noexcept {
      return dm_.data();
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
      return dm_.begin( i );
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
      return dm_.cbegin( i );
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
      return dm_.cbegin( i );
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
      return dm_.end( i );
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
      return dm_.cend( i );
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
      return dm_.cend( i );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return dm_.columns();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return dm_.rows();
   }
   //**********************************************************************************************

   //**Spacing function****************************************************************************
   /*!\brief Returns the spacing between the beginning of two rows.
   //
   // \return The spacing between the beginning of two rows.
   */
   inline size_t spacing() const noexcept {
      return dm_.spacing();
   }
   //**********************************************************************************************

   //**Reset function******************************************************************************
   /*!\brief Resets the matrix elements.
   //
   // \return void
   */
   inline void reset() {
      return dm_.reset();
   }
   //**********************************************************************************************

   //**Resize function*****************************************************************************
   /*!\brief Changing the size of the matrix.
   //
   // \param m The new number of rows of the matrix.
   // \param n The new number of columns of the matrix.
   // \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
   // \return void
   */
   inline void resize( size_t m, size_t n, bool preserve=true ) {
      dm_.resize( m, n, preserve );
   }
   //**********************************************************************************************

   //**IsIntact function***************************************************************************
   /*!\brief Returns whether the invariants of the matrix are intact.
   //
   // \return \a true in case the matrix's invariants are intact, \a false otherwise.
   */
   inline bool isIntact() const noexcept {
      using blaze::isIntact;
      return isIntact( dm_ );
   }
   //**********************************************************************************************

   //**CanAliased function*************************************************************************
   /*!\brief Returns whether the matrix can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this matrix, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool canAlias( const Other* alias ) const noexcept
   {
      return dm_.canAlias( alias );
   }
   //**********************************************************************************************

   //**IsAliased function**************************************************************************
   /*!\brief Returns whether the matrix is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this matrix, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool isAliased( const Other* alias ) const noexcept
   {
      return dm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**IsAligned function**************************************************************************
   /*!\brief Returns whether the matrix is properly aligned in memory.
   //
   // \return \a true in case the matrix is aligned, \a false if not.
   */
   inline bool isAligned() const noexcept
   {
      return dm_.isAligned();
   }
   //**********************************************************************************************

   //**CanSMPAssign function***********************************************************************
   /*!\brief Returns whether the matrix can be used in SMP assignments.
   //
   // \return \a true in case the matrix can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept
   {
      return dm_.canSMPAssign();
   }
   //**********************************************************************************************

   //**Load function*******************************************************************************
   /*!\brief Load of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \return The loaded SIMD element.
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE SIMDType load( size_t i, size_t j ) const noexcept
   {
      return dm_.load( j, i );
   }
   //**********************************************************************************************

   //**Loada function******************************************************************************
   /*!\brief Aligned load of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \return The loaded SIMD element.
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept
   {
      return dm_.loada( j, i );
   }
   //**********************************************************************************************

   //**Loadu function******************************************************************************
   /*!\brief Unaligned load of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \return The loaded SIMD element.
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept
   {
      return dm_.loadu( j, i );
   }
   //**********************************************************************************************

   //**Store function******************************************************************************
   /*!\brief Store of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void store( size_t i, size_t j, const SIMDType& value ) noexcept
   {
      dm_.store( j, i, value );
   }
   //**********************************************************************************************

   //**Storea function******************************************************************************
   /*!\brief Aligned store of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept
   {
      dm_.storea( j, i, value );
   }
   //**********************************************************************************************

   //**Storeu function*****************************************************************************
   /*!\brief Unaligned store of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept
   {
      dm_.storeu( j, i, value );
   }
   //**********************************************************************************************

   //**Stream function*****************************************************************************
   /*!\brief Aligned, non-temporal store of a SIMD element of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range [0..M-1].
   // \param j Access index for the column. The index has to be in the range [0..N-1].
   // \param value The SIMD element to be stored.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors.
   */
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept
   {
      dm_.stream( j, i, value );
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
      dm_.assign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose addition assignment of matrices***************************************************
   /*!\brief Implementation of the transpose addition assignment of a matrix.
   //
   // \param rhs The right-hand side matrix to be added.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename MT2  // Type of the right-hand side matrix
           , bool SO2 >    // Storage order of the right-hand side matrix
   inline void addAssign( const Matrix<MT2,SO2>& rhs )
   {
      dm_.addAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose subtraction assignment of matrices************************************************
   /*!\brief Implementation of the transpose subtraction assignment of a matrix.
   //
   // \param rhs The right-hand side matrix to be subtracted.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename MT2  // Type of the right-hand side matrix
           , bool SO2 >    // Storage order of the right-hand side matrix
   inline void subAssign( const Matrix<MT2,SO2>& rhs )
   {
      dm_.subAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose Schur product assignment of matrices**********************************************
   /*!\brief Implementation of the transpose Schur product assignment of a matrix.
   //
   // \param rhs The right-hand side matrix for the Schur product.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename MT2  // Type of the right-hand side matrix
           , bool SO2 >    // Storage order of the right-hand side matrix
   inline void schurAssign( const Matrix<MT2,SO2>& rhs )
   {
      dm_.schurAssign( trans( *rhs ) );
   }
   //**********************************************************************************************

   //**Transpose multiplication assignment of matrices*********************************************
   // No special implementation for the transpose multiplication assignment of matrices.
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   MT& dm_;  //!< The dense matrix operand.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT );
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
/*!\brief Resetting the dense matrix contained in a DMatTransposer.
// \ingroup dense_matrix_expression
//
// \param m The dense matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline void reset( DMatTransposer<MT,SO>& m )
{
   m.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given DMatTransposer are intact.
// \ingroup dense_matrix_expression
//
// \param m The dense matrix to be tested.
// \return \a true in caes the given matrix's invariants are intact, \a false otherwise.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline bool isIntact( const DMatTransposer<MT,SO>& m ) noexcept
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
struct Size< DMatTransposer<MT,SO>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO >
struct Size< DMatTransposer<MT,SO>, 1UL >
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
struct MaxSize< DMatTransposer<MT,SO>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO >
struct MaxSize< DMatTransposer<MT,SO>, 1UL >
   : public MaxSize<MT,1UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct HasConstDataAccess< DMatTransposer<MT,SO> >
   : public HasConstDataAccess<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct HasMutableDataAccess< DMatTransposer<MT,SO> >
   : public HasMutableDataAccess<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsAligned< DMatTransposer<MT,SO> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGUOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsContiguous< DMatTransposer<MT,SO> >
   : public IsContiguous<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsPadded< DMatTransposer<MT,SO> >
   : public IsPadded<MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
