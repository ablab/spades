//=================================================================================================
/*!
//  \file blaze/math/adaptors/SymmetricMatrix.h
//  \brief Header file for the implementation of a symmetric matrix adaptor
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/adaptors/symmetricmatrix/BaseTemplate.h>
#include <blaze/math/adaptors/symmetricmatrix/DenseNonScalar.h>
#include <blaze/math/adaptors/symmetricmatrix/DenseScalar.h>
#include <blaze/math/adaptors/symmetricmatrix/SparseNonScalar.h>
#include <blaze/math/adaptors/symmetricmatrix/SparseScalar.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/Exception.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/DeclDiagTrait.h>
#include <blaze/math/traits/DeclHermTrait.h>
#include <blaze/math/traits/DeclLowTrait.h>
#include <blaze/math/traits/DeclStrLowTrait.h>
#include <blaze/math/traits/DeclStrUppTrait.h>
#include <blaze/math/traits/DeclSymTrait.h>
#include <blaze/math/traits/DeclUniLowTrait.h>
#include <blaze/math/traits/DeclUniUppTrait.h>
#include <blaze/math/traits/DeclUppTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/YieldsDiagonal.h>
#include <blaze/math/typetraits/YieldsHermitian.h>
#include <blaze/math/typetraits/YieldsIdentity.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/typetraits/IsBuiltin.h>


namespace blaze {

//=================================================================================================
//
//  SYMMETRICMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SymmetricMatrix operators */
//@{
template< typename MT, bool SO, bool DF, bool SF >
void reset( SymmetricMatrix<MT,SO,DF,SF>& m );

template< typename MT, bool SO, bool DF, bool SF >
void reset( SymmetricMatrix<MT,SO,DF,SF>& m, size_t i );

template< typename MT, bool SO, bool DF, bool SF >
void clear( SymmetricMatrix<MT,SO,DF,SF>& m );

template< RelaxationFlag RF, typename MT, bool SO, bool DF, bool SF >
bool isDefault( const SymmetricMatrix<MT,SO,DF,SF>& m );

template< typename MT, bool SO, bool DF, bool SF >
bool isIntact( const SymmetricMatrix<MT,SO,DF,SF>& m );

template< typename MT, bool SO, bool DF, bool SF >
void swap( SymmetricMatrix<MT,SO,DF,SF>& a, SymmetricMatrix<MT,SO,DF,SF>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given symmetric matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool SF >    // Scalar flag
inline void reset( SymmetricMatrix<MT,SO,DF,SF>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the specified row/column of the given symmetric matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given symmetric matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool SF >    // Scalar flag
inline void reset( SymmetricMatrix<MT,SO,DF,SF>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given symmetric matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be cleared.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool SF >    // Scalar flag
inline void clear( SymmetricMatrix<MT,SO,DF,SF>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given symmetric matrix is in default state.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
//
// This function checks whether the matrix is in default state. For instance, in case the
// matrix is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all matrix elements are 0 and \a false in case any matrix element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::SymmetricMatrix<int> A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( A ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the adapted matrix
        , bool SO            // Storage order of the adapted matrix
        , bool DF            // Density flag
        , bool SF >          // Scalar flag
inline bool isDefault( const SymmetricMatrix<MT,SO,DF,SF>& m )
{
   return isDefault<RF>( m.matrix_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given symmetric matrix are intact.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the symmetric matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   SymmetricMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool SF >    // Scalar flag
inline bool isIntact( const SymmetricMatrix<MT,SO,DF,SF>& m )
{
   return m.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup symmetric_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool SF >    // Scalar flag
inline void swap( SymmetricMatrix<MT,SO,DF,SF>& a, SymmetricMatrix<MT,SO,DF,SF>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given symmetric dense matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense matrix by means of the specified matrix
// inversion algorithm \c IF. The The inversion fails if the given matrix is singular and not
// invertible. In this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a m may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert( SymmetricMatrix<MT,SO,true,true>& m )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   if( IF == asUniLower || IF == asUniUpper ) {
      BLAZE_INTERNAL_ASSERT( isIdentity( m ), "Violation of preconditions detected" );
      return;
   }

   constexpr InversionFlag flag( ( IF == byLU || IF == byLDLT || IF == byLDLH ||
                                   IF == asGeneral || IF == asSymmetric || IF == asHermitian )
                                 ? ( byLDLT )
                                 : ( ( IF == byLLH )
                                     ?( byLLH )
                                     :( asDiagonal ) ) );

   MT tmp( m.matrix_ );
   invert<flag>( tmp );
   m.matrix_ = std::move( tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( m ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a matrix to a symmetric matrix.
// \ingroup symmetric_matrix
//
// \param lhs The target left-hand side symmetric matrix.
// \param rhs The right-hand side matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , bool SF       // Scalar flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryAssign( const SymmetricMatrix<MT1,SO1,DF,SF>& lhs,
                       const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (*rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (*rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const size_t M( (*rhs).rows()    );
   const size_t N( (*rhs).columns() );

   if( ( row + M <= column ) || ( column + N <= row ) )
      return true;

   const bool   lower( row > column );
   const size_t size ( min( row + M, column + N ) - ( lower ? row : column ) );

   if( size < 2UL )
      return true;

   const size_t subrow( lower ? 0UL : column - row );
   const size_t subcol( lower ? row - column : 0UL );

   return isSymmetric( submatrix( *rhs, subrow, subcol, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a symmetric matrix.
// \ingroup symmetric_matrix
//
// \param lhs The target left-hand side symmetric matrix.
// \param rhs The right-hand side matrix to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , bool SF       // Scalar flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryAddAssign( const SymmetricMatrix<MT1,SO1,DF,SF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a symmetric
//        matrix.
// \ingroup symmetric_matrix
//
// \param lhs The target left-hand side symmetric matrix.
// \param rhs The right-hand side matrix to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , bool SF       // Scalar flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool trySubAssign( const SymmetricMatrix<MT1,SO1,DF,SF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the Schur product assignment of a matrix to a symmetric
//        matrix.
// \ingroup symmetric_matrix
//
// \param lhs The target left-hand side symmetric matrix.
// \param rhs The right-hand side matrix for the Schur product.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , bool SF       // Scalar flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool trySchurAssign( const SymmetricMatrix<MT1,SO1,DF,SF>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
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
template< typename MT, bool SO, bool DF, bool SF >
struct Size< SymmetricMatrix<MT,SO,DF,SF>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO, bool DF, bool SF >
struct Size< SymmetricMatrix<MT,SO,DF,SF>, 1UL >
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
template< typename MT, bool SO, bool DF, bool SF >
struct MaxSize< SymmetricMatrix<MT,SO,DF,SF>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO, bool DF, bool SF >
struct MaxSize< SymmetricMatrix<MT,SO,DF,SF>, 1UL >
   : public MaxSize<MT,1UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSQUARE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsSquare< SymmetricMatrix<MT,SO,DF,SF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsUniform< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsUniform<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsSymmetric< SymmetricMatrix<MT,SO,DF,SF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsHermitian< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsBuiltin< ElementType_t<MT> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsStrictlyLower< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsZero<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsStrictlyUpper< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsZero<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsAdaptor< SymmetricMatrix<MT,SO,DF,SF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsRestricted< SymmetricMatrix<MT,SO,DF,SF> >
   : public TrueType
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
template< typename MT, bool SO, bool SF >
struct HasConstDataAccess< SymmetricMatrix<MT,SO,true,SF> >
   : public TrueType
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
template< typename MT, bool SO, bool DF, bool SF >
struct IsAligned< SymmetricMatrix<MT,SO,DF,SF> >
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
template< typename MT, bool SO, bool DF, bool SF >
struct IsContiguous< SymmetricMatrix<MT,SO,DF,SF> >
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
template< typename MT, bool SO, bool DF, bool SF >
struct IsPadded< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsPadded<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsResizable< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsResizable<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSHRINKABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct IsShrinkable< SymmetricMatrix<MT,SO,DF,SF> >
   : public IsShrinkable<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REMOVEADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct RemoveAdaptor< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = MT;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct AddTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( ( IsSymmetric_v<T1> && IsSymmetric_v<T2> ) ||
                                    ( IsSymmetric_v<T1> && IsDiagonal_v<T2> ) ||
                                    ( IsDiagonal_v<T1> && IsSymmetric_v<T2> ) ) &&
                                  !( IsDiagonal_v<T1> && IsDiagonal_v<T2> ) &&
                                  !( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = SymmetricMatrix< typename AddTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SubTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( ( IsSymmetric_v<T1> && IsSymmetric_v<T2> ) ||
                                    ( IsSymmetric_v<T1> && IsDiagonal_v<T2> ) ||
                                    ( IsDiagonal_v<T1> && IsSymmetric_v<T2> ) ) &&
                                  !( IsDiagonal_v<T1> && IsDiagonal_v<T2> ) &&
                                  !( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = SymmetricMatrix< typename SubTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHURTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SchurTraitEval1< T1, T2
                      , EnableIf_t< IsMatrix_v<T1> &&
                                    IsMatrix_v<T2> &&
                                    ( IsSymmetric_v<T1> && IsSymmetric_v<T2> ) &&
                                    !( IsDiagonal_v<T1> || IsDiagonal_v<T2> ) &&
                                    !( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                    !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = SymmetricMatrix< typename SchurTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsScalar_v<T2> &&
                                   ( IsSymmetric_v<T1> && !IsDiagonal_v<T1> && !IsUniform_v<T1> ) > >
{
   using Type = SymmetricMatrix< typename MultTraitEval2<T1,T2>::Type >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsSymmetric_v<T2> && !IsDiagonal_v<T2> && !IsUniform_v<T2> ) > >
{
   using Type = SymmetricMatrix< typename MultTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  KRONTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct KronTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsSymmetric_v<T1> && IsSymmetric_v<T2> ) &&
                                   !( IsDiagonal_v<T1> && IsDiagonal_v<T2> ) &&
                                   !( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = SymmetricMatrix< typename KronTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct DivTraitEval1< T1, T2
                    , EnableIf_t< IsSymmetric_v<T1> && !IsDiagonal_v<T1> && IsScalar_v<T2> > >
{
   using Type = SymmetricMatrix< typename DivTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP >
struct UnaryMapTraitEval1< T, OP
                         , EnableIf_t< YieldsSymmetric_v<OP,T> &&
                                       !YieldsHermitian_v<OP,T> &&
                                       !YieldsDiagonal_v<OP,T> &&
                                       !YieldsIdentity_v<OP,T> > >
{
   using Type = SymmetricMatrix< typename UnaryMapTraitEval2<T,OP>::Type, StorageOrder_v<T> >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< YieldsSymmetric_v<OP,T1,T2> &&
                                        !YieldsHermitian_v<OP,T1,T2> &&
                                        !YieldsDiagonal_v<OP,T1,T2> &&
                                        !YieldsIdentity_v<OP,T1,T2> > >
{
   using Type = SymmetricMatrix< typename BinaryMapTraitEval2<T1,T2,OP>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REPEATTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, size_t R0, size_t R1 >
struct RepeatTraitEval1< T, R0, R1, inf
                       , EnableIf_t< R0 != inf && R1 != inf && R0 == R1 &&
                                     IsSymmetric_v<T> > >
{
   using Type = SymmetricMatrix< typename RepeatTraitEval2<T,R0,R1,inf>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSYMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclSymTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = SymmetricMatrix<MT,SO,DF,SF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLHERMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclHermTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = HermitianMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclLowTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = DiagonalMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUNILOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclUniLowTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = IdentityMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSTRLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclStrLowTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = ZeroMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclUppTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = DiagonalMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUNIUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclUniUppTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = IdentityMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSTRUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclStrUppTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = ZeroMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLDIAGTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF >
struct DeclDiagTrait< SymmetricMatrix<MT,SO,DF,SF> >
{
   using Type = DiagonalMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, bool SF1, typename MT2, bool SO2, bool DF2, bool SF2 >
struct HighType< SymmetricMatrix<MT1,SO1,DF1,SF1>, SymmetricMatrix<MT2,SO2,DF2,SF2> >
{
   using Type = SymmetricMatrix< typename HighType<MT1,MT2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, bool SF1, typename MT2, bool SO2, bool DF2, bool SF2 >
struct LowType< SymmetricMatrix<MT1,SO1,DF1,SF1>, SymmetricMatrix<MT2,SO2,DF2,SF2> >
{
   using Type = SymmetricMatrix< typename LowType<MT1,MT2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I, size_t N >
struct SubmatrixTraitEval1< MT, I, I, N, N
                          , EnableIf_t< I != inf && N != inf &&
                                        IsSymmetric_v<MT> &&
                                        !IsDiagonal_v<MT> &&
                                        !IsUniform_v<MT> &&
                                        !IsZero_v<MT> > >
{
   using Type = SymmetricMatrix< typename SubmatrixTraitEval2<MT,I,I,N,N>::Type >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
