//=================================================================================================
/*!
//  \file blaze/math/Matrix.h
//  \brief Header file for all basic Matrix functionality
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

#ifndef _BLAZE_MATH_MATRIX_H_
#define _BLAZE_MATH_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <iosfwd>
#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/views/Band.h>
#include <blaze/math/views/Elements.h>
#include <blaze/util/EnableIf.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Matrix functions */
//@{
template< typename MT, bool SO >
bool isSymmetric( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isHermitian( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isUniform( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isLower( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isUniLower( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isStrictlyLower( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isUpper( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isUniUpper( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isStrictlyUpper( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isDiagonal( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
bool isIdentity( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
auto trace( const Matrix<MT,SO>& m );

template< bool RF, typename MT >
decltype(auto) reverse( MT&& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is symmetric.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is symmetric, \a false if not.
//
// This function checks if the given dense or sparse matrix is symmetric. The matrix is considered
// to be symmetric if it is a square matrix whose transpose is equal to itself (\f$ A = A^T \f$).
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isSymmetric( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isSymmetric<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a symmetric matrix:

   \code
   if( isSymmetric( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isSymmetric( const Matrix<MT,SO>& m )
{
   return isSymmetric<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is Hermitian.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is Hermitian, \a false if not.
//
// This function checks if the given dense or sparse matrix is an Hermitian matrix. The matrix
// is considered to be an Hermitian matrix if it is a square matrix whose conjugate transpose is
// equal to itself (\f$ A = \overline{A^T} \f$), i.e. each matrix element \f$ a_{ij} \f$ is equal
// to the complex conjugate of the element \f$ a_{ji} \f$. The following code example demonstrates
// the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isHermitian( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isHermitian<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an Hermitian matrix:

   \code
   if( isHermitian( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isHermitian( const Matrix<MT,SO>& m )
{
   return isHermitian<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is a uniform matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is a uniform matrix. The matrix
// is considered to be uniform if all its elements are identical. The following code example
// demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniform( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniform<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a uniform matrix:

   \code
   if( isUniform( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isUniform( const Matrix<MT,SO>& m )
{
   return isUniform<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is a lower triangular matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is a lower triangular matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is a lower triangular matrix. The
// matrix is considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        l_{1,0} & l_{1,1} & 0       & \cdots & 0       \\
                        l_{2,0} & l_{2,1} & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially lower triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a lower triangular matrix:

   \code
   if( isLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isLower( const Matrix<MT,SO>& m )
{
   return isLower<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is a lower unitriangular matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is a lower unitriangular matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is a lower unitriangular matrix. The
// matrix is considered to be lower unitriangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 1       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 1       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 1      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a lower unitriangular matrix:

   \code
   if( isUniLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isUniLower( const Matrix<MT,SO>& m )
{
   return isUniLower<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is a strictly lower triangular matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is a strictly lower triangular matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is a strictly lower triangular matrix.
// The matrix is considered to be strictly lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 0       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 0       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 0      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isStrictlyLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isStrictlyLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a strictly lower triangular
// matrix:

   \code
   if( isStrictlyLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isStrictlyLower( const Matrix<MT,SO>& m )
{
   return isStrictlyLower<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is an upper triangular matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is an upper triangular matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is an upper triangular matrix. The
// matrix is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        u_{0,0} & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & u_{1,1} & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & u_{2,2} & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & u_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially upper triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an upper triangular matrix:

   \code
   if( isUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isUpper( const Matrix<MT,SO>& m )
{
   return isUpper<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is an upper unitriangular matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is an upper unitriangular matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is an upper unitriangular matrix. The
// matrix is considered to be upper unitriangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 1       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 1       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 1       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an upper unitriangular matrix:

   \code
   if( isUniUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isUniUpper( const Matrix<MT,SO>& m )
{
   return isUniUpper<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given matrix is a strictly upper triangular matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is a strictly upper triangular matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is a strictly upper triangular matrix.
// The matrix is considered to be strictly upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 0       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 0       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 0       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isStrictlyUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isStrictlyUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a strictly upper triangular
// matrix:

   \code
   if( isStrictlyUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isStrictlyUpper( const Matrix<MT,SO>& m )
{
   return isStrictlyUpper<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the give matrix is diagonal.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is diagonal, \a false if not.
//
// This function checks if the given dense or sparse matrix is diagonal, i.e. if the non-diagonal
// elements are default elements. In case of integral or floating point data types, a diagonal
// matrix has the form

                        \f[\left(\begin{array}{*{5}{c}}
                        aa     & 0      & 0      & \cdots & 0  \\
                        0      & bb     & 0      & \cdots & 0  \\
                        0      & 0      & cc     & \cdots & 0  \\
                        \vdots & \vdots & \vdots & \ddots & 0  \\
                        0      & 0      & 0      & 0      & xx \\
                        \end{array}\right)\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially diagonal. The
// following example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isDiagonal( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDiagonal<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a diagonal matrix:

   \code
   if( isDiagonal( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isDiagonal( const Matrix<MT,SO>& m )
{
   return isDiagonal<relaxed>( *m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the give matrix is an identity matrix.
// \ingroup matrix
//
// \param m The matrix to be checked.
// \return \a true if the matrix is an identity matrix, \a false if not.
//
// This function checks if the given dense or sparse matrix is an identity matrix, i.e. if the
// diagonal elements are 1 and the non-diagonal elements are 0. In case of integral or floating
// point data types, an identity matrix has the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & 0      & 0      & \cdots & 0 \\
                        0      & 1      & 0      & \cdots & 0 \\
                        0      & 0      & 1      & \cdots & 0 \\
                        \vdots & \vdots & \vdots & \ddots & 0 \\
                        0      & 0      & 0      & 0      & 1 \\
                        \end{array}\right)\f]

// The following example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isIdentity( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isIdentity<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an identity matrix:

   \code
   if( isIdentity( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline bool isIdentity( const Matrix<MT,SO>& m )
{
   return isIdentity<relaxed>( *m );
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Computes the trace of the given square matrix.
// \ingroup matrix
//
// \param m Reference to a constant matrix object.
// \return The trace of the matrix.
// \exception std::invalid_argument Invalid input matrix for trace computation.
//
// This function computes the trace of the given square matrix, i.e. sums the elements on its
// diagonal:

            \f[ trace(A) = a_{11} + a_{22} + ... + a_{nn} = \sum_{i=1}^{n} a_{ii} \f]

// In case the given matrix is not a square matrix a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline auto trace( const Matrix<MT,SO>& m )
{
   if( !isSquare( *m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix for trace computation" );
   }

   return sum( diagonal( *m ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c reverse() function for reversing the rows of a matrix.
// \ingroup matrix
//
// \param m The matrix to be reversed.
// \return The reversed matrix.
*/
template< bool RF      // Reverse flag
        , typename MT  // Type of the matrix
        , EnableIf_t< RF == rowwise >* = nullptr >
inline decltype(auto) reverse_backend( MT&& m )
{
   return rows( std::forward<MT>( m ), [max=m.rows()-1UL]( size_t i ){ return max - i; }, m.rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c reverse() function for reversing the columns of a matrix.
// \ingroup matrix
//
// \param m The matrix to be reversed.
// \return The reversed matrix.
*/
template< bool RF      // Reverse flag
        , typename MT  // Type of the matrix
        , EnableIf_t< RF == columnwise >* = nullptr >
inline decltype(auto) reverse_backend( MT&& m )
{
   return columns( std::forward<MT>( m ), [max=m.columns()-1UL]( size_t i ){ return max - i; }, m.columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reverse the rows or columns of a matrix.
// \ingroup matrix
//
// \param m The matrix to be reversed.
// \return The reversed matrix.
//
// This function reverses the rows or matrices of a dense or sparse matrix. In case the compile
// time flag \a RF is set to \a blaze::rowwise, the rows of the matrix are reversed, in case \a RF
// is set to \a blaze::columnwise, the columns of the matrix are reversed. The following examples
// gives an impression of both alternatives:

   \code
   blaze::DynamicMatrix<int,rowMajor> A{ { 1, 0, 2, 3 },
                                         { 2, 4, 0, 1 },
                                         { 0, 3, 1, 0 } };
   blaze::DynamicMatrix<int> B;

   // Reversing the rows result in the matrix
   //
   //    ( 0 3 1 0 )
   //    ( 2 4 0 1 )
   //    ( 1 0 2 3 )
   //
   B = reverse<rowwise>( A );

   // Reversing the columns result in the matrix
   //
   //    ( 3 2 0 1 )
   //    ( 1 0 4 2 )
   //    ( 0 1 3 0 )
   //
   B = reverse<columnwise>( A );
   \endcode
*/
template< bool RF        // Reverse flag
        , typename MT >  // Type of the matrix
inline decltype(auto) reverse( MT&& m )
{
   return reverse_backend<RF>( std::forward<MT>( m ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Matrix operators */
//@{
template< typename MT, bool SO >
std::ostream& operator<<( std::ostream& os, const Matrix<MT,SO>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for dense and sparse matrices.
// \ingroup matrix
//
// \param os Reference to the output stream.
// \param m Reference to a constant matrix object.
// \return Reference to the output stream.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline std::ostream& operator<<( std::ostream& os, const Matrix<MT,SO>& m )
{
   CompositeType_t<MT> tmp( *m );

   for( size_t i=0UL; i<tmp.rows(); ++i ) {
      os << "( ";
      for( size_t j=0UL; j<tmp.columns(); ++j ) {
         os << std::setw(12) << tmp(i,j) << " ";
      }
      os << ")\n";
   }

   return os;
}
//*************************************************************************************************

} // namespace blaze

#endif
