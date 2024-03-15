//=================================================================================================
/*!
//  \file blaze/math/views/rows/Rows.h
//  \brief Rows documentation
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

#ifndef _BLAZE_MATH_VIEWS_ROWS_ROWS_H_
#define _BLAZE_MATH_VIEWS_ROWS_ROWS_H_


//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup rows Rows
// \ingroup views
//
// Row selections provide views on arbitrary compositions of rows of dense and sparse matrices.
// These views act as a reference to the selected rows and represent them as another dense or
// sparse matrix. This reference is valid and can be used in every way any other dense or sparse
// matrix can be used as long as the matrix containing the rows is not resized or entirely
// destroyed. The row selection also acts as an alias to the matrix elements in the specified
// range: Changes made to the rows (e.g. modifying values, inserting or erasing elements) are
// immediately visible in the matrix and changes made via the matrix are immediately visible
// in the rows.
//
//
// \n \section rows_setup Setup of Row Selections
//
// A row selection can be created very conveniently via the \c rows() function. It can be included
// via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Rows.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The indices of the rows to be selected can be specified either at compile time or at runtime
// (by means of an initializer list, array or vector):

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Selecting the rows 4, 6, 8, and 10 (compile time arguments)
   auto rs1 = rows<4UL,6UL,8UL,10UL>( A );

   // Selecting the rows 3, 2, and 1 (runtime arguments via an initializer list)
   const std::initializer_list<size_t> list{ 3UL, 2UL, 1UL };
   auto rs2 = rows( A, { 3UL, 2UL, 1UL } );
   auto rs3 = rows( A, list );

   // Selecting the rows 1, 2, 3, 3, 2, and 1 (runtime arguments via a std::array)
   const std::array<size_t> array{ 1UL, 2UL, 3UL, 3UL, 2UL, 1UL };
   auto rs4 = rows( A, array );
   auto rs5 = rows( A, array.data(), array.size() );

   // Selecting the row 4 fives times (runtime arguments via a std::vector)
   const std::vector<size_t> vector{ 4UL, 4UL, 4UL, 4UL, 4UL };
   auto rs6 = rows( A, vector );
   auto rs7 = rows( A, vector.data(), vector.size() );
   \endcode

// Note that it is possible to alias the rows of the underlying matrix in any order. Also note
// that it is possible to use the same index multiple times.
//
// Alternatively it is possible to pass a callable such as a lambda or functor that produces the
// indices:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 18UL );

   // Selecting all even rows of the matrix, i.e. selecting the rows 0, 2, 4, 6, and 8
   auto rs1 = rows( A, []( size_t i ){ return i*2UL; }, 5UL );

   // Selecting all odd rows of the matrix, i.e. selecting the rows 1, 3, 5, and 7
   auto rs2 = rows( x, []( size_t i ){ return i*2UL+1UL; }, 4UL );

   // Reversing the rows of the matrix, i.e. selecting the rows 8, 7, 6, 5, 4, 3, 2, 1, and 0
   auto rs3 = rows( v, [max=A.rows()-1UL]( size_t i ){ return max-i; }, 9UL );
   \endcode

// The \c rows() function returns an expression representing the view on the selected rows. The
// type of this expression depends on the given arguments, primarily the type of the matrix and
// the compile time arguments. If the type is required, it can be determined via the \c decltype
// specifier:

   \code
   using MatrixType = blaze::DynamicMatrix<int>;
   using RowsType = decltype( blaze::rows<3UL,0UL,4UL,8UL>( std::declval<MatrixType>() ) );
   \endcode

// The resulting view can be treated as any other dense or sparse matrix, i.e. it can be assigned
// to, it can be copied from, and it can be used in arithmetic operations. Note, however, that a
// row selection will always be treated as a row-major matrix, regardless of the storage order of
// the matrix containing the rows. The view can also be used on both sides of an assignment: It
// can either be used as an alias to grant write access to specific rows of a matrix primitive
// on the left-hand side of an assignment or to grant read-access to specific rows of a matrix
// primitive or expression on the right-hand side of an assignment. The following example
// demonstrates this in detail:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   blaze::DynamicMatrix<double,blaze::columnMajor> B;
   blaze::CompressedMatrix<double,blaze::rowMajor> C;
   // ... Resizing and initialization

   // Selecting the rows 1, 3, 5, and 7 of A
   auto rs = rows( A, { 1UL, 3UL, 5UL, 7UL } );

   // Setting rows 1, 3, 5, and 7 of A to row 4 of B
   rs = rows( B, { 4UL, 4UL, 4UL, 4UL } );

   // Setting the rows 2, 4, 6, and 8 of A to C
   rows( A, { 2UL, 4UL, 6UL, 8UL } ) = C;

   // Setting the first 4 rows of A to the rows 5, 4, 3, and 2 of C
   submatrix( A, 0UL, 0UL, 4UL, A.columns() ) = rows( C, { 5UL, 4UL, 3UL, 2UL } );

   // Rotating the result of the addition between rows 1, 3, 5, and 7 of A and C
   B = rows( rs + C, { 2UL, 3UL, 0UL, 1UL } );
   \endcode

// \n \section rows_element_access Element access
//
// The elements of a row selection can be directly accessed via the function call operator:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a view on the first four rows of A in reverse order
   auto rs = rows( A, { 3UL, 2UL, 1UL, 0UL } );

   // Setting the element (0,0) of the row selection, which corresponds
   // to the element at position (3,0) in matrix A
   rs(0,0) = 2.0;
   \endcode

// Alternatively, the elements of a row selection can be traversed via (const) iterators. Just as
// with matrices, in case of non-const row selection, \c begin() and \c end() return an iterator,
// which allows to manipuate the elements, in case of constant row selection an iterator to
// immutable elements is returned:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a selection of rows of matrix A
   auto rs = rows( A, { 16UL, 32UL, 64UL, 128UL } );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( auto it=rs.begin(0); it!=rs.end(0); ++it ) {
      *it = ...;  // OK: Write access to the dense value.
      ... = *it;  // OK: Read access to the dense value.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( auto it=rs.cbegin(1); it!=rs.cend(1); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense value.
   }
   \endcode

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a selection of rows of matrix A
   auto rs = rows( A, { 16UL, 32UL, 64UL, 128UL } );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( auto it=rs.begin(0); it!=rs.end(0); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( auto it=rs.cbegin(1); it!=rs.cend(1); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section rows_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse row selection can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 256UL, 512UL );  // Non-initialized matrix of size 256x512

   auto rs = rows( A, { 10UL, 20UL, 30UL, 40UL } );  // View on the rows 10, 20, 30, and 40 of A

   // The function call operator provides access to all possible elements of the sparse row
   // selection, including the zero elements. In case the function call operator is used to
   // access an element that is currently not stored in the sparse row selection, the element
   // is inserted into the row selection.
   rs(2,4) = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the row selection it is inserted into the row selection, if it is already
   // contained in the row selection its value is modified.
   rs.set( 2UL, 5UL, -1.2 );

   // An alternative for inserting elements into the row selection is the \c insert() function.
   // However, it inserts the element only in case the element is not already contained in the
   // row selection.
   rs.insert( 2UL, 6UL, 3.7 );

   // Just as in the case of sparse matrices, elements can also be inserted via the \c append()
   // function. In case of row selections, \c append() also requires that the appended element's
   // index is strictly larger than the currently largest non-zero index in the according row
   // of the row selection and that the according row's capacity is large enough to hold the new
   // element. Note however that due to the nature of a row selection, which may be an alias to
   // an arbitrary collection of rows, the \c append() function does not work as efficiently for
   // a row selection as it does for a matrix.
   rs.reserve( 2UL, 10UL );
   rs.append( 2UL, 10UL, -2.1 );
   \endcode

// \n \section rows_common_operations Common Operations
//
// A view on specific rows of a matrix can be used like any other dense or sparse matrix. For
// instance, the current size of the matrix, i.e. the number of rows or columns can be obtained
// via the \c rows() and \c columns() functions, the current total capacity via the \c capacity()
// function, and the number of non-zero elements via the \c nonZeros() function. However, since
// row selections are views on specific rows of a matrix, several operations are not possible,
// such as resizing and swapping:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a view on the rows 8, 16, 24, and 32 of matrix A
   auto rs = rows( A, { 8UL, 16UL, 24UL, 32UL } );

   rs.rows();      // Returns the number of rows of the row selection
   rs.columns();   // Returns the number of columns of the row selection
   rs.capacity();  // Returns the capacity of the row selection
   rs.nonZeros();  // Returns the number of non-zero elements contained in the row selection

   rs.resize( 10UL, 8UL );  // Compilation error: Cannot resize a row selection

   auto rs2 = rows( A, 9UL, 17UL, 25UL, 33UL );
   swap( rs, rs2 );  // Compilation error: Swap operation not allowed
   \endcode

// \n \section rows_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse row selections can be used in all arithmetic operations that any other
// dense or sparse matrix can be used in. The following example gives an impression of the use
// of dense row selctions within arithmetic operations. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and
// sparse matrices with fitting element types:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D1, D2, D3;
   blaze::CompressedMatrix<double,blaze::rowMajor> S1, S2;

   blaze::CompressedVector<double,blaze::columnVector> a, b;

   // ... Resizing and initialization

   std::initializer_list<size_t> indices1{ 0UL, 3UL, 6UL,  9UL, 12UL, 15UL, 18UL, 21UL };
   std::initializer_list<size_t> indices2{ 1UL, 4UL, 7UL, 10UL, 13UL, 16UL, 19UL, 22UL };
   std::initializer_list<size_t> indices3{ 2UL, 5UL, 8UL, 11UL, 14UL, 17UL, 20UL, 23UL };

   auto rs = rows( D1, indices1 );  // Selecting the every third row of D1 in the range [0..21]

   rs = D2;                    // Dense matrix assignment to the selected rows
   rows( D1, indices2 ) = S1;  // Sparse matrix assignment to the selected rows

   D3 = rs + D2;                    // Dense matrix/dense matrix addition
   S2 = S1 - rows( D1, indices2 );  // Sparse matrix/dense matrix subtraction
   D2 = rs % rows( D1, indices3 );  // Dense matrix/dense matrix Schur product
   D2 = rows( D1, indices2 ) * D1;  // Dense matrix/dense matrix multiplication

   rows( D1, indices2 ) *= 2.0;      // In-place scaling of the second selection of rows
   D2 = rows( D1, indices3 ) * 2.0;  // Scaling of the elements in the third selection of rows
   D2 = 2.0 * rows( D1, indices3 );  // Scaling of the elements in the third selection of rows

   rows( D1, indices1 ) += D2;  // Addition assignment
   rows( D1, indices2 ) -= S1;  // Subtraction assignment
   rows( D1, indices3 ) %= rs;  // Schur product assignment

   a = rows( D1, indices1 ) * b;  // Dense matrix/sparse vector multiplication
   \endcode

// \n \section rows_on_column_major_matrix Row Selections on Column-Major Matrices
//
// Especially noteworthy is that row selections can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// possible to traverse a row of a column-major matrix or a column of a row-major matrix. For
// instance:

   \code
   blaze::DynamicMatrix<int,blaze::columnMajor> A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 1st and 3rd row of a column-major matrix A
   auto rs = rows( A, { 1UL, 3UL } );

   // Traversing row 0 of the selection, which corresponds to the 1st row of matrix A
   for( auto it=rs.begin( 0UL ); it!=rs.end( 0UL ); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a row selection on a matrix stored in a column-major fashion
// can result in a considerable performance decrease in comparison to a row selection on a matrix
// with row-major storage format. This is due to the non-contiguous storage of the matrix elements.
// Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two column-major matrices
   blaze::DynamicMatrix<double,blaze::columnMajor> A( 128UL, 128UL );
   blaze::DynamicMatrix<double,blaze::columnMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th, 30th, and 45th row of the multiplication between A and B ...
   blaze::DynamicMatrix<double,blaze::rowMajor> x = rows( A * B, { 15UL, 30UL, 45UL } );

   // ... is essentially the same as the following computation, which multiplies
   // the 15th, 30th, and 45th row of the column-major matrix A with B.
   blaze::DynamicMatrix<double,blaze::rowMajor> x = rows( A, { 15UL, 30UL, 45UL } ) * B;
   \endcode

// Although Blaze performs the resulting matrix/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix \c A would result in a more efficient evaluation.
*/
//*************************************************************************************************

} // namespace blaze

#endif
