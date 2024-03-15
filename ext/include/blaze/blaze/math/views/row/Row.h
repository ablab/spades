//=================================================================================================
/*!
//  \file blaze/math/views/row/Row.h
//  \brief Row documentation
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

#ifndef _BLAZE_MATH_VIEWS_ROW_ROW_H_
#define _BLAZE_MATH_VIEWS_ROW_ROW_H_


//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup row Row
// \ingroup views
//
// Rows provide views on a specific row of a dense or sparse matrix. As such, rows act as a
// reference to a specific row. This reference is valid and can be used in every way any other
// row vector can be used as long as the matrix containing the row is not resized or entirely
// destroyed. The row also acts as an alias to the row elements: Changes made to the elements
// (e.g. modifying values, inserting or erasing elements) are immediately visible in the matrix
// and changes made via the matrix are immediately visible in the row.
//
//
// \n \section row_setup Setup of Rows
//
// \image html row.png
// \image latex row.eps "Row view" width=250pt
//
// A reference to a dense or sparse row can be created very conveniently via the \c row() function.
// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Row.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The row index must be in the range from \f$[0..M-1]\f$, where \c M is the total number of rows
// of the matrix, and can be specified both at compile time or at runtime:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a reference to the 1st row of matrix A (compile time index)
   auto row1 = row<1UL>( A );

   // Creating a reference to the 2nd row of matrix A (runtime index)
   auto row2 = row( A, 2UL );
   \endcode

// The \c row() function returns an expression representing the row view. The type of this
// expression depends on the given row arguments, primarily the type of the matrix and the compile
// time arguments. If the type is required, it can be determined via \c decltype specifier:

   \code
   using MatrixType = blaze::DynamicMatrix<int>;
   using RowType = decltype( blaze::row<1UL>( std::declval<MatrixType>() ) );
   \endcode

// The resulting view can be treated as any other row vector, i.e. it can be assigned to, it can
// be copied from, and it can be used in arithmetic operations. The reference can also be used on
// both sides of an assignment: The row can either be used as an alias to grant write access to a
// specific row of a matrix primitive on the left-hand side of an assignment or to grant read-access
// to a specific row of a matrix primitive or expression on the right-hand side of an assignment.
// The following example demonstrates this in detail:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   blaze::CompressedVector<double,blaze::rowVector> y;
   blaze::DynamicMatrix<double,blaze::rowMajor> A, B;
   blaze::CompressedMatrix<double,blaze::rowMajor> C, D;
   // ... Resizing and initialization

   // Setting the 2nd row of matrix A to x
   auto row2 = row( A, 2UL );
   row2 = x;

   // Setting the 3rd row of matrix B to y
   row( B, 3UL ) = y;

   // Setting x to the 4th row of the result of the matrix multiplication
   x = row( A * B, 4UL );

   // Setting y to the 2nd row of the result of the sparse matrix multiplication
   y = row( C * D, 2UL );
   \endcode

// \n \section row_element_access Element access
//
// The elements of a row can be directly accessed with the subscript operator:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a view on the 4th row of matrix A
   auto row4 = row( A, 4UL );

   // Setting the 1st element of the dense row, which corresponds
   // to the 1st element in the 4th row of matrix A
   row4[1] = 2.0;
   \endcode

// The numbering of the row elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of columns of the referenced matrix. Alternatively, the elements of a
// row can be traversed via iterators. Just as with vectors, in case of non-const rows, \c begin()
// and \c end() return an iterator, which allows to manipulate the elements, in case of constant
// rows an iterator to immutable elements is returned:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of matrix A
   auto row31 = row( A, 31UL );

   // Traversing the elements via iterators to non-const elements
   for( auto it=row31.begin(); it!=row31.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense row value
      ... = *it;  // OK: Read access to the dense row value.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=row31.cbegin(); it!=row31.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense row value.
   }
   \endcode

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of matrix A
   auto row31 = row( A, 31UL );

   // Traversing the elements via iterators to non-const elements
   for( auto it=row31.begin(); it!=row31.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=row31.cbegin(); it!=row31.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_row_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse row can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 10UL, 100UL );  // Non-initialized 10x100 matrix

   auto row0( row( A, 0UL ) );  // Reference to the 0th row of A

   // The subscript operator provides access to all possible elements of the sparse row,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse row, the element is inserted into the row.
   row0[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the row it is inserted into the row, if it is already contained in
   // the row its value is modified.
   row0.set( 45UL, -1.2 );

   // An alternative for inserting elements into the row is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the row.
   row0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse row is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the row and that the row's capacity is large
   // enough to hold the new element.
   row0.reserve( 10UL );
   row0.append( 51UL, -2.1 );
   \endcode

// \n \section row_common_operations Common Operations
//
// A row view can be used like any other row vector. For instance, the current number of row
// elements can be obtained via the \c size() function, the current capacity via the \c capacity()
// function, and the number of non-zero elements via the \c nonZeros() function. However, since
// rows are references to specific rows of a matrix, several operations are not possible, such as
// resizing and swapping. The following example shows this by means of a dense row view:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd row of matrix A
   auto row2 = row( A, 2UL );

   row2.size();          // Returns the number of elements in the row
   row2.capacity();      // Returns the capacity of the row
   row2.nonZeros();      // Returns the number of non-zero elements contained in the row

   row2.resize( 84UL );  // Compilation error: Cannot resize a single row of a matrix

   auto row3 = row( A, 3UL );
   swap( row2, row3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section row_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse rows can be used in all arithmetic operations that any other dense or
// sparse row vector can be used in. The following example gives an impression of the use of
// dense rows within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse rows with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::rowVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::rowVector> c( 2UL );
   c[1] = 3.0;

   blaze::DynamicMatrix<double,blaze::rowMajor> A( 4UL, 2UL );  // Non-initialized 4x2 matrix

   auto row0( row( A, 0UL ) );  // Reference to the 0th row of A

   row0[0] = 0.0;        // Manual initialization of the 0th row of A
   row0[1] = 0.0;
   row( A, 1UL ) = 1.0;  // Homogeneous initialization of the 1st row of A
   row( A, 2UL ) = a;    // Dense vector initialization of the 2nd row of A
   row( A, 3UL ) = c;    // Sparse vector initialization of the 3rd row of A

   b = row0 + a;              // Dense vector/dense vector addition
   b = c + row( A, 1UL );     // Sparse vector/dense vector addition
   b = row0 * row( A, 2UL );  // Component-wise vector multiplication

   row( A, 1UL ) *= 2.0;     // In-place scaling of the 1st row
   b = row( A, 1UL ) * 2.0;  // Scaling of the 1st row
   b = 2.0 * row( A, 1UL );  // Scaling of the 1st row

   row( A, 2UL ) += a;              // Addition assignment
   row( A, 2UL ) -= c;              // Subtraction assignment
   row( A, 2UL ) *= row( A, 0UL );  // Multiplication assignment

   double scalar = row( A, 1UL ) * trans( c );  // Scalar/dot/inner product between two vectors

   A = trans( c ) * row( A, 1UL );  // Outer product between two vectors
   \endcode

// \n \section row_on_column_major_matrix Rows on Column-Major Matrices
//
// Especially noteworthy is that row views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// possible to traverse a row of a column-major matrix or a column of a row-major matrix. For
// instance:

   \code
   blaze::DynamicMatrix<int,blaze::columnMajor> A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 1st row of a column-major matrix A
   auto row1 = row( A, 1UL );

   for( auto it=row1.begin(); it!=row1.end(); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a row view on a matrix stored in a column-major fashion
// can result in a considerable performance decrease in comparison to a row view on a matrix
// with row-major storage format. This is due to the non-contiguous storage of the matrix
// elements. Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two column-major matrices
   blaze::DynamicMatrix<double,blaze::columnMajor> A( 128UL, 128UL );
   blaze::DynamicMatrix<double,blaze::columnMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th row of the multiplication between A and B ...
   blaze::DynamicVector<double,blaze::rowVector> x = row( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // the 15th row of the column-major matrix A with B.
   blaze::DynamicVector<double,blaze::rowVector> x = row( A, 15UL ) * B;
   \endcode

// Although Blaze performs the resulting vector/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix \c A would result in a more efficient evaluation.
*/
//*************************************************************************************************

#endif
