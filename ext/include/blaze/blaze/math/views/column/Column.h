//=================================================================================================
/*!
//  \file blaze/math/views/column/Column.h
//  \brief Column documentation
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_COLUMN_H_
#define _BLAZE_MATH_VIEWS_COLUMN_COLUMN_H_


//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup column Column
// \ingroup views
//
// Just as rows provide a view on a specific row of a matrix, columns provide views on a specific
// column of a dense or sparse matrix. As such, columns act as a reference to a specific column.
// This reference is valid and can be used in every way any other column vector can be used as
// long as the matrix containing the column is not resized or entirely destroyed. Changes made to
// the elements (e.g. modifying values, inserting or erasing elements) are immediately visible in
// the matrix and changes made via the matrix are immediately visible in the column.
//
//
// \n \section column_setup Setup of Columns
//
// \image html column.png
// \image latex column.eps "Column view" width=250pt
//
// A reference to a dense or sparse column can be created very conveniently via the \c column()
// function. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Column.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The column index must be in the range from \f$[0..N-1]\f$, where \c N is the total number of
// columns of the matrix, and can be specified both at compile time or at runtime:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   // Creating a reference to the 1st column of matrix A (compile time index)
   auto col1 = column<1UL>( A );

   // Creating a reference to the 2nd column of matrix A (runtime index)
   auto col2 = column( A, 2UL );
   \endcode

// The \c column() function returns an expression representing the column view. The type of this
// expression depends on the given column arguments, primarily the type of the matrix and the
// compile time arguments. If the type is required, it can be determined via \c decltype specifier:

   \code
   using MatrixType = blaze::DynamicMatrix<int>;
   using ColumnType = decltype( blaze::column<1UL>( std::declval<MatrixType>() ) );
   \endcode

// The resulting view can be treated as any other column vector, i.e. it can be assigned to, it
// can be copied from, and it can be used in arithmetic operations. The reference can also be used
// on both sides of an assignment: The column can either be used as an alias to grant write access
// to a specific column of a matrix primitive on the left-hand side of an assignment or to grant
// read-access to a specific column of a matrix primitive or expression on the right-hand side
// of an assignment. The following example demonstrates this in detail:

   \code
   blaze::DynamicVector<double,blaze::columnVector> x;
   blaze::CompressedVector<double,blaze::columnVector> y;
   blaze::DynamicMatrix<double,blaze::columnMajor> A, B;
   blaze::CompressedMatrix<double,blaze::columnMajor> C, D;
   // ... Resizing and initialization

   // Setting the 1st column of matrix A to x
   auto col1 = column( A, 1UL );
   col1 = x;

   // Setting the 4th column of matrix B to y
   column( B, 4UL ) = y;

   // Setting x to the 2nd column of the result of the matrix multiplication
   x = column( A * B, 2UL );

   // Setting y to the 2nd column of the result of the sparse matrix multiplication
   y = column( C * D, 2UL );
   \endcode

// \n \section column_element_access Element access
//
// The elements of the dense column can be directly accessed with the subscript operator.

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   // Creating a view on the 4th column of matrix A
   auto col4 = column( A, 4UL );

   // Setting the 1st element of the dense column, which corresponds
   // to the 1st element in the 4th column of matrix A
   col4[1] = 2.0;
   \endcode

// The numbering of the column elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of rows of the referenced matrix. Alternatively, the elements of a column
// can be traversed via iterators. Just as with vectors, in case of non-const columns, \c begin()
// and \c end() return an iterator, which allows to manipulate the elements, in case of constant
// columns an iterator to immutable elements is returned:

   \code
   blaze::DynamicMatrix<int,blaze::columnMajor> A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st column of matrix A
   auto col31 = column( A, 31UL );

   // Traversing the elements via iterators to non-const elements
   for( auto it=col31.begin(); it!=col31.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense column value
      ... = *it;  // OK: Read access to the dense column value.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=col31.cbegin(); it!=col31.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense column value.
   }
   \endcode

   \code
   blaze::CompressedMatrix<int,blaze::columnMajor> A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st column of matrix A
   auto col31 = column( A, 31UL );

   // Traversing the elements via iterators to non-const elements
   for( auto it=col31.begin(); it!=col31.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=col31.cbegin(); it!=col31.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_column_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse column can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A( 100UL, 10UL );  // Non-initialized 100x10 matrix

   auto col0( column( A, 0UL ) );  // Reference to the 0th column of A

   // The subscript operator provides access to all possible elements of the sparse column,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse column, the element is inserted into the column.
   col0[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the column it is inserted into the column, if it is already contained
   // in the column its value is modified.
   col0.set( 45UL, -1.2 );

   // An alternative for inserting elements into the column is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the column.
   col0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse column is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the column and that the column's capacity is
   // large enough to hold the new element.
   col0.reserve( 10UL );
   col0.append( 51UL, -2.1 );
   \endcode

// \n \section column_common_operations Common Operations
//
// A column view can be used like any other column vector. For instance, the current number of
// column elements can be obtained via the \c size() function, the current capacity via the
// \c capacity() function, and the number of non-zero elements via the \c nonZeros() function.
// However, since columns are references to specific columns of a matrix, several operations
// are not possible, such as resizing and swapping:

   \code
   blaze::DynamicMatrix<int,blaze::columnMajor> A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd column of matrix A
   auto col2 = column( A, 2UL );

   col2.size();          // Returns the number of elements in the column
   col2.capacity();      // Returns the capacity of the column
   col2.nonZeros();      // Returns the number of non-zero elements contained in the column

   col2.resize( 84UL );  // Compilation error: Cannot resize a single column of a matrix

   auto col3 = column( A, 3UL );
   swap( col2, col3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section column_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse columns can be used in all arithmetic operations that any other dense or
// sparse column vector can be used in. The following example gives an impression of the use of
// dense columns within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse columns with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::columnVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::columnVector> c( 2UL );
   c[1] = 3.0;

   blaze::DynamicMatrix<double,blaze::columnMajor> A( 2UL, 4UL );  // Non-initialized 2x4 matrix

   auto col0( column( A, 0UL ) );  // Reference to the 0th column of A

   col0[0] = 0.0;           // Manual initialization of the 0th column of A
   col0[1] = 0.0;
   column( A, 1UL ) = 1.0;  // Homogeneous initialization of the 1st column of A
   column( A, 2UL ) = a;    // Dense vector initialization of the 2nd column of A
   column( A, 3UL ) = c;    // Sparse vector initialization of the 3rd column of A

   b = col0 + a;                 // Dense vector/dense vector addition
   b = c + column( A, 1UL );     // Sparse vector/dense vector addition
   b = col0 * column( A, 2UL );  // Component-wise vector multiplication

   column( A, 1UL ) *= 2.0;     // In-place scaling of the 1st column
   b = column( A, 1UL ) * 2.0;  // Scaling of the 1st column
   b = 2.0 * column( A, 1UL );  // Scaling of the 1st column

   column( A, 2UL ) += a;                 // Addition assignment
   column( A, 2UL ) -= c;                 // Subtraction assignment
   column( A, 2UL ) *= column( A, 0UL );  // Multiplication assignment

   double scalar = trans( c ) * column( A, 1UL );  // Scalar/dot/inner product between two vectors

   A = column( A, 1UL ) * trans( c );  // Outer product between two vectors
   \endcode

// \n \section column_on_row_major_matrix Columns on a Row-Major Matrix
//
// Especially noteworthy is that column views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// possible to traverse a row of a column-major matrix or a column of a row-major matrix. For
// instance:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 1st column of a column-major matrix A
   auto col1 = column( A, 1UL );

   for( auto it=col1.begin(); it!=col1.end(); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a column view on a matrix stored in a row-major fashion
// can result in a considerable performance decrease in comparison to a column view on a matrix
// with column-major storage format. This is due to the non-contiguous storage of the matrix
// elements. Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two row-major matrices
   blaze::DynamicMatrix<double,blaze::rowMajor> A( 128UL, 128UL );
   blaze::DynamicMatrix<double,blaze::rowMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th column of the multiplication between A and B ...
   blaze::DynamicVector<double,blaze::columnVector> x = column( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // A with the 15th column of the row-major matrix B.
   blaze::DynamicVector<double,blaze::columnVector> x = A * column( B, 15UL );
   \endcode

// Although Blaze performs the resulting matrix/vector multiplication as efficiently as possible
// using a column-major storage order for matrix \c A would result in a more efficient evaluation.
*/
//*************************************************************************************************

} // namespace blaze

#endif
