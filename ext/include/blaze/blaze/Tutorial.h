//=================================================================================================
/*!
//  \file blaze/Tutorial.h
//  \brief Tutorial of the Blaze library
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

#ifndef _BLAZE_TUTORIAL_H_
#define _BLAZE_TUTORIAL_H_


//=================================================================================================
//
//  BLAZE TUTORIAL
//
//=================================================================================================

//**Mainpage***************************************************************************************
/*!\mainpage
//
// \image html blaze300x150.jpg
//
// This is the API for the \b Blaze high performance C++ math library. It gives a complete
// overview of the individual features and sublibraries of \b Blaze. To get a first impression
// on \b Blaze, the short \ref getting_started tutorial is a good place to start. Afterwards,
// the following long tutorial covers the most important aspects of the \b Blaze math library.
// The tabs at the top of the page allow a direct access to the individual modules, namespaces,
// classes, and files of the \b Blaze library.\n\n
//
// \section table_of_content Table of Contents
//
// <ul>
//    <li> \ref configuration_and_installation </li>
//    <li> \ref getting_started </li>
//    <li> \ref vectors
//       <ul>
//          <li> \ref vector_types
//             <ul>
//                <li> \ref vector_types_dense_vectors </li>
//                <li> \ref vector_types_sparse_vectors </li>
//             </ul>
//          </li>
//          <li> \ref vector_operations
//             <ul>
//                <li> \ref vector_operations_constructors </li>
//                <li> \ref vector_operations_assignment </li>
//                <li> \ref vector_operations_element_access </li>
//                <li> \ref vector_operations_element_insertion </li>
//                <li> \ref vector_operations_element_removal </li>
//                <li> \ref vector_operations_element_lookup </li>
//                <li> \ref vector_operations_non_modifying_operations </li>
//                <li> \ref vector_operations_modifying_operations </li>
//                <li> \ref vector_operations_arithmetic_operations </li>
//                <li> \ref vector_operations_reduction_operations </li>
//                <li> \ref vector_operations_norms </li>
//                <li> \ref vector_operations_scalar_expansion </li>
//                <li> \ref vector_operations_vector_expansion </li>
//                <li> \ref vector_operations_statistic_operations </li>
//                <li> \ref vector_operations_declaration_operations </li>
//                <li> \ref vector_operations_vector_generators </li>
//             </ul>
//          </li>
//       </ul>
//    </li>
//    <li> \ref matrices
//       <ul>
//          <li> \ref matrix_types
//             <ul>
//                <li> \ref matrix_types_dense_matrices </li>
//                <li> \ref matrix_types_sparse_matrices </li>
//             </ul>
//          </li>
//          <li> \ref matrix_operations
//             <ul>
//                <li> \ref matrix_operations_constructors </li>
//                <li> \ref matrix_operations_assignment </li>
//                <li> \ref matrix_operations_element_access </li>
//                <li> \ref matrix_operations_element_insertion </li>
//                <li> \ref matrix_operations_element_removal </li>
//                <li> \ref matrix_operations_element_lookup </li>
//                <li> \ref matrix_operations_non_modifying_operations </li>
//                <li> \ref matrix_operations_modifying_operations </li>
//                <li> \ref matrix_operations_arithmetic_operations </li>
//                <li> \ref matrix_operations_reduction_operations </li>
//                <li> \ref matrix_operations_norms </li>
//                <li> \ref matrix_operations_scalar_expansion </li>
//                <li> \ref matrix_operations_statistic_operations </li>
//                <li> \ref matrix_operations_declaration_operations </li>
//                <li> \ref matrix_operations_matrix_generators </li>
//                <li> \ref matrix_operations_matrix_inversion </li>
//                <li> \ref matrix_operations_matrix_exponential </li>
//                <li> \ref matrix_operations_decomposition </li>
//                <li> \ref matrix_operations_linear_systems </li>
//                <li> \ref matrix_operations_eigenvalues </li>
//                <li> \ref matrix_operations_singularvalues </li>
//             </ul>
//          </li>
//       </ul>
//    </li>
//    <li> \ref adaptors
//       <ul>
//          <li> \ref adaptors_symmetric_matrices </li>
//          <li> \ref adaptors_hermitian_matrices </li>
//          <li> \ref adaptors_triangular_matrices </li>
//       </ul>
//    </li>
//    <li> \ref views
//       <ul>
//          <li> \ref views_subvectors </li>
//          <li> \ref views_element_selections </li>
//          <li> \ref views_submatrices </li>
//          <li> \ref views_rows </li>
//          <li> \ref views_row_selections </li>
//          <li> \ref views_columns </li>
//          <li> \ref views_column_selections </li>
//          <li> \ref views_bands </li>
//       </ul>
//    </li>
//    <li> \ref arithmetic_operations
//       <ul>
//          <li> \ref addition </li>
//          <li> \ref subtraction </li>
//          <li> \ref scalar_multiplication </li>
//          <li> \ref vector_vector_multiplication
//             <ul>
//                <li> \ref componentwise_multiplication </li>
//                <li> \ref inner_product </li>
//                <li> \ref outer_product </li>
//                <li> \ref cross_product </li>
//                <li> \ref vector_kronecker_product </li>
//             </ul>
//          </li>
//          <li> \ref vector_vector_division </li>
//          <li> \ref matrix_vector_multiplication </li>
//          <li> \ref matrix_matrix_multiplication
//             <ul>
//                <li> \ref schur_product </li>
//                <li> \ref matrix_product </li>
//                <li> \ref matrix_kronecker_product </li>
//             </ul>
//          </li>
//       </ul>
//    </li>
//    <li> \ref bitwise_operations
//       <ul>
//          <li> \ref bitwise_shift </li>
//          <li> \ref bitwise_and </li>
//          <li> \ref bitwise_or </li>
//          <li> \ref bitwise_xor </li>
//       </ul>
//    </li>
//    <li> \ref logical_operations
//       <ul>
//          <li> \ref logical_not </li>
//          <li> \ref logical_and </li>
//          <li> \ref logical_or </li>
//       </ul>
//    </li>
//    <li> \ref shared_memory_parallelization
//       <ul>
//          <li> \ref hpx_parallelization </li>
//          <li> \ref cpp_threads_parallelization </li>
//          <li> \ref boost_threads_parallelization </li>
//          <li> \ref openmp_parallelization </li>
//          <li> \ref serial_execution </li>
//       </ul>
//    </li>
//    <li> \ref serialization
//       <ul>
//          <li> \ref vector_serialization </li>
//          <li> \ref matrix_serialization </li>
//       </ul>
//    </li>
//    <li> \ref customization
//       <ul>
//          <li> \ref configuration_files </li>
//          <li> \ref vector_and_matrix_customization
//             <ul>
//                <li> \ref custom_data_members </li>
//                <li> \ref custom_operations </li>
//                <li> \ref custom_data_types </li>
//             </ul>
//          </li>
//          <li> \ref grouping_tagging </li>
//          <li> \ref error_reporting_customization </li>
//       </ul>
//    </li>
//    <li> \ref blas_functions </li>
//    <li> \ref lapack_functions </li>
//    <li> \ref block_vectors_and_matrices </li>
//    <li> \ref intra_statement_optimization </li>
//    <li> \ref faq </li>
//    <li> \ref issue_creation_guidelines </li>
//    <li> \ref blaze_references </li>
// </ul>
*/
//*************************************************************************************************


//**Configuration and Installation*****************************************************************
/*!\page configuration_and_installation Configuration and Installation
//
// \tableofcontents
//
//
// Since \b Blaze is a header-only library, setting up the \b Blaze library on a particular system
// is a fairly easy two step process. In the following, this two step process is explained in
// detail, preceded only by a short summary of the requirements.
//
//
// \n \section requirements Requirements
// <hr>
//
// For maximum performance the \b Blaze library expects you to have a BLAS library installed
// (<a href="http://software.intel.com/en-us/articles/intel-mkl/">Intel MKL</a>,
// <a href="http://developer.amd.com/libraries/acml/">ACML</a>,
// <a href="http://math-atlas.sourceforge.net">Atlas</a>,
// <a href="http://www.tacc.utexas.edu/tacc-projects/gotoblas2">Goto</a>, ...). If you don't
// have a BLAS library installed on your system, \b Blaze will still work and will not be reduced
// in functionality, but performance may be limited. Thus it is strongly recommended to install a
// BLAS library.
//
// Additionally, for computing the determinant of a dense matrix, for the decomposition of dense
// matrices, for the dense matrix inversion, and for the computation of eigenvalues and singular
// values \b Blaze requires <a href="https://en.wikipedia.org/wiki/LAPACK">LAPACK</a>. When either
// of these features is used it is necessary to link the LAPACK library to the final executable.
// If no LAPACK library is available the use of these features will result in a linker error.
//
// Furthermore, it is possible to use Boost threads to run numeric operations in parallel. In this
// case the Boost library is required to be installed on your system. It is recommended to use the
// newest Boost library available, but \b Blaze requires at minimum the Boost version 1.54.0. If
// you don't have Boost installed on your system, you can download it for free from
// <a href="http://www.boost.org">www.boost.org</a>.
//
//
// \n \section step_1_installation Step 1: Installation
// <hr>
//
// \subsection step_1_cmake Installation via CMake
//
// The first step is the installation of the \b Blaze header files. The most convenient way
// to do this is via <a href="https://cmake.org">CMake</a>. Linux and macOS users can use the
// following two lines to copy the \b Blaze headers in the <tt>./blaze</tt> subdirectory to
// the directory \c ${CMAKE_INSTALL_PREFIX}/include and the package configuration files to
// \c ${CMAKE_INSTALL_PREFIX}/share/blaze/cmake.

   \code
   cmake -DCMAKE_INSTALL_PREFIX=/usr/local/
   sudo make install
   \endcode

// Windows users can do the same via the cmake-gui. Alternatively, it is possible to include
// \b Blaze by adding the following lines in any \c CMakeLists.txt file:

   \code
   find_package( blaze )
   if( blaze_FOUND )
      add_library( blaze_target INTERFACE )
      target_link_libraries( blaze_target INTERFACE blaze::blaze )
   endif()
   \endcode

// Alternatively \b Blaze provides the <tt>./cmake/Blaze_Import</tt> CMake function to import
// the \b Blaze library into CMake based projects. This approach includes the configuration
// step (see \ref step_2_configuration). To do so you need to import the function file like
// any other module/function into your CMake project:

   \code
   list(APPEND CMAKE_MODULE_PATH ${BLAZE_LIBRARY_PATH}/cmake)
   include(Blaze_Import)
   \endcode

// After importing the function script you can import and use the \b Blaze library:

   \code
   Blaze_Import(ARGUMENTS)
   target_link_libraries(TARGET Blaze)
   \endcode

// In this example, \c TARGET is the executable/library using \b Blaze and \c ARGUMENTS is the
// configuration you want for building \b Blaze. To configure \b Blaze using the import function
// you can set the input arguments like this example:

   \code
   Blaze_Import(
      QUIET
      BLAS on
      LAPACK on
      THREADING Boost
      CACHE_SIZE auto
      VECTORIZATION on
      STORAGE_ORDER rowMajor
      THRESHOLD_DMATDVECMULT 100000UL
      THRESHOLD_SMP_DVECDVECADD 1000000UL
   )
   \endcode

// For more details about available configuration options please have a look at
// \ref configuration_files and the <tt>Blaze_Import.cmake</tt> function script.
//
// \n \subsection step_1_vcpkg Installation via the VC++ Packaging Tool
//
// An alternate way to install \b Blaze for Windows users is Microsoft's
// <a href="https://github.com/Microsoft/vcpkg">VC++ Packaging Tool (vcpkg)</a>. \b Blaze can
// be installed via the command line:

   \code
   C:\src\vcpkg> .\vcpkg install blaze
   \endcode

// The tool automatically downloads the latest \b Blaze release and copies the header files to
// the common include directory. Please note that since \b Blaze is a header-only library the
// attempt to install any static or dynamic library will fail!
//
// \n \subsection step_1_installation_unix Manual Installation on Linux/macOS
//
// Since \b Blaze only consists of header files, the <tt>./blaze</tt> subdirectory can be simply
// copied to a standard include directory (note that this requires root privileges):

   \code
   cp -r ./blaze /usr/local/include
   \endcode

// Alternatively, on Unix-based machines (which includes Linux and Mac OS X) the
// \c CPLUS_INCLUDE_PATH environment variable can be set. The specified directory will be
// searched after any directories specified on the command line with the option \c -I and
// before the standard default directories (such as \c /usr/local/include and \c /usr/include).
// Assuming a user named 'Jon', the environment variable can be set as follows:

   \code
   CPLUS_INCLUDE_PATH=/usr/home/jon/blaze
   export CPLUS_INCLUDE_PATH
   \endcode

// Last but not least, the <tt>./blaze</tt> subdirectory can be explicitly specified on the
// command line. The following example demonstrates this by means of the GNU C++ compiler:

   \code
   g++ -I/usr/home/jon/blaze -o BlazeTest BlazeTest.cpp
   \endcode

// \n \subsection step_1_installation_windows Manual Installation on Windows
//
// Windows doesn't have a standard include directory. Therefore the \b Blaze header files can be
// copied to any other directory or simply left in the default \b Blaze directory. However, the
// chosen include directory has to be explicitly specified as include path. In Visual Studio,
// this is done via the project property pages, configuration properties, C/C++, General settings.
// Here the additional include directories can be specified.
//
//
// \n \section step_2_configuration Step 2: Configuration
// <hr>
//
// The second step is the configuration and customization of the \b Blaze library. Many aspects
// of \b Blaze can be adapted to specific requirements, environments and architectures. The most
// convenient way to configure \b Blaze is to modify the headers in the <tt>./blaze/config/</tt>
// subdirectory by means of <a href="https://cmake.org">CMake</a>. Alternatively these header
// files can be customized manually. In both cases, however, the files are modified. If this is
// not an option it is possible to configure \b Blaze via the command line (see the tutorial
// section \ref configuration_files or the documentation in the configuration files).
//
// Since the default settings are reasonable for most systems this step can also be skipped.
// However, in order to achieve maximum performance a customization of at least the following
// configuration files is required:
//
//  - <b><tt><blaze/config/BLAS.h></tt></b>: Via this configuration file \b Blaze can be enabled
//    to use a third-party BLAS library for several basic linear algebra functions (such as for
//    instance dense matrix multiplications). In case no BLAS library is used, all linear algebra
//    functions use the default implementations of the \b Blaze library and therefore BLAS is not a
//    requirement for the compilation process. However, please note that performance may be limited.
//  - <b><tt><blaze/config/CacheSize.h></tt></b>: This file contains the hardware specific cache
//    settings. \b Blaze uses this information to optimize its cache usage. For maximum performance
//    it is recommended to adapt these setting to a specific target architecture.
//  - <b><tt><blaze/config/Thresholds.h></tt></b>: This file contains all thresholds for the
//    customization of the \b Blaze compute kernels. In order to tune the kernels for a specific
//    architecture and to maximize performance it can be necessary to adjust the thresholds,
//    especially for a parallel execution (see \ref shared_memory_parallelization).
//
// For an overview of other customization options and more details, please see the section
// \ref configuration_files.
//
//
// \n \section blaze_version Blaze Version
// <hr>
//
// The current major and minor version number of the \b Blaze library can be found in the
// <b><tt><blaze/system/Version.h></tt></b> header file. It is automatically included via the
// <b><tt><blaze/Blaze.h></tt></b> header file. The file contains the two following macros,
// which can for instance be used for conditional compilation:

   \code
   #define BLAZE_MAJOR_VERSION 3
   #define BLAZE_MINOR_VERSION 8
   #define BLAZE_PATCH_VERSION 0
   \endcode

// \n Next: \ref getting_started
*/
//*************************************************************************************************


//**Getting Started********************************************************************************
/*!\page getting_started Getting Started
//
// This short tutorial serves the purpose to give a quick overview of the way mathematical
// expressions have to be formulated in \b Blaze. Starting with \ref vector_types, the following
// long tutorial covers the most important aspects of the \b Blaze math library.
//
//
// \n \section getting_started_vector_example A First Example
//
// \b Blaze is written such that using mathematical expressions is as close to mathematical
// textbooks as possible and therefore as intuitive as possible. In nearly all cases the seemingly
// easiest solution is the right solution and most users experience no problems when trying to
// use \b Blaze in the most natural way. The following example gives a first impression of the
// formulation of a vector addition in \b Blaze:

   \code
   #include <iostream>
   #include <blaze/Math.h>

   using blaze::StaticVector;
   using blaze::DynamicVector;

   int main()
   {
      // Instantiation of a static 3D column vector. The vector is directly initialized as
      //   ( 4 -2  5 )
      StaticVector<int,3UL> a{ 4, -2, 5 };

      // Instantiation of a dynamic 3D column vector. Via the subscript operator the values are set to
      //   ( 2  5 -3 )
      DynamicVector<int> b( 3UL );
      b[0] = 2;
      b[1] = 5;
      b[2] = -3;

      // Adding the vectors a and b
      DynamicVector<int> c = a + b;

      // Printing the result of the vector addition
      std::cout << "c =\n" << c << "\n";
   }
   \endcode

// Note that the entire \b Blaze math library can be included via the \c blaze/Math.h header
// file. Alternatively, the entire \b Blaze library, including both the math and the entire
// utility module, can be included via the \c blaze/Blaze.h header file. Also note that all
// classes and functions of \b Blaze are contained in the blaze namespace.\n\n
//
// Assuming that this program resides in a source file called \c FirstExample.cpp, it can be
// compiled for instance via the GNU C++ compiler:

   \code
   g++ -std=c++14 -O3 -DNDEBUG -mavx -o FirstExample FirstExample.cpp
   \endcode

// Note the definition of the \c NDEBUG preprocessor symbol. In order to achieve maximum
// performance, it is necessary to compile the program in release mode, which deactivates
// all debugging functionality inside \b Blaze. It is also strongly recommended to specify
// the available architecture specific instruction set (as for instance the AVX instruction
// set, which if available can be activated via the \c -mavx flag). This allows \b Blaze
// to optimize computations via vectorization.\n\n
//
// When running the resulting executable \c FirstExample, the output of the last line of
// this small program is

   \code
   c =
   (           6 )
   (           3 )
   (           2 )
   \endcode

// \n \section getting_started_matrix_example An Example Involving Matrices
//
// Similarly easy and intuitive are expressions involving matrices:

   \code
   #include <iostream>
   #include <blaze/Math.h>

   using namespace blaze;

   int main()
   {
      // Instantiating a dynamic 3D column vector
      DynamicVector<int> x{ 4, -1, 3 };

      // Instantiating a dynamic 2x3 row-major matrix, preinitialized with 0. Via the function call
      // operator three values of the matrix are explicitly set to get the matrix
      //   ( 1  0  4 )
      //   ( 0 -2  0 )
      DynamicMatrix<int> A( 2UL, 3UL, 0 );
      A(0,0) =  1;
      A(0,2) =  4;
      A(1,1) = -2;

      // Performing a matrix/vector multiplication
      DynamicVector<int> y = A * x;

      // Printing the resulting vector
      std::cout << "y =\n" << y << "\n";

      // Instantiating a static column-major matrix. The matrix is directly initialized as
      //   (  3 -1 )
      //   (  0  2 )
      //   ( -1  0 )
      StaticMatrix<int,3UL,2UL,columnMajor> B{ { 3, -1 }, { 0, 2 }, { -1, 0 } };

      // Performing a matrix/matrix multiplication
      DynamicMatrix<int> C = A * B;

      // Printing the resulting matrix
      std::cout << "C =\n" << C << "\n";
   }
   \endcode

// The output of this program is

   \code
   y =
   (          16 )
   (           2 )

   C =
   (           -1           -1 )
   (            0           -4 )
   \endcode

// \n \section getting_started_complex_example A Complex Example
//
// The following example is much more sophisticated. It shows the implementation of the Conjugate
// Gradient (CG) algorithm (http://en.wikipedia.org/wiki/Conjugate_gradient) by means of the
// \b Blaze library:
//
// \image html cg.jpg
//
// In this example it is not important to understand the CG algorithm itself, but to see the
// advantage of the API of the \b Blaze library. In the \b Blaze implementation we will use a
// sparse matrix/dense vector multiplication for a 2D Poisson equation using \f$ N \times N \f$
// unknowns. It becomes apparent that the core of the algorithm is very close to the mathematical
// formulation and therefore has huge advantages in terms of readability and maintainability,
// while the performance of the code is close to the expected theoretical peak performance:

   \code
   #include <blaze/Math.h>

   int main()
   {
      const size_t N ( 1000UL );
      const size_t iterations( 10UL );

      const size_t NN( N*N );

      blaze::CompressedMatrix<double,rowMajor> A( NN, NN );
      blaze::DynamicVector<double,columnVector> x( NN, 1.0 ), b( NN, 0.0 ), r( NN ), p( NN ), Ap( NN );
      double alpha, beta, delta;

      // ... Initializing the sparse matrix A

      // Performing the CG algorithm
      r = b - A * x;
      p = r;
      delta = (r,r);

      for( size_t iteration=0UL; iteration<iterations; ++iteration )
      {
         Ap = A * p;
         alpha = delta / (p,Ap);
         x += alpha * p;
         r -= alpha * Ap;
         beta = (r,r);
         if( std::sqrt( beta ) < 1E-8 ) break;
         p = r + ( beta / delta ) * p;
         delta = beta;
      }
   }
   \endcode

// \n Hopefully this short tutorial gives a good first impression of how mathematical expressions
// are formulated with \b Blaze. The following long tutorial, starting with \ref vector_types,
// will cover all aspects of the \b Blaze math library, i.e. it will introduce all vector and
// matrix types, all possible operations on vectors and matrices, and of course all possible
// mathematical expressions.
//
// \n Previous: \ref configuration_and_installation &nbsp; &nbsp; Next: \ref vectors
*/
//*************************************************************************************************


//**Vectors****************************************************************************************
/*!\page vectors Vectors
//
// \tableofcontents
//
//
// \n \section vectors_general General Concepts
// <hr>
//
// The \b Blaze library currently offers five dense vector types (\ref vector_types_static_vector,
// \ref vector_types_dynamic_vector, \ref vector_types_hybrid_vector, \ref vector_types_custom_vector,
// and \ref vector_types_uniform_vector) and two sparse vector types (\ref vector_types_compressed_vector
// and \ref vector_types_zero_vector). All vectors can be specified as either column vectors or row
// vectors:

   \code
   using blaze::DynamicVector;
   using blaze::columnVector;
   using blaze::rowVector;

   // Setup of the 3-dimensional dense column vector
   //
   //    ( 1 )
   //    ( 2 )
   //    ( 3 )
   //
   DynamicVector<int,columnVector> a{ 1, 2, 3 };

   // Setup of the 3-dimensional dense row vector
   //
   //    ( 4  5  6 )
   //
   DynamicVector<int,rowVector> b{ 4, 5, 6 };
   \endcode

// Per default, all vectors in \b Blaze are column vectors:

   \code
   // Instantiation of a 3-dimensional column vector
   blaze::DynamicVector<int> c( 3UL );
   \endcode

// \n \section vectors_details Vector Details
// <hr>
//
//  - \ref vector_types
//  - \ref vector_operations
//
//
// \n \section vectors_examples Examples
// <hr>

   \code
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::rowVector;
   using blaze::columnVector;

   StaticVector<int,6UL> a;            // Instantiation of a 6-dimensional static column vector
   CompressedVector<int,rowVector> b;  // Instantiation of a compressed row vector
   DynamicVector<int,columnVector> c;  // Instantiation of a dynamic column vector

   // ... Resizing and initialization

   c = a + trans( b );
   \endcode

// \n Previous: \ref getting_started &nbsp; &nbsp; Next: \ref vector_types
*/
//*************************************************************************************************


//**Vector Types***********************************************************************************
/*!\page vector_types Vector Types
//
// \tableofcontents
//
//
// \n \section vector_types_dense_vectors Dense Vectors
// <hr>
//
// \subsection vector_types_static_vector StaticVector
//
// The blaze::StaticVector class template is the representation of a fixed size vector with
// statically allocated elements of arbitrary type. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/StaticVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the number of elements, the transpose flag, the alignment, the
// padding, and the group tag of the vector can be specified via the six template parameters:

   \code
   namespace blaze {

   template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
   class StaticVector;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the vector elements. StaticVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c N   : specifies the total number of vector elements. It is expected that StaticVector is
//             only used for tiny and small vectors.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c AF  : specifies whether the first element of the vector is properly aligned with
//             respect to the available instruction set (SSE, AVX, ...). Possible values are
//             \c blaze::aligned and \c blaze::unaligned. The default value is
//             \c blaze::defaultAlignmentFlag.
//  - \c PF  : specifies whether the vector should be padded to maximize the efficiency of
//             vectorized operations. Possible values are \c blaze::padded and \c blaze::unpadded.
//             The default value is \c blaze::defaultPaddingFlag.
//  - \c Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::StaticVector is perfectly suited for small to medium vectors whose size is known at
// compile time:

   \code
   // Definition of a 3-dimensional integral column vector
   blaze::StaticVector<int,3UL> a;

   // Definition of a 4-dimensional single precision column vector
   blaze::StaticVector<float,4UL,blaze::columnVector> b;

   // Definition of an unaligned, unpadded 6-dimensional double precision row vector
   blaze::StaticVector<double,6UL,blaze::rowVector,blaze::unaligned,blaze::unpadded> c;
   \endcode

// \subsubsection vector_types_static_vector_alignment Alignment
//
// In case \c AF is set to \c blaze::aligned, the elements of a blaze::StaticVector are possibly
// over-aligned to meet the alignment requirements of the available instruction set (SSE, AVX,
// AVX-512, ...). The alignment for fundamental types (\c short, \c int, \c float, \c double, ...)
// and complex types (\c complex<float>, \c complex<double>, ...) is 16 bytes for SSE, 32 bytes
// for AVX, and 64 bytes for AVX-512. All other types are aligned according to their intrinsic
// alignment:

   \code
   struct Int { int i; };

   using VT1 = blaze::StaticVector<double,3UL>;
   using VT2 = blaze::StaticVector<complex<float>,2UL>;
   using VT3 = blaze::StaticVector<Int,5UL>;

   alignof( VT1 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( VT2 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( VT3 );  // Evaluates to 'alignof( Int )'
   \endcode

// Note that an aligned blaze::StaticVector instance may be bigger than the sum of its data
// elements:

   \code
   sizeof( VT1 );  // Evaluates to 32 for both SSE and AVX
   sizeof( VT2 );  // Evaluates to 16 for SSE and 32 for AVX
   sizeof( VT3 );  // Evaluates to 20; no special alignment requirements
   \endcode

// Please note that for this reason an aligned blaze::StaticVector cannot be used in containers
// using dynamic memory such as \c std::vector without additionally providing an allocator that
// can provide over-aligned memory:

   \code
   using Type = blaze::StaticVector<double,3UL>;
   using Allocator = blaze::AlignedAllocator<Type>;

   std::vector<Type> v1;  // Might be misaligned for AVX or AVX-512
   std::vector<Type,Allocator> v2;  // Properly aligned for AVX or AVX-512
   \endcode

// \subsubsection vector_types_static_vector_padding Padding
//
// Adding padding elements to the end of a blaze::StaticVector can have a significant impact on
// the performance. For instance, assuming that AVX is available, then two padded 3-dimensional
// vectors of double precision values can be added via a single SIMD addition operation:

   \code
   using blaze::StaticVector;
   using blaze::columnVector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   StaticVector<double,3UL,columnVector,aligned,padded> a1, b1, c1;
   StaticVector<double,3UL,columnVector,unaligned,unpadded> a2, b2, c2;

   // ... Initialization

   c1 = a1 + b1;  // AVX-based vector addition; maximum performance
   c2 = a2 + b2;  // Scalar vector addition; limited performance

   sizeof( a1 );  // Evaluates to 32 for SSE and AVX, and 64 for AVX-512
   sizeof( a2 );  // Evaluates to 24 for SSE, AVX, and AVX-512 (minimum size)
   \endcode

// Due to padding, the first addition will run at maximum performance. On the flip side, the size
// of each vector instance is increased due to the padding elements. The total size of an instance
// depends on the number of elements and width of the available instruction set (16 bytes for
// SSE, 32 bytes for AVX, and 64 bytes for AVX-512).
//
// The second addition will be limited in performance since due to the number of elements some of
// the elements need to be handled in a scalar operation. However, the size of an \c unaligned,
// \c unpadded blaze::StaticVector instance is guaranteed to be the sum of its elements.
//
// Please also note that \b Blaze will zero initialize the padding elements in order to achieve
// maximum performance!
//
//
// \n \subsection vector_types_dynamic_vector DynamicVector
//
// The blaze::DynamicVector class template is the representation of an arbitrary sized vector
// with dynamically allocated elements of arbitrary type. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/DynamicVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the transpose flag, the type of the allocator, and the group tag of
// the vector can be specified via the four template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Alloc, typename Tag >
   class DynamicVector;

   } // namespace blaze
   \endcode

//  - \c Type : specifies the type of the vector elements. DynamicVector can be used with any
//              non-cv-qualified, non-reference, non-pointer element type.
//  - \c TF   : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//              vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c Alloc: specifies the type of allocator used to allocate dynamic memory. The default type
//              of allocator is \c blaze::AlignedAllocator.
//  - \c Tag  : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//              See \ref grouping_tagging for details.
//
// The blaze::DynamicVector is the default choice for all kinds of dense vectors and the best
// choice for medium to large vectors. Its size can be modified at runtime:

   \code
   // Definition of a 3-dimensional integral column vector
   blaze::DynamicVector<int> a( 3UL );

   // Definition of a 4-dimensional single precision column vector
   blaze::DynamicVector<float,blaze::columnVector> b( 4UL );

   // Definition of a double precision row vector with size 0
   blaze::DynamicVector<double,blaze::rowVector> c;
   \endcode

// \subsubsection vector_types_dynamic_vector_allocators Allocators
//
// Via the third template parameter it is possible to customize the memory allocation of a
// \c blaze::DynamicVector. The provided allocator is expected to represent an implementation of
// the allocator concept of the standard library (see for instance
// <a href="https://en.cppreference.com/w/cpp/container/vector">std::vector</a> and
// <a href="https://en.cppreference.com/w/cpp/memory/allocator">std::allocator</a>). In
// addition, the provided allocator is also required to provide properly (over-)aligned memory
// for fundamental and complex numbers. For instance, in case SSE vectorization is possible, the
// returned memory must be at least 16-byte aligned. In case AVX is active, the memory must be at
// least 32-byte aligned, and in case of AVX-512 the memory must be even 64-byte aligned.
//
//
// \n \subsection vector_types_hybrid_vector HybridVector
//
// The blaze::HybridVector class template combines the advantages of the blaze::StaticVector and
// the blaze::DynamicVector class templates. It represents a fixed size vector with statically
// allocated elements, but still can be dynamically resized (within the bounds of the available
// memory). It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/HybridVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the maximum number of elements, the transpose flag, the alignment,
// the padding, and the group tag of the vector can be specified via the six template parameters:

   \code
   namespace blaze {

   template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
   class HybridVector;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the vector elements. HybridVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c N   : specifies the maximum number of vector elements. It is expected that HybridVector
//             is only used for tiny and small vectors.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c AF  : specifies whether the first element of the vector is properly aligned with
//             respect to the available instruction set (SSE, AVX, ...). Possible values are
//             \c blaze::aligned and \c blaze::unaligned. The default value is
//             \c blaze::defaultAlignmentFlag.
//  - \c PF  : specifies whether the vector should be padded to maximize the efficiency of
//             vectorized operations. Possible values are \c blaze::padded and \c blaze::unpadded.
//             The default value is \c blaze::defaultPaddingFlag.
//  - \c Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::HybridVector is a suitable choice for small to medium vectors, whose size is not
// known at compile time or not fixed at runtime, but whose maximum size is known at compile
// time:

   \code
   // Definition of a 3-dimensional integral column vector with a maximum size of 6
   blaze::HybridVector<int,6UL> a( 3UL );

   // Definition of a 4-dimensional single precision column vector with a maximum size of 16
   blaze::HybridVector<float,16UL,blaze::columnVector> b( 4UL );

   // Definition of a unaligned, unpadded double precision row vector with size 0 and a maximum size of 6
   blaze::HybridVector<double,6UL,blaze::rowVector,blaze::unaligned,blaze::unpadded> c;
   \endcode

// \subsubsection vector_types_hybrid_vector_alignment Alignment
//
// In case \c AF is set to \c blaze::aligned, the elements of a blaze::HybridVector are possibly
// over-aligned to meet the alignment requirements of the available instruction set (SSE, AVX,
// AVX-512, ...). The alignment for fundamental types (\c short, \c int, \c float, \c double, ...)
// and complex types (\c complex<float>, \c complex<double>, ...) is 16 bytes for SSE, 32 bytes
// for AVX, and 64 bytes for AVX-512. All other types are aligned according to their intrinsic
// alignment:

   \code
   struct Int { int i; };

   using VT1 = blaze::HybridVector<double,3UL>;
   using VT2 = blaze::HybridVector<complex<float>,2UL>;
   using VT3 = blaze::HybridVector<Int,5UL>;

   alignof( VT1 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( VT2 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( VT3 );  // Evaluates to 'alignof( Int )'
   \endcode

// Note that an aligned blaze::HybridVector instance may be bigger than an according unaligned
// blaze::HybridVector:

   \code
   sizeof( VT1 );  // Evaluates to 32 for both SSE and AVX
   sizeof( VT2 );  // Evaluates to 16 for SSE and 32 for AVX
   sizeof( VT3 );  // Evaluates to 20; no special alignment requirements
   \endcode

// Please note that for this reason an aligned blaze::HybridVector cannot be used in containers
// using dynamic memory such as \c std::vector without additionally providing an allocator that
// can provide over-aligned memory:

   \code
   using Type = blaze::HybridVector<double,3UL>;
   using Allocator = blaze::AlignedAllocator<Type>;

   std::vector<Type> v1;  // Might be misaligned for AVX or AVX-512
   std::vector<Type,Allocator> v2;  // Properly aligned for AVX or AVX-512
   \endcode

// \subsubsection vector_types_hybrid_vector_padding Padding
//
// Adding padding elements to the end of a blaze::HybridVector can have a significant impact on
// the performance. For instance, assuming that AVX is available, then two padded 3-dimensional
// vectors of double precision values can be added via a single SIMD addition operation:

   \code
   using blaze::HybridVector;
   using blaze::columnVector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   HybridVector<double,3UL,columnVector,aligned,padded> a1, b1, c1;
   HybridVector<double,3UL,columnVector,unaligned,unpadded> a2, b2, c2;

   // ... Resizing and initialization

   c1 = a1 + b1;  // AVX-based vector addition; maximum performance
   c2 = a2 + b2;  // Scalar vector addition; limited performance

   sizeof( a1 );  // Evaluates to 48 for SSE,  64 and AVX, and 128 for AVX-512
   sizeof( a2 );  // Evaluates to 32 for SSE, AVX, and AVX-512 (minimum size)
   \endcode

// Due to padding, the first addition will run at maximum performance. On the flip side, the size
// of each vector instance is increased due to the padding elements. The total size of an instance
// depends on the number of elements and width of the available instruction set (16 bytes for
// SSE, 32 bytes for AVX, and 64 bytes for AVX-512).
//
// The second addition will be limited in performance since due to the number of elements some of
// the elements need to be handled in a scalar operation. However, the size of an \c unaligned,
// \c unpadded blaze::HybridVector instance is guaranteed to be the sum of its elements plus the
// necessary data members to store the current size.
//
// Please also note that \b Blaze will zero initialize the padding elements in order to achieve
// maximum performance!
//
//
// \n \subsection vector_types_custom_vector CustomVector
//
// The blaze::CustomVector class template provides the functionality to represent an external
// array of elements of arbitrary type and a fixed size as a native \b Blaze dense vector data
// structure. Thus in contrast to all other dense vector types a custom vector does not perform
// any kind of memory allocation by itself, but it is provided with an existing array of element
// during construction. A custom vector can therefore be considered an alias to the existing
// array. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/CustomVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the properties of the given array of elements, the transpose flag,
// and the group tag of the vector can be specified via the following five template parameters:

   \code
   namespace blaze {

   template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag >
   class CustomVector;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the vector elements. blaze::CustomVector can be used with
//             any non-cv-qualified, non-reference, non-pointer element type.
//  - \c AF  : specifies whether the represented, external arrays are properly aligned with
//             respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::aligned
//             or \c blaze::unaligned).
//  - \c PF  : specified whether the represented, external arrays are properly padded with
//             respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::padded
//             or \c blaze::unpadded).
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::CustomVector is the right choice if any external array needs to be represented as
// a \b Blaze dense vector data structure or if a custom memory allocation strategy needs to be
// realized:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of an unmanaged custom column vector for unaligned, unpadded integer arrays
   using UnalignedUnpadded = CustomVector<int,unaligned,unpadded,columnVector>;
   std::vector<int> vec( 7UL );
   UnalignedUnpadded a( &vec[0], 7UL );

   // Definition of a managed custom column vector for unaligned but padded 'float' arrays
   using UnalignedPadded = CustomVector<float,unaligned,padded,columnVector>;
   std::unique_ptr<float[]> memory1( new float[16] );
   UnalignedPadded b( memory1.get(), 9UL, 16UL );

   // Definition of a managed custom row vector for aligned, unpadded 'double' arrays
   using AlignedUnpadded = CustomVector<double,aligned,unpadded,rowVector>;
   std::unique_ptr<double[],Deallocate> memory2( blaze::allocate<double>( 7UL ) );
   AlignedUnpadded c( memory2.get(), 7UL );

   // Definition of a managed custom row vector for aligned, padded 'complex<double>' arrays
   using cplx = complex<double>;
   using AlignedPadded = CustomVector<cplx,aligned,padded,columnVector>;
   std::unique_ptr<cplx[],Deallocate> memory3( allocate<cplx>( 8UL ) );
   AlignedPadded d( memory3.get(), 5UL, 8UL );
   \endcode

// In comparison with the remaining \b Blaze dense vector types blaze::CustomVector has several
// special characteristics. All of these result from the fact that a custom vector is not
// performing any kind of memory allocation, but instead is given an existing array of elements.
// The following sections discuss all of these characteristics:
//
//  -# <b>\ref vector_types_custom_vector_memory_management</b>
//  -# <b>\ref vector_types_custom_vector_copy_operations</b>
//  -# <b>\ref vector_types_custom_vector_alignment</b>
//  -# <b>\ref vector_types_custom_vector_padding</b>
//
// \subsubsection vector_types_custom_vector_memory_management Memory Management
//
// The blaze::CustomVector class template acts as an adaptor for an existing array of elements. As
// such it provides everything that is required to use the array just like a native \b Blaze dense
// vector data structure. However, this flexibility comes with the price that the user of a custom
// vector is responsible for the resource management.
//
// The following examples give an impression of several possible types of custom vectors:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of a 3-dimensional custom vector with unaligned, unpadded and externally
   // managed integer array. Note that the std::vector must be guaranteed to outlive the
   // custom vector!
   std::vector<int> vec( 3UL );
   CustomVector<int,unaligned,unpadded> a( &vec[0], 3UL );

   // Definition of a custom vector with size 3 and capacity 16 with aligned, padded and
   // externally managed integer array. Note that the std::unique_ptr must be guaranteed
   // to outlive the custom vector!
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 16UL ) );
   CustomVector<int,aligned,padded> b( memory.get(), 3UL, 16UL );
   \endcode

// \subsubsection vector_types_custom_vector_copy_operations Copy Operations
//
// As with all dense vectors it is possible to copy construct a custom vector:

   \code
   using blaze::CustomVector;
   using blaze::unaligned;
   using blaze::unpadded;

   using CustomType = CustomVector<int,unaligned,unpadded>;

   std::vector<int> vec( 5UL, 10 );  // Vector of 5 integers of the value 10
   CustomType a( &vec[0], 5UL );     // Represent the std::vector as Blaze dense vector
   a[1] = 20;                        // Also modifies the std::vector

   CustomType b( a );  // Creating a copy of vector a
   b[2] = 20;          // Also affects vector a and the std::vector
   \endcode

// It is important to note that a custom vector acts as a reference to the specified array. Thus
// the result of the copy constructor is a new custom vector that is referencing and representing
// the same array as the original custom vector.
//
// In contrast to copy construction, just as with references, copy assignment does not change
// which array is referenced by the custom vector, but modifies the values of the array:

   \code
   std::vector<int> vec2( 5UL, 4 );  // Vector of 5 integers of the value 4
   CustomType c( &vec2[0], 5UL );    // Represent the std::vector as Blaze dense vector

   a = c;  // Copy assignment: Set all values of vector a and b to 4.
   \endcode

// \subsubsection vector_types_custom_vector_alignment Alignment
//
// In case the custom vector is specified as \c aligned the passed array must be guaranteed to
// be aligned according to the requirements of the used instruction set (SSE, AVX, ...). For
// instance, if AVX is active an array of integers must be 32-bit aligned:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unpadded;

   // Allocation of 32-bit aligned memory
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 5UL ) );

   CustomVector<int,aligned,unpadded> a( memory.get(), 5UL );
   \endcode

// In case the alignment requirements are violated, a \c std::invalid_argument exception is
// thrown.
//
// \subsubsection vector_types_custom_vector_padding Padding
//
// Adding padding elements to the end of an array can have a significant impact on the performance.
// For instance, assuming that AVX is available, then two aligned, padded, 3-dimensional vectors
// of double precision values can be added via a single SIMD addition operation:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::padded;

   using CustomType = CustomVector<double,aligned,padded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 4UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 4UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 4UL ) );

   // Creating padded custom vectors of size 3 and a capacity of 4
   CustomType a( memory1.get(), 3UL, 4UL );
   CustomType b( memory2.get(), 3UL, 4UL );
   CustomType c( memory3.get(), 3UL, 4UL );

   // ... Initialization

   c = a + b;  // AVX-based vector addition
   \endcode

// In this example, maximum performance is possible. However, in case no padding elements are
// inserted, a scalar addition has to be used:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unpadded;

   using CustomType = CustomVector<double,aligned,unpadded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 3UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 3UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 3UL ) );

   // Creating unpadded custom vector of size 3
   CustomType a( allocate<double>( 3UL ), 3UL );
   CustomType b( allocate<double>( 3UL ), 3UL );
   CustomType c( allocate<double>( 3UL ), 3UL );

   // ... Initialization

   c = a + b;  // Scalar vector addition
   \endcode

// Note the different number of constructor parameters for unpadded and padded custom vectors:
// In contrast to unpadded vectors, where during the construction only the size of the array
// has to be specified, during the construction of a padded custom vector it is additionally
// necessary to explicitly specify the capacity of the array.
//
// The number of padding elements is required to be sufficient with respect to the available
// instruction set: In case of an aligned padded custom vector the added padding elements must
// guarantee that the capacity is greater or equal than the size and a multiple of the SIMD vector
// width. In case of unaligned padded vectors the number of padding elements can be greater or
// equal the number of padding elements of an aligned padded custom vector. In case the padding
// is insufficient with respect to the available instruction set, a \c std::invalid_argument
// exception is thrown.
//
// Please also note that \b Blaze will zero initialize the padding elements in order to achieve
// maximum performance!
//
//
// \n \subsection vector_types_uniform_vector UniformVector
//
// The blaze::UniformVector class template is the representation of an arbitrary sized uniform
// vector with elements of arbitrary type. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/UniformVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the transpose flag, and the group tag of the vector can be specified
// via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class UniformVector;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the vector elements. UniformVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::UniformVector is the best choice for uniform vectors of any size. Its size can be
// modified at runtime:

   \code
   // Definition of a 3-dimensional integral column vector
   blaze::UniformVector<int> a( 3UL );

   // Definition of a 4-dimensional single precision column vector
   blaze::UniformVector<float,blaze::columnVector> b( 4UL );

   // Definition of a double precision row vector with size 0
   blaze::UniformVector<double,blaze::rowVector> c;
   \endcode

// \n \section vector_types_sparse_vectors Sparse Vectors
// <hr>
//
// \subsection vector_types_compressed_vector CompressedVector
//
// The blaze::CompressedVector class is the representation of an arbitrarily sized sparse
// vector, which stores only non-zero elements of arbitrary type. It can be included via the
// header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/CompressedVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the transpose flag, and the group tag of the vector can be specified
// via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class CompressedVector;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the vector elements. CompressedVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::CompressedVector is the right choice for all kinds of sparse vectors:

   \code
   // Definition of a 3-dimensional integral column vector
   blaze::CompressedVector<int> a( 3UL );

   // Definition of a 4-dimensional single precision column vector with capacity for 3 non-zero elements
   blaze::CompressedVector<float,blaze::columnVector> b( 4UL, 3UL );

   // Definition of a double precision row vector with size 0
   blaze::CompressedVector<double,blaze::rowVector> c;
   \endcode

// \n \subsection vector_types_zero_vector ZeroVector
//
// The blaze::ZeroVector class template is the representation of an immutable, arbitrary sized
// zero vector with elements of arbitrary type. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/ZeroVector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the transpose flag, and the group tag of the vector can be specified
// via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class ZeroVector;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the vector elements. ZeroVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - \c Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::ZeroVector is the perfect choice to represent a zero vector:

   \code
   // Definition of a 3-dimensional integral zero column vector
   blaze::ZeroVector<int> a( 3UL );

   // Definition of a 6-dimensional single precision zero column vector
   blaze::ZeroVector<float,blaze::columnVector> b( 6UL );

   // Definition of a double precision row vector with size 0
   blaze::ZeroVector<double,blaze::rowVector> c;
   \endcode

// \n Previous: \ref vectors &nbsp; &nbsp; Next: \ref vector_operations
*/
//*************************************************************************************************


//**Vector Operations******************************************************************************
/*!\page vector_operations Vector Operations
//
// \tableofcontents
//
//
// \n \section vector_operations_constructors Constructors
// <hr>
//
// Instantiating and setting up a vector is very easy and intuitive. However, there are a few
// rules to take care of:
//  - In case the last template parameter (the transpose flag) is omitted, the vector is per
//    default a column vector.
//  - The elements of a \c StaticVector or \c HybridVector are default initialized (i.e. built-in
//    data types are initialized to 0, class types are initialized via the default constructor).
//  - Newly allocated elements of a \c DynamicVector or \c CompressedVector remain uninitialized
//    if they are of built-in type and are default constructed if they are of class type.
//
// \n \subsection vector_operations_default_construction Default Construction

   \code
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   // All vectors can be default constructed. Whereas the size
   // of StaticVectors is fixed via the second template parameter,
   // the initial size of a default constructed DynamicVector or
   // CompressedVector is 0.
   StaticVector<int,2UL> v1;                // Instantiation of a 2D integer column vector.
                                            // All elements are initialized to 0.
   StaticVector<long,3UL,columnVector> v2;  // Instantiation of a 3D long integer column vector.
                                            // Again, all elements are initialized to 0L.
   DynamicVector<float> v3;                 // Instantiation of a dynamic single precision column
                                            // vector of size 0.
   DynamicVector<double,rowVector> v4;      // Instantiation of a dynamic double precision row
                                            // vector of size 0.
   CompressedVector<int> v5;                // Instantiation of a compressed integer column
                                            // vector of size 0.
   CompressedVector<double,rowVector> v6;   // Instantiation of a compressed double precision row
                                            // vector of size 0.
   \endcode

// \n \subsection vector_operations_size_construction Construction with Specific Size
//
// The \c DynamicVector, \c HybridVector and \c CompressedVector classes offer a constructor that
// allows to immediately give the vector the required size. Whereas both dense vectors (i.e.
// \c DynamicVector and \c HybridVector) use this information to allocate memory for all vector
// elements, \c CompressedVector merely acquires the size but remains empty.

   \code
   DynamicVector<int,columnVector> v7( 9UL );      // Instantiation of an integer dynamic column vector
                                                   // of size 9. The elements are NOT initialized!
   HybridVector< complex<float>, 5UL > v8( 2UL );  // Instantiation of a column vector with two single
                                                   // precision complex values. The elements are
                                                   // default constructed.
   CompressedVector<int,rowVector> v9( 10UL );     // Instantiation of a compressed row vector with
                                                   // size 10. Initially, the vector provides no
                                                   // capacity for non-zero elements.
   \endcode

// \n \subsection vector_operations_initialization_constructors Initialization Constructors
//
// All dense vector classes offer a constructor that allows for a direct, homogeneous initialization
// of all vector elements. In contrast, for sparse vectors the predicted number of non-zero elements
// can be specified

   \code
   StaticVector<int,3UL,rowVector> v10( 2 );            // Instantiation of a 3D integer row vector.
                                                        // All elements are initialized to 2.
   DynamicVector<float> v11( 3UL, 7.0F );               // Instantiation of a dynamic single precision
                                                        // column vector of size 3. All elements are
                                                        // set to 7.0F.
   CompressedVector<float,rowVector> v12( 15UL, 3UL );  // Instantiation of a single precision column
                                                        // vector of size 15, which provides enough
                                                        // space for at least 3 non-zero elements.
   \endcode

// \n \subsection vector_operations_array_construction Array Construction
//
// Alternatively, all dense vector classes offer a constructor for an initialization with a dynamic
// or static array, or with a \c std::array. If the vector is initialized from a dynamic array, the
// constructor expects the actual size of the array as first argument, the array as second argument.
// In case of a static array or \c std::array, the fixed size of the array is used:

   \code
   const unique_ptr<double[]> array1( new double[2] );
   // ... Initialization of the dynamic array
   blaze::StaticVector<double,2UL> v13( 2UL, array1.get() );

   const int array2[4] = { 4, -5, -6, 7 };
   blaze::StaticVector<int,4UL> v14( array2 );

   const std::array<float,3UL> array3{ 1.1F, 2.2F, 3.3F };
   blaze::StaticVector<float,3UL> v15( array3 );
   \endcode

// \n \subsection vector_operations_initializer_list_construction Initializer List Construction
//
// In addition, all dense and sparse vector classes can be directly initialized by means of an
// initializer list:

   \code
   blaze::DynamicVector<float> v16{ 1.0F, 2.0F, 3.0F, 4.0F };
   blaze::CompressedVector<int> v17{ 0, 2, 0, 0, 5, 0, 7, 0 };
   \endcode

// Dynamically sized vectors (such as e.g. \ref vector_types_hybrid_vector,
// \ref vector_types_dynamic_vector or \ref vector_types_compressed_vector) are sized according
// to the size of the initializer list and all their elements are (copy) assigned the values of
// the list. For fixed size vectors (such as e.g. \ref vector_types_static_vector) missing values
// are initialized as default and in case the size of the initializer list exceeds the size
// of the vector a \c std::invalid_argument exception is thrown. In case of sparse vectors, only
// the non-zero elements are used to initialize the vector.
//
// \n \subsection vector_operations_copy_construction Copy Construction
//
// All dense and sparse vectors can be created as the copy of any other dense or sparse vector
// with the same transpose flag (i.e. blaze::rowVector or blaze::columnVector).

   \code
   StaticVector<int,9UL,columnVector> v18( v7 );  // Instantiation of the dense column vector v17
                                                  // as copy of the dense column vector v7.
   DynamicVector<int,rowVector> v19( v9 );        // Instantiation of the dense row vector v18 as
                                                  // copy of the sparse row vector v9.
   CompressedVector<int,columnVector> v20( v1 );  // Instantiation of the sparse column vector v19
                                                  // as copy of the dense column vector v1.
   CompressedVector<float,rowVector> v21( v12 );  // Instantiation of the sparse row vector v20 as
                                                  // copy of the row vector v12.
   \endcode

// Note that it is not possible to create a \c StaticVector as a copy of a vector with a different
// size:

   \code
   StaticVector<int,5UL,columnVector> v22( v7 );  // Runtime error: Size does not match!
   StaticVector<int,4UL,rowVector> v23( v10 );    // Compile time error: Size does not match!
   \endcode

// \n \section vector_operations_assignment Assignment
// <hr>
//
// There are several types of assignment to dense and sparse vectors:
// \ref vector_operations_homogeneous_assignment, \ref vector_operations_array_assignment,
// \ref vector_operations_copy_assignment, and \ref vector_operations_compound_assignment.
//
// \n \subsection vector_operations_homogeneous_assignment Homogeneous Assignment
//
// Sometimes it may be necessary to assign the same value to all elements of a dense vector.
// For this purpose, the assignment operator can be used:

   \code
   blaze::StaticVector<int,3UL> v1;
   blaze::DynamicVector<double> v2;

   // Setting all integer elements of the StaticVector to 2
   v1 = 2;

   // Setting all double precision elements of the DynamicVector to 5.0
   v2 = 5.0;
   \endcode

// \n \subsection vector_operations_array_assignment Array Assignment
//
// Dense vectors can also be assigned a static array or \c std::array:

   \code
   blaze::StaticVector<float,2UL> v1;
   blaze::DynamicVector<double,rowVector> v2;

   const float array1[2] = { 1.0F, 2.0F };
   const std::array<double,5UL> array2{ 2.1, 4.0, -1.7, 8.6, -7.2 };

   v1 = array1;
   v2 = array2;
   \endcode

// \n \subsection vector_operations_initializer_list_assignment Initializer List Assignment
//
// Alternatively, it is possible to directly assign an initializer list to a dense or sparse
// vector:

   \code
   blaze::DynamicVector<float> v1;
   blaze::CompressedVector<double,rowVector> v2;

   v1 = { 1.0F, 2.0F };
   v2 = { 2.1, 0.0, -1.7, 0.0, -7.2 };
   \endcode

// Dynamically sized vectors (such as e.g. \ref vector_types_hybrid_vector,
// \ref vector_types_dynamic_vector or \ref vector_types_compressed_vector) are resized according
// to the size of the initializer list and all their elements are (copy) assigned the values of
// the list. For fixed size vectors (such as e.g. \ref vector_types_static_vector) missing values
// are reset to their default value and in case the size of the initializer list exceeds the size
// of the vector a \c std::invalid_argument exception is thrown. In case of sparse vectors, only
// the non-zero elements are considered.
//
// \n \subsection vector_operations_copy_assignment Copy Assignment
//
// For all vector types it is generally possible to assign another vector with the same transpose
// flag (i.e. blaze::columnVector or blaze::rowVector). Note that in case of \c StaticVectors, the
// assigned vector is required to have the same size as the \c StaticVector since the size of a
// \c StaticVector cannot be adapted!

   \code
   blaze::StaticVector<int,3UL,columnVector> v1;
   blaze::DynamicVector<int,columnVector>    v2( 3UL );
   blaze::DynamicVector<float,columnVector>  v3( 5UL );
   blaze::CompressedVector<int,columnVector> v4( 3UL );
   blaze::CompressedVector<float,rowVector>  v5( 3UL );

   // ... Initialization of the vectors

   v1 = v2;  // OK: Assignment of a 3D dense column vector to another 3D dense column vector
   v1 = v4;  // OK: Assignment of a 3D sparse column vector to a 3D dense column vector
   v1 = v3;  // Runtime error: Cannot assign a 5D vector to a 3D static vector
   v1 = v5;  // Compilation error: Cannot assign a row vector to a column vector
   \endcode

// \n \subsection vector_operations_compound_assignment Compound Assignment
//
// Next to plain assignment, it is also possible to use addition assignment, subtraction
// assignment, and multiplication assignment. Note however, that in contrast to plain assignment
// the size and the transpose flag of the vectors has be to equal in order to able to perform a
// compound assignment.

   \code
   blaze::StaticVector<int,5UL,columnVector>   v1;
   blaze::DynamicVector<int,columnVector>      v2( 5UL );
   blaze::CompressedVector<float,columnVector> v3( 7UL );
   blaze::DynamicVector<float,rowVector>       v4( 7UL );
   blaze::CompressedVector<float,rowVector>    v5( 7UL );

   // ... Initialization of the vectors

   v1 += v2;  // OK: Addition assignment between two column vectors of the same size
   v1 += v3;  // Runtime error: No compound assignment between vectors of different size
   v1 -= v4;  // Compilation error: No compound assignment between vectors of different transpose flag
   v4 *= v5;  // OK: Multiplication assignment between two row vectors of the same size
   \endcode

// \n \section vector_operations_element_access Element Access
// <hr>
//
// \subsection vector_operations_subscript_operator_1 Subscript Operator
//
// The easiest and most intuitive way to access a dense or sparse vector is via the subscript
// operator. The indices to access a vector are zero-based:

   \code
   blaze::DynamicVector<int> v1( 5UL );
   v1[0] = 1;
   v1[1] = 3;
   // ...

   blaze::CompressedVector<float> v2( 5UL );
   v2[2] = 7.3F;
   v2[4] = -1.4F;
   \endcode

// Whereas using the subscript operator on a dense vector only accesses the already existing
// element, accessing an element of a sparse vector via the subscript operator potentially
// inserts the element into the vector and may therefore be more expensive. Consider the
// following example:

   \code
   blaze::CompressedVector<int> v1( 10UL );

   for( size_t i=0UL; i<v1.size(); ++i ) {
      ... = v1[i];
   }
   \endcode

// Although the compressed vector is only used for read access within the for loop, using the
// subscript operator temporarily inserts 10 non-zero elements into the vector. Therefore the
// preferred way to traverse the non-zero elements of a sparse vector is to use iterators.
//
// \n \subsection vector_operations_iterators Iterators
//
// An alternate way to traverse the elements contained in a dense or sparse vector is by means
// of iterators. For that purpose, all vectors provide the \c begin(), \c cbegin(), \c end(),
// and \c cend() members functions. In case of non-const vectors, \c begin() and \c end() return
// an \c Iterator, which allows a manipulation of the (non-zero) value. In case of a constant
// vector or in case \c cbegin() or \c cend() are used a \c ConstIterator is returned. Iterators
// on dense vectors traverse all elements of the vector, including the zero elements. Iterators
// on sparse vectors only traverse the non-zero elements.
//
// The following two examples demonstrate how to traverse the elements of a dense and sparse
// vector, respectively:

   \code
   using blaze::DynamicVector;

   DynamicVector<int> v1( 10UL );

   // Traversing all elements contained in the vector by Iterator
   for( DynamicVector<int>::Iterator it=v1.begin(); it!=v1.end(); ++it ) {
      *it = ...;  // OK: Write access to the value of the element.
      ... = *it;  // OK: Read access to the value of the element.
   }

   // Traversing all elements contained in the vector by ConstIterator
   for( DynamicVector<int>::ConstIterator it=v1.cbegin(); it!=v1.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the value of the element.
   }

   // Traversing the vector elements by means of a range-based for loop
   for( int& i : v1 ) {
      i = ...;  // OK: Write access to the value of the element.
      ... = i;  // OK: Read access to the value of the element.
   }
   \endcode

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v2( 10UL );

   // ... Initialization of the vector

   // Traversing the non-zero elements contained in the vector by Iterator
   for( CompressedVector<int>::Iterator it=v2.begin(); it!=v2.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the non-zero element.
   }

   // Traversing the non-zero elements contained in the vector by ConstIterator
   for( CompressedVector<int>::ConstIterator it=v2.cbegin(); it!=v2.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the non-zero element.
   }
   \endcode

// Note that \c begin(), \c cbegin(), \c end(), and \c cend() are also available as free functions:

   \code
   for( CompressedVector<int>::Iterator it=begin( v2 ); it!=end( v2 ); ++it ) {
      // ...
   }

   for( CompressedVector<int>::ConstIterator it=cbegin( v2 ); it!=cend( v2 ); ++it ) {
      // ...
   }
   \endcode

// \n \subsection vector_operations_data .data() / data()
//
// Sometimes it is necessary to acquire a pointer to the first element of the underlying array
// of a dense vector. For that purpose the \c data() member function or the free \c data() function
// can be used:

   \code
   // Instantiating a dynamic vector with 10 elements
   blaze::DynamicVector<int> v( 10UL );
   v.data();   // Returns a pointer to the first element of the dynamic vector
   data( v );  // Same effect as the member function
   \endcode

// \n \section vector_operations_element_insertion Element Insertion
// <hr>
//
// In contrast to dense vectors, that store all elements independent of their value and that
// offer direct access to all elements, sparse vectors only store the non-zero elements contained
// in the vector. Therefore it is necessary to explicitly add elements to the vector.
//
// \n \subsection vector_operations_subscript_operator_2 Subscript Operator
//
// The first option to add elements to a sparse vector is the subscript operator:

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v1( 3UL );
   v1[1] = 2;
   \endcode

// In case the element at the given index is not yet contained in the vector, it is automatically
// inserted. Otherwise the old value is replaced by the new value 2. The operator returns a
// reference to the sparse vector element.
//
// \n \subsection vector_operations_set .set()
//
// An alternative to the subscript operator is the \c set() function: In case the element is not
// yet contained in the vector the element is inserted, else the element's value is modified:

   \code
   // Insert or modify the value at index 3
   v1.set( 3, 1 );
   \endcode

// \n \subsection vector_operations_insert .insert()
//
// The insertion of elements can be better controlled via the \c insert() function. In contrast to
// the subscript operator and the \c set() function it emits an exception in case the element is
// already contained in the vector. In order to check for this case, the \c find() function can be
// used:

   \code
   // In case the element at index 4 is not yet contained in the matrix it is inserted
   // with a value of 6.
   if( v1.find( 4 ) == v1.end() )
      v1.insert( 4, 6 );
   \endcode

// \n \subsection vector_operations_append .append()
//
// Although the \c insert() function is very flexible, due to performance reasons it is not suited
// for the setup of large sparse vectors. A very efficient, yet also very low-level way to fill
// a sparse vector is the \c append() function. It requires the sparse vector to provide enough
// capacity to insert a new element. Additionally, the index of the new element must be larger
// than the index of the previous element. Violating these conditions results in undefined
// behavior!

   \code
   v1.reserve( 10 );    // Reserving space for 10 non-zero elements
   v1.append( 5, -2 );  // Appending the element -2 at index 5
   v1.append( 6,  4 );  // Appending the element 4 at index 6
   // ...
   \endcode

// \n \section vector_operations_element_removal Element Removal
// <hr>
//
// \subsection vector_operations_erase .erase()
//
// The \c erase() member functions can be used to remove elements from a sparse vector. The
// following example gives an impression of the five different flavors of \c erase():

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v( 42 );
   // ... Initialization of the vector

   // Erasing the element at index 21
   v.erase( 21 );

   // Erasing a single element via iterator
   v.erase( v.find( 4 ) );

   // Erasing all non-zero elements in the range [7..24]
   v.erase( v.lowerBound( 7 ), v.upperBound( 24 ) );

   // Erasing all non-zero elements with a value larger than 9 by passing a unary predicate
   v.erase( []( int i ){ return i > 9; } );

   // Erasing all non-zero elements in the range [30..40] with a value larger than 5
   v.erase( v.lowerBound( 30 ), v.upperBound( 40 ), []( int i ){ return i > 5; } );
   \endcode

// \n \section vector_operations_element_lookup Element Lookup
// <hr>
//
// A sparse vector only stores the non-zero elements contained in the vector. Therefore, whenever
// accessing a vector element at a specific index a lookup operation is required. Whereas the
// subscript operator is performing this lookup automatically, it is also possible to use the
// \c find(), \c lowerBound(), and \c upperBound() member functions for a manual lookup.
//
// \n \subsection vector_operations_find .find() / find()
//
// The \c find() function can be used to check whether a specific element is contained in a sparse
// vector. It specifically searches for the element at the given index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past the
// last non-zero element of the compressed vector (the \c end() iterator) is returned. Note that
// the returned iterator is subject to invalidation due to inserting operations via the subscript
// operator, the \c set() function or the \c insert() function!

   \code
   using blaze::CompressedVector;

   CompressedVector<int> a( 42 );
   // ... Initialization of the vector

   // Searching the element at index 7. In case the element is not
   // contained in the vector, the end() iterator is returned.
   CompressedVector<int>::Iterator pos( a.find( 7 ) );

   if( pos != a.end( 7 ) ) {
      // ...
   }
   \endcode

// Alternatively, the free function \c find() can be used to find a specific element in a sparse
// vector:

   \code
   find( a, 7 );  // Searching the element at index 7; same effect as the member function
   \endcode

// \n \subsection vector_operations_lowerbound .lowerBound() / lowerBound()
//
// The \c lowerBound() function returns an iterator to the first element with an index not less
// then the given index. In combination with the \c upperBound() function this function can be
// used to create a pair of iterators specifying a range of indices. Note that the returned
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the \c set() function or the \c insert() function!

   \code
   using blaze::CompressedVector;

   CompressedVector<int> a( 42 );
   // ... Initialization of the vector

   // Searching the lower bound of index 17.
   CompressedVector<int>::Iterator pos1( a.lowerBound( 17 ) );

   // Searching the upper bound of index 28
   CompressedVector<int>::Iterator pos2( a.upperBound( 28 ) );

   // Erasing all elements in the specified range
   a.erase( pos1, pos2 );
   \endcode

// Alternatively, the free function \c lowerBound() can be used to:

   \code
   lowerBound( a, 17 );  // Searching the lower bound of index 17; same effect as the member function
   \endcode

// \n \subsection vector_operations_upperbound .upperBound() / upperBound()
//
// The \c upperBound() function returns an iterator to the first element with an index greater then
// the given index. In combination with the \c lowerBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned iterator is
// subject to invalidation due to inserting operations via the subscript operator, the \c set()
// function or the \c insert() function!

   \code
   using blaze::CompressedVector;

   CompressedVector<int> a( 42 );
   // ... Initialization of the vector

   // Searching the lower bound of index 17.
   CompressedVector<int>::Iterator pos1( a.lowerBound( 17 ) );

   // Searching the upper bound of index 28
   CompressedVector<int>::Iterator pos2( a.upperBound( 28 ) );

   // Erasing all elements in the specified range
   a.erase( pos1, pos2 );
   \endcode

// Alternatively, the free function \c upperBound() can be used to:

   \code
   upperBound( a, 28 );  // Searching the upper bound of index 28; same effect as the member function
   \endcode

// \n \section vector_operations_non_modifying_operations Non-Modifying Operations
// <hr>
//
// \subsection vector_operations_size .size() / size()
//
// Via the \c size() member function, the current size of a dense or sparse vector can be queried:

   \code
   // Instantiating a dynamic vector with size 10
   blaze::DynamicVector<int> v1( 10UL );
   v1.size();  // Returns 10

   // Instantiating a compressed vector with size 12 and capacity for 3 non-zero elements
   blaze::CompressedVector<double> v2( 12UL, 3UL );
   v2.size();  // Returns 12
   \endcode

// Alternatively, the free function \c size() can be used to query to current size of a vector.
// In contrast to the member function, the free function can also be used to query the size of
// vector expressions:

   \code
   size( v1 );  // Returns 10, i.e. has the same effect as the member function
   size( v2 );  // Returns 12, i.e. has the same effect as the member function

   blaze::DynamicMatrix<int> A( 15UL, 12UL );
   size( A * v2 );  // Returns 15, i.e. the size of the resulting vector
   \endcode

// \n \subsection vector_operations_capacity .capacity() / capacity()
//
// Via the \c capacity() (member) function the internal capacity of a dense or sparse vector
// can be queried. Note that the capacity of a vector doesn't have to be equal to the size
// of a vector. In case of a dense vector the capacity will always be greater or equal than
// the size of the vector, in case of a sparse vector the capacity may even be less than
// the size.

   \code
   v1.capacity();   // Returns at least 10
   \endcode

// For symmetry reasons, there is also a free function /c capacity() available that can be used
// to query the capacity:

   \code
   capacity( v1 );  // Returns at least 10, i.e. has the same effect as the member function
   \endcode

// Note, however, that it is not possible to query the capacity of a vector expression:

   \code
   capacity( A * v1 );  // Compilation error!
   \endcode

// \n \subsection vector_operations_nonzeros .nonZeros() / nonZeros()
//
// For both dense and sparse vectors the number of non-zero elements can be determined via the
// \c nonZeros() member function. Sparse vectors directly return their number of non-zero
// elements, dense vectors traverse their elements and count the number of non-zero elements.

   \code
   v1.nonZeros();  // Returns the number of non-zero elements in the dense vector
   v2.nonZeros();  // Returns the number of non-zero elements in the sparse vector
   \endcode

// There is also a free function \c nonZeros() available to query the current number of non-zero
// elements:

   \code
   nonZeros( v1 );  // Returns the number of non-zero elements in the dense vector
   nonZeros( v2 );  // Returns the number of non-zero elements in the sparse vector
   \endcode

// The free \c nonZeros() function can also be used to query the number of non-zero elements in
// a vector expression. However, the result is not the exact number of non-zero elements, but
// may be a rough estimation:

   \code
   nonZeros( A * v1 );  // Estimates the number of non-zero elements in the vector expression
   \endcode

// \n \subsection vector_operations_isempty isEmpty()
//
// The \c isEmpty() function returns whether the total number of elements of the vector is zero:

   \code
   blaze::DynamicVector<int> a;  // Create an empty vector
   isEmpty( a );                 // Returns true
   a.resize( 10 );               // Resize to 10 elements
   isEmpty( a );                 // Returns false
   \endcode

// \n \subsection vector_operations_isnan isnan()
//
// The \c isnan() function provides the means to check a dense or sparse vector for non-a-number
// elements:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   if( isnan( a ) ) { ... }
   \endcode

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isnan( a ) ) { ... }
   \endcode

// If at least one element of the vector is not-a-number, the function returns \c true, otherwise
// it returns \c false.
//
//
// \n \subsection vector_operations_isinf isinf()
//
// The \c isinf() function checks the given dense or sparse vector for infinite (\c inf) elements:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   if( isinf( a ) ) { ... }
   \endcode

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isinf( a ) ) { ... }
   \endcode

// If at least one element of the vector is infinite, the function returns \c true, otherwise it
// returns \c false.
//
//
// \n \subsection vector_operations_isfinite isfinite()
//
// The \c isfinite() function checks if all elements of the given dense or sparse vector are
// finite elements (i.e. normal, subnormal or zero elements, but not infinite or NaN):

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   if( isfinite( a ) ) { ... }
   \endcode

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isfinite( a ) ) { ... }
   \endcode

// If all elements of the vector are finite, the function returns \c true, otherwise it returns
// \c false.
//
//
// \n \subsection vector_operations_isdefault isDefault()
//
// The \c isDefault() function returns whether the given dense or sparse vector is in default state:

   \code
   blaze::HybridVector<int,20UL> a;
   // ... Resizing and initialization
   if( isDefault( a ) ) { ... }
   \endcode

// A vector is in default state if it appears to just have been default constructed. All resizable
// vectors (\c HybridVector, \c DynamicVector, or \c CompressedVector) and \c CustomVector are
// in default state if its size is equal to zero. A non-resizable vector (\c StaticVector, all
// subvectors, element selections, rows, and columns) is in default state if all its elements are
// in default state. For instance, in case the vector is instantiated for a built-in integral or
// floating point data type, the function returns \c true in case all vector elements are 0 and
// \c false in case any vector element is not 0.
//
//
// \n \subsection vector_operations_isUniform isUniform()
//
// In order to check if all vector elements are identical, the \c isUniform() function can be used:

   \code
   blaze::DynamicVector<int> a;
   // ... Resizing and initialization
   if( isUniform( a ) ) { ... }
   \endcode

// Note that in case of sparse vectors the zero elements are also taken into account!
//
//
// \n \subsection vector_operations_isZero isZero()
//
// In order to check if all vector elements are zero, the \c isZero() function can be used:

   \code
   blaze::DynamicVector<int> a;
   // ... Resizing and initialization
   if( isZero( a ) ) { ... }
   \endcode

// \n \subsection vector_operations_length length() / sqrLength()
//
// In order to calculate the length (magnitude) of a dense or sparse vector, both the \c length()
// and \c sqrLength() function can be used:

   \code
   blaze::StaticVector<float,3UL,rowVector> v{ -1.2F, 2.7F, -2.3F };

   const float len    = length   ( v );  // Computes the current length of the vector
   const float sqrlen = sqrLength( v );  // Computes the square length of the vector
   \endcode

// Note that both functions can only be used for vectors with built-in or complex element type!
//
//
// \n \subsection vector_operations_vector_trans trans()
//
// As already mentioned, vectors can either be column vectors (blaze::columnVector) or row vectors
// (blaze::rowVector). A column vector cannot be assigned to a row vector and vice versa. However,
// vectors can be transposed via the \c trans() function:

   \code
   blaze::DynamicVector<int,columnVector> v1( 4UL );
   blaze::CompressedVector<int,rowVector> v2( 4UL );

   v1 = v2;            // Compilation error: Cannot assign a row vector to a column vector
   v1 = trans( v2 );   // OK: Transposing the row vector to a column vector and assigning it
                       //     to the column vector v1
   v2 = trans( v1 );   // OK: Transposing the column vector v1 and assigning it to the row vector v2
   v1 += trans( v2 );  // OK: Addition assignment of two column vectors
   \endcode

// \n \subsection vector_operations_ctrans ctrans()
//
// It is also possible to compute the conjugate transpose of a vector. This operation is available
// via the \c ctrans() function:

   \code
   blaze::CompressedVector< complex<float>, rowVector > v1( 4UL );
   blaze::DynamicVector< complex<float>, columnVector > v2( 4UL );

   v1 = ctrans( v2 );  // Compute the conjugate transpose vector
   \endcode

// Note that the \c ctrans() function has the same effect as manually applying the \c conj() and
// \c trans() function in any order:

   \code
   v1 = trans( conj( v2 ) );  // Computing the conjugate transpose vector
   v1 = conj( trans( v2 ) );  // Computing the conjugate transpose vector
   \endcode

// \n \subsection vector_operations_reverse reverse()
//
// Via the \c reverse() function is is possible to reverse the elements of a dense or sparse
// vector. The following examples demonstrates this by means of a dense vector:

   \code
   blaze::DynamicVector<int> a{ 1, 2, 3, 4, 5 };
   blaze::DynamicVector<int> b;

   b = reverse( a );  // Results in ( 5 4 3 2 1 )
   \endcode

// \n \subsection vector_operations_evaluate eval() / evaluate()
//
// The \c evaluate() function forces an evaluation of the given vector expression and enables
// an automatic deduction of the correct result type of an operation. The following code example
// demonstrates its intended use for the multiplication of a dense and a sparse vector:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   auto c = evaluate( a * b );
   \endcode

// In this scenario, the \c evaluate() function assists in deducing the exact result type of
// the operation via the \c auto keyword. Please note that if \c evaluate() is used in this
// way, no temporary vector is created and no copy operation is performed. Instead, the result
// is directly written to the target vector due to the return value optimization (RVO). However,
// if \c evaluate() is used in combination with an explicit target type, a temporary will be
// created and a copy operation will be performed if the used type differs from the type
// returned from the function:

   \code
   CompressedVector<double> d( a * b );  // No temporary & no copy operation
   DynamicVector<double> e( a * b );     // Temporary & copy operation
   d = evaluate( a * b );                // Temporary & copy operation
   \endcode

// Sometimes it might be desirable to explicitly evaluate a sub-expression within a larger
// expression. However, please note that \c evaluate() is not intended to be used for this
// purpose. This task is more elegantly and efficiently handled by the \c eval() function:

   \code
   blaze::DynamicVector<double> a, b, c, d;

   d = a + evaluate( b * c );  // Unnecessary creation of a temporary vector
   d = a + eval( b * c );      // No creation of a temporary vector
   \endcode

// In contrast to the \c evaluate() function, \c eval() can take the complete expression
// into account and therefore can guarantee the most efficient way to evaluate it (see also
// \ref intra_statement_optimization).
//
// \n \subsection vector_operations_noalias noalias()
//
// The \b Blaze library is able to reliably detect aliasing during the assignment of vectors.
// In case the aliasing would lead to an incorrect result, \b Blaze introduces an intermediate
// temporary of the appropriate type to break the aliasing. For instance, in the following
// example \b Blaze performs an alias detection in both assignments, but only, in the second
// assignment it detects a problematic aliasing and uses an intermediate temporary in order
// to be able to compute the correct result:

   \code
   blaze::DynamicVector<double> x, y;
   blaze::DynamicMatrix<double> A;

   x = x + y;  // No problematic aliasing of x, no intermediate temporary is required.
   x = A * x;  // Problematic aliasing of x; intermediate temporary required!
   \endcode

// The detection of aliasing effects, however, takes a small runtime effort. In order to disable
// the aliasing detection, the \c noalias() function can be used:

   \code
   blaze::DynamicVector<double> x, y;
   blaze::DynamicMatrix<double> A;

   x = noalias( x + y );  // No alias detection performed, no intermediate temporary.
   x = noalias( A * x );  // No alias detection performed, no intermediate temporary.
                          // Note that the final result will be incorrect!
   \endcode

// \warning The \c noalias() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Using \c noalias() in a situation
// where an aliasing effect occurs leads to undefined behavior (which can be violated invariants
// or wrong computation results)!
//
// \n \subsection vector_operations_nosimd nosimd()
//
// By default, \b Blaze attempts to vectorize all operations by means of SSE, AVX, etc. in order
// to achieve maximum performance. However, via the \c nosimd() operation it is possible to disable
// the SIMD evaluation of any operation:

   \code
   blaze::DynamicVector<double> x, y;
   blaze::DynamicMatrix<double> A;

   x = nosimd( x + y );  // Disables SIMD for the vector/vector addition
   x = nosimd( A * x );  // Disables SIMD for the matrix/vector multiplication
   \endcode

// Please note that the main purpose of the \c nosimd() operation is to enable an easy performance
// comparison between the vectorized and non-vectorized evaluation. Using the \c nosimd() operation
// will likely result in significantly reduced performance!
//
//
// \n \section vector_operations_modifying_operations Modifying Operations
// <hr>
//
// \subsection vector_operations_resize_reserve .resize() / .reserve()
//
// The size of a \c StaticVector is fixed by the second template parameter and a \c CustomVector
// cannot be resized. In contrast, the size of \c DynamicVectors, \c HybridVectors as well as
// \c CompressedVectors can be changed via the \c resize() function:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   DynamicVector<int,columnVector> v1;
   CompressedVector<int,rowVector> v2( 4 );
   v2[1] = -2;
   v2[3] = 11;

   // Adapting the size of the dynamic and compressed vectors. The (optional) second parameter
   // specifies whether the existing elements should be preserved. Per default, the existing
   // elements are preserved.
   v1.resize( 5UL );         // Resizing vector v1 to 5 elements. Elements of built-in type remain
                             // uninitialized, elements of class type are default constructed.
   v1.resize( 3UL, false );  // Resizing vector v1 to 3 elements. The old elements are lost, the
                             // new elements are NOT initialized!
   v2.resize( 8UL, true );   // Resizing vector v2 to 8 elements. The old elements are preserved.
   v2.resize( 5UL, false );  // Resizing vector v2 to 5 elements. The old elements are lost.
   \endcode

// Note that resizing a vector invalidates all existing views (see e.g. \ref views_subvectors)
// on the vector:

   \code
   blaze::DynamicVector<int,rowVector> v1( 10UL );  // Creating a dynamic vector of size 10
   auto sv = subvector( v1, 2UL, 5UL );             // Creating a view on the range [2..6]
   v1.resize( 6UL );                                // Resizing the vector invalidates the view
   \endcode

// When the internal capacity of a vector is no longer sufficient, the allocation of a larger
// junk of memory is triggered. In order to avoid frequent reallocations, the \c reserve()
// function can be used up front to set the internal capacity:

   \code
   blaze::DynamicVector<int> v1;
   v1.reserve( 100 );
   v1.size();      // Returns 0
   v1.capacity();  // Returns at least 100
   \endcode

// Note that the size of the vector remains unchanged, but only the internal capacity is set
// according to the specified value!
//
// \n \subsection vector_operations_shrinkToFit .shrinkToFit()
//
// The internal capacity of vectors with dynamic memory is preserved in order to minimize the
// number of reallocations. For that reason, the \c resize() and \c reserve() functions can lead
// to memory overhead. The \c shrinkToFit() member function can be used to minimize the internal
// capacity:

   \code
   blaze::DynamicVector<int> v1( 1000UL );  // Create a vector of 1000 integers
   v1.resize( 10UL );                       // Resize to 10, but the capacity is preserved
   v1.shrinkToFit();                        // Remove the unused capacity
   \endcode

// Please note that due to padding the capacity might not be reduced exactly to \c size(). Please
// also note that in case a reallocation occurs, all iterators (including \c end() iterators), all
// pointers and references to elements of the vector are invalidated.
//
// \subsection vector_operations_reset_clear reset() / clear()
//
// In order to reset all elements of a vector, the \c reset() function can be used:

   \code
   // Setup of a single precision column vector, whose elements are initialized with 2.0F.
   blaze::DynamicVector<float> v1( 3UL, 2.0F );

   // Resetting all elements to 0.0F. Only the elements are reset, the size of the vector is unchanged.
   reset( v1 );  // Resetting all elements
   v1.size();    // Returns 3: size and capacity remain unchanged
   \endcode

// In order to return a vector to its default state (i.e. the state of a default constructed
// vector), the \c clear() function can be used:

   \code
   // Setup of a single precision column vector, whose elements are initialized with -1.0F.
   blaze::DynamicVector<float> v1( 5, -1.0F );

   // Resetting the entire vector.
   clear( v1 );  // Resetting the entire vector
   v1.size();    // Returns 0: size is reset, but capacity remains unchanged
   \endcode

// Note that resetting or clearing both dense and sparse vectors does not change the capacity
// of the vectors.
//
//
// \n \subsection vector_operations_swap swap()
//
// Via the \c swap() function it is possible to completely swap the contents of two vectors of
// the same type:

   \code
   blaze::DynamicVector<int,columnVector> v1( 10UL );
   blaze::DynamicVector<int,columnVector> v2( 20UL );

   swap( v1, v2 );  // Swapping the contents of v1 and v2
   \endcode

// \n \section vector_operations_arithmetic_operations Arithmetic Operations
// <hr>
//
// \subsection vector_operations_normalize normalize()
//
// The \c normalize() function can be used to scale any non-zero vector to a length of 1. In
// case the vector does not contain a single non-zero element (i.e. is a zero vector), the
// \c normalize() function returns a zero vector.

   \code
   blaze::DynamicVector<float,columnVector>     v1( 10UL );
   blaze::CompressedVector<double,columnVector> v2( 12UL );

   v1 = normalize( v1 );  // Normalizing the dense vector v1
   length( v1 );          // Returns 1 (or 0 in case of a zero vector)
   v1 = normalize( v2 );  // Assigning v1 the normalized vector v2
   length( v1 );          // Returns 1 (or 0 in case of a zero vector)
   \endcode

// Note that the \c normalize() function only works for floating point vectors. The attempt to
// use it for an integral vector results in a compile time error.
//
//
// \n \subsection vector_operations_min_max min() / max()
//
// The \c min() and \c max() functions can be used for a single vector, multiple vectors, and
// a vector and a scalar.
//
// <b>Single Vector</b>
//
// If passed a single vector, the functions return the smallest and largest element of the given
// dense vector or the smallest and largest non-zero element of the given sparse vector,
// respectively:

   \code
   blaze::StaticVector<int,4UL,rowVector> a{ -5, 2,  7, -4 };

   min( a );  // Returns -5
   max( a );  // Returns 7
   \endcode

   \code
   blaze::CompressedVector<int> b{ 1, 0, 3, 0 };

   min( b );  // Returns 1
   max( b );  // Returns 3
   \endcode

// For more information on the unary \c min() and \c max() reduction operations see the
// \ref vector_operations_reduction_operations section.
//
// <b>Multiple Vectors</b>
//
// If passed two or more dense vectors, the \c min() and \c max() functions compute the
// componentwise minimum or maximum of the given vectors, respectively:

   \code
   blaze::StaticVector<int,4UL,rowVector> c{ -5, 1, -7, 4 };
   blaze::StaticVector<int,4UL,rowVector> d{ -5, 3,  0, 2 };

   min( a, c );     // Results in the vector ( -5, 1, -7, -4 )
   max( a, c, d );  // Results in the vector ( -5, 3,  7,  4 )
   \endcode

// Please note that sparse vectors can only be used in the unary \c min() and \c max() functions.
// Also note that all forms of the \c min() and \c max() functions can be used to compute the
// smallest and largest element of a vector expression:

   \code
   min( a + b + c );  // Returns -9, i.e. the smallest value of the resulting vector
   max( a - b - c );  // Returns 11, i.e. the largest value of the resulting vector

   min( a + c, c - d );  // Results in ( -10 -2 -7 0 )
   max( a - c, c + d );  // Results in ( 0 4 14 6 )
   \endcode

// <b>Vector and Scalar</b>
//
// If passed a dense vector and a scalar, the \c min() and \c max() functions compute the
// componentwise minimum or maximum between the given vector and a uniform vector represented by
// the scalar value:

   \code
   min( a, 0 );  // Results in ( -5, 0, 0, -4 )
   min( 0, a );  // Results in ( -5, 0, 0, -4 )
   max( a, 0 );  // Results in ( 0, 2, 7, 0 )
   max( 0, a );  // Results in ( 0, 2, 7, 0 )
   \endcode

// \n \subsection vector_operators_softmax softmax()
//
// The <a href="https://en.wikipedia.org/wiki/Softmax_function">softmax function</a>, also called
// the normalized exponential function, of a given dense vector can be computed via \c softmax().
// The resulting dense vector consists of real values in the range (0..1], which add up to 1.

   \code
   blaze::StaticVector<double,7UL,rowVector> x{ 1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0 };
   blaze::StaticVector<double,7UL,rowVector> y;

   // Evaluating the softmax function
   y = softmax( x );     // Results in ( 0.024 0.064 0.175 0.475 0.024 0.064 0.175 )
   double s = sum( y );  // Results in 1
   \endcode

// \n \subsection vector_operators_abs abs()
//
// The \c abs() function can be used to compute the absolute values of each element of a vector.
// For instance, the following computation

   \code
   blaze::StaticVector<int,3UL,rowVector> a{ -1, 2, -3 };
   blaze::StaticVector<int,3UL,rowVector> b( abs( a ) );
   \endcode

// results in the vector

                          \f$ b = \left(\begin{array}{*{1}{c}}
                          1 \\
                          2 \\
                          3 \\
                          \end{array}\right)\f$

// \n \subsection vector_operators_sign sign()
//
// The \c sign() function can be used to evaluate the sign of each element of a vector \a a. For
// each element \c i the corresponding result is 1 if \a a[i] is greater than zero, 0 if \a a[i]
// is zero, and -1 if \a a[i] is less than zero. For instance, the following use of the \c sign()
// function

   \code
   blaze::StaticVector<int,3UL,rowVector> a{ -1, 2, 0 };
   blaze::StaticVector<int,3UL,rowVector> b( sign( a ) );
   \endcode

// results in the vector

                          \f$ b = \left(\begin{array}{*{1}{c}}
                          -1 \\
                           1 \\
                           0 \\
                          \end{array}\right)\f$

// \n \subsection vector_operations_rounding_functions floor() / ceil() / trunc() / round()
//
// The \c floor(), \c ceil(), \c trunc(), and \c round() functions can be used to round down/up
// each element of a vector, respectively:

   \code
   blaze::StaticVector<double,3UL,rowVector> a, b;

   b = floor( a );  // Rounding down each element of the vector
   b = ceil ( a );  // Rounding up each element of the vector
   b = trunc( a );  // Truncating each element of the vector
   b = round( a );  // Rounding each element of the vector
   \endcode

// \n \subsection vector_operators_conj conj()
//
// The \c conj() function can be applied on a dense or sparse vector to compute the complex
// conjugate of each element of the vector:

   \code
   using blaze::StaticVector;

   using cplx = std::complex<double>;

   // Creating the vector
   //    ( (-2,-1) )
   //    ( ( 1, 1) )
   StaticVector<cplx,2UL> a{ cplx(-2.0,-1.0), cplx(1.0,1.0) };

   // Computing the vector of complex conjugates
   //    ( (-2, 1) )
   //    ( ( 1,-1) )
   StaticVector<cplx,2UL> b;
   b = conj( a );
   \endcode

// Additionally, vectors can be conjugated in-place via the \c conjugate() function:

   \code
   blaze::DynamicVector<cplx> c( 5UL );

   conjugate( c );  // In-place conjugate operation.
   c = conj( c );   // Same as above
   \endcode

// \n \subsection vector_operators_real real()
//
// The \c real() function can be used on a dense or sparse vector to extract the real part of
// each element of the vector:

   \code
   using blaze::StaticVector;

   using cplx = std::complex<double>;

   // Creating the vector
   //    ( (-2,-1) )
   //    ( ( 1, 1) )
   StaticVector<cplx,2UL> a{ cplx(-2.0,-1.0), cplx(1.0,1.0) };

   // Extracting the real part of each vector element
   //    ( -2 )
   //    (  1 )
   StaticVector<double,2UL> b;
   b = real( a );
   \endcode

// \n \subsection vector_operators_imag imag()
//
// The \c imag() function can be used on a dense or sparse vector to extract the imaginary part
// of each element of the vector:

   \code
   using blaze::StaticVector;

   using cplx = std::complex<double>;

   // Creating the vector
   //    ( (-2,-1) )
   //    ( ( 1, 1) )
   StaticVector<cplx,2UL> a{ cplx(-2.0,-1.0), cplx(1.0,1.0) };

   // Extracting the imaginary part of each vector element
   //    ( -1 )
   //    (  1 )
   StaticVector<double,2UL> b;
   b = imag( a );
   \endcode

// \n \subsection vector_operators_arg arg()
//
// The \c arg() function can be used on a dense or sparse vector to compute the phase angle for
// each element of the vector:

   \code
   using blaze::StaticVector;

   using cplx = std::complex<double>;

   // Creating the vector
   //    ( (-2,-1) )
   //    ( ( 1, 1) )
   StaticVector<cplx,2UL> a{ cplx(-2.0,-1.0), cplx(1.0,1.0) };

   // Compute the phase angle of each vector element
   //    ( -2.67795  )
   //    (  0.785398 )
   StaticVector<double,2UL> b;
   b = arg( a );
   \endcode

// \n \subsection vector_operations_sqrt sqrt() / invsqrt()
//
// Via the \c sqrt() and \c invsqrt() functions the (inverse) square root of each element of a
// vector can be computed:

   \code
   blaze::DynamicVector<double> a, b, c;

   b = sqrt( a );     // Computes the square root of each element
   c = invsqrt( a );  // Computes the inverse square root of each element
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_cbrt cbrt() / invcbrt()
//
// The \c cbrt() and \c invcbrt() functions can be used to compute the the (inverse) cubic root
// of each element of a vector:

   \code
   blaze::HybridVector<double,3UL> a, b, c;

   b = cbrt( a );     // Computes the cubic root of each element
   c = invcbrt( a );  // Computes the inverse cubic root of each element
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_hypot hypot()
//
// The \c hypot() function can be used to compute the componentwise hypotenous for a pair of
// dense vectors:

   \code
   blaze::StaticVector<double,3UL> a, b, c;

   c = hypot( a, b );  // Computes the componentwise hypotenuous
   \endcode

// \n \subsection vector_operations_clamp clamp()
//
// The \c clamp() function can be used to restrict all elements of a vector to a specific range:

   \code
   blaze::DynamicVector<double> a, b

   b = clamp( a, -1.0, 1.0 );  // Restrict all elements to the range [-1..1]
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_pow pow()
//
// The \c pow() function can be used to compute the exponential value of each element of a vector.
// If passed a vector and a numeric exponent, the function computes the exponential value of each
// element of the vector using the same exponent. If passed a second vector, the function computes
// the componentwise exponential value:

   \code
   blaze::StaticVector<double,3UL> a, b, c;

   c = pow( a, 1.2 );  // Computes the exponential value of each element
   c = pow( a, b );    // Computes the componentwise exponential value
   \endcode

// \n \subsection vector_operations_exp exp() / exp2() / exp10()
//
// \c exp(), \c exp2() and \c exp10() compute the base e/2/10 exponential of each element of a
// vector, respectively:

   \code
   blaze::DynamicVector<double> a, b;

   b = exp( a );    // Computes the base e exponential of each element
   b = exp2( a );   // Computes the base 2 exponential of each element
   b = exp10( a );  // Computes the base 10 exponential of each element
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_log log() / log2() / log10() / log1p() / lgamma()
//
// The \c log(), \c log2(), \c log10(), \c log1p() and \c lgamma() functions can be used to
// compute the natural, binary and common logarithm of each element of a vector:

   \code
   blaze::StaticVector<double,3UL> a, b;

   b = log( a );     // Computes the natural logarithm of each element
   b = log2( a );    // Computes the binary logarithm of each element
   b = log10( a );   // Computes the common logarithm of each element
   b = log1p( a );   // Computes the natural logarithm of x+1 of each element
   b = lgamma( a );  // Computes the natural logarithm of the absolute value of the gamma function
   \endcode

// \n \subsection vector_operations_trigonometric_functions sin() / cos() / tan() / asin() / acos() / atan()
//
// The following trigonometric functions are available for both dense and sparse vectors:

   \code
   blaze::DynamicVector<double> a, b;

   b = sin( a );  // Computes the sine of each element of the vector
   b = cos( a );  // Computes the cosine of each element of the vector
   b = tan( a );  // Computes the tangent of each element of the vector

   b = asin( a );  // Computes the inverse sine of each element of the vector
   b = acos( a );  // Computes the inverse cosine of each element of the vector
   b = atan( a );  // Computes the inverse tangent of each element of the vector
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_hyperbolic_functions sinh() / cosh() / tanh() / asinh() / acosh() / atanh()
//
// The following hyperbolic functions are available for both dense and sparse vectors:

   \code
   blaze::DynamicVector<double> a, b;

   b = sinh( a );  // Computes the hyperbolic sine of each element of the vector
   b = cosh( a );  // Computes the hyperbolic cosine of each element of the vector
   b = tanh( a );  // Computes the hyperbolic tangent of each element of the vector

   b = asinh( a );  // Computes the inverse hyperbolic sine of each element of the vector
   b = acosh( a );  // Computes the inverse hyperbolic cosine of each element of the vector
   b = atanh( a );  // Computes the inverse hyperbolic tangent of each element of the vector
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_atan2 atan2()
//
// The multi-valued inverse tangent is available for a pair of dense vectors:

   \code
   blaze::DynamicVector<double> a, b, c;

   c = atan2( a, b );  // Computes the componentwise multi-valued inverse tangent
   \endcode

// \n \subsection vector_operations_erf erf() / erfc()
//
// The \c erf() and \c erfc() functions compute the (complementary) error function of each
// element of a vector:

   \code
   blaze::StaticVector<double,3UL,rowVector> a, b;

   b = erf( a );   // Computes the error function of each element
   b = erfc( a );  // Computes the complementary error function of each element
   \endcode

// Note that in case of sparse vectors only the non-zero elements are taken into account!
//
//
// \n \subsection vector_operations_map map() / forEach()
//
// Via the \c map() functions it is possible to execute componentwise custom operations on vectors.
// The unary \c map() function can be used to apply a custom operation on each element of a dense
// or sparse vector. For instance, the following example demonstrates a custom square root
// computation via a lambda:

   \code
   blaze::DynamicVector<double> a, b;

   b = map( a, []( double d ) { return std::sqrt( d ); } );
   \endcode

// The N-ary \c map() functions can be used to apply an operation componentwise to the elements
// of N dense vectors (where \f$ N <= 6 \f$). The following example demonstrates the merging of
// two column vectors of double precision values into a vector of double precision complex numbers:

   \code
   blaze::DynamicVector<double> real{ 2.1, -4.2,  1.0,  0.6 };
   blaze::DynamicVector<double> imag{ 0.3,  1.4,  2.9, -3.4 };

   blaze::DynamicVector< complex<double> > cplx;

   // Creating the vector
   //    ( ( 2.1,  0.3) )
   //    ( (-4.2,  1.4) )
   //    ( ( 1.0,  2.9) )
   //    ( ( 0.6, -3.4) )
   cplx = map( real, imag, []( double r, double i ){ return complex<double>( r, i ); } );
   \endcode

// Applying the map() function to a column vector and a row vector results in the outer map of
// the two vectors. The following example demonstrates the outer sum of a column vector and a
// row vector:

   \code
   blaze::DynamicVector<int,columnVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   // Results in the matrix
   //
   //       (  1  5  0  6 )
   //   A = (  4  8  3  9 )
   //       ( -2  2 -3  3 )
   //
   blaze::StaticMatrix<int,3UL,4UL> M1 = map( v1, v2, []( int a, int b ){ return a + b; } );
   \endcode

// Although the computation in the two previous examples can be parallelized it is not vectorized
// and thus cannot perform at peak performance. However, it is also possible to create vectorized
// custom operations. See \ref custom_operations for a detailed overview of the possibilities of
// custom operations.
//
// Please note that unary custom operations on vectors have been introduced in \b Blaze 3.0 in
// form of the \c forEach() function. With the introduction of binary custom functions, the
// \c forEach() function has been renamed to \c map(). The \c forEach() function can still be
// used, but the function might be deprecated in future releases of \b Blaze.
//
//
// \n \subsection vector_operations_select select()
//
// The \c select() function performs a componentwise, conditional selection of elements. Given
// the three dense vectors \c cond, \c a, and \c b, in case an element in the \c cond vector
// evaluates to \c true, the according element of \a a is selected, in case the \a cond element
// evaluates to \c false, the according element of \a b is selected. The following example
// demonstrates the use of the \a select() function:

   \code
   blaze::DynamicVector<bool> cond{ true, false, true false };
   blaze::DynamicVector<int> a{ 1, -1, 1, -1 };
   blaze::DynamicVector<int> b{ -2, 2, -2, 2 };
   blaze::DynamicVector<int> c;
   // ... Resizing and initialization

   c = select( cond, a, b );  // Results in ( 1, 2, 1, 2 )
   \endcode

// \n \section vector_operations_reduction_operations Reduction Operations
// <hr>
//
// \subsection vector_operations_reduction_operations_reduce reduce()
//
// The \c reduce() function performs a total reduction of the elements of the given dense vector
// or the non-zero elements of the given sparse vector. The following examples demonstrate the
// total reduction of a dense and sparse vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization

   const double totalsum1 = reduce( a, blaze::Add() );
   const double totalsum2 = reduce( a, []( double a, double b ){ return a + b; } );
   \endcode

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization

   const double totalmin1 = reduce( a, blaze::Min() );
   const double totalmin2 = reduce( a, []( double a, double b ){ return blaze::min( a, b ); } );
   \endcode

// As demonstrated in the examples it is possible to pass any binary callable as custom reduction
// operation. However, for instance in the case of lambdas the vectorization of the reduction
// operation is compiler dependent and might not perform at peak performance. However, it is also
// possible to create vectorized custom operations. See \ref custom_operations for a detailed
// overview of the possibilities of custom operations.
//
// Please note that the evaluation order of the \c reduce() function is unspecified. Thus the
// behavior is non-deterministic if the given reduction operation is not associative or not
// commutative. Also, the operation is undefined if the given reduction operation modifies the
// values.
//
// \n \subsection vector_operations_reduction_operations_sum sum()
//
// The \c sum() function reduces the elements of the given dense vector or the non-zero elements
// of the given sparse vector by means of addition:

   \code
   blaze::DynamicVector<int> a{ 1, 2, 3, 4 };

   const int totalsum = sum( a );  // Results in 10
   \endcode

   \code
   blaze::CompressedVector<int> a{ 1, 2, 3, 4 };

   const int totalsum = sum( a );  // Results in 10
   \endcode

// Please note that the evaluation order of the \c sum() function is unspecified.
//
// \n \subsection vector_operations_reduction_operations_prod prod()
//
// The \c prod() function reduces the elements of the given dense vector or the non-zero elements
// of the given sparse vector by means of multiplication:

   \code
   blaze::DynamicVector<int> a{ 1, 2, 3, 4 };

   const int totalprod = prod( a );  // Results in 24
   \endcode

   \code
   blaze::CompressedVector<int> a{ 1, 2, 3, 4 };

   const int totalprod = prod( a );  // Results in 24
   \endcode

// \n \subsection vector_operations_reduction_operations_min min()
//
// The unary \c min() function returns the smallest element of the given dense vector or the
// smallest non-zero element of the given sparse vector. It can only be used for element types
// that support the smaller-than relationship. In case the given vector currently has a size
// of 0, the returned value is the default value (e.g. 0 in case of fundamental data types).

   \code
   blaze::DynamicVector<int> a{ 1, -2, 3, 0 };

   const int totalmin = min( a );  // Results in -2
   \endcode

   \code
   blaze::CompressedVector<int> a{ 1, 0, 3, 0 };

   const int totalmin = min( a );  // Results in 1
   \endcode

// \note In case the sparse vector is not completely filled, the implicit zero elements are NOT
// taken into account. In the previous example the compressed vector has only 2 non-zero elements.
// However, the minimum of the vector is 1.
//
// \n \subsection vector_operations_reduction_operations_max max()
//
// The unary \c max() function returns the largest element of the given dense vector or the
// largest non-zero element of the given sparse vector. It can only be used for element types
// that support the smaller-than relationship. In case the given vector currently has a size
// of 0, the returned value is the default value (e.g. 0 in case of fundamental data types).

   \code
   blaze::DynamicVector<int> a{ 1, -2, 3, 0 };

   const int totalmax = max( a );  // Results in 3
   \endcode

   \code
   blaze::CompressedVector<int> a{ -1, 0, -3, 0 };

   const int totalmin = max( a );  // Results in -1
   \endcode

// \note In case the sparse vector is not completely filled, the implicit zero elements are NOT
// taken into account. In the previous example the compressed vector has only 2 non-zero elements.
// However, the maximum of the vector is -1.
//
// \n \subsection vector_operations_reduction_operations_argmin argmin()
//
// The \c argmin() function returns the index of the first smallest element of the given dense
// vector. This function can only be used for element types that support the smaller-than
// relationship. In case the given vector currently has a size of 0, the returned index is 0.

   \code
   blaze::DynamicVector<int> a{ 1, -2, 3, 0 };
   const size_t minindex = argmin( a );  // Results in 1
   \endcode

// \n \subsection vector_operations_reduction_operations_argmax argmax()
//
// The \c argmax() function returns the index of the first largest element of the given dense
// vector. This function can only be used for element types that support the smaller-than
// relationship. In case the given vector currently has a size of 0, the returned index is 0.

   \code
   blaze::DynamicVector<int> a{ 1, -2, 3, 0 };
   const size_t maxindex = argmax( a );  // Results in 2
   \endcode

// \n \section vector_operations_norms Norms
// <hr>
//
// \subsection vector_operations_norms_norm norm()
//
// The \c norm() function computes the L2 norm of the given dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = norm( a );
   const double norm2 = norm( b );
   \endcode

// \n \subsection vector_operations_norms_sqrnorm sqrNorm()
//
// The \c sqrNorm() function computes the squared L2 norm of the given dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = sqrNorm( a );
   const double norm2 = sqrNorm( b );
   \endcode

// \n \subsection vector_operations_norms_l1norm l1Norm()
//
// The \c l1Norm() function computes the squared L1 norm of the given dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = l1Norm( a );
   const double norm2 = l1Norm( b );
   \endcode

// \n \subsection vector_operations_norms_l2norm l2Norm()
//
// The \c l2Norm() function computes the squared L2 norm of the given dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = l2Norm( a );
   const double norm2 = l2Norm( b );
   \endcode

// \n \subsection vector_operations_norms_l3norm l3Norm()
//
// The \c l3Norm() function computes the squared L3 norm of the given dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = l3Norm( a );
   const double norm2 = l3Norm( b );
   \endcode

// \n \subsection vector_operations_norms_l4norm l4Norm()
//
// The \c l4Norm() function computes the squared L4 norm of the given dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = l4Norm( a );
   const double norm2 = l4Norm( b );
   \endcode

// \n \subsection vector_operations_norms_lpnorm lpNorm()
//
// The \c lpNorm() function computes the general Lp norm of the given dense or sparse vector,
// where the norm is specified by either a compile time or a runtime argument:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = lpNorm<2>( a );    // Compile time argument
   const double norm2 = lpNorm( b, 2.3 );  // Runtime argument
   \endcode

// \n \subsection vector_operations_norms_maxnorm linfNorm() / maxNorm()
//
// The \c linfNorm() and \c maxNorm() functions compute the infinity/maximum norm of the given
// dense or sparse vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   const double norm1 = linfNorm( a );
   const double norm2 = maxNorm( b );
   \endcode

// \n \section vector_operations_scalar_expansion Scalar Expansion
// <hr>
//
// By means of the \c uniform() function it is possible to expand a scalar value into a dense,
// uniform vector. By default, the resulting uniform vector is a column vector, but it is possible
// to specify the transpose flag explicitly:

   \code
   using blaze::columnVector;

   int scalar = 5;

   blaze::DynamicVector<int,columnVector> v;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3-dimensional uniform column vector
   //
   //    ( 5 )
   //    ( 5 )
   //    ( 5 )
   //
   v = uniform( 3UL, scalar );
   v = uniform<columnVector>( 3UL, scalar );
   \endcode

// \n \section vector_operations_vector_expansion Vector Expansion
// <hr>
//
// Via the \c expand() function it is possible to convert a dense or sparse vector into a matrix.
// A column vector is expanded into a column-major matrix, a row vector is expanded into a
// row-major matrix. As demonstrated by the following examples, \c expand() can be used with both
// runtime and compile time parameters:

   \code
   blaze::DynamicVector<int,columnVector> a{ 1, 2, 3 };
   blaze::CompressedVector<int,rowVector> b{ 1, 0, 3, 0, 5 };

   // Expand the dense column vector ( 1 2 3 ) into a dense 3x5 column-major matrix
   //
   //   ( 1 1 1 1 1 )
   //   ( 2 2 2 2 2 )
   //   ( 3 3 3 3 3 )
   //
   expand( a, 5 );  // Runtime parameter
   expand<5>( a );  // Compile time parameter

   // Expand the sparse row vector ( 1 0 3 0 5 ) into a sparse 3x5 row-major matrix
   //
   //   ( 1 0 3 0 5 )
   //   ( 1 0 3 0 5 )
   //   ( 1 0 3 0 5 )
   //
   expand( b, 3 );  // Runtime parameter
   expand<3>( b );  // Compile time parameter
   \endcode

// \n \section vector_operations_vector_repetition Vector Repetition
// <hr>
//
// Via the \c repeat() function it is possible to repeat a dense or sparse vector multiple times
// to represent a larger vector. Repeating a column vector results in a column vector, repeating
// a row vector results in a row vector. As demonstrated by the following examples, \c repeat()
// can be used with both runtime and compile time parameters:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<int,columnVector> a1{ 1, 0, -2 };
   blaze::CompressedVector<int,rowVector> b1{ 0, -1, 7 };

   blaze::DynamicVector<int,columnVector> a2;
   blaze::CompressedVector<int,rowVector> b2;

   // ... Resizing and initialization

   // Repeating the dense column vector ( 1  0 -2 ) three times results in
   //
   //   ( 1  0 -2  1  0 -2  1  0 -2 )
   //
   a2 = repeat( a1, 3UL );
   a2 = repeat<3UL>( a1 );

   // Repeating the sparse row vector ( 0 -1  7 ) three times results in
   //
   //   ( 0 -1  7  0 -1  7  0 -1  7 )
   //
   b2 = repeat( b1, 3UL );
   b2 = repeat<3UL>( b1 );
   \endcode

// \n \section vector_operations_statistic_operations Statistic Operations
// <hr>
//
// \subsection vector_operations_mean mean()
//
// The <a href="https://en.wikipedia.org/wiki/Arithmetic_mean">(arithmetic) mean</a> of a dense or
// sparse vector can be computed via the \c mean() function. In case of a sparse vector, both the
// non-zero and zero elements are taken into account. The following example demonstrates the
// computation of the mean of a dense vector:

   \code
   blaze::DynamicVector<int> v{ 1, 4, 3, 6, 7 };

   const double m = mean( v );  // Results in 4.2 (i.e. 21/5)
   \endcode

// In case the size of the given vector is 0, a \c std::invalid_argument is thrown.
//
// \n \subsection vector_operations_var var()
//
// The <a href="https://en.wikipedia.org/wiki/Variance">variance</a> of a dense or sparse vector
// can be computed via the \c var() function. In case of a sparse vector, both the non-zero and
// zero elements are taken into account. The following example demonstrates the computation of
// the variance of a dense vector:

   \code
   blaze::DynamicVector<int> v{ 1, 4, 3, 6, 7 };

   const double v = var( v );  // Results in 5.7
   \endcode

// In case the size of the given vector is smaller than 2, a \c std::invalid_argument is thrown.
//
// \n \subsection vector_operations_stddev stddev()
//
// The <a href="https://en.wikipedia.org/wiki/Standard_deviation">standard deviation</a> of a
// dense or sparse vector can be computed via the \c stddev() function. In case of a sparse
// vector, both the non-zero and zero elements are taken into account. The following example
// demonstrates the computation of the standard deviation of a dense vector:

   \code
   blaze::DynamicVector<int> v{ 1, 4, 3, 6, 7 };

   const double s = stddev( v );  // Results in 2.38747
   \endcode

// In case the size of the given vector is smaller than 2, a \c std::invalid_argument is thrown.
//
//
// \n \section vector_operations_declaration_operations Declaration Operations
// <hr>
//
// \subsection vector_operations_declzero declzero()
//
// The \c declzero() operation can be used to explicitly declare any vector or vector expression
// as zero vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization

   b = declzero( a );
   \endcode

// Any vector or vector expression that has been declared as zero vector via \c declzero() will
// gain all the benefits of a zero vector, which range from reduced runtime checking to a
// considerable speed-up in computations:

   \code
   using blaze::DynamicVector;

   DynamicVector<double> a, b, c;
   // ... Resizing and initialization

   isZero( declzero( a ) );  // Will always return true without runtime effort

   c = declzero( a ) + b;  // Declare the left operand of the vector addition as a
                           // zero vector, i.e. no addition needs to be performed
   \endcode

// \warning The \c declzero() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-zero vector or
// vector expression as zero vector via the \c declzero() operation leads to undefined behavior
// (which can be violated invariants or wrong computation results)!
//
//
// \n \section vector_operations_vector_generators Vector Generators
// <hr>
//
// \subsection vector_operations_generate generate()
//
// The \c generate() function returns a dense vector filled elementwise via the given custom
// operation. By default, the returned vector is a column vector, but this setting can be changed
// via the \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch (see \ref transpose_flag). Alternatively it is
// possible to specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \c generate() function:

   \code
   using blaze::generate;
   using blaze::columnVector;
   using blaze::rowVector;

   // Generates the homogeneous integer vector ( 2, 2, 2, 2, 2 )
   blaze::DynamicVector<int,columnVector> a;
   a = generate( 5UL, []( size_t index ){ return 2; } );

   // Generates the linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   blaze::DynamicVector<float,columnVector> b;
   b = generate( 4UL, []( size_t index ){ return 2.1F + 1.1F*index; } );

   // Generates the logarithmically spaced double vector ( 1.0, 10.0, 100.0, 1000.0 )
   blaze::DynamicVector<double,columnVector> c;
   c = generate<columnVector>( 4UL, []( size_t index ){ return blaze::exp10( 1.0 + 1.0*index ); } );

   // Generates the vector of integer vectors ( ( 1, 2 ), ( 2, 3 ), ( 3, 4 ), ( 4, 5 ) )
   using VT = blaze::StaticVector<int,2UL>;
   blaze::StaticVector<VT,4UL,rowVector> d;
   d = generate<rowVector>( []( size_t index ) { return evaluate( VT{ 1, 2 } + index ); } );
   \endcode

// \n \subsection vector_operations_linspace linspace()
//
// The \c linspace() function returns a dense vector filled with linearly spaced elements. By
// default, the returned vector is a column vector, but this setting can be changed via the
// \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch (see \ref transpose_flag). Alternatively it is possible
// to specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \c linspace() function:

   \code
   using blaze::linspace;
   using blaze::columnVector;
   using blaze::rowVector;

   // Generates the linearly spaced integer vector ( 2, 3, 4, 5, 6 )
   blaze::DynamicVector<int,columnVector> a;
   a = linspace( 5UL, 2, 6 );

   // Generates the linearly spaced integer vector ( 6, 5, 4, 3, 2 )
   blaze::DynamicVector<int,columnVector> b;
   b = linspace<columnVector>( 5UL, 6, 2 );

   // Generates the linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   blaze::DynamicVector<float,rowVector> c;
   c = linspace<rowVector>( 4UL, 2.1F, 5.4F );
   \endcode

// \n \subsection vector_operations_logspace logspace()
//
// The \c logspace() function returns a dense vector filled with logarithmically spaced elements.
// By default, the returned vector is a column vector, but this setting can be changed via the
// \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch (see \ref transpose_flag). Alternatively it is possible
// to specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \c logspace() function:

   \code
   using blaze::logspace;
   using blaze::columnVector;
   using blaze::rowVector;

   // Generates the logarithmically spaced double vector ( 1, 10, 100, 1000 )
   blaze::DynamicVector<int,columnVector> a;
   a = logspace( 4UL, 0, 3 );

   // Generates the logarithmically spaced double vector ( 1000.0, 100.0, 10.0, 1.0 )
   blaze::DynamicVector<double,rowVector> b;
   b = logspace<rowVector>( 4UL, 3.0, 0.0 );
   \endcode

// \n \subsection vector_operations_uniform uniform()
//
// The \c uniform() function creates a uniform vector of the given size. By default, the
// resulting uniform vector is a column vector, but this setting can be changed via the
// \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch (see \ref transpose_flag). Alternatively it is
// possible to specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \c uniform() function:

   \code
   using blaze::uniform;
   using blaze::columnVector;
   using blaze::rowVector;

   // Creates the uniform column vector ( 1, 1, 1, 1, 1 )
   auto u1 = uniform( 5UL, 1 );

   // Creates the uniform column vector ( 1.2, 1.2, 1.2 )
   auto u2 = uniform<columnVector>( 3UL, 1.2 );

   // Creates the uniform row vector ( 5U, 5U, 5U, 5U )
   auto u3 = uniform<rowVector>( 4UL, 5U );
   \endcode

// \n \subsection vector_operations_zero zero()
//
// The \c zero() function creates a zero vector of the given element type and size. By default,
// the resulting zero vector is a column vector, but this setting can be changed via the
// \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch (see \ref transpose_flag). Alternatively it is
// possible to specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \c zero() function:

   \code
   using blaze::zero;
   using blaze::columnVector;
   using blaze::rowVector;

   // Creates the zero column vector ( 0, 0, 0, 0, 0 )
   auto z1 = zero<int>( 5UL );

   // Creates the zero column vector ( 0.0, 0.0, 0.0 )
   auto z2 = zero<double,columnVector>( 3UL );

   // Creates the zero row vector ( 0U, 0U, 0U, 0U )
   auto z3 = zero<unsigned int,rowVector>( 4UL );
   \endcode

// \n Previous: \ref vector_types &nbsp; &nbsp; Next: \ref matrices
*/
//*************************************************************************************************


//**Matrices***************************************************************************************
/*!\page matrices Matrices
//
// \tableofcontents
//
//
// \n \section matrices_general General Concepts
// <hr>
//
// The \b Blaze library currently offers five dense matrix types (\ref matrix_types_static_matrix,
// \ref matrix_types_dynamic_matrix, \ref matrix_types_hybrid_matrix, \ref matrix_types_custom_matrix,
// and \ref matrix_types_uniform_matrix) and three sparse matrix types (\ref matrix_types_compressed_matrix,
// \ref matrix_types_identity_matrix, and \ref matrix_types_zero_matrix). All matrices can either
// be stored as row-major matrices or column-major matrices:

   \code
   using blaze::DynamicMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Setup of the 2x3 row-major dense matrix
   //
   //    ( 1  2  3 )
   //    ( 4  5  6 )
   //
   DynamicMatrix<int,rowMajor> A{ { 1, 2, 3 },
                                  { 4, 5, 6 } };

   // Setup of the 3x2 column-major dense matrix
   //
   //    ( 1  4 )
   //    ( 2  5 )
   //    ( 3  6 )
   //
   DynamicMatrix<int,columnMajor> B{ { 1, 4 },
                                     { 2, 5 },
                                     { 3, 6 } };
   \endcode

// Per default, all matrices in \b Blaze are row-major matrices:

   \code
   // Instantiation of a 3x3 row-major matrix
   blaze::DynamicMatrix<int> C( 3UL, 3UL );
   \endcode

// \n \section matrices_details Matrix Details
// <hr>
//
//  - \ref matrix_types
//  - \ref matrix_operations
//
//
// \n \section matrices_examples Examples
// <hr>

   \code
   using blaze::StaticMatrix;
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   StaticMatrix<double,6UL,20UL> A;      // Instantiation of a 6x20 row-major static matrix
   CompressedMatrix<double,rowMajor> B;  // Instantiation of a row-major compressed matrix
   DynamicMatrix<double,columnMajor> C;  // Instantiation of a column-major dynamic matrix

   // ... Resizing and initialization

   C = A * B;
   \endcode

// \n Previous: \ref vector_operations &nbsp; &nbsp; Next: \ref matrix_types
*/
//*************************************************************************************************


//**Matrix Types***********************************************************************************
/*!\page matrix_types Matrix Types
//
// \tableofcontents
//
//
// \n \section matrix_types_dense_matrices Dense Matrices
// <hr>
//
// \subsection matrix_types_static_matrix StaticMatrix
//
// The blaze::StaticMatrix class template is the representation of a fixed size matrix with
// statically allocated elements of arbitrary type. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/StaticMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the number of rows and columns, the storage order of the matrix,
// the alignment, the padding, and the group tag of the matrix can be specified via the seven
// template parameters:

   \code
   namespace blaze {

   template< typename Type, size_t M, size_t N, bool SO, AlignmentFlag AF, PaddingFlag PF, typename Tag >
   class StaticMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. StaticMatrix can be used with any
//             non-cv-qualified, non-reference element type.
//  - \c M   : specifies the total number of rows of the matrix.
//  - \c N   : specifies the total number of columns of the matrix. Note that it is expected
//             that StaticMatrix is only used for tiny and small matrices.
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c AF  : specifies whether the first element of every row/column is properly aligned with
//             respect to the available instruction set (SSE, AVX, ...). Possible values are
//             \c blaze::aligned and \c blaze::unaligned. The default value is
//             \c blaze::defaultAlignmentFlag.
//  - \c PF  : specifies whether every row/column of the matrix should be padded to maximize the
//             efficiency of vectorized operations. Possible values are \c blaze::padded and
//             \c blaze::unpadded. The default value is \c blaze::defaultPaddingFlag.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::StaticMatrix is perfectly suited for small to medium matrices whose dimensions are
// known at compile time:

   \code
   // Definition of a 3x4 integral row-major matrix
   blaze::StaticMatrix<int,3UL,4UL> A;

   // Definition of a 4x6 single precision row-major matrix
   blaze::StaticMatrix<float,4UL,6UL,blaze::rowMajor> B;

   // Definition of an unaligned, unpadded 6x4 double precision column-major matrix
   blaze::StaticMatrix<double,6UL,4UL,blaze::columnMajor,blaze::unaligned,blaze::unpadded> C;
   \endcode

// \subsubsection matrix_types_static_matrix_alignment Alignment
//
// In case \c AF is set to \c blaze::aligned, the elements of a blaze::StaticMatrix are possibly
// over-aligned to meet the alignment requirements of the available instruction set (SSE, AVX,
// AVX-512, ...). The alignment for fundamental types (\c short, \c int, \c float, \c double, ...)
// and complex types (\c complex<float>, \c complex<double>, ...) is 16 bytes for SSE, 32 bytes
// for AVX, and 64 bytes for AVX-512. All other types are aligned according to their intrinsic
// alignment:

   \code
   struct Int { int i; };

   using MT1 = blaze::StaticMatrix<double,3UL,5UL>;
   using MT2 = blaze::StaticMatrix<complex<float>,2UL,3UL>;
   using MT3 = blaze::StaticMatrix<Int,5UL,4UL>;

   alignof( MT1 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( MT2 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( MT3 );  // Evaluates to 'alignof( Int )'
   \endcode

// Note that an aligned blaze::StaticMatrix instance may be bigger than the sum of its data
// elements:

   \code
   sizeof( MT1 );  // Evaluates to 160 for SSE, and 192 for AVX and AVX-512
   sizeof( MT2 );  // Evaluates to 64 for SSE and AVX and 128 for AVX-512
   sizeof( MT3 );  // Evaluates to 80; no special alignment requirements
   \endcode

// Please note that for this reason a blaze::StaticMatrix cannot be used in containers using
// dynamic memory such as \c std::vector without additionally providing an allocator that can
// provide over-aligned memory:

   \code
   using Type = blaze::StaticMatrix<double,3UL,5UL>;
   using Allocator = blaze::AlignedAllocator<Type>;

   std::vector<Type> v1;  // Might be misaligned for AVX or AVX-512
   std::vector<Type,Allocator> v2;  // Properly aligned for AVX or AVX-512
   \endcode

// \subsubsection matrix_types_static_matrix_padding Padding
//
// Adding padding elements to the end of every row or column of a blaze::StaticMatrix can have a
// significant impact on the performance. For instance, assuming that AVX is available, then two
// padded 3x3 matrices of double precision values can be added with three SIMD addition operations:

   \code
   using blaze::StaticMatrix;
   using blaze::rowMajor;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   StaticMatrix<double,3UL,3UL,rowMajor,aligned,padded> A1, B1, C1;
   StaticMatrix<double,3UL,3UL,rowMajor,unaligned,unpadded> A2, B2, C2;

   // ... Initialization

   C1 = A1 + B1;  // AVX-based matrix addition; maximum performance
   C2 = A2 + B2;  // Scalar matrix addition; limited performance

   sizeof( A1 );  // Evaluates to 96 for SSE and AVX, and 192 for AVX-512
   sizeof( A2 );  // Evaluates to 72 for SSE, AVX, and AVX-512 (minimum size)
   \endcode

// Due to padding, the first addition will run at maximum performance. On the flip side, the size
// of each matrix instance is increased due to the padding elements. The total size of an instance
// depends on the number of elements and width of the available instruction set (16 bytes for
// SSE, 32 bytes for AVX, and 64 bytes for AVX-512).
//
// The second addition will be limited in performance since due to the number of elements some of
// the elements need to be handled in a scalar operation. However, the size of an \c unaligned,
// \c unpadded blaze::StaticMatrix instance is guaranteed to be the sum of its elements.
//
// Please also note that \b Blaze will zero initialize the padding elements in order to achieve
// maximum performance!
//
//
// \n \subsection matrix_types_dynamic_matrix DynamicMatrix
//
// The blaze::DynamicMatrix class template is the representation of an arbitrary sized matrix
// with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. It can be included
// via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/DynamicMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the storage order, the type of the allocator, and the group tag of
// the matrix can be specified via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Alloc, typename Tag >
   class DynamicMatrix;

   } // namespace blaze
   \endcode

//  - \c Type : specifies the type of the matrix elements. DynamicMatrix can be used with any
//              non-cv-qualified, non-reference element type.
//  - \c SO   : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//              matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c Alloc: specifies the type of allocator used to allocate dynamic memory. The default type
//              of allocator is \c blaze::AlignedAllocator.
//  - \c Tag  : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//              See \ref grouping_tagging for details.
//
// The blaze::DynamicMatrix is the default choice for all kinds of dense matrices and the best
// choice for medium to large matrices. The number of rows and columns can be modified at runtime:

   \code
   // Definition of a 3x4 integral row-major matrix
   blaze::DynamicMatrix<int> A( 3UL, 4UL );

   // Definition of a 4x6 single precision row-major matrix
   blaze::DynamicMatrix<float,blaze::rowMajor> B( 4UL, 6UL );

   // Definition of a double precision column-major matrix with 0 rows and columns
   blaze::DynamicMatrix<double,blaze::columnMajor> C;
   \endcode

// \subsubsection matrix_types_dynamic_matrix_allocators Allocators
//
// Via the third template parameter it is possible to customize the memory allocation of a
// \c blaze::DynamicMatrix. The provided allocator is expected to represent an implementation of
// the allocator concept of the standard library (see for instance
// <a href="https://en.cppreference.com/w/cpp/container/vector">std::vector</a> and
// <a href="https://en.cppreference.com/w/cpp/memory/allocator">std::allocator</a>). In
// addition, the provided allocator is also required to provide properly (over-)aligned memory
// for fundamental and complex numbers. For instance, in case SSE vectorization is possible, the
// returned memory must be at least 16-byte aligned. In case AVX is active, the memory must be at
// least 32-byte aligned, and in case of AVX-512 the memory must be even 64-byte aligned.
//
//
// \n \subsection matrix_types_hybrid_matrix HybridMatrix
//
// The HybridMatrix class template combines the flexibility of a dynamically sized matrix with
// the efficiency and performance of a fixed size matrix. It is implemented as a crossing between
// the blaze::StaticMatrix and the blaze::DynamicMatrix class templates: Similar to the static
// matrix it uses static stack memory instead of dynamically allocated memory and similar to the
// dynamic matrix it can be resized (within the extend of the static memory). It can be included
// via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/HybridMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the maximum number of rows and columns, the storage order of the
// matrix, the alignment, the padding, and the group tag of the matrix can be specified via the
// seven template parameters:

   \code
   namespace blaze {

   template< typename Type, size_t M, size_t N, bool SO, AlignmentFlag AF, PaddingFlag PF, typename Tag >
   class HybridMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. HybridMatrix can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c M   : specifies the maximum number of rows of the matrix.
//  - \c N   : specifies the maximum number of columns of the matrix. Note that it is expected
//             that HybridMatrix is only used for tiny and small matrices.
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c AF  : specifies whether the first element of every row/column is properly aligned with
//             respect to the available instruction set (SSE, AVX, ...). Possible values are
//             \c blaze::aligned and \c blaze::unaligned. The default value is
//             \c blaze::defaultAlignmentFlag.
//  - \c PF  : specifies whether every row/column of the matrix should be padded to maximize the
//             efficiency of vectorized operations. Possible values are \c blaze::padded and
//             \c blaze::unpadded. The default value is \c blaze::defaultPaddingFlag.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::HybridMatrix is a suitable choice for small to medium matrices, whose dimensions
// are not known at compile time or not fixed at runtime, but whose maximum dimensions are known
// at compile time:

   \code
   // Definition of a 3x4 integral row-major matrix with maximum dimensions of 6x8
   blaze::HybridMatrix<int,6UL,8UL> A( 3UL, 4UL );

   // Definition of a 4x6 single precision row-major matrix with maximum dimensions of 12x16
   blaze::HybridMatrix<float,12UL,16UL,blaze::rowMajor> B( 4UL, 6UL );

   // Definition of an unaligned, unpadded 0x0 double precision column-major matrix and maximum dimensions of 6x6
   blaze::HybridMatrix<double,6UL,6UL,blaze::columnMajor,blaze::unaligned,blaze::unpadded> C;
   \endcode

// \subsubsection matrix_types_hybrid_matrix_alignment Alignment
//
// In case \c AF is set to \c blaze::aligned, the elements of a blaze::HybridMatrix are possibly
// over-aligned to meet the alignment requirements of the available instruction set (SSE, AVX,
// AVX-512, ...). The alignment for fundamental types (\c short, \c int, \c float, \c double, ...)
// and complex types (\c complex<float>, \c complex<double>, ...) is 16 bytes for SSE, 32 bytes
// for AVX, and 64 bytes for AVX-512. All other types are aligned according to their intrinsic
// alignment:

   \code
   struct Int { int i; };

   using MT1 = blaze::HybridMatrix<double,3UL,5UL>;
   using MT2 = blaze::HybridMatrix<complex<float>,2UL,3UL>;
   using MT3 = blaze::HybridMatrix<Int,5UL,4UL>;

   alignof( MT1 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( MT2 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( MT3 );  // Evaluates to 'alignof( Int )'
   \endcode

// Note that an aligned blaze::HybridMatrix instance may be bigger than an according unaligned
// blaze::HybridMatrix:

   \code
   sizeof( MT1 );  // Evaluates to 160 for SSE, 224 for AVX, and 256 for AVX-512
   sizeof( MT2 );  // Evaluates to 80 for SSE, 96 for AVX, and 192 for AVX-512
   sizeof( MT3 );  // Evaluates to 96; no special alignment requirements
   \endcode

// Please note that for this reason a blaze::HybridMatrix cannot be used in containers using
// dynamic memory such as \c std::vector without additionally providing an allocator that can
// provide over-aligned memory:

   \code
   using Type = blaze::HybridMatrix<double,3UL,5UL>;
   using Allocator = blaze::AlignedAllocator<Type>;

   std::vector<Type> v1;  // Might be misaligned for AVX or AVX-512
   std::vector<Type,Allocator> v2;  // Properly aligned for AVX or AVX-512
   \endcode

// \subsubsection matrix_types_hybrid_matrix_padding Padding
//
// Adding padding elements to the end of every row or column of a blaze::HybridMatrix can have a
// significant impact on the performance. For instance, assuming that AVX is available, then two
// padded 3x3 matrices of double precision values can be added with three SIMD addition operations:

   \code
   using blaze::HybridMatrix;
   using blaze::rowMajor;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   HybridMatrix<double,3UL,3UL,rowMajor,aligned,padded> A1, B1, C1;
   HybridMatrix<double,3UL,3UL,rowMajor,unaligned,unpadded> A2, B2, C2;

   // ... Initialization

   C1 = A1 + B1;  // AVX-based matrix addition; maximum performance
   C2 = A2 + B2;  // Scalar matrix addition; limited performance

   sizeof( A1 );  // Evaluates to 112 for SSE, 128 for AVX, and 256 for AVX-512
   sizeof( A2 );  // Evaluates to 88 for SSE, AVX, and AVX-512 (minimum size)
   \endcode

// Due to padding, the first addition will run at maximum performance. On the flip side, the size
// of each matrix instance is increased due to the padding elements. The total size of an instance
// depends on the number of elements and width of the available instruction set (16 bytes for
// SSE, 32 bytes for AVX, and 64 bytes for AVX-512).
//
// The second addition will be limited in performance since due to the number of elements some of
// the elements need to be handled in a scalar operation. However, the size of an \c unaligned,
// \c unpadded blaze::HybridMatrix instance is guaranteed to be the sum of its elements plus the.
// necessary data members to store the current number of rows and columns.
//
// Please also note that \b Blaze will zero initialize the padding elements in order to achieve
// maximum performance!
//
//
// \n \subsection matrix_types_custom_matrix CustomMatrix
//
// The blaze::CustomMatrix class template provides the functionality to represent an external
// array of elements of arbitrary type and a fixed size as a native \b Blaze dense matrix data
// structure. Thus in contrast to all other dense matrix types a custom matrix does not perform
// any kind of memory allocation by itself, but it is provided with an existing array of element
// during construction. A custom matrix can therefore be considered an alias to the existing
// array. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/CustomMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the properties of the given array of elements, the storage order,
// and the group tag of the matrix can be specified via the following five template parameters:

   \code
   namespace blaze {

   template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag >
   class CustomMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. blaze::CustomMatrix can be used with
//             any non-cv-qualified, non-reference, non-pointer element type.
//  - \c AF  : specifies whether the represented, external arrays are properly aligned with
//             respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::aligned
//             or \c blaze::unaligned).
//  - \c PF  : specified whether the represented, external arrays are properly padded with
//             respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::padded
//             or \c blaze::unpadded).
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::CustomMatrix is the right choice if any external array needs to be represented as
// a \b Blaze dense matrix data structure or if a custom memory allocation strategy needs to be
// realized:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of an unmanaged 3x4 custom matrix for unaligned, unpadded integer arrays
   using UnalignedUnpadded = CustomMatrix<int,unaligned,unpadded,rowMajor>;
   std::vector<int> vec( 12UL )
   UnalignedUnpadded A( &vec[0], 3UL, 4UL );

   // Definition of a managed 5x6 custom matrix for unaligned but padded 'float' arrays
   using UnalignedPadded = CustomMatrix<float,unaligned,padded,columnMajor>;
   std::unique_ptr<float[]> memory1( new float[40] );
   UnalignedPadded B( memory1.get(), 5UL, 6UL, 8UL );

   // Definition of a managed 12x13 custom matrix for aligned, unpadded 'double' arrays
   using AlignedUnpadded = CustomMatrix<double,aligned,unpadded,rowMajor>;
   std::unique_ptr<double[],Deallocate> memory2( blaze::allocate<double>( 192UL ) );
   AlignedUnpadded C( memory2.get(), 12UL, 13UL, 16UL );

   // Definition of a 7x14 custom matrix for aligned, padded 'complex<double>' arrays
   using cplx = complex<double>;
   using AlignedPadded = CustomMatrix<cplx,aligned,padded,columnMajor>;
   std::unique_ptr<cplx[],Deallocate> memory3( blaze::allocate<cplx>( 112UL ) );
   AlignedPadded D( memory3.get(), 7UL, 14UL, 16UL );
   \endcode

// In comparison with the remaining \b Blaze dense matrix types blaze::CustomMatrix has several
// special characteristics. All of these result from the fact that a custom matrix is not
// performing any kind of memory allocation, but instead is given an existing array of elements.
// The following sections discuss all of these characteristics:
//
//  -# <b>\ref matrix_types_custom_matrix_memory_management</b>
//  -# <b>\ref matrix_types_custom_matrix_copy_operations</b>
//  -# <b>\ref matrix_types_custom_matrix_alignment</b>
//  -# <b>\ref matrix_types_custom_matrix_padding</b>
//
// \subsubsection matrix_types_custom_matrix_memory_management Memory Management
//
// The blaze::CustomMatrix class template acts as an adaptor for an existing array of elements. As
// such it provides everything that is required to use the array just like a native \b Blaze dense
// matrix data structure. However, this flexibility comes with the price that the user of a custom
// matrix is responsible for the resource management.
//
// The following examples give an impression of several possible types of custom matrices:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of a 3x4 custom row-major matrix with unaligned, unpadded and externally
   // managed integer array. Note that the std::vector must be guaranteed to outlive the
   // custom matrix!
   std::vector<int> vec( 12UL );
   CustomMatrix<int,unaligned,unpadded> A( &vec[0], 3UL, 4UL );

   // Definition of a custom 8x12 matrix for an aligned and padded integer array of
   // capacity 128 (including 8 padding elements per row). Note that the std::unique_ptr
   // must be guaranteed to outlive the custom matrix!
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 128UL ) );
   CustomMatrix<int,aligned,padded> B( memory.get(), 8UL, 12UL, 16UL );
   \endcode

// \subsubsection matrix_types_custom_matrix_copy_operations Copy Operations
//
// As with all dense matrices it is possible to copy construct a custom matrix:

   \code
   using blaze::CustomMatrix;
   using blaze::unaligned;
   using blaze::unpadded;

   using CustomType = CustomMatrix<int,unaligned,unpadded>;

   std::vector<int> vec( 6UL, 10 );    // Vector of 6 integers of the value 10
   CustomType A( &vec[0], 2UL, 3UL );  // Represent the std::vector as Blaze dense matrix
   a[1] = 20;                          // Also modifies the std::vector

   CustomType B( a );  // Creating a copy of vector a
   b[2] = 20;          // Also affects matrix A and the std::vector
   \endcode

// It is important to note that a custom matrix acts as a reference to the specified array. Thus
// the result of the copy constructor is a new custom matrix that is referencing and representing
// the same array as the original custom matrix.
//
// In contrast to copy construction, just as with references, copy assignment does not change
// which array is referenced by the custom matrices, but modifies the values of the array:

   \code
   std::vector<int> vec2( 6UL, 4 );     // Vector of 6 integers of the value 4
   CustomType C( &vec2[0], 2UL, 3UL );  // Represent the std::vector as Blaze dense matrix

   A = C;  // Copy assignment: Set all values of matrix A and B to 4.
   \endcode

// \subsubsection matrix_types_custom_matrix_alignment Alignment
//
// In case the custom matrix is specified as \c aligned the passed array must adhere to some
// alignment restrictions based on the alignment requirements of the used data type and the
// used instruction set (SSE, AVX, ...). The restriction applies to the first element of each
// row/column: In case of a row-major matrix the first element of each row must be properly
// aligned, in case of a column-major matrix the first element of each column must be properly
// aligned. For instance, if a row-major matrix is used and AVX is active the first element of
// each row must be 32-bit aligned:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::padded;
   using blaze::rowMajor;

   // Allocation of 32-bit aligned memory
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 40UL ) );

   CustomMatrix<int,aligned,padded,rowMajor> A( memory.get(), 5UL, 6UL, 8UL );
   \endcode

// In the example, the row-major matrix has six columns. However, since with AVX eight integer
// values are loaded together the matrix is padded with two additional elements. This guarantees
// that the first element of each row is 32-bit aligned. In case the alignment requirements are
// violated, a \c std::invalid_argument exception is thrown.
//
// \subsubsection matrix_types_custom_matrix_padding Padding
//
// Adding padding elements to the end of each row/column can have a significant impact on the
// performance. For instance, assuming that AVX is available, then two aligned, padded, 3x3 double
// precision matrices can be added via three SIMD addition operations:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::padded;

   using CustomType = CustomMatrix<double,aligned,padded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 12UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 12UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 12UL ) );

   // Creating padded custom 3x3 matrix with an additional padding element in each row
   CustomType A( memory1.get(), 3UL, 3UL, 4UL );
   CustomType B( memory2.get(), 3UL, 3UL, 4UL );
   CustomType C( memory3.get(), 3UL, 3UL, 4UL );

   // ... Initialization

   C = A + B;  // AVX-based matrix addition
   \endcode

// In this example, maximum performance is possible. However, in case no padding elements are
// inserted a scalar addition has to be used:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unpadded;

   using CustomType = CustomMatrix<double,aligned,unpadded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 9UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 9UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 9UL ) );

   // Creating unpadded custom 3x3 matrix
   CustomType A( memory1.get(), 3UL, 3UL );
   CustomType B( memory2.get(), 3UL, 3UL );
   CustomType C( memory3.get(), 3UL, 3UL );

   // ... Initialization

   C = A + B;  // Scalar matrix addition
   \endcode

// Note that the construction of padded and unpadded aligned matrices looks identical. However,
// in case of padded matrices, \b Blaze will zero initialize the padding element and use them
// in all computations in order to achieve maximum performance. In case of an unpadded matrix
// \b Blaze will ignore the elements with the downside that it is not possible to load a complete
// row to an AVX register, which makes it necessary to fall back to a scalar addition.
//
// The number of padding elements is required to be sufficient with respect to the available
// instruction set: In case of an aligned padded custom matrix the added padding elements must
// guarantee that the total number of elements in each row/column is a multiple of the SIMD
// vector width. In case of an unaligned padded matrix the number of padding elements can be
// greater or equal the number of padding elements of an aligned padded custom matrix. In case
// the padding is insufficient with respect to the available instruction set, a
// \c std::invalid_argument exception is thrown.
//
//
// \n \subsection matrix_types_uniform_matrix UniformMatrix
//
// The blaze::UniformMatrix class template is the representation of an arbitrary sized uniform
// matrix with elements of arbitrary type. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/UniformMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the storage order, and the group tag of the matrix can be specified
// via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Tag >
   class UniformMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. UniformMatrix can be used with any
//             non-cv-qualified, non-reference element type.
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::UniformVector is the best choice for uniform matrices of any size. The number of
// rows and columns can be modified at runtime:

   \code
   // Definition of a 3x4 integral row-major matrix
   blaze::UniformMatrix<int> A( 3UL, 4UL );

   // Definition of a 4x6 single precision row-major matrix
   blaze::UniformMatrix<float,blaze::rowMajor> B( 4UL, 6UL );

   // Definition of a double precision column-major matrix with 0 rows and columns
   blaze::UniformMatrix<double,blaze::columnMajor> C;
   \endcode

// \n \section matrix_types_sparse_matrices Sparse Matrices
// <hr>
//
// \subsection matrix_types_compressed_matrix CompressedMatrix
//
// The blaze::CompressedMatrix class template is the representation of an arbitrary sized sparse
// matrix with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. It can be
// included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/CompressedMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the storage order, and the group tag of the matrix can be specified
// via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Tag >
   class CompressedMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. CompressedMatrix can be used with
//             any non-cv-qualified, non-reference, non-pointer element type.
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::CompressedMatrix is the right choice for all kinds of sparse matrices:

   \code
   // Definition of a 3x4 integral row-major matrix
   blaze::CompressedMatrix<int> A( 3UL, 4UL );

   // Definition of a 4x6 single precision row-major matrix
   blaze::CompressedMatrix<float,blaze::rowMajor> B( 4UL, 6UL );

   // Definition of a double precision column-major matrix with 0 rows and columns
   blaze::CompressedMatrix<double,blaze::columnMajor> C;
   \endcode

// \n \subsection matrix_types_identity_matrix IdentityMatrix
//
// The blaze::IdentityMatrix class template is the representation of an immutable, arbitrary
// sized identity matrix with \f$ N \cdot N \f$ elements of arbitrary type. It can be included
// via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/IdentityMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements and the storage order of the matrix can be specified via the three
// template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO >
   class IdentityMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. IdentityMatrix can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::IdentityMatrix is the perfect choice to represent an identity matrix:

   \code
   // Definition of a 3x3 integral row-major identity matrix
   blaze::IdentityMatrix<int> A( 3UL );

   // Definition of a 6x6 single precision row-major identity matrix
   blaze::IdentityMatrix<float,blaze::rowMajor> B( 6UL );

   // Definition of a double precision column-major identity matrix with 0 rows and columns
   blaze::IdentityMatrix<double,blaze::columnMajor> C;
   \endcode

// \n \subsection matrix_types_zero_matrix ZeroMatrix
//
// The blaze::ZeroMatrix class template is the representation of an immutable, arbitrary sized
// zero matrix with \f$ M \cdot N \f$ elements of arbitrary type. It can be included via the
// header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/ZeroMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the elements, the storage order, and the group tag of the matrix can be specified
// via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Tag >
   class ZeroMatrix;

   } // namespace blaze
   \endcode

//  - \c Type: specifies the type of the matrix elements. ZeroMatrix can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//             matrix. The default value is \c blaze::defaultStorageOrder.
//  - \c Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//             See \ref grouping_tagging for details.
//
// The blaze::ZeroMatrix is the perfect choice to represent a zero matrix:

   \code
   // Definition of a 3x5 integral row-major zero matrix
   blaze::ZeroMatrix<int> A( 3UL, 5UL );

   // Definition of a 6x4 single precision row-major zero matrix
   blaze::ZeroMatrix<float,blaze::rowMajor> B( 6UL, 4UL );

   // Definition of a double precision column-major zero matrix with 0 rows and columns
   blaze::ZeroMatrix<double,blaze::columnMajor> C;
   \endcode

// \n Previous: \ref matrices &nbsp; &nbsp; Next: \ref matrix_operations
*/
//*************************************************************************************************


//**Matrix Operations******************************************************************************
/*!\page matrix_operations Matrix Operations
//
// \tableofcontents
//
//
// \n \section matrix_operations_constructors Constructors
// <hr>
//
// Matrices are just as easy and intuitive to create as vectors. Still, there are a few rules
// to be aware of:
//  - In case the last template parameter (the storage order) is omitted, the matrix is per
//    default stored in row-major order.
//  - The elements of a \c StaticMatrix or \c HybridMatrix are default initialized (i.e. built-in
//    data types are initialized to 0, class types are initialized via the default constructor).
//  - Newly allocated elements of a \c DynamicMatrix or \c CompressedMatrix remain uninitialized
//    if they are of built-in type and are default constructed if they are of class type.
//
// \n \subsection matrix_operations_default_construction Default Construction

   \code
   using blaze::StaticMatrix;
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   // All matrices can be default constructed. Whereas the size of
   // a StaticMatrix is fixed via the second and third template
   // parameter, the initial size of a constructed DynamicMatrix
   // or CompressedMatrix is 0.
   StaticMatrix<int,2UL,2UL> M1;             // Instantiation of a 2x2 integer row-major
                                             // matrix. All elements are initialized to 0.
   DynamicMatrix<float> M2;                  // Instantiation of a single precision dynamic
                                             // row-major matrix with 0 rows and 0 columns.
   DynamicMatrix<double,columnMajor> M3;     // Instantiation of a double precision dynamic
                                             // column-major matrix with 0 rows and 0 columns.
   CompressedMatrix<int> M4;                 // Instantiation of a compressed integer
                                             // row-major matrix of size 0x0.
   CompressedMatrix<double,columnMajor> M5;  // Instantiation of a compressed double precision
                                             // column-major matrix of size 0x0.
   \endcode

// \n \subsection matrix_operations_size_construction Construction with Specific Size
//
// The \c DynamicMatrix, \c HybridMatrix, and \c CompressedMatrix classes offer a constructor
// that allows to immediately give the matrices a specific number of rows and columns:

   \code
   DynamicMatrix<int> M6( 5UL, 4UL );                   // Instantiation of a 5x4 dynamic row-major
                                                        // matrix. The elements are not initialized.
   HybridMatrix<double,5UL,9UL> M7( 3UL, 7UL );         // Instantiation of a 3x7 hybrid row-major
                                                        // matrix. The elements are not initialized.
   CompressedMatrix<float,columnMajor> M8( 8UL, 6UL );  // Instantiation of an empty 8x6 compressed
                                                        // column-major matrix.
   \endcode

// Note that dense matrices (in this case \c DynamicMatrix and \c HybridMatrix) immediately
// allocate enough capacity for all matrix elements. Sparse matrices on the other hand (in this
// example \c CompressedMatrix) merely acquire the size, but don't necessarily allocate memory.
//
//
// \n \subsection matrix_operations_initialization_constructors Initialization Constructors
//
// All dense matrix classes offer a constructor for a direct, homogeneous initialization of all
// matrix elements. In contrast, for sparse matrices the predicted number of non-zero elements
// can be specified.

   \code
   StaticMatrix<int,4UL,3UL,columnMajor> M9( 7 );  // Instantiation of a 4x3 integer column-major
                                                   // matrix. All elements are initialized to 7.
   DynamicMatrix<float> M10( 2UL, 5UL, 2.0F );     // Instantiation of a 2x5 single precision row-major
                                                   // matrix. All elements are initialized to 2.0F.
   CompressedMatrix<int> M11( 3UL, 4UL, 4 );       // Instantiation of a 3x4 integer row-major
                                                   // matrix with capacity for 4 non-zero elements.
   \endcode

// \n \subsection matrix_operations_array_construction Array Construction
//
// Alternatively, all dense matrix classes offer a constructor for an initialization with a dynamic
// or static array, or with a \c std::array. If the matrix is initialized from a dynamic array, the
// constructor expects the dimensions of values provided by the array as first and second argument,
// the array as third argument. In case of a static array or \c std::array, the fixed size of the
// array is used:

   \code
   const std::unique_ptr<double[]> array1( new double[6] );
   // ... Initialization of the dynamic array
   blaze::StaticMatrix<double,2UL,3UL> M12( 2UL, 3UL, array1.get() );

   int array2[2][2] = { { 4, -5 }, { -6, 7 } };
   blaze::StaticMatrix<int,2UL,2UL,rowMajor> M13( array2 );

   const std::array<std::array<float,3UL>,2UL> array3{ { { 1, 2, 3 }, { 4, 5, 6 } } };
   blaze::StaticMatrix<int,2UL,3UL> M14( array3 );
   \endcode

// \n \subsection matrix_operations_initializer_list_construction
//
// In addition, all dense and sparse matrix classes can be directly initialized by means of an
// initializer list:

   \code
   blaze::DynamicMatrix<float,columnMajor> M15{ {  3.1F,  6.4F },
                                                { -0.9F, -1.2F },
                                                {  4.8F,  0.6F } };
   blaze::CompressedMatrix<int,rowMajor> M16{ { 3 },
                                              { 1 },
                                              { 0, 2 } };
   \endcode

// Dynamically sized matrices (such as e.g. \ref matrix_types_hybrid_matrix,
// \ref matrix_types_dynamic_matrix or \ref matrix_types_compressed_matrix) are sized according
// to the size of the initializer list and all their elements are (copy) assigned the values of
// the list. For fixed size matrices (such as e.g. \ref matrix_types_static_matrix) missing values
// are initialized as default and in case the size of the top-level initializer list does not
// match the number of rows of the matrix or the size of any nested list exceeds the number of
// columns, a \c std::invalid_argument exception is thrown. In case of sparse matrices, only
// the non-zero elements are used to initialize the matrix.
//
// \n \subsection matrix_operations_copy_construction Copy Construction
//
// All dense and sparse matrices can be created as a copy of another dense or sparse matrix.

   \code
   StaticMatrix<int,5UL,4UL,rowMajor> M17( M6 );    // Instantiation of the dense row-major matrix M16
                                                    // as copy of the dense row-major matrix M6.
   DynamicMatrix<float,columnMajor> M18( M8 );      // Instantiation of the dense column-major matrix M17
                                                    // as copy of the sparse column-major matrix M8.
   CompressedMatrix<double,columnMajor> M19( M7 );  // Instantiation of the compressed column-major matrix
                                                    // M18 as copy of the dense row-major matrix M7.
   CompressedMatrix<float,rowMajor> M20( M8 );      // Instantiation of the compressed row-major matrix
                                                    // M19 as copy of the compressed column-major matrix M8.
   \endcode

// Note that it is not possible to create a \c StaticMatrix as a copy of a matrix with a different
// number of rows and/or columns:

   \code
   StaticMatrix<int,4UL,5UL,rowMajor> M21( M6 );     // Runtime error: Number of rows and columns
                                                     // does not match!
   StaticMatrix<int,4UL,4UL,columnMajor> M22( M9 );  // Compile time error: Number of columns does
                                                     // not match!
   \endcode

// \n \section matrix_operations_assignment Assignment
// <hr>
//
// There are several types of assignment to dense and sparse matrices:
// \ref matrix_operations_homogeneous_assignment, \ref matrix_operations_array_assignment,
// \ref matrix_operations_copy_assignment, and \ref matrix_operations_compound_assignment.
//
//
// \n \subsection matrix_operations_homogeneous_assignment Homogeneous Assignment
//
// It is possible to assign the same value to all elements of a dense matrix. All dense matrix
// classes provide an according assignment operator:

   \code
   blaze::StaticMatrix<int,3UL,2UL> M1;
   blaze::DynamicMatrix<double> M2;

   // Setting all integer elements of the StaticMatrix to 4
   M1 = 4;

   // Setting all double precision elements of the DynamicMatrix to 3.5
   M2 = 3.5
   \endcode

// \n \subsection matrix_operations_array_assignment Array Assignment
//
// Dense matrices can also be assigned a static array:

   \code
   blaze::StaticMatrix<int,2UL,2UL,rowMajor> M1;
   blaze::StaticMatrix<int,2UL,2UL,columnMajor> M2;
   blaze::DynamicMatrix<double> M3;

   int array1[2][2] = { { 1, 2 }, { 3, 4 } };
   double array2[3][2] = { { 3.1, 6.4 }, { -0.9, -1.2 }, { 4.8, 0.6 } };

   M1 = array1;
   M2 = array1;
   M3 = array2;
   \endcode

// Note that the dimensions of the static array have to match the size of a \c StaticMatrix,
// whereas a \c DynamicMatrix is resized according to the array dimensions:

                          \f$ M3 = \left(\begin{array}{*{2}{c}}
                           3.1 &  6.4 \\
                          -0.9 & -1.2 \\
                           4.8 &  0.6 \\
                          \end{array}\right)\f$

// \n \subsection matrix_operations_initializer_list_assignment Initializer List Assignment
//
// Alternatively, it is possible to directly assign an initializer list to a dense or sparse
// matrix:

   \code
   blaze::DynamicMatrix<double> M1;
   blaze::CompressedMatrix<int> M2;

   M1 = { { 3.1, 6.4 }, { -0.9, -1.2 }, { 4.8, 0.6 } };
   M2 = { { 1, 0 }, {}, { 0, 1 }, { 2 } };
   \endcode

// Dynamically sized matrices (such as e.g. \ref matrix_types_hybrid_matrix,
// \ref matrix_types_dynamic_matrix or \ref matrix_types_compressed_matrix) are resized according
// to the size of the initializer list and all their elements are (copy) assigned the values of
// the list. For fixed size matrices (such as e.g. \ref matrix_types_static_matrix) missing values
// are reset to their default value and in case the size of the top-level initializer list does
// not match the number of rows of the matrix or the size of any nested list exceeds the number
// of columns, a \c std::invalid_argument exception is thrown. In case of sparse matrices, only
// the non-zero elements are considered.
//
// \n \subsection matrix_operations_copy_assignment Copy Assignment
//
// All kinds of matrices can be assigned to each other. The only restriction is that since a
// \c StaticMatrix cannot change its size, the assigned matrix must match both in the number of
// rows and in the number of columns.

   \code
   blaze::StaticMatrix<int,3UL,2UL,rowMajor>  M1;
   blaze::DynamicMatrix<int,rowMajor>         M2( 3UL, 2UL );
   blaze::DynamicMatrix<float,rowMajor>       M3( 5UL, 2UL );
   blaze::CompressedMatrix<int,rowMajor>      M4( 3UL, 2UL );
   blaze::CompressedMatrix<float,columnMajor> M5( 3UL, 2UL );

   // ... Initialization of the matrices

   M1 = M2;  // OK: Assignment of a 3x2 dense row-major matrix to another 3x2 dense row-major matrix
   M1 = M4;  // OK: Assignment of a 3x2 sparse row-major matrix to a 3x2 dense row-major matrix
   M1 = M3;  // Runtime error: Cannot assign a 5x2 matrix to a 3x2 static matrix
   M1 = M5;  // OK: Assignment of a 3x2 sparse column-major matrix to a 3x2 dense row-major matrix
   \endcode

// \n \subsection matrix_operations_compound_assignment Compound Assignment
//
// Compound assignment is also available for matrices: addition assignment, subtraction assignment,
// and multiplication assignment. In contrast to plain assignment, however, the number of rows
// and columns of the two operands have to match according to the arithmetic operation.

   \code
   blaze::StaticMatrix<int,2UL,3UL,rowMajor>    M1;
   blaze::DynamicMatrix<int,rowMajor>           M2( 2UL, 3UL );
   blaze::CompressedMatrix<float,columnMajor>   M3( 2UL, 3UL );
   blaze::CompressedMatrix<float,rowMajor>      M4( 2UL, 4UL );
   blaze::StaticMatrix<float,2UL,4UL,rowMajor>  M5;
   blaze::CompressedMatrix<float,rowMajor>      M6( 3UL, 2UL );

   // ... Initialization of the matrices

   M1 += M2;  // OK: Addition assignment between two row-major matrices of the same dimensions
   M1 -= M3;  // OK: Subtraction assignment between between a row-major and a column-major matrix
   M1 += M4;  // Runtime error: No compound assignment between matrices of different size
   M1 -= M5;  // Compilation error: No compound assignment between matrices of different size
   M2 *= M6;  // OK: Multiplication assignment between two row-major matrices
   \endcode

// Note that the multiplication assignment potentially changes the number of columns of the
// target matrix:

                          \f$\left(\begin{array}{*{3}{c}}
                          2 & 0 & 1 \\
                          0 & 3 & 2 \\
                          \end{array}\right) \times
                          \left(\begin{array}{*{2}{c}}
                          4 & 0 \\
                          1 & 0 \\
                          0 & 3 \\
                          \end{array}\right) =
                          \left(\begin{array}{*{2}{c}}
                          8 & 3 \\
                          3 & 6 \\
                          \end{array}\right)\f$

// Since a \c StaticMatrix cannot change its size, only a square StaticMatrix can be used in a
// multiplication assignment with other square matrices of the same dimensions.
//
//
// \n \section matrix_operations_element_access Element Access
// <hr>
//
// \subsection matrix_operations_function_call_operator_1 Function Call Operator
//
// The easiest way to access a specific dense or sparse matrix element is via the function call
// operator. The indices to access a matrix are zero-based:

   \code
   blaze::DynamicMatrix<int> M1( 4UL, 6UL );
   M1(0,0) = 1;
   M1(0,1) = 3;
   // ...

   blaze::CompressedMatrix<double> M2( 5UL, 3UL );
   M2(0,2) =  4.1;
   M2(1,1) = -6.3;
   \endcode

// Since dense matrices allocate enough memory for all contained elements, using the function
// call operator on a dense matrix directly returns a reference to the accessed value. In case
// of a sparse matrix, if the accessed value is currently not contained in the matrix, the
// value is inserted into the matrix prior to returning a reference to the value, which can
// be much more expensive than the direct access to a dense matrix. Consider the following
// example:

   \code
   blaze::CompressedMatrix<int> M1( 4UL, 4UL );

   for( size_t i=0UL; i<M1.rows(); ++i ) {
      for( size_t j=0UL; j<M1.columns(); ++j ) {
         ... = M1(i,j);
      }
   }
   \endcode

// Although the compressed matrix is only used for read access within the for loop, using the
// function call operator temporarily inserts 16 non-zero elements into the matrix. Therefore
// the preferred way to traverse the non-zero elements of a sparse matrix is to use iterators.
//
// \n \subsection matrix_operations_iterators Iterators
//
// An alternate way to traverse the elements contained in a dense or sparse matrix is by means
// of iterators. For that purpose, all matrices provide the \c begin(), \c cbegin(), \c end(),
// and \c cend() members functions. Note that it is not possible to traverse all elements of the
// matrix, but that it is only possible to traverse elements in a row-wise fashion (in case of
// a row-major matrix) or in a column-wise fashion (in case of a column-major matrix). In case of
// non-const matrices, \c begin() and \c end() return an \c Iterator, which allows a manipulation
// of the (non-zero) value. In case of a constant matrix or in case \c cbegin() or \c cend() are
// used a \c ConstIterator is returned. Iterators on dense matrices traverse all elements of the
// matrix, including the zero elements. Iterators on sparse matrices only traverse the non-zero
// elements.
//
// The following two examples demonstrate how to traverse the elements of a dense and sparse
// matrix, respectively:

   \code
   using blaze::DynamicMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   DynamicMatrix<int,rowMajor> M1( 4UL, 6UL );
   DynamicMatrix<int,columnMajor> M2( 4UL, 6UL );

   // Traversing all elements contained in the row-major matrix by Iterator
   for( size_t i=0UL; i<M1.rows(); ++i ) {
      for( DynamicMatrix<int,rowMajor>::Iterator it=M1.begin(i); it!=M1.end(i); ++it ) {
         *it = ...;  // OK: Write access to the value of the element.
         ... = *it;  // OK: Read access to the value of the element.
      }
   }

   // Traversing all elements contained in the column-major matrix by ConstIterator
   for( size_t j=0UL; j<M2.columns(); ++j ) {
      for( DynamicMatrix<int,columnMajor>::ConstIterator it=M2.cbegin(j); it!=M2.cend(j); ++it ) {
         *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
         ... = *it;  // OK: Read access to the value of the element.
      }
   }
   \endcode

   \code
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   CompressedMatrix<int,rowMajor> M3( 4UL, 6UL );
   CompressedMatrix<int,columnMajor> M4( 4UL, 6UL );

   // Traversing the non-zero elements contained in the row-major matrix by Iterator
   for( size_t i=0UL; i<M3.rows(); ++i ) {
      for( CompressedMatrix<int,rowMajor>::Iterator it=M3.begin(i); it!=M3.end(i); ++it ) {
         it->value() = ...;  // OK: Write access to the value of the non-zero element.
         ... = it->value();  // OK: Read access to the value of the non-zero element.
         it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
         ... = it->index();  // OK: Read access to the index of the non-zero element.
      }
   }

   // Traversing the non-zero elements contained in the column-major matrix by ConstIterator
   for( size_t j=0UL; j<M4.columns(); ++j ) {
      for( CompressedMatrix<int,columnMajor>::ConstIterator it=M4.cbegin(j); it!=M4.cend(j); ++it ) {
         it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
         ... = it->value();  // OK: Read access to the value of the non-zero element.
         it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
         ... = it->index();  // OK: Read access to the index of the non-zero element.
      }
   }
   \endcode

// Note that \c begin(), \c cbegin(), \c end(), and \c cend() are also available as free functions:

   \code
   for( size_t i=0UL; i<M3.rows(); ++i ) {
      for( CompressedMatrix<int,rowMajor>::Iterator it=begin( M3, i ); it!=end( M3, i ); ++it ) {
         // ...
      }
   }

   for( size_t j=0UL; j<M4.columns(); ++j ) {
      for( CompressedMatrix<int,columnMajor>::ConstIterator it=cbegin( M4, j ); it!=cend( M4, j ); ++it ) {
         // ...
      }
   }
   \endcode

// \n \subsection matrix_operations_data .data() / data()
//
// Sometimes it is necessary to acquire a pointer to the first element of the underlying array
// of a dense matrix. For that purpose the \c data() member function or the free \c data() function
// can be used:

   \code
   // Instantiating a dynamic vector with 10 elements
   blaze::DynamicMatrix<int> A( 5UL, 7UL );
   A.data();   // Returns a pointer to the first element of the dynamic matrix
   data( A );  // Same effect as the member function
   \endcode

// Note that you can NOT assume that all matrix elements lie adjacent to each other! The dense
// matrix may use techniques such as padding to improve the alignment of the data. Whereas the
// number of elements within a row/column are given by the \ref matrix_operations_rows "rows()" and
// \ref matrix_operations_columns "columns()" functions, respectively, the total number of elements including
// padding is given by the \ref matrix_operations_spacing "spacing()" function.
//
//
// \n \section matrix_operations_element_insertion Element Insertion
// <hr>
//
// Whereas a dense matrix always provides enough capacity to store all matrix elements, a sparse
// matrix only stores the non-zero elements. Therefore it is necessary to explicitly add elements
// to the matrix.
//
// \n \subsection matrix_operations_function_call_operator_2 Function Call Operator
//
// The first possibility to add elements to a sparse matrix is the function call operator:

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int> M1( 3UL, 4UL );
   M1(1,2) = 9;
   \endcode

// In case the element at the given position is not yet contained in the sparse matrix, it is
// automatically inserted. Otherwise the old value is replaced by the new value 2. The operator
// returns a reference to the sparse vector element.
//
// \n \subsection matrix_operations_set .set()
//
// An alternative to the function call operator is the \c set() function: In case the element is
// not yet contained in the matrix the element is inserted, else the element's value is modified:

   \code
   // Insert or modify the value at position (2,0)
   M1.set( 2, 0, 1 );
   \endcode

// \n \subsection matrix_operations_insert .insert()

// The insertion of elements can be better controlled via the \c insert() function. In contrast
// to the function call operator and the \c set() function it emits an exception in case the
// element is already contained in the matrix. In order to check for this case, the \c find()
// function can be used:

   \code
   // In case the element at position (2,3) is not yet contained in the matrix it is inserted
   // with a value of 4.
   if( M1.find( 2, 3 ) == M1.end( 2 ) )
      M1.insert( 2, 3, 4 );
   \endcode

// \n \subsection matrix_operations_append .append()
//
// Although the \c insert() function is very flexible, due to performance reasons it is not
// suited for the setup of large sparse matrices. A very efficient, yet also very low-level
// way to fill a sparse matrix is the \c append() function. It requires the sparse matrix to
// provide enough capacity to insert a new element in the specified row/column. Additionally,
// the index of the new element must be larger than the index of the previous element in the
// same row/column. Violating these conditions results in undefined behavior!

   \code
   M1.reserve( 0, 3 );     // Reserving space for three non-zero elements in row 0
   M1.append( 0, 1,  2 );  // Appending the element 2 in row 0 at column index 1
   M1.append( 0, 2, -4 );  // Appending the element -4 in row 0 at column index 2
   // ...
   \endcode

// The most efficient way to fill a sparse matrix with elements, however, is a combination of
// \c reserve(), \c append(), and the \c finalize() function:

   \code
   // Setup of the compressed row-major matrix
   //
   //       ( 0 1 0 2 0 )
   //   A = ( 0 0 0 0 0 )
   //       ( 3 0 0 0 0 )
   //
   blaze::CompressedMatrix<int> M1( 3UL, 5UL );
   M1.reserve( 3 );       // Reserving enough space for 3 non-zero elements
   M1.append( 0, 1, 1 );  // Appending the value 1 in row 0 with column index 1
   M1.append( 0, 3, 2 );  // Appending the value 2 in row 0 with column index 3
   M1.finalize( 0 );      // Finalizing row 0
   M1.finalize( 1 );      // Finalizing the empty row 1 to prepare row 2
   M1.append( 2, 0, 3 );  // Appending the value 3 in row 2 with column index 0
   M1.finalize( 2 );      // Finalizing row 2
   \endcode

// \note The \c finalize() function has to be explicitly called for each row or column, even
// for empty ones!
// \note Although \c append() does not allocate new memory, it still invalidates all iterators
// returned by the \c end() functions!
//
//
// \n \section matrix_operations_element_removal Element Removal
// <hr>
//
// \subsection matrix_operations_erase .erase()
//
// The \c erase() member functions can be used to remove elements from a sparse matrix. The
// following example gives an impression of the five different flavors of \c erase():

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int,rowMajor> A( 42, 53 );
   // ... Initialization of the matrix

   // Erasing the element at position (21,23)
   A.erase( 21, 23 );

   // Erasing a single element in row 17 via iterator
   A.erase( 17, A.find( 4 ) );

   // Erasing all non-zero elements in the range [7..24] of row 33
   A.erase( 33, A.lowerBound( 33, 7 ), A.upperBound( 33, 24 ) );

   // Erasing all non-zero elements with a value larger than 9 by passing a unary predicate
   A.erase( []( int i ){ return i > 9; } );

   // Erasing all non-zero elements in the range [30..40] of row 37 with a value larger than 5
   CompressedMatrix<int,rowMajor>::Iterator pos1( A.lowerBound( 37, 30 ) );
   CompressedMatrix<int,rowMajor>::Iterator pos2( A.upperBound( 37, 40 ) );
   A.erase( 37, pos1, pos2, []( int i ){ return i > 5; } );
   \endcode

// \n \section matrix_operations_element_lookup Element Lookup
// <hr>
//
// A sparse matrix only stores the non-zero elements contained in the matrix. Therefore, whenever
// accessing a matrix element at a specific position a lookup operation is required. Whereas the
// function call operator is performing this lookup automatically, it is also possible to use the
// \c find(), \c lowerBound(), and \c upperBound() member functions for a manual lookup.
//
// \n \subsection matrix_operations_find .find() / find()
//
// The \c find() function can be used to check whether a specific element is contained in the
// sparse matrix. It specifically searches for the element at the specified position. In case
// the element is found, the function returns an iterator to the element. Otherwise an iterator
// just past the last non-zero element of the according row or column (the \c end() iterator)
// is returned. Note that the returned iterator is subject to invalidation due to inserting
// operations via the function call operator, the \c set() function or the \c insert() function!

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int,rowMajor> A( 42, 53 );
   // ... Initialization of the matrix

   // Searching the element at position (7,17). In case the element is not
   // contained in the vector, the end() iterator of row 7 is returned.
   CompressedMatrix<int,rowMajor>::Iterator pos( A.find( 7, 17 ) );

   if( pos != A.end( 7 ) ) {
      // ...
   }
   \endcode

// Alternatively, the free function \c find() can be used to find a specific element in a sparse
// matrix:

   \code
   find( A, 7, 17 );  // Searching the element at position (7,17); same effect as the member function
   \endcode

// \n \subsection matrix_operations_lowerbound .lowerBound() / lowerBound()
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index not less then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index not less then the given row
// index. In combination with the \c upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the \c set()
// function or the \c insert() function!

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int,rowMajor> A( 42, 53 );
   // ... Initialization of the matrix

   // Searching the lower bound of column index 17 in row 7.
   CompressedMatrix<int,rowMajor>::Iterator pos1( A.lowerBound( 7, 17 ) );

   // Searching the upper bound of column index 28 in row 7
   CompressedMatrix<int,rowMajor>::Iterator pos2( A.upperBound( 7, 28 ) );

   // Erasing all elements in the specified range
   A.erase( 7, pos1, pos2 );
   \endcode

// Alternatively, the free function \c lowerBound() can be used to:

   \code
   lowerBound( A, 7, 17 );  // Searching the lower bound of (7,17); same effect as the member function
   \endcode

// \n \subsection matrix_operations_upperbound .upperBound() / upperBound()
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index greater then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index greater then the given row
// index. In combination with the \c lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the \c set()
// function or the \c insert() function!

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int,columnMajor> A( 42, 53 );
   // ... Initialization of the matrix

   // Searching the lower bound of row index 17 in column 9.
   CompressedMatrix<int,columnMajor>::Iterator pos1( A.lowerBound( 17, 9 ) );

   // Searching the upper bound of row index 28 in column 9
   CompressedMatrix<int,columnMajor>::Iterator pos2( A.upperBound( 28, 9 ) );

   // Erasing all elements in the specified range
   A.erase( 9, pos1, pos2 );
   \endcode

// Alternatively, the free function \c upperBound() can be used to:

   \code
   upperBound( A, 28, 9 );  // Searching the upper bound of (28,9); same effect as the member function
   \endcode

// \n \section matrix_operations_non_modifying_operations Non-Modifying Operations
// <hr>
//
// \subsection matrix_operations_rows .rows() / rows()
//
// The current number of rows of a matrix can be acquired via the \c rows() member function:

   \code
   // Instantiating a dynamic matrix with 10 rows and 8 columns
   blaze::DynamicMatrix<int> M1( 10UL, 8UL );
   M1.rows();  // Returns 10

   // Instantiating a compressed matrix with 8 rows and 12 columns
   blaze::CompressedMatrix<double> M2( 8UL, 12UL );
   M2.rows();  // Returns 8
   \endcode

// Alternatively, the free functions \c rows() can be used to query the current number of rows of
// a matrix. In contrast to the member function, the free function can also be used to query the
// number of rows of a matrix expression:

   \code
   rows( M1 );  // Returns 10, i.e. has the same effect as the member function
   rows( M2 );  // Returns 8, i.e. has the same effect as the member function

   rows( M1 * M2 );  // Returns 10, i.e. the number of rows of the resulting matrix
   \endcode

// \n \subsection matrix_operations_columns .columns() / columns()
//
// The current number of columns of a matrix can be acquired via the \c columns() member function:

   \code
   // Instantiating a dynamic matrix with 6 rows and 8 columns
   blaze::DynamicMatrix<int> M1( 6UL, 8UL );
   M1.columns();   // Returns 8

   // Instantiating a compressed matrix with 8 rows and 7 columns
   blaze::CompressedMatrix<double> M2( 8UL, 7UL );
   M2.columns();   // Returns 7
   \endcode

// There is also a free function \c columns() available, which can also be used to query the number
// of columns of a matrix expression:

   \code
   columns( M1 );  // Returns 8, i.e. has the same effect as the member function
   columns( M2 );  // Returns 7, i.e. has the same effect as the member function

   columns( M1 * M2 );  // Returns 7, i.e. the number of columns of the resulting matrix
   \endcode

// \subsection matrix_operations_size size()
//
// The \c size() function returns the total number of elements of a matrix:

   \code
   // Instantiating a dynamic matrix with 6 rows and 8 columns
   blaze::DynamicMatrix<int> M1( 6UL, 8UL );
   size( M1 );   // Returns 48

   // Instantiating a compressed matrix with 8 rows and 7 columns
   blaze::CompressedMatrix<double> M2( 8UL, 7UL );
   size( M2 );  // Returns 56
   \endcode

// \subsection matrix_operations_spacing .spacing() / spacing()
//
// The total number of elements of a row or column of a dense matrix, including potential padding
// elements, can be acquired via the \c spacing member function. In case of a row-major matrix
// (i.e. in case the storage order is set to blaze::rowMajor) the function returns the spacing
// between two rows, in case of a column-major matrix (i.e. in case the storage flag is set to
// blaze::columnMajor) the function returns the spacing between two columns:

   \code
   // Instantiating a row-major dynamic matrix with 7 rows and 8 columns
   blaze::DynamicMatrix<int,blaze::rowMajor> M1( 7UL, 8UL );
   M1.spacing();  // Returns the total number of elements in a row

   // Instantiating a column-major dynamic matrix with 8 rows and 12 columns
   blaze::CompressedMatrix<double> M2( 8UL, 12UL );
   M2.spacing();  // Returns the total number of element in a column
   \endcode

// Alternatively, the free functions \c spacing() can be used to query the current number of
// elements in a row/column.

   \code
   spacing( M1 );  // Returns the total number of elements in a row
   spacing( M2 );  // Returns the total number of elements in a column
   \endcode

// \n \subsection matrix_operations_capacity .capacity() / capacity()
//
// The \c capacity() member function returns the internal capacity of a dense or sparse matrix.
// Note that the capacity of a matrix doesn't have to be equal to the size of a matrix. In case of
// a dense matrix the capacity will always be greater or equal than the total number of elements
// of the matrix. In case of a sparse matrix, the capacity will usually be much less than the
// total number of elements.

   \code
   blaze::DynamicMatrix<float> M1( 5UL, 7UL );
   blaze::StaticMatrix<float,7UL,4UL> M2;
   M1.capacity();  // Returns at least 35
   M2.capacity();  // Returns at least 28
   \endcode

// There is also a free function \c capacity() available to query the capacity. However, please
// note that this function cannot be used to query the capacity of a matrix expression:

   \code
   capacity( M1 );  // Returns at least 35, i.e. has the same effect as the member function
   capacity( M2 );  // Returns at least 28, i.e. has the same effect as the member function

   capacity( M1 * M2 );  // Compilation error!
   \endcode

// \n \subsection matrix_operations_nonzeros .nonZeros() / nonZeros()
//
// For both dense and sparse matrices the current number of non-zero elements can be queried
// via the \c nonZeros() member function. In case of matrices there are two flavors of the
// \c nonZeros() function: One returns the total number of non-zero elements in the matrix,
// the second returns the number of non-zero elements in a specific row (in case of a row-major
// matrix) or column (in case of a column-major matrix). Sparse matrices directly return their
// number of non-zero elements, dense matrices traverse their elements and count the number of
// non-zero elements.

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 3UL, 5UL );

   // ... Initializing the dense matrix

   M1.nonZeros();     // Returns the total number of non-zero elements in the dense matrix
   M1.nonZeros( 2 );  // Returns the number of non-zero elements in row 2
   \endcode

   \code
   blaze::CompressedMatrix<double,columnMajor> M2( 4UL, 7UL );

   // ... Initializing the sparse matrix

   M2.nonZeros();     // Returns the total number of non-zero elements in the sparse matrix
   M2.nonZeros( 3 );  // Returns the number of non-zero elements in column 3
   \endcode

// The free \c nonZeros() function can also be used to query the number of non-zero elements in a
// matrix expression. However, the result is not the exact number of non-zero elements, but may be
// a rough estimation:

   \code
   nonZeros( M1 );     // Has the same effect as the member function
   nonZeros( M1, 2 );  // Has the same effect as the member function

   nonZeros( M2 );     // Has the same effect as the member function
   nonZeros( M2, 3 );  // Has the same effect as the member function

   nonZeros( M1 * M2 );  // Estimates the number of non-zero elements in the matrix expression
   \endcode

// \n \subsection matrix_operations_isempty isEmpty()
//
// The \c isEmpty() function returns whether the total number of elements of the matrix is zero:

   \code
   blaze::DynamicMatrix<int> A;  // Create an empty matrix
   isEmpty( A );                 // Returns true
   A.resize( 5, 0 );             // Resize to a 5x0 matrix
   isEmpty( A );                 // Returns true
   A.resize( 5, 3 );             // Resize to a 5x3 matrix
   isEmpty( A );                 // Returns false
   \endcode

// \n \subsection matrix_operations_isnan isnan()
//
// The \c isnan() function provides the means to check a dense or sparse matrix for non-a-number
// elements:

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isnan( A ) ) { ... }
   \endcode

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isnan( A ) ) { ... }
   \endcode

// If at least one element of the matrix is not-a-number, the function returns \c true, otherwise
// it returns \c false.
//
//
// \n \subsection matrix_operations_isinf isinf()
//
// The \c isinf() function checks the given dense or sparse matrix for infinite (\c inf) elements:

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isinf( A ) ) { ... }
   \endcode

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isinf( A ) ) { ... }
   \endcode

// If at least one element of the matrix is infinite, the function returns \c true, otherwise it
// returns \c false.
//
//
// \n \subsection matrix_operations_isfinite isfinite()
//
// The \c isfinite() function checks if all elements of the given dense or sparse matrix are
// finite elements (i.e. normal, subnormal or zero elements, but not infinite or NaN):

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isfinite( A ) ) { ... }
   \endcode

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isfinite( A ) ) { ... }
   \endcode

// If all elements of the matrix are finite, the function returns \c true, otherwise it returns
// \c false.
//
//
// \n \subsection matrix_operations_isdefault isDefault()
//
// The \c isDefault() function returns whether the given dense or sparse matrix is in default state:

   \code
   blaze::HybridMatrix<int,5UL,4UL> A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode

// A matrix is in default state if it appears to just have been default constructed. All resizable
// matrices (\c HybridMatrix, \c DynamicMatrix, or \c CompressedMatrix) and \c CustomMatrix are in
// default state if its size is equal to zero. A non-resizable matrix (\c StaticMatrix and all
// submatrices) is in default state if all its elements are in default state. For instance, in case
// the matrix is instantiated for a built-in integral or floating point data type, the function
// returns \c true in case all matrix elements are 0 and \c false in case any matrix element is
// not 0.
//
//
// \n \subsection matrix_operations_isSquare isSquare()
//
// Whether a dense or sparse matrix is a square matrix (i.e. if the number of rows is equal to the
// number of columns) can be checked via the \c isSquare() function:

   \code
   blaze::DynamicMatrix<double> A;
   // ... Resizing and initialization
   if( isSquare( A ) ) { ... }
   \endcode

// \n \subsection matrix_operations_issymmetric isSymmetric()
//
// Via the \c isSymmetric() function it is possible to check whether a dense or sparse matrix
// is symmetric:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isSymmetric( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be symmetric!
//
//
// \n \subsection matrix_operations_isUniform isUniform()
//
// In order to check if all matrix elements are identical, the \c isUniform() function can be used:

   \code
   blaze::DynamicMatrix<int> A;
   // ... Resizing and initialization
   if( isUniform( A ) ) { ... }
   \endcode

// Note that in case of a sparse matrix also the zero elements are also taken into account!
//
//
// \n \subsection matrix_operations_isZero isZero()
//
// In order to check if all matrix elements are zero, the \c isZero() function can be used:

   \code
   blaze::DynamicMatrix<int> A;
   // ... Resizing and initialization
   if( isZero( A ) ) { ... }
   \endcode

// \n \subsection matrix_operations_islower isLower()
//
// Via the \c isLower() function it is possible to check whether a dense or sparse matrix is
// lower triangular:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isLower( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be lower triangular!
//
//
// \n \subsection matrix_operations_isunilower isUniLower()
//
// Via the \c isUniLower() function it is possible to check whether a dense or sparse matrix is
// lower unitriangular:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isUniLower( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be lower unitriangular!
//
//
// \n \subsection matrix_operations_isstrictlylower isStrictlyLower()
//
// Via the \c isStrictlyLower() function it is possible to check whether a dense or sparse matrix
// is strictly lower triangular:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isStrictlyLower( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be strictly lower triangular!
//
//
// \n \subsection matrix_operations_isUpper isUpper()
//
// Via the \c isUpper() function it is possible to check whether a dense or sparse matrix is
// upper triangular:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isUpper( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be upper triangular!
//
//
// \n \subsection matrix_operations_isuniupper isUniUpper()
//
// Via the \c isUniUpper() function it is possible to check whether a dense or sparse matrix is
// upper unitriangular:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isUniUpper( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be upper unitriangular!
//
//
// \n \subsection matrix_operations_isstrictlyupper isStrictlyUpper()
//
// Via the \c isStrictlyUpper() function it is possible to check whether a dense or sparse matrix
// is strictly upper triangular:

   \code
   blaze::DynamicMatrix<float> A;
   // ... Resizing and initialization
   if( isStrictlyUpper( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be strictly upper triangular!
//
//
// \n \subsection matrix_operations_isdiagonal isDiagonal()
//
// The \c isDiagonal() function checks if the given dense or sparse matrix is a diagonal matrix,
// i.e. if it has only elements on its diagonal and if the non-diagonal elements are default
// elements:

   \code
   blaze::CompressedMatrix<float> A;
   // ... Resizing and initialization
   if( isDiagonal( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be diagonal!
//
//
// \n \subsection matrix_operations_isidentity isIdentity()
//
// The \c isIdentity() function checks if the given dense or sparse matrix is an identity matrix,
// i.e. if all diagonal elements are 1 and all non-diagonal elements are 0:

   \code
   blaze::CompressedMatrix<float> A;
   // ... Resizing and initialization
   if( isIdentity( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be identity matrices!
//
//
// \n \subsection matrix_operations_ispositivedefinite isPositiveDefinite()
//
// The \c isPositiveDefinite() function checks if the given dense matrix is positive definite.

   \code
   blaze::DynamicMatrix<double> A;
   // ... Initialization
   if( isPositiveDefinite( A ) ) { ... }
   \endcode

// Note that non-square matrices are never considered to be positive definite!
//
// \note The \c isPositiveDefinite() function can only be used for dense matrices with \c float,
// \c double, \c complex<float> or \c complex<double> element type. The attempt to call the
// function with matrices of any other element type or with a sparse matrix results in a compile
// time error!
//
// \note The function is depending on LAPACK kernels. Thus the function can only be used if a
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error
// will be created.
//
//
// \n \subsection matrix_operations_matrix_trans trans()
//
// Matrices can be transposed via the \c trans() function. Row-major matrices are transposed into
// a column-major matrix and vice versa:

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 5UL, 2UL );
   blaze::CompressedMatrix<int,columnMajor> M2( 3UL, 7UL );

   M1 = M2;            // Assigning a column-major matrix to a row-major matrix
   M1 = trans( M2 );   // Assigning the transpose of M2 (i.e. a row-major matrix) to M1
   M1 += trans( M2 );  // Addition assignment of two row-major matrices
   \endcode

// \n \subsection matrix_operations_ctrans ctrans()
//
// The conjugate transpose of a dense or sparse matrix (also called adjoint matrix, Hermitian
// conjugate, or transjugate) can be computed via the \c ctrans() function:

   \code
   blaze::DynamicMatrix< complex<float>, rowMajor > M1( 5UL, 2UL );
   blaze::CompressedMatrix< complex<float>, columnMajor > M2( 2UL, 5UL );

   M1 = ctrans( M2 );  // Compute the conjugate transpose matrix
   \endcode

// Note that the \c ctrans() function has the same effect as manually applying the \c conj() and
// \c trans() function in any order:

   \code
   M1 = trans( conj( M2 ) );  // Computing the conjugate transpose matrix
   M1 = conj( trans( M2 ) );  // Computing the conjugate transpose matrix
   \endcode

// \n \subsection matrix_operations_reverse reverse()
//
// Via the \c reverse() function is is possible to reverse the rows or columns of a dense or sparse
// matrix. The following examples gives an impression of both alternatives:

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

// \n \subsection matrix_operations_evaluate eval() / evaluate()
//
// The \c evaluate() function forces an evaluation of the given matrix expression and enables
// an automatic deduction of the correct result type of an operation. The following code example
// demonstrates its intended use for the multiplication of a lower and a strictly lower dense
// matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::StrictlyLowerMatrix;

   LowerMatrix< DynamicMatrix<double> > A;
   StrictlyLowerMatrix< DynamicMatrix<double> > B;
   // ... Resizing and initialization

   auto C = evaluate( A * B );
   \endcode

// In this scenario, the \c evaluate() function assists in deducing the exact result type of
// the operation via the \c auto keyword. Please note that if \c evaluate() is used in this
// way, no temporary matrix is created and no copy operation is performed. Instead, the result
// is directly written to the target matrix due to the return value optimization (RVO). However,
// if \c evaluate() is used in combination with an explicit target type, a temporary will be
// created and a copy operation will be performed if the used type differs from the type
// returned from the function:

   \code
   StrictlyLowerMatrix< DynamicMatrix<double> > D( A * B );  // No temporary & no copy operation
   LowerMatrix< DynamicMatrix<double> > E( A * B );          // Temporary & copy operation
   DynamicMatrix<double> F( A * B );                         // Temporary & copy operation
   D = evaluate( A * B );                                    // Temporary & copy operation
   \endcode

// Sometimes it might be desirable to explicitly evaluate a sub-expression within a larger
// expression. However, please note that \c evaluate() is not intended to be used for this
// purpose. This task is more elegantly and efficiently handled by the \c eval() function:

   \code
   blaze::DynamicMatrix<double> A, B, C, D;

   D = A + evaluate( B * C );  // Unnecessary creation of a temporary matrix
   D = A + eval( B * C );      // No creation of a temporary matrix
   \endcode

// In contrast to the \c evaluate() function, \c eval() can take the complete expression
// into account and therefore can guarantee the most efficient way to evaluate it (see also
// \ref intra_statement_optimization).
//
// \n \subsection matrix_operations_noalias noalias()
//
// The \b Blaze library is able to reliably detect aliasing during the assignment of matrices.
// In case the aliasing would lead to an incorrect result, \b Blaze introduces an intermediate
// temporary of the appropriate type to break the aliasing. For instance, in the following
// example \b Blaze performs an alias detection in both assignments, but only, in the second
// assignment it detects a problematic aliasing and uses an intermediate temporary in order
// to be able to compute the correct result:

   \code
   blaze::DynamicMatrix<double> A, B;

   A = A + B;  // No problematic aliasing of A, no intermediate temporary is required.
   A = A * B;  // Problematic aliasing of A; intermediate temporary required!
   \endcode

// The detection of aliasing effects, however, takes a small runtime effort. In order to disable
// the aliasing detection, the \c noalias() function can be used:

   \code
   blaze::DynamicMatrix<double> A, B;

   A = noalias( A + B );  // No alias detection performed, no intermediate temporary.
   A = noalias( A * B );  // No alias detection performed, no intermediate temporary.
                          // Note that the final result will be incorrect!
   \endcode

// \warning The \c noalias() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Using \c noalias() in a situation
// where an aliasing effect occurs leads to undefined behavior (which can be violated invariants
// or wrong computation results)!
//
// \n \subsection matrix_operations_nosimd nosimd()
//
// By default, \b Blaze attempts to vectorize all operations by means of SSE, AVX, etc. in order
// to achieve maximum performance. However, via the \c nosimd() operation it is possible to disable
// the SIMD evaluation of any operation:

   \code
   blaze::DynamicMatrix<double> A, B;

   A = nosimd( A + B );  // Disables SIMD for the matrix/matrix addition
   A = nosimd( A * B );  // Disables SIMD for the matrix/matrix multiplication
   \endcode

// Please note that the main purpose of the \c nosimd() operation is to enable an easy performance
// comparison between the vectorized and non-vectorized evaluation. Using the \c nosimd() operation
// will likely result in significantly reduced performance!
//
//
// \n \section matrix_operations_modifying_operations Modifying Operations
// <hr>
//
// \subsection matrix_operations_resize_reserve .resize() / .reserve()
//
// The dimensions of a \c StaticMatrix are fixed at compile time by the second and third template
// parameter and a \c CustomMatrix cannot be resized. In contrast, the number or rows and columns
// of \c DynamicMatrix, \c HybridMatrix, and \c CompressedMatrix can be changed at runtime:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   DynamicMatrix<int,rowMajor> M1;
   CompressedMatrix<int,columnMajor> M2( 3UL, 2UL );

   // Adapting the number of rows and columns via the resize() function. The (optional)
   // third parameter specifies whether the existing elements should be preserved. Per
   // default, the existing elements are preserved.
   M1.resize( 2UL, 2UL );         // Resizing matrix M1 to 2x2 elements. Elements of built-in type
                                  // remain uninitialized, elements of class type are default
                                  // constructed.
   M1.resize( 3UL, 1UL, false );  // Resizing M1 to 3x1 elements. The old elements are lost, the
                                  // new elements are NOT initialized!
   M2.resize( 5UL, 7UL, true );   // Resizing M2 to 5x7 elements. The old elements are preserved.
   M2.resize( 3UL, 2UL, false );  // Resizing M2 to 3x2 elements. The old elements are lost.
   \endcode

// Note that resizing a matrix invalidates all existing views (see e.g. \ref views_submatrices)
// on the matrix:

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 10UL, 20UL );  // Creating a 10x20 matrix
   auto row8 = row( M1, 8UL );  // Creating a view on the 8th row of the matrix
   M1.resize( 6UL, 20UL );      // Resizing the matrix invalidates the view
   \endcode

// When the internal capacity of a matrix is no longer sufficient, the allocation of a larger
// junk of memory is triggered. In order to avoid frequent reallocations, the \c reserve()
// function can be used up front to set the internal capacity:

   \code
   blaze::DynamicMatrix<int> M1;
   M1.reserve( 100 );
   M1.rows();      // Returns 0
   M1.capacity();  // Returns at least 100
   \endcode

// Additionally it is possible to reserve memory in a specific row (for a row-major matrix) or
// column (for a column-major matrix):

   \code
   blaze::CompressedMatrix<int> M1( 4UL, 6UL );
   M1.reserve( 1, 4 );  // Reserving enough space for four non-zero elements in row 1
   \endcode

// \n \subsection matrix_operations_shrinkToFit .shrinkToFit()
//
// The internal capacity of matrices with dynamic memory is preserved in order to minimize the
// number of reallocations. For that reason, the \c resize() and \c reserve() functions can lead
// to memory overhead. The \c shrinkToFit() member function can be used to minimize the internal
// capacity:

   \code
   blaze::DynamicMatrix<int> M1( 100UL, 100UL );  // Create a 100x100 integer matrix
   M1.resize( 10UL, 10UL );                       // Resize to 10x10, but the capacity is preserved
   M1.shrinkToFit();                              // Remove the unused capacity
   \endcode

// Please note that due to padding the capacity might not be reduced exactly to \c rows() times
// \c columns(). Please also note that in case a reallocation occurs, all iterators (including
// \c end() iterators), all pointers and references to elements of this matrix are invalidated.
//
//
// \subsection matrix_operations_reset_clear reset() / clear
//
// In order to reset all elements of a dense or sparse matrix, the \c reset() function can be
// used. The number of rows and columns of the matrix are preserved:

   \code
   // Setting up a single precision row-major matrix, whose elements are initialized with 2.0F.
   blaze::DynamicMatrix<float> M1( 4UL, 5UL, 2.0F );

   // Resetting all elements to 0.0F.
   reset( M1 );  // Resetting all elements
   M1.rows();    // Returns 4: size and capacity remain unchanged
   \endcode

// Alternatively, only a single row or column of the matrix can be resetted:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor>    M1( 7UL, 6UL, 5 );  // Setup of a row-major matrix
   blaze::DynamicMatrix<int,blaze::columnMajor> M2( 4UL, 5UL, 4 );  // Setup of a column-major matrix

   reset( M1, 2UL );  // Resetting the 2nd row of the row-major matrix
   reset( M2, 3UL );  // Resetting the 3rd column of the column-major matrix
   \endcode

// In order to reset a row of a column-major matrix or a column of a row-major matrix, use a
// row or column view (see \ref views_rows and views_colums).
//
// In order to return a matrix to its default state (i.e. the state of a default constructed
// matrix), the \c clear() function can be used:

   \code
   // Setting up a single precision row-major matrix, whose elements are initialized with 2.0F.
   blaze::DynamicMatrix<float> M1( 4UL, 5UL, 2.0F );

   // Resetting all elements to 0.0F.
   clear( M1 );  // Resetting the entire matrix
   M1.rows();    // Returns 0: size is reset, but capacity remains unchanged
   \endcode

// \n \subsection matrix_operations_matrix_transpose transpose()
//
// In addition to the non-modifying \c trans() function, matrices can be transposed in-place via
// the \c transpose() function:

   \code
   blaze::DynamicMatrix<int,rowMajor> M( 5UL, 2UL );

   transpose( M );  // In-place transpose operation.
   M = trans( M );  // Same as above
   \endcode

// Note however that the transpose operation fails if ...
//
//  - ... the given matrix has a fixed size and is non-square;
//  - ... the given matrix is a triangular matrix;
//  - ... the given submatrix affects the restricted parts of a triangular matrix;
//  - ... the given submatrix would cause non-deterministic results in a symmetric/Hermitian matrix.
//
//
// \n \subsection matrix_operations_ctranspose ctranspose()
//
// The \c ctranspose() function can be used to perform an in-place conjugate transpose operation:

   \code
   blaze::DynamicMatrix<int,rowMajor> M( 5UL, 2UL );

   ctranspose( M );  // In-place conjugate transpose operation.
   M = ctrans( M );  // Same as above
   \endcode

// Note however that the conjugate transpose operation fails if ...
//
//  - ... the given matrix has a fixed size and is non-square;
//  - ... the given matrix is a triangular matrix;
//  - ... the given submatrix affects the restricted parts of a triangular matrix;
//  - ... the given submatrix would cause non-deterministic results in a symmetric/Hermitian matrix.
//
//
// \n \subsection matrix_operations_swap swap()
//
// Via the \c \c swap() function it is possible to completely swap the contents of two matrices
// of the same type:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> M1( 10UL, 15UL );
   blaze::DynamicMatrix<int,blaze::rowMajor> M2( 20UL, 10UL );

   swap( M1, M2 );  // Swapping the contents of M1 and M2
   \endcode

// \n \section matrix_operations_arithmetic_operations Arithmetic Operations
// <hr>
//
// \subsection matrix_operations_min_max min() / max()
//
// The \c min() and \c max() functions can be used for a single matrix, multiple matrices, and
// a matrix and a scalar.
//
// <b>Single Matrix</b>
//
// If passed a single matrix, the functions return the smallest and largest element of the given
// dense matrix or the smallest and largest non-zero element of the given sparse matrix,
// respectively:

   \code
   blaze::StaticMatrix<int,2UL,3UL> A{ { -5, 2, 7 },
                                       { -4, 0, 1 } };

   min( A );  // Returns -5
   max( A );  // Returns 7
   \endcode

   \code
   blaze::CompressedMatrix<int> B{ { 1, 0, 3 },
                                   { 0, 0, 0 } };

   min( B );  // Returns 1
   max( B );  // Returns 3
   \endcode

// For more information on the unary \c min() and \c max() reduction operations see the
// \ref matrix_operations_reduction_operations section.
//
// <b>Multiple Matrices</b>
//
// If passed two or more dense matrices, the \c min() and \c max() functions compute the
// componentwise minimum or maximum of the given matrices, respectively:

   \code
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> C{ { -5, 1, -7 }, { 4, 1, 0 } };
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> D{ { -5, 3,  0 }, { 2, 2, -2 } };

   min( A, C );     // Results in the matrix ( -5, 1, -7 ) ( -4, 0, 0 )
   max( A, C, D );  // Results in the matrix ( -5, 3, 7 ) (  4, 2, 1 )
   \endcode

// Please note that sparse matrices can only be used in the unary \c min() and \c max() functions.
// Also note that all forms of the \c min() and \c max() functions can be used to compute the
// smallest and largest element of a matrix expression:

   \code
   min( A + B + C );  // Returns -9, i.e. the smallest value of the resulting matrix
   max( A - B - C );  // Returns 11, i.e. the largest value of the resulting matrix
   \endcode

// <b>Matrix and Scalar</b>
//
// If passed a dense matrix and a scalar, the \c min() and \c max() functions compute the
// componentwise minimum or maximum between the given matrix and a uniform matrix represented by
// the scalar value:

   \code
   min( A, 0 );  // Results in the matrix ( 0, 2, 7 ) ( 0, 0, 1 )
   min( 0, A );  // Results in the matrix ( 0, 2, 7 ) ( 0, 0, 1 )
   max( A, 0 );  // Results in the matrix ( -5, 0, 0 ) ( -4, 0, 0 )
   max( 0, A );  // Results in the matrix ( -5, 0, 0 ) ( -4, 0, 0 )
   \endcode

// \n \subsection matrix_operators_softmax softmax()
//
// The <a href="https://en.wikipedia.org/wiki/Softmax_function">softmax function</a>, also called
// the normalized exponential function, of a given dense matrix can be computed via \c softmax().
// The resulting dense matrix consists of real values in the range (0..1], which add up to 1.

   \code
   blaze::StaticMatrix<double,3UL,3UL> A{ { 1.0, 2.0, 3.0 }
                                        , { 4.0, 1.0, 2.0 }
                                        , { 3.0, 4.0, 1.0 } };
   blaze::StaticMatrix<double,3UL,3UL> B;

   // Evaluating the softmax function
   B = softmax( A );  // Results in ( 0.0157764  0.0428847  0.116573  )
                      //            ( 0.316878   0.0157764  0.0428847 )
                      //            ( 0.116573   0.316878   0.0157764 )

   double b = sum( B );  // Results in 1
   \endcode

// Alternatively it is possible to compute a row- or columnwise \c softmax() function. The
// resulting dense matrix consists of real values in the range (0..1], which add up to the number
// of rows or columns, respectively.

   \code
   using blaze::rowwise;
   using blaze::columnwise;

   blaze::StaticMatrix<double,3UL,3UL> C, D;

   // Evaluating the rowwise softmax function
   C = softmax<rowwise>( A );  // Results in ( 0.0900306  0.244728   0.665241 )
                               //            ( 0.843795   0.0420101  0.114195 )
                               //            ( 0.259496   0.705385   0.035119 )

   double c = sum( C );  // Results in 3 (the number of rows of A)

   // Evaluating the columnwise softmax function
   D = softmax<columnwise>( A );  // Results in ( 0.035119  0.114195   0.665241  )
                                  //            ( 0.705385  0.0420101  0.244728  )
                                  //            ( 0.259496  0.843795   0.0900306 )

   double d = sum( D );  // Results in 3 (the number of columns of A)
   \endcode

// \n \subsection matrix_operators_trace trace()
//
// The \c trace() function sums the diagonal elements of a square dense or sparse matrix:

   \code
   blaze::StaticMatrix<int,3UL,3UL> A{ { -1,  2, -3 }
                                     , { -4, -5,  6 }
                                     , {  7, -8, -9 } };

   trace( A );  // Returns the sum of the diagonal elements, i.e. -15
   \endcode

// In case the given matrix is not a square matrix, a \c std::invalid_argument exception is
// thrown.
//
//
// \n \subsection matrix_operations_matrix_determinant det()
//
// The determinant of a square dense matrix can be computed by means of the \c det() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization
   double d = det( A );  // Compute the determinant of A
   \endcode

// In case the given dense matrix is not a square matrix, a \c std::invalid_argument exception is
// thrown.
//
// \note The \c det() function can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The function is depending on LAPACK kernels. Thus the function can only be used if a
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error
// will be created.
//
//
// \n \subsection matrix_operators_rank rank()
//
// The \c rank() function computes the rank of a given dense matrix:

   \code
   blaze::DynamicMatrix<double> A( 5UL, 8UL );
   // ... Initialization
   rank( A );
   \endcode

// The rank is determined as the number of singular values greater than a given tolerance. This
// tolerance is computed as

   \code
   tolerance = max(m,n) * max(s) * epsilon,
   \endcode

// where \c m is the number of rows of the dense matrix, \c n is the number of columns of the
// dense matrix, \c max(s) is the maximum singular value of the dense matrix and \c epsilon is
// the difference between 1 and the least value greater than 1 that is representable by the
// floating point type of the singular values.
//
// \note The \c rank() function can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The function is depending on LAPACK kernels. Thus the function can only be used if a
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error
// will be created.
//
//
// \n \subsection matrix_operators_abs abs()
//
// The \c abs() function can be used to compute the absolute values of each element of a matrix.
// For instance, the following computation

   \code
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> A{ { -1,  2, -3 },
                                                {  4, -5,  6 } };
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> B( abs( A ) );
   \endcode

// results in the matrix

                          \f$ B = \left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          4 & 5 & 6 \\
                          \end{array}\right)\f$

// \n \subsection matrix_operators_sign sign()
//
// The \c sign() function can be used to evaluate the sign of each element of a matrix \a A. For
// each element \c (i,j) the corresponding result is 1 if \a A(i,j) is greater than zero, 0 if
// \a A(i,j) is zero, and -1 if \a A(i,j) is less than zero. For instance, the following use of
// the \c sign() function

   \code
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> A{ { -1,  2,  0 },
                                                {  4,  0, -6 } };
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> B( sign( A ) );
   \endcode

// results in the matrix

                          \f$ B = \left(\begin{array}{*{3}{c}}
                          -1 & 1 &  0 \\
                           1 & 0 & -1 \\
                          \end{array}\right)\f$

// \n \subsection matrix_operators_rounding_functions floor() / ceil() / trunc() / round()
//
// The \c floor(), \c ceil(), \c trunc(), and \c round() functions can be used to round down/up
// each element of a matrix, respectively:

   \code
   blaze::StaticMatrix<double,3UL,3UL> A, B;

   B = floor( A );  // Rounding down each element of the matrix
   B = ceil ( A );  // Rounding up each element of the matrix
   B = trunc( A );  // Truncating each element of the matrix
   B = round( A );  // Rounding each element of the matrix
   \endcode

// \n \subsection matrix_operators_conj conj()
//
// The \c conj() function can be applied on a dense or sparse matrix to compute the complex
// conjugate of each element of the matrix:

   \code
   using blaze::StaticMatrix;

   using cplx = std::complex<double>;

   // Creating the matrix
   //    ( (1,0)  (-2,-1) )
   //    ( (1,1)  ( 0, 1) )
   StaticMatrix<cplx,2UL,2UL> A{ { cplx( 1.0, 0.0 ), cplx( -2.0, -1.0 ) },
                                 { cplx( 1.0, 1.0 ), cplx(  0.0,  1.0 ) } };

   // Computing the matrix of conjugate values
   //    ( (1, 0)  (-2, 1) )
   //    ( (1,-1)  ( 0,-1) )
   StaticMatrix<cplx,2UL,2UL> B;
   B = conj( A );
   \endcode

// Additionally, matrices can be conjugated in-place via the \c conjugate() function:

   \code
   blaze::DynamicMatrix<cplx> C( 5UL, 2UL );

   conjugate( C );  // In-place conjugate operation.
   C = conj( C );   // Same as above
   \endcode

// \n \subsection matrix_operators_real real()
//
// The \c real() function can be used on a dense or sparse matrix to extract the real part of
// each element of the matrix:

   \code
   using blaze::StaticMatrix;

   using cplx = std::complex<double>;

   // Creating the matrix
   //    ( (1,0)  (-2,-1) )
   //    ( (1,1)  ( 0, 1) )
   StaticMatrix<cplx,2UL,2UL> A{ { cplx( 1.0, 0.0 ), cplx( -2.0, -1.0 ) },
                                 { cplx( 1.0, 1.0 ), cplx(  0.0,  1.0 ) } };

   // Extracting the real part of each matrix element
   //    ( 1 -2 )
   //    ( 1  0 )
   StaticMatrix<double,2UL,2UL> B;
   B = real( A );
   \endcode

// \n \subsection matrix_operators_imag imag()
//
// The \c imag() function can be used on a dense or sparse matrix to extract the imaginary part
// of each element of the matrix:

   \code
   using blaze::StaticMatrix;

   using cplx = std::complex<double>;

   // Creating the matrix
   //    ( (1,0)  (-2,-1) )
   //    ( (1,1)  ( 0, 1) )
   StaticMatrix<cplx,2UL,2UL> A{ { cplx( 1.0, 0.0 ), cplx( -2.0, -1.0 ) },
                                 { cplx( 1.0, 1.0 ), cplx(  0.0,  1.0 ) } };

   // Extracting the imaginary part of each matrix element
   //    ( 0 -1 )
   //    ( 1  1 )
   StaticMatrix<double,2UL,2UL> B;
   B = imag( A );
   \endcode

// \n \subsection matrix_operators_arg arg()
//
// The \c arg() function can be used on a dense or sparse matrix to compute the phase angle for
// each element of the matrix:

   \code
   using blaze::StaticMatrix;

   using cplx = std::complex<double>;

   // Creating the matrix
   //    ( (1,0)  (-2,-1) )
   //    ( (1,1)  ( 0, 1) )
   StaticMatrix<cplx,2UL,2UL> A{ { cplx( 1.0, 0.0 ), cplx( -2.0, -1.0 ) },
                                 { cplx( 1.0, 1.0 ), cplx(  0.0,  1.0 ) } };

   // Computing the phase angle of each matrix element
   //    ( 0.0      -2.67795 )
   //    ( 0.785398  1.5708  )
   StaticMatrix<double,2UL,2UL> B;
   B = arg( A );
   \endcode

// \n \subsection matrix_operators_sqrt sqrt() / invsqrt()
//
// Via the \c sqrt() and \c invsqrt() functions the (inverse) square root of each element of a
// matrix can be computed:

   \code
   blaze::StaticMatrix<double,3UL,3UL> A, B, C;

   B = sqrt( A );     // Computes the square root of each element
   C = invsqrt( A );  // Computes the inverse square root of each element
   \endcode

// Note that in case of sparse matrices only the non-zero elements are taken into account!
//
//
// \n \subsection matrix_operators_cbrt cbrt() / invcbrt()
//
// The \c cbrt() and \c invcbrt() functions can be used to compute the the (inverse) cubic root
// of each element of a matrix:

   \code
   blaze::DynamicMatrix<double> A, B, C;

   B = cbrt( A );     // Computes the cubic root of each element
   C = invcbrt( A );  // Computes the inverse cubic root of each element
   \endcode

// Note that in case of sparse matrices only the non-zero elements are taken into account!
//
//
// \n \subsection matrix_operations_hypot hypot()
//
// The \c hypot() function can be used to compute the componentwise hypotenous for a pair of
// dense matrices:

   \code
   blaze::StaticMatrix<double,3UL,3UL> A, B, C;

   C = hypot( A, B );  // Computes the componentwise hypotenuous
   \endcode

// \n \subsection matrix_operators_clamp clamp()
//
// The \c clamp() function can be used to restrict all elements of a matrix to a specific range:

   \code
   blaze::DynamicMatrix<double> A, B;

   B = clamp( A, -1.0, 1.0 );  // Restrict all elements to the range [-1..1]
   \endcode

// Note that in case of sparse matrices only the non-zero elements are taken into account!
//
//
// \n \subsection matrix_operators_pow pow()
//
// The \c pow() function can be used to compute the exponential value of each element of a matrix.
// If passed a matrix and a numeric exponent, the function computes the exponential value of each
// element of the matrix using the same exponent. If passed a second matrix, the function computes
// the componentwise exponential value:

   \code
   blaze::StaticMatrix<double,3UL,3UL> A, B, C;

   C = pow( A, 1.2 );  // Computes the exponential value of each element
   C = pow( A, B );    // Computes the componentwise exponential value
   \endcode

// \n \subsection matrix_operators_exp exp() / exp2() / exp10()
//
// \c exp(), \c exp2() and \c exp10() compute the base e/2/10 exponential of each element of a
// matrix, respectively:

   \code
   blaze::HybridMatrix<double,3UL,3UL> A, B;

   B = exp( A );    // Computes the base e exponential of each element
   B = exp2( A );   // Computes the base 2 exponential of each element
   B = exp10( A );  // Computes the base 10 exponential of each element
   \endcode

// Note that in case of sparse matrices only the non-zero elements are taken into account!
//
//
// \n \subsection matrix_operators_log log() / log2() / log10() / log1p() / lgamma()
//
// The \c log(), \c log2(), \c log10(), \c log1p() and \c lgamma() functions can be used to
// compute the natural, binary and common logarithm of each element of a matrix:

   \code
   blaze::StaticMatrix<double,3UL,3UL> A, B;

   B = log( A );     // Computes the natural logarithm of each element
   B = log2( A );    // Computes the binary logarithm of each element
   B = log10( A );   // Computes the common logarithm of each element
   B = log1p( A );   // Computes the natural logarithm of x+1 of each element
   B = lgamma( A );  // Computes the natural logarithm of the absolute value of the gamma function
   \endcode

// \n \subsection matrix_operators_trigonometric_functions sin() / cos() / tan() / asin() / acos() / atan()
//
// The following trigonometric functions are available for both dense and sparse matrices:

   \code
   blaze::DynamicMatrix<double> A, B;

   B = sin( A );  // Computes the sine of each element of the matrix
   B = cos( A );  // Computes the cosine of each element of the matrix
   B = tan( A );  // Computes the tangent of each element of the matrix

   B = asin( A );  // Computes the inverse sine of each element of the matrix
   B = acos( A );  // Computes the inverse cosine of each element of the matrix
   B = atan( A );  // Computes the inverse tangent of each element of the matrix
   \endcode

// Note that in case of sparse matrices only the non-zero elements are taken into account!
//
//
// \n \subsection matrix_operators_hyperbolic_functions sinh() / cosh() / tanh() / asinh() / acosh() / atanh()
//
// The following hyperbolic functions are available for both dense and sparse matrices:

   \code
   blaze::DynamicMatrix<double> A, B;

   B = sinh( A );  // Computes the hyperbolic sine of each element of the matrix
   B = cosh( A );  // Computes the hyperbolic cosine of each element of the matrix
   B = tanh( A );  // Computes the hyperbolic tangent of each element of the matrix

   B = asinh( A );  // Computes the inverse hyperbolic sine of each element of the matrix
   B = acosh( A );  // Computes the inverse hyperbolic cosine of each element of the matrix
   B = atanh( A );  // Computes the inverse hyperbolic tangent of each element of the matrix
   \endcode

// \n \subsection matrix_operations_atan2 atan2()
//
// The multi-valued inverse tangent is available for a pair of dense matrices:

   \code
   blaze::DynamicMatrix<double> A, B, C;

   C = atan2( A, B );  // Computes the componentwise multi-valued inverse tangent
   \endcode

// \n \subsection matrix_operators_erf erf() / erfc()
//
// The \c erf() and \c erfc() functions compute the (complementary) error function of each
// element of a matrix:

   \code
   blaze::StaticMatrix<double,3UL,3UL> A, B;

   B = erf( A );   // Computes the error function of each element
   B = erfc( A );  // Computes the complementary error function of each element
   \endcode

// Note that in case of sparse matrices only the non-zero elements are taken into account!
//
//
// \n \subsection matrix_operations_map map() / forEach()
//
// Via the \c map() functions it is possible to execute componentwise custom operations on matrices.
// The unary \c map() function can be used to apply a custom operation on each element of a
// dense or sparse matrix. For instance, the following example demonstrates a custom square root
// computation via a lambda:

   \code
   blaze::DynamicMatrix<double> A, B;

   B = map( A, []( double d ) { return std::sqrt( d ); } );
   \endcode

// The N-ary \c map() functions can be used to apply an operation componentwise to the elements
// of N dense matrices (where \f$ N <= 6 \f$). The following example demonstrates the merging of
// two matrices of double precision values into a matrix of double precision complex numbers:

   \code
   blaze::DynamicMatrix<double> real{ { 2.1, -4.2 }, { 1.0,  0.6 } };
   blaze::DynamicMatrix<double> imag{ { 0.3,  1.4 }, { 2.9, -3.4 } };

   blaze::DynamicMatrix< complex<double> > cplx;

   // Creating the matrix
   //    ( ( 2.1,  0.3) (-4.2,  1.4) )
   //    ( ( 1.0,  2.9) ( 0.6, -3.4) )
   cplx = map( real, imag, []( double r, double i ){ return complex<double>( r, i ); } );
   \endcode

// Although the computation can be parallelized it is not vectorized and thus cannot perform at
// peak performance. However, it is also possible to create vectorized custom operations. See
// \ref custom_operations for a detailed overview of the possibilities of custom operations.
//
// Please note that unary custom operations on vectors have been introduced in \b Blaze 3.0 in
// form of the \c forEach() function. With the introduction of binary custom functions, the
// \c forEach() function has been renamed to \c map(). The \c forEach() function can still be
// used, but the function might be deprecated in future releases of \b Blaze.
//
//
// \n \subsection matrix_operations_select select()
//
// The \c select() function performs a componentwise, conditional selection of elements. Given
// the three dense matrices \c cond, \c A, and \c B, in case an element in the \c cond vector
// evaluates to \c true, the according element of \a A is selected, in case the \a cond element
// evaluates to \c false, the according element of \a B is selected. The following example
// demonstrates the use of the \a select() function:

   \code
   blaze::DynamicMatrix<bool> cond{ { true, false }, { true false } };
   blaze::DynamicMatrix<int> A{ { 1, -1 }, { 1, -1 } };
   blaze::DynamicMatrix<int> B{ { -2, 2 }, { -2, 2 } };
   blaze::DynamicMatrix<int> C;
   // ... Resizing and initialization

   C = select( cond, A, B );  // Results in ( 1, 2 ) ( 1, 2 )
   \endcode

// \n \section matrix_operations_reduction_operations Reduction Operations
// <hr>
//
// \subsection matrix_operations_reduction_operations_reduce reduce()
//
// The \c reduce() function performs either a total reduction, a rowwise reduction or a columnwise
// reduction of the elements of the given dense matrix or the non-zero elements of the given sparse
// matrix. The following examples demonstrate the total reduction of a dense and sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   // ... Resizing and initialization

   const double totalsum1 = reduce( A, blaze::Add() );
   const double totalsum2 = reduce( A, []( double a, double b ){ return a + b; } );
   \endcode

   \code
   blaze::CompressedMatrix<double> A;
   // ... Resizing and initialization

   const double totalsum1 = reduce( A, blaze::Add() );
   const double totalsum2 = reduce( A, []( double a, double b ){ return a + b; } );
   \endcode

// By specifying \c blaze::columnwise or \c blaze::rowwise the \c reduce() function performs a
// column-wise or row-wise reduction, respectively. In case \c blaze::columnwise is specified, the
// (non-zero) elements of the matrix are reduced column-wise and the result is a row vector. In
// case \c blaze::rowwise is specified, the (non-zero) elements of the matrix are reduced row-wise
// and the result is a column vector:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   blaze::DynamicVector<double,rowVector> colsum1, colsum2;
   // ... Resizing and initialization

   colsum1 = reduce<columnwise>( A, blaze::Add() );
   colsum2 = reduce<columnwise>( B, []( double a, double b ){ return a + b; } );
   \endcode

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   blaze::DynamicVector<double,columnVector> rowsum1, rowsum2;
   // ... Resizing and initialization

   rowsum1 = reduce<rowwise>( A, blaze::Add() );
   rowsum2 = reduce<rowwise>( B, []( double a, double b ){ return a + b; } );
   \endcode

// As demonstrated in the examples it is possible to pass any binary callable as custom reduction
// operation. However, for instance in the case of lambdas the vectorization of the reduction
// operation is compiler dependent and might not perform at peak performance. However, it is also
// possible to create vectorized custom operations. See \ref custom_operations for a detailed
// overview of the possibilities of custom operations.
//
// Please note that the evaluation order of the \c reduce() function is unspecified. Thus the
// behavior is non-deterministic if the given reduction operation is not associative or not
// commutative. Also, the operation is undefined if the given reduction operation modifies the
// values.
//
// \n \subsection matrix_operations_reduction_operations_sum sum()
//
// The \c sum() function reduces the elements of the given dense vector or the non-zero elements
// of the given sparse vector by means of addition:

   \code
   blaze::DynamicMatrix<int> A{ { 1, 2 }, { 3, 4 } };

   const int totalsum = sum( A );  // Results in 10
   \endcode

   \code
   blaze::CompressedMatrix<int> a{ { 1, 2 }, { 3, 4 } };

   const int totalsum = sum( A );  // Results in 10
   \endcode

// By specifying \c blaze::columnwise or \c blaze::rowwise the \c sum() function performs a
// column-wise or row-wise summation, respectively. In case \c blaze::columnwise is specified,
// the (non-zero) elements of the matrix are summed up column-wise and the result is a row vector.
// In case \c blaze::rowwise is specified, the (non-zero) elements of the matrix are summed up
// row-wise and the result is a column vector:

   \code
   using blaze::columnwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::DynamicVector<int,rowVector> colsum1, colsum2;

   colsum1 = sum<columnwise>( A );  // Results in ( 2, 3, 6 )
   colsum2 = sum<columnwise>( B );  // Same result
   \endcode

   \code
   using blaze::rowwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::DynamicVector<int,columnVector> rowsum1, rowsum2;

   rowsum1 = sum<rowwise>( A );  // Results in ( 3, 8 )
   rowsum2 = sum<rowwise>( B );  // Same result
   \endcode

// Please note that the evaluation order of the \c sum() function is unspecified.
//
// \n \subsection matrix_operations_reduction_operations_prod prod()
//
// The \c prod() function reduces the elements of the given dense vector or the non-zero elements
// of the given sparse vector by means of multiplication:

   \code
   blaze::DynamicMatrix<int> A{ { 1, 2 }, { 3, 4 } };

   const int totalprod = prod( A );  // Results in 24
   \endcode

   \code
   blaze::CompressedMatrix<int> A{ { 1, 2 }, { 3, 4 } };

   const int totalprod = prod( A );  // Results in 24
   \endcode

// By specifying \c blaze::columnwise or \c blaze::rowwise the \c prod() function performs a
// column-wise or row-wise multiplication, respectively. In case \c blaze::columnwise is specified,
// the (non-zero) elements of the matrix are multiplied column-wise and the result is a row vector.
// In case \c blaze::rowwise is specified, the (non-zero) elements of the matrix are multiplied
// row-wise and the result is a column vector:

   \code
   using blaze::columnwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::DynamicVector<int,rowVector> colprod1, colprod2;

   colprod1 = prod<columnwise>( A );  // Results in ( 1, 0, 8 )
   colprod2 = prod<columnwise>( A );  // Results in ( 1, 3, 8 )
   \endcode

   \code
   using blaze::rowwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::DynamicVector<int,columnVector> rowprod1, rowprod2;

   rowprod1 = prod<rowwise>( A );  // Results in ( 0, 12 )
   rowprod2 = prod<rowwise>( A );  // Results in ( 2, 12 )
   \endcode

// Please note that the evaluation order of the \c prod() function is unspecified.
//
// \n \subsection matrix_operations_reduction_operations_min min()
//
// The unary \c min() function returns the smallest element of the given dense matrix or the
// smallest non-zero element of the given sparse matrix. This function can only be used for
// element types that support the smaller-than relationship. In case the given matrix currently
// has either 0 rows or 0 columns, the returned value is the default value (e.g. 0 in case of
// fundamental data types).

   \code
   blaze::DynamicMatrix<int> A{ { 1, 2 }, { 3, 4 } };

   const int totalmin = min( A );  // Results in 1
   \endcode

   \code
   blaze::CompressedMatrix<int> A{ { 1, 0 }, { 3, 0 } };

   const int totalmin = min( A );  // Results in 1
   \endcode

// \note In case the sparse matrix is not completely filled, the implicit zero elements are NOT
// taken into account. In the previous example the compressed matrix has only 2 non-zero elements.
// However, the minimum of this matrix is 1.
//
// By specifying \c blaze::columnwise or \c blaze::rowwise the \c min() function determines the
// smallest (non-zero) element in each row or column, respectively. In case \c blaze::columnwise
// is specified, the smallest (non-zero) element of each column is determined and the result is
// a row vector. In case \c blaze::rowwise is specified, the smallest (non-zero) element of each
// row is determined and the result is a column vector.

   \code
   using blaze::columnwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::DynamicVector<int,rowVector> colmin1, colmin2;

   colmin1 = min<columnwise>( A );  // Results in ( 1, 0, 2 )
   colmin2 = min<columnwise>( B );  // Results in ( 1, 3, 2 )
   \endcode

   \code
   using blaze::rowwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::DynamicVector<int,columnVector> rowmin1, rowmin2;

   rowmin1 = min<rowwise>( A );  // Results in ( 0, 1 )
   rowmin2 = min<rowwise>( B );  // Results in ( 1, 1 )
   \endcode

// \note In case the sparse matrix is not completely filled, the implicit zero elements are NOT
// taken into account.
//
// \n \subsection matrix_operations_reduction_operations_max max()
//
// The unary \c max() function returns the largest element of the given dense matrix or the
// largest non-zero element of the given sparse matrix. This function can only be used for
// element types that support the smaller-than relationship. In case the given matrix currently
// has either 0 rows or 0 columns, the returned value is the default value (e.g. 0 in case of
// fundamental data types).

   \code
   blaze::DynamicMatrix<int> A{ { 1, 2 }, { 3, 4 } };

   const int totalmax = max( A );  // Results in 4
   \endcode

   \code
   blaze::CompressedMatrix<int> A{ { -1, 0 }, { -3, 0 } };

   const int totalmax = max( A );  // Results in -1
   \endcode

// \note In case the sparse matrix is not completely filled, the implicit zero elements are NOT
// taken into account. In the previous example the compressed matrix has only 2 non-zero elements.
// However, the maximum of this matrix is -1.
//
// By specifying \c blaze::columnwise or \c blaze::rowwise the \c max() function determines the
// largest (non-zero) element in each row or column, respectively. In case \c blaze::columnwise
// is specified, the largest (non-zero) element of each column is determined and the result is
// a row vector. In case \c blaze::rowwise is specified, the largest (non-zero) element of each
// row is determined and the result is a column vector.

   \code
   using blaze::columnwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { -1, 0, -2 }, { -1, -3, -4 } };
   blaze::DynamicVector<int,rowVector> colmax1, colmax2;

   colmax1 = max<columnwise>( A );  // Results in ( 1, 3, 4 )
   colmax2 = max<columnwise>( B );  // Results in ( -1, -3, -2 )
   \endcode

   \code
   using blaze::rowwise;

   blaze::DynamicMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   blaze::CompressedMatrix<int> B{ { -1, 0, -2 }, { -1, -3, -4 } };
   blaze::DynamicVector<int,columnVector> rowmax1, rowmax2;

   rowmax1 = max<rowwise>( A );  // Results in ( 2, 4 )
   rowmax2 = max<rowwise>( B );  // Results in ( -1, -1 )
   \endcode

// \note In case the sparse matrix is not completely filled, the implicit zero elements are NOT
// taken into account.
//
//
// \n \section matrix_operations_norms Norms
// <hr>
//
// \subsection matrix_operations_norms_norm norm()
//
// The \c norm() function computes the L2 norm of the given dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = norm( A );
   const double norm2 = norm( B );
   \endcode

// \n \subsection matrix_operations_norms_sqrnorm sqrNorm()
//
// The \c sqrNorm() function computes the squared L2 norm of the given dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = sqrNorm( A );
   const double norm2 = sqrNorm( B );
   \endcode

// \n \subsection matrix_operations_norms_l1norm l1Norm()
//
// The \c l1Norm() function computes the squared L1 norm of the given dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = l1Norm( A );
   const double norm2 = l1Norm( B );
   \endcode

// \n \subsection matrix_operations_norms_l2norm l2Norm()
//
// The \c l2Norm() function computes the squared L2 norm of the given dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = l2Norm( A );
   const double norm2 = l2Norm( B );
   \endcode

// \n \subsection matrix_operations_norms_l3norm l3Norm()
//
// The \c l3Norm() function computes the squared L3 norm of the given dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = l3Norm( A );
   const double norm2 = l3Norm( B );
   \endcode

// \n \subsection matrix_operations_norms_l4norm l4Norm()
//
// The \c l4Norm() function computes the squared L4 norm of the given dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = l4Norm( A );
   const double norm2 = l4Norm( B );
   \endcode

// \n \subsection matrix_operations_norms_lpnorm lpNorm()
//
// The \c lpNorm() function computes the general Lp norm of the given dense or sparse matrix,
// where the norm is specified by either a compile time or a runtime argument:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = lpNorm<2>( A );    // Compile time argument
   const double norm2 = lpNorm( B, 2.3 );  // Runtime argument
   \endcode

// \n \subsection matrix_operations_norms_maxnorm linfNorm() / maxNorm()
//
// The \c linfNorm() and \c maxNorm() functions compute the infinity/maximum norm of the given
// dense or sparse matrix:

   \code
   blaze::DynamicMatrix<double> A;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization

   const double norm1 = linfNorm( A );
   const double norm2 = maxNorm( B );
   \endcode

// \n \section matrix_operations_scalar_expansion Scalar Expansion
// <hr>
//
// By means of the \c uniform() function it is possible to expand a scalar value into a dense,
// uniform matrix. By default, the resulting uniform matrix is a row-major matrix, but it is
// possible to specify the storage order explicitly:

   \code
   using blaze::rowMajor;

   int scalar = 5;

   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3x5 row-major matrix
   //
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //
   A = uniform( 3UL, 5UL, scalar );
   A = uniform<columnMajor>( 3UL, 5UL, scalar );
   \endcode

// \n \section matrix_operations_matrix_repetition Matrix Repetition
// <hr>
//
// Via the \c repeat() function it is possible to repeat a dense or sparse matrix multiple times
// to represent a larger matrix. Repeating a row-major matrix results in a row-major matrix,
// repeating a column-major matrix results in a column-major matrix. As demonstrated by the
// following examples, \c repeat() can be used with both runtime and compile time parameters:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<int,rowMajor> A1{ { 1, 0, -2 }, { 0, 5, 0 } };
   blaze::CompressedMatrix<int,columnMajor> B1{ { 0, -1 }, { 0, 4 }, { 7, 0 } };

   blaze::DynamicMatrix<int,rowMajor> A2;
   blaze::CompressedMatrix<int,columnMajor> B2;

   // ... Resizing and initialization

   // Repeating the 2x3 dense row-major matrix 'A1' 2x rowwise and 3x columnwise results in
   //
   //    (  1  0 -2  1  0 -2  1  0 -2 )
   //    (  0  5  0  0  5  0  0  5  0 )
   //    (  1  0 -2  1  0 -2  1  0 -2 )
   //    (  0  5  0  0  5  0  0  5  0 )
   //
   A2 = repeat( A1, 2UL, 3UL );
   A2 = repeat<2UL,3UL>( A1 );

   // Repeating the 3x2 sparse column-major matrix 'B1' 2x rowwise and 3x columnwise results in
   //
   //    (  0 -1  0 -1  0 -1 )
   //    (  0  4  0  4  0  4 )
   //    (  7  0  7  0  7  0 )
   //    (  0 -1  0 -1  0 -1 )
   //    (  0  4  0  4  0  4 )
   //    (  7  0  7  0  7  0 )
   //
   B2 = repeat( B1, 2UL, 3UL );
   B2 = repeat<2UL,3UL>( B1 );
   \endcode

// \n \section matrix_operations_statistic_operations Statistic Operations
// <hr>
//
// \subsection matrix_operations_mean mean()
//
// The <a href="https://en.wikipedia.org/wiki/Arithmetic_mean">(arithmetic) mean</a> of a dense or
// sparse matrix can be computed via the \c mean() function. In case of a sparse matrix, both the
// non-zero and zero elements are taken into account. The following example demonstrates the
// computation of the mean of a dense matrix:

   \code
   blaze::DynamicMatrix<int> A{ { 1, 4, 3, 6, 7 }
                              , { 2, 6, 3, 1, 0 } };

   const double m = mean( A );  // Results in 3.3 (i.e. 33/10)
   \endcode

// In case the number of rows or columns of the given matrix is 0, a \c std::invalid_argument is
// thrown.
//
// Alternatively it is possible to compute the row- or columnwise mean:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicMatrix<int> A{ { 1, 4, 3, 6, 7 }
                              , { 2, 6, 3, 1, 0 } };

   blaze::DynamicVector<double,columnVector> rm;
   blaze::DynamicVector<double,rowVector> cm;

   rm = mean<rowwise>( A );     // Results in ( 4.2  2.4 )
   cm = mean<columnwise>( A );  // Results in ( 1.5  5.0  3.0  3.5  3.5 )
   \endcode

// In case the rowwise mean is computed and the number of columns of the given matrix is 0 or
// in case the columnwise mean is computed and the number of rows of the given matrix is 0, a
// \c std::invalid_argument is thrown.
//
// \n \subsection matrix_operations_var var()
//
// The <a href="https://en.wikipedia.org/wiki/Variance">variance</a> of a dense or sparse matrix
// can be computed via the \c var() function. In case of a sparse vector, both the non-zero and
// zero elements are taken into account. The following example demonstrates the computation of
// the variance of a dense matrix:

   \code
   blaze::DynamicMatrix<int> A{ { 1, 3, 2 }
                              , { 2, 6, 4 }
                              , { 9, 6, 3 } };

   const double v = var( A );  // Results in 6.5
   \endcode

// In case the size of the given matrix is smaller than 2, a \c std::invalid_argument is thrown.
//
// Alternatively it is possible to compute the row- or columnwise variance:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicMatrix<int> A{ { 1, 3, 2 }
                              , { 2, 6, 4 }
                              , { 9, 6, 3 } };

   blaze::DynamicVector<double,columnVector> rv;
   blaze::DynamicVector<double,rowVector> cv;

   rv = var<rowwise>( A );     // Results in ( 1  4  9 )
   cv = var<columnwise>( A );  // Results in ( 19  3  1 )
   \endcode

// In case the rowwise varoamce is computed and the number of columns of the given matrix is
// smaller than 2 or in case the columnwise mean is computed and the number of rows of the given
// matrix is smaller than 2, a \c std::invalid_argument is thrown.
//
// \n \subsection matrix_operations_stddev stddev()
//
// The <a href="https://en.wikipedia.org/wiki/Standard_deviation">standard deviation</a> of a
// dense or sparse matrix can be computed via the \c stddev() function. In case of a sparse
// vector, both the non-zero and zero elements are taken into account. The following example
// demonstrates the computation of the standard deviation of a dense matrix:

   \code
   blaze::DynamicMatrix<int> A{ { 1, 3, 2 }
                              , { 2, 6, 4 }
                              , { 9, 6, 3 } };

   const double s = stddev( A );  // Results in sqrt(6.5)
   \endcode

// In case the size of the given matrix is smaller than 2, a \c std::invalid_argument is thrown.
//
// Alternatively it is possible to compute the row- or columnwise standard deviation:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicMatrix<int> A{ { 1, 3, 2 }
                              , { 2, 6, 4 }
                              , { 9, 6, 3 } };

   blaze::DynamicVector<double,columnVector> rs;
   blaze::DynamicVector<double,rowVector> cs;

   rs = stddev<rowwise>( A );     // Results in ( 1  2  3 )
   cs = stddev<columnwise>( A );  // Results in ( sqrt(19)  sqrt(3)  1 )
   \endcode

// In case the rowwise standard deviation is computed and the number of columns of the given
// matrix is smaller than 2 or in case the columnwise mean is computed and the number of rows of
// the given matrix is smaller than 2, a \c std::invalid_argument is thrown.
//
//
// \n \section matrix_operations_declaration_operations Declaration Operations
// <hr>
//
// \subsection matrix_operations_declsym declsym()
//
// The \c declsym() operation can be used to explicitly declare any matrix or matrix expression
// as symmetric:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declsym( A );
   \endcode

// Any matrix or matrix expression that has been declared as symmetric via \c declsym() will
// gain all the benefits of a symmetric matrix, which range from reduced runtime checking to
// a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   DynamicMatrix<double> A, B, C;
   SymmetricMatrix< DynamicMatrix<double> > S;
   // ... Resizing and initialization

   isSymmetric( declsym( A ) );  // Will always return true without runtime effort

   S = declsym( A );  // Omit any runtime check for symmetry

   C = declsym( A * B );  // Declare the result of the matrix multiplication as symmetric,
                          // i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declsym() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-symmetric matrix or
// matrix expression as symmetric via the \c declsym() operation leads to undefined behavior
// (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declherm declherm()
//
// The \c declherm() operation can be used to explicitly declare any matrix or matrix expression
// as Hermitian:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declherm( A );
   \endcode

// Any matrix or matrix expression that has been declared as Hermitian via \c declherm() will
// gain all the benefits of an Hermitian matrix, which range from reduced runtime checking to
// a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;

   DynamicMatrix<double> A, B, C;
   HermitianMatrix< DynamicMatrix<double> > S;
   // ... Resizing and initialization

   isHermitian( declherm( A ) );  // Will always return true without runtime effort

   S = declherm( A );  // Omit any runtime check for Hermitian symmetry

   C = declherm( A * B );  // Declare the result of the matrix multiplication as Hermitian,
                           // i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declherm() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-Hermitian matrix or
// matrix expression as Hermitian via the \c declherm() operation leads to undefined behavior
// (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_decllow decllow()
//
// The \c decllow() operation can be used to explicitly declare any matrix or matrix expression
// as lower triangular:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = decllow( A );
   \endcode

// Any matrix or matrix expression that has been declared as lower triangular via \c decllow()
// will gain all the benefits of a lower triangular matrix, which range from reduced runtime
// checking to a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;

   DynamicMatrix<double> A, B, C;
   LowerMatrix< DynamicMatrix<double> > L;
   // ... Resizing and initialization

   isLower( decllow( A ) );  // Will always return true without runtime effort

   L = decllow( A );  // Omit any runtime check for A being a lower matrix

   C = decllow( A * B );  // Declare the result of the matrix multiplication as lower triangular,
                          // i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c decllow() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-lower matrix or
// matrix expression as lower triangular via the \c decllow() operation leads to undefined
// behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declunilow declunilow()
//
// The \c declunilow() operation can be used to explicitly declare any matrix or matrix expression
// as lower unitriangular:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declunilow( A );
   \endcode

// Any matrix or matrix expression that has been declared as lower unitriangular via \c declunilow()
// will gain all the benefits of a lower unitriangular matrix, which range from reduced runtime
// checking to a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::UniLowerMatrix;

   DynamicMatrix<double> A, B, C;
   UniLowerMatrix< DynamicMatrix<double> > L;
   // ... Resizing and initialization

   isUniLower( declunilow( A ) );  // Will always return true without runtime effort

   L = declunilow( A );  // Omit any runtime check for A being an unilower matrix

   C = declunilow( A * B );  // Declare the result of the matrix multiplication as lower
                             // unitriangular, i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declunilow() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-unilower matrix or
// matrix expression as lower unitriangular via the \c declunilow() operation leads to undefined
// behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declstrlow declstrlow()
//
// The \c declstrlow() operation can be used to explicitly declare any matrix or matrix expression
// as strictly lower triangular:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declstrlow( A );
   \endcode

// Any matrix or matrix expression that has been declared as strictly lower triangular via
// \c declstrlow() will gain all the benefits of a strictly lower triangular matrix, which range
// from reduced runtime checking to a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::StrictlyLowerMatrix;

   DynamicMatrix<double> A, B, C;
   StrictlyLowerMatrix< DynamicMatrix<double> > L;
   // ... Resizing and initialization

   isStrictlyLower( declstrlow( A ) );  // Will always return true without runtime effort

   L = declstrlow( A );  // Omit any runtime check for A being a strictly lower matrix

   C = declstrlow( A * B );  // Declare the result of the matrix multiplication as strictly lower
                             // triangular, i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declstrlow() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-strictly-lower matrix
// or matrix expression as strictly lower triangular via the \c declstrlow() operation leads to
// undefined behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declupp declupp()
//
// The \c declupp() operation can be used to explicitly declare any matrix or matrix expression
// as upper triangular:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declupp( A );
   \endcode

// Any matrix or matrix expression that has been declared as upper triangular via \c declupp()
// will gain all the benefits of an upper triangular matrix, which range from reduced runtime
// checking to a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::UpperMatrix;

   DynamicMatrix<double> A, B, C;
   UpperMatrix< DynamicMatrix<double> > U;
   // ... Resizing and initialization

   isUpper( declupp( A ) );  // Will always return true without runtime effort

   U = declupp( A );  // Omit any runtime check for A being an upper matrix

   C = declupp( A * B );  // Declare the result of the matrix multiplication as upper triangular,
                          // i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declupp() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-upper matrix or
// matrix expression as upper triangular via the \c declupp() operation leads to undefined
// behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_decluniupp decluniupp()
//
// The \c decluniupp() operation can be used to explicitly declare any matrix or matrix expression
// as upper unitriangular:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = decluniupp( A );
   \endcode

// Any matrix or matrix expression that has been declared as upper unitriangular via \c decluniupp()
// will gain all the benefits of a upper unitriangular matrix, which range from reduced runtime
// checking to a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::UniUpperMatrix;

   DynamicMatrix<double> A, B, C;
   UniUpperMatrix< DynamicMatrix<double> > L;
   // ... Resizing and initialization

   isUniUpper( decluniupp( A ) );  // Will always return true without runtime effort

   L = decluniupp( A );  // Omit any runtime check for A being an uniupper matrix

   C = decluniupp( A * B );  // Declare the result of the matrix multiplication as upper
                             // unitriangular, i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c decluniupp() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-uniupper matrix or
// matrix expression as upper unitriangular via the \c decluniupp() operation leads to undefined
// behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declstrupp declstrupp()
//
// The \c declstrupp() operation can be used to explicitly declare any matrix or matrix expression
// as strictly upper triangular:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declstrupp( A );
   \endcode

// Any matrix or matrix expression that has been declared as strictly upper triangular via
// \c declstrupp() will gain all the benefits of a strictly upper triangular matrix, which range
// from reduced runtime checking to a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::StrictlyUpperMatrix;

   DynamicMatrix<double> A, B, C;
   StrictlyUpperMatrix< DynamicMatrix<double> > L;
   // ... Resizing and initialization

   isStrictlyUpper( declstrupp( A ) );  // Will always return true without runtime effort

   L = declstrupp( A );  // Omit any runtime check for A being a strictly upper matrix

   C = declstrupp( A * B );  // Declare the result of the matrix multiplication as strictly upper
                             // triangular, i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declstrupp() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-strictly-upper matrix
// or matrix expression as strictly upper triangular via the \c declstrupp() operation leads to
// undefined behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_decldiag decldiag()
//
// The \c decldiag() operation can be used to explicitly declare any matrix or matrix expression
// as diagonal:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = decldiag( A );
   \endcode

// Any matrix or matrix expression that has been declared as diagonal via \c decldiag() will
// gain all the benefits of a diagonal matrix, which range from reduced runtime checking to
// a considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::DiagonalMatrix;

   DynamicMatrix<double> A, B, C;
   DiagonalMatrix< DynamicMatrix<double> > D;
   // ... Resizing and initialization

   isDiagonal( decldiag( A ) );  // Will always return true without runtime effort

   D = decldiag( A );  // Omit any runtime check for A being a diagonal matrix

   C = decldiag( A * B );  // Declare the result of the matrix multiplication as diagonal,
                           // i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c decldiag() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-diagonal matrix
// or matrix expression as diagonal via the \c decldiag() operation leads to undefined
// behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declid declid()
//
// The \c declid() operation can be used to explicitly declare any matrix or matrix expression
// as identity matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declid( A );
   \endcode

// Any matrix or matrix expression that has been declared as identity matrix via \c declid() will
// gain all the benefits of an identity matrix, which range from reduced runtime checking to a
// considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;
   using blaze::DiagonalMatrix;

   DynamicMatrix<double> A, B, C;
   DiagonalMatrix< DynamicMatrix<double> > D;
   // ... Resizing and initialization

   isIdentity( declid( A ) );  // Will always return true without runtime effort

   D = declid( A );  // Omit any runtime check for A being a diagonal matrix

   C = declid( A ) * B;  // Declare the left operand of the matrix multiplication as an
                         // identity matrix, i.e. perform an optimized matrix multiplication
   \endcode

// \warning The \c declid() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-identity matrix
// or matrix expression as identity matrix via the \c declid() operation leads to undefined
// behavior (which can be violated invariants or wrong computation results)!
//
//
// \n \subsection matrix_operations_declzero declzero()
//
// The \c declzero() operation can be used to explicitly declare any matrix or matrix expression
// as zero matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization

   B = declzero( A );
   \endcode

// Any matrix or matrix expression that has been declared as zero matrix via \c declzero() will
// gain all the benefits of a zero matrix, which range from reduced runtime checking to a
// considerable speed-up in computations:

   \code
   using blaze::DynamicMatrix;

   DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization

   isZero( declzero( A ) );  // Will always return true without runtime effort

   C = declzero( A ) + B;  // Declare the left operand of the matrix addition as a
                           // zero matrix, i.e. no addition needs to be performed
   \endcode

// \warning The \c declzero() operation has the semantics of a cast: The caller is completely
// responsible and the system trusts the given information. Declaring a non-zero matrix or
// matrix expression as zero matrix via the \c declzero() operation leads to undefined behavior
// (which can be violated invariants or wrong computation results)!
//
//
// \n \section matrix_operations_matrix_generators Matrix Generators
// <hr>
//
// \subsection matrix_operations_generate generate()
//
// The \c generate() function returns a dense matrix filled elementwise via the given custom
// binary operation. By default, the returned matrix is a row-major matrix, but this setting can
// be changed via the \c BLAZE_DEFAULT_STORAGE_ORDER switch (see \ref storage_order). Alternatively
// it is possible to specify the storage order explicitly.\n
// The following example demonstrates the use of the \c generate() function:

   \code
   using blaze::generate;
   using blaze::rowMajor;
   using blaze::columnMajor>

   // Generates the uniform integer matrix ( ( 2, 2, 2 ), ( 2, 2, 2 ) )
   blaze::DynamicMatrix<int,rowMajor> A;
   A = generate( 2UL, 3UL, []( size_t i, size_t j ){ return 2; } );

   // Generates the linearly spaced float matrix ( ( 2.1, 3.2, 4.3 ), ( 5.4, 6.5, 7.6 ) )
   blaze::DynamicMatrix<float,rowMajor> B;
   B = generate( 2UL, 3UL, []( size_t i, size_t j ){ return 2.1F + 1.1F*(i*3UL+j); } );

   // Generates the logarithmically spaced double vector ( ( 1.0, 10.0 ), ( 100.0, 1000.0 ) )
   blaze::DynamicMatrix<double,rowMajor> C;
   C = generate<rowMajor>( 2UL, 2UL, []( size_t i, size_t j ) { return blaze::exp10( 1.0 + 1.0*(i*2UL+j) ); } );

   // Generates the vector of integer vectors ( ( 1, 2 ), ( 2, 3 ), ( 3, 4 ), ( 4, 5 ) )
   using VT = StaticVector<int,2UL>;
   blaze::DynamicMatrix<VT,columnMajor> D;
   D = generate<columnMajor>( 2UL, 2UL, []( size_t i, size_t j ) { return evaluate( VT{ 1, 2 } + (i*2UL+j) ); } );
   \endcode

// \n \subsection matrix_operations_uniform uniform()
//
// The \c uniform() function creates a uniform matrix of the given size. By default, the
// resulting uniform matrix is a row-major matrix, but this setting can be changed via the
// \c BLAZE_DEFAULT_STORAGE_ORDER switch (see \ref storage_order). Alternatively it is
// possible to specify the storage order explicitly.\n
// The following example demonstrates the use of the \c uniform() function:

   \code
   using blaze::uniform;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Creates the uniform row-major matrix
   //    ( 1, 1, 1, 1, 1 )
   //    ( 1, 1, 1, 1, 1 )
   auto U1 = uniform( 2UL, 5UL, 1 );

   // Creates the uniform row-major matrix
   //    ( 1.2, 1.2 )
   //    ( 1.2, 1.2 )
   //    ( 1.2, 1.2 )
   auto U2 = uniform<rowMajor>( 3UL, 2UL, 1.2 );

   // Creates the uniform column-major matrix
   //   ( 5U, 5U, 5U, 5U, 5U, 5U, 5U )
   //   ( 5U, 5U, 5U, 5U, 5U, 5U, 5U )
   auto U3 = uniform<columnMajor>( 2UL, 7UL, 5U );
   \endcode

// \n \subsection matrix_operations_zero zero()
//
// The \c zero() function creates a zero matrix of the given element type and size. By default,
// the resulting zero matrix is a row-major matrix, but this setting can be changed via the
// \c BLAZE_DEFAULT_STORAGE_ORDER switch (see \ref storage_order). Alternatively it is possible
// to specify the storage order explicitly.\n
// The following example demonstrates the use of the \c zero() function:

   \code
   using blaze::zero;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Creates the row-major zero matrix
   //    ( 0, 0, 0, 0, 0 )
   //    ( 0, 0, 0, 0, 0 )
   auto Z1 = zero<int>( 2UL, 5UL );

   // Creates the row-major zero matrix
   //    ( 0.0, 0.0 )
   //    ( 0.0, 0.0 )
   //    ( 0.0, 0.0 )
   auto Z2 = zero<double,rowMajor>( 3UL, 2UL );

   // Creates the column-major zero matrix
   //    ( 0U, 0U, 0U, 0U, 0U, 0U, 0U )
   //    ( 0U, 0U, 0U, 0U, 0U, 0U, 0U )
   auto Z3 = zero<unsigned int,columnMajor>( 2UL, 7UL );
   \endcode

// \n \section matrix_operations_matrix_inversion Matrix Inversion
// <hr>
//
// The inverse of a square dense matrix can be computed via the \c inv() function:

   \code
   blaze::DynamicMatrix<float,blaze::rowMajor> A, B;
   // ... Resizing and initialization
   B = inv( A );  // Compute the inverse of A
   \endcode

// Alternatively, an in-place inversion of a dense matrix can be performed via the \c invert()
// function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization
   invert( A );  // In-place matrix inversion
   \endcode

// Both the \c inv() and the \c invert() functions will automatically select the most suited matrix
// inversion algorithm depending on the size and type of the given matrix. For small matrices of
// up to 6x6, both functions use manually optimized kernels for maximum performance. For matrices
// larger than 6x6 the inversion is performed by means of the most suited matrix decomposition
// method: In case of a general matrix the LU decomposition is used, for symmetric matrices the
// LDLT decomposition is applied, for Hermitian matrices the LDLH decomposition is performed, and
// for triangular matrices the inverse is computed via a forward or back substitution.
//
// In case the type of the matrix does not provide additional compile time information about its
// structure (symmetric, lower, upper, diagonal, ...), the information can be provided manually
// by means of \ref matrix_operations_declaration_operations when calling the \c invert() function:

   \code
   invert( declsym( A ) );     // In-place inversion of a symmetric matrix
   invert( declherm( A ) );    // In-place inversion of an Hermitian matrix
   invert( decllow( A ) );     // In-place inversion of a lower triangular matrix
   invert( declunilow( A ) );  // In-place inversion of a lower unitriangular matrix
   invert( declupp( A ) );     // In-place inversion of an upper triangular matrix
   invert( decluniupp( A ) );  // In-place inversion of an upper unitriangular matrix
   invert( decldiag( A ) );    // In-place inversion of a diagonal matrix
   \endcode

// Alternatively, via the \c invert() function it is possible to explicitly specify the inversion
// algorithm:

   \code
   using blaze::byLU;
   using blaze::byLDLT;
   using blaze::byLDLH;
   using blaze::byLLH;

   // In-place inversion of a general matrix by means of an LU decomposition
   invert<byLU>( A );

   // In-place inversion of a symmetric indefinite matrix by means of a Bunch-Kaufman decomposition
   invert<byLDLT>( A );

   // In-place inversion of an Hermitian indefinite matrix by means of a Bunch-Kaufman decomposition
   invert<byLDLH>( A );

   // In-place inversion of a positive definite matrix by means of a Cholesky decomposition
   invert<byLLH>( A );
   \endcode

// Whereas the inversion by means of an LU decomposition works for every general square matrix,
// the inversion by LDLT only works for symmetric indefinite matrices, the inversion by LDLH is
// restricted to Hermitian indefinite matrices and the Cholesky decomposition (LLH) only works
// for Hermitian positive definite matrices. Please note that it is in the responsibility of the
// function caller to guarantee that the selected algorithm is suited for the given matrix. In
// case this precondition is violated the result can be wrong and might not represent the inverse
// of the given matrix!
//
// For both the \c inv() and \c invert() function the matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \c std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The functions invert the dense matrix by means of LAPACK kernels. Thus the functions can
// only be used if a fitting LAPACK library is available and linked to the executable. Otherwise
// a linker error will be created.
//
// \note It is not possible to use any kind of view on the expression object returned by the
// \c inv() function. Also, it is not possible to access individual elements via the function call
// operator on the expression object:

   \code
   row( inv( A ), 2UL );  // Compilation error: Views cannot be used on an inv() expression!
   inv( A )(1,2);         // Compilation error: It is not possible to access individual elements!
   \endcode

// \note The inversion functions do not provide any exception safety guarantee, i.e. in case an
// exception is thrown the matrix may already have been modified.
//
//
// \n \section matrix_operations_matrix_exponential Matrix Exponential
// <hr>
//
// The matrix exponential of a \f$N \times N\f$ matrix \f$ X \f$ is defined as

                  \f[ e^X = \sum\limits_{k=0}^\infty \frac{1}{k!} X^k. \f]

// In order to compute the matrix exponential of a square dense matrix, the \c matexp() function
// can be used:

   \code
   blaze::DynamicMatrix<float,blaze::rowMajor> A, B;
   // ... Resizing and initialization
   B = matexp( A );  // Compute the exponential of A
   \endcode

// \note The matrix exponential can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note It is not possible to use any kind of view on the expression object returned by the
// \c matexp() function. Also, it is not possible to access individual elements via the function
// call operator on the expression object:

   \code
   row( matexp( A ), 2UL );  // Compilation error: Views cannot be used on an matexp() expression!
   matexp( A )(1,2);         // Compilation error: It is not possible to access individual elements!
   \endcode

// \n \section matrix_operations_decomposition Matrix Decomposition
// <hr>
//
// \note All decomposition functions can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The functions decompose a dense matrix by means of LAPACK kernels. Thus the functions can
// only be used if a fitting LAPACK library is available and linked to the executable. Otherwise
// a linker error will be created.
//
// \subsection matrix_operations_decomposition_lu LU Decomposition
//
// The LU decomposition of a dense matrix can be computed via the \c lu() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> L, U, P;

   lu( A, L, U, P );  // LU decomposition of a row-major matrix

   assert( A == L * U * P );
   \endcode

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::columnMajor> L, U, P;

   lu( A, L, U, P );  // LU decomposition of a column-major matrix

   assert( A == P * L * U );
   \endcode

// The function works for both \c rowMajor and \c columnMajor matrices. Note, however, that the
// three matrices \c A, \c L and \c U are required to have the same storage order. Also, please
// note that the way the permutation matrix \c P needs to be applied differs between row-major and
// column-major matrices, since the algorithm uses column interchanges for row-major matrices and
// row interchanges for column-major matrices.
//
// Furthermore, \c lu() can be used with adaptors. For instance, the following example demonstrates
// the LU decomposition of a symmetric matrix into a lower and upper triangular matrix:

   \code
   blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   blaze::LowerMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > L;
   blaze::UpperMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > U;
   blaze::DynamicMatrix<double,blaze::columnMajor> P;

   lu( A, L, U, P );  // LU decomposition of A
   \endcode

// \n \subsection matrix_operations_decomposition_llh Cholesky Decomposition
//
// The Cholesky (LLH) decomposition of a dense matrix can be computed via the \c llh() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> L;

   llh( A, L );  // LLH decomposition of a row-major matrix

   assert( A == L * ctrans( L ) );
   \endcode

// The function works for both \c rowMajor and \c columnMajor matrices and the two matrices \c A
// and \c L can have any storage order.
//
// Furthermore, \c llh() can be used with adaptors. For instance, the following example demonstrates
// the LLH decomposition of a symmetric matrix into a lower triangular matrix:

   \code
   blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   blaze::LowerMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > L;

   llh( A, L );  // Cholesky decomposition of A
   \endcode

// \n \subsection matrix_operations_decomposition_qr QR Decomposition
//
// The QR decomposition of a dense matrix can be computed via the \c qr() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::columnMajor> Q;
   blaze::DynamicMatrix<double,blaze::rowMajor> R;

   qr( A, Q, R );  // QR decomposition of a row-major matrix

   assert( A == Q * R );
   \endcode

// The function works for both \c rowMajor and \c columnMajor matrices and the three matrices
// \c A, \c Q and \c R can have any storage order.
//
// Furthermore, \c qr() can be used with adaptors. For instance, the following example demonstrates
// the QR decomposition of a symmetric matrix into a general matrix and an upper triangular matrix:

   \code
   blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> Q;
   blaze::UpperMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > R;

   qr( A, Q, R );  // QR decomposition of A
   \endcode

// \n \subsection matrix_operations_decomposition_rq RQ Decomposition
//
// Similar to the QR decomposition, the RQ decomposition of a dense matrix can be computed via
// the \c rq() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> R;
   blaze::DynamicMatrix<double,blaze::columnMajor> Q;

   rq( A, R, Q );  // RQ decomposition of a row-major matrix

   assert( A == R * Q );
   \endcode

// The function works for both \c rowMajor and \c columnMajor matrices and the three matrices
// \c A, \c R and \c Q can have any storage order.
//
// Also the \c rq() function can be used in combination with matrix adaptors. For instance, the
// following example demonstrates the RQ decomposition of an Hermitian matrix into a general
// matrix and an upper triangular matrix:

   \code
   blaze::HermitianMatrix< blaze::DynamicMatrix<complex<double>,blaze::columnMajor> > A;
   // ... Resizing and initialization

   blaze::UpperMatrix< blaze::DynamicMatrix<complex<double>,blaze::columnMajor> > R;
   blaze::DynamicMatrix<complex<double>,blaze::rowMajor> Q;

   rq( A, R, Q );  // RQ decomposition of A
   \endcode

// \n \subsection matrix_operations_decomposition_ql QL Decomposition
//
// The QL decomposition of a dense matrix can be computed via the \c ql() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> Q;
   blaze::DynamicMatrix<double,blaze::columnMajor> L;

   ql( A, Q, L );  // QL decomposition of a row-major matrix

   assert( A == Q * L );
   \endcode

// The function works for both \c rowMajor and \c columnMajor matrices and the three matrices
// \c A, \c Q and \c L can have any storage order.
//
// Also the \c ql() function can be used in combination with matrix adaptors. For instance, the
// following example demonstrates the QL decomposition of a symmetric matrix into a general
// matrix and a lower triangular matrix:

   \code
   blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> Q;
   blaze::LowerMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > L;

   ql( A, Q, L );  // QL decomposition of A
   \endcode

// \n \subsection matrix_operations_decomposition_lq LQ Decomposition
//
// The LQ decomposition of a dense matrix can be computed via the \c lq() function:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> L;
   blaze::DynamicMatrix<double,blaze::columnMajor> Q;

   lq( A, L, Q );  // LQ decomposition of a row-major matrix

   assert( A == L * Q );
   \endcode

// The function works for both \c rowMajor and \c columnMajor matrices and the three matrices
// \c A, \c L and \c Q can have any storage order.
//
// Furthermore, \c lq() can be used with adaptors. For instance, the following example demonstrates
// the LQ decomposition of an Hermitian matrix into a lower triangular matrix and a general matrix:

   \code
   blaze::HermitianMatrix< blaze::DynamicMatrix<complex<double>,blaze::columnMajor> > A;
   // ... Resizing and initialization

   blaze::LowerMatrix< blaze::DynamicMatrix<complex<double>,blaze::columnMajor> > L;
   blaze::DynamicMatrix<complex<double>,blaze::rowMajor> Q;

   lq( A, L, Q );  // LQ decomposition of A
   \endcode

// \n \section matrix_operations_linear_systems Linear Systems
// <hr>
//
// The \c solve() function computes a solution for the given dense linear system of equations (LSE)
// \f$ A*x=b \f$, where \c A is the given system matrix, \c x is the solution vector, and \c b is
// the given dense right-hand side vector:

   \code
   blaze::DynamicMatrix<double> A;  // The square general system matrix
   blaze::DynamicVector<double> b;  // The right-hand side vector
   // ... Resizing and initialization

   blaze::DynamicVector<double> x;  // The solution vector

   solve( A, x, b );   // Computing the solution x
   x = solve( A, b );  // Alternative syntax
   \endcode

// Alternatively, \c solve() computes a solution for the given dense LSE \f$ A*X=B \f$, where \c A
// is the given dense system matrix, the columns of \c X are the solution vectors, and the columns
// of \c B are the given right-hand side vectors:

   \code
   blaze::DynamicMatrix<double> A;  // The square general system matrix
   blaze::DynamicMatrix<double> B;  // The right-hand side matrix
   // ... Resizing and initialization

   blaze::DynamicMatrix<double> X;  // The solution matrix

   solve( A, X, B );   // Computing the solutions X
   X = solve( A, B );  // Alternative syntax
   \endcode

// Both \c solve() functions will automatically select the most suited direct solver algorithm
// depending on the size and type of the given system matrix. For small matrices of up to 6x6,
// both functions use manually optimized kernels for maximum performance. For matrices larger
// than 6x6 the computation is performed by means of the most suited LAPACK solver method (see
// \ref lapack_linear_system_solver).
//
// In case the type of the matrix does not provide additional compile time information about
// its structure (symmetric, lower, upper, diagonal, ...), the information can be provided
// manually by means of \ref matrix_operations_declaration_operations when calling the \c solve()
// functions:

   \code
   blaze::DynamicMatrix<double> A;  // The square lower system matrix
   blaze::DynamicVector<double> b;  // The right-hand side vector
   // ... Resizing and initialization

   blaze::DynamicVector<double> x;  // The solution vector

   solve( declsym( A ), x, b );     // Solving the LSE with a symmetric system matrix
   solve( declherm( A ), x, b );    // Solving the LSE with an Hermitian system matrix
   solve( decllow( A ), x, b );     // Solving the LSE with a lower system matrix
   solve( declunilow( A ), x, b );  // Solving the LSE with an unilower system matrix
   solve( declupp( A ), x, b );     // Solving the LSE with an upper system matrix
   solve( decluniupp( A ), x, b );  // Solving the LSE with an uniupper system matrix
   solve( decldiag( A ), x, b );    // Solving the LSE with a diagonal system matrix
   \endcode

// For both \c solve() functions the computation fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the size of the right-hand side vector doesn't match the dimensions of the system matrix;
//  - ... the number of rows of the right-hand side matrix doesn't match the dimensions of the system matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \c std::invalid_argument exception is thrown.
//
// \note The \c solve() functions can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The functions may make use of LAPACK kernels. Thus the functions can only be used if a
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error will
// be created.
//
// \note It is not possible to use any kind of view on the expression object returned by the
// two-argument \c solve() function. Also, it is not possible to access individual elements via
// the function call operator on the expression object:

   \code
   row( solve( A, b ), 2UL );  // Compilation error: Views cannot be used on an solve() expression!
   solve( A, b )[2];           // Compilation error: It is not possible to access individual elements!

   rows( solve( A, B ), { 2UL, 4UL } );  // Compilation error: Views cannot be used on an solve() expression!
   solve( A, B )(1,2);                   // Compilation error: It is not possible to access individual elements!
   \endcode

// \note The \c solve() functions do not provide any exception safety guarantee, i.e. in case an
// exception is thrown the solution vector or matrix may already have been modified.
//
//
// \n \section matrix_operations_eigenvalues Eigenvalues/Eigenvectors
// <hr>
//
// The eigenvalues and eigenvectors of a dense matrix can be computed via the \c eigen() functions.
// The following examples give an impression of the computation of eigenvalues and eigenvectors
// for a general, a symmetric, and an Hermitian matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The general matrix A
   // ... Initialization

   DynamicVector<complex<double>,columnVector> w( 5UL );   // The vector for the complex eigenvalues
   DynamicMatrix<complex<double>,rowMajor> V( 5UL, 5UL );  // The matrix for the left eigenvectors

   w = eigen( A );    // Computing only the eigenvalues of A (one argument)
   eigen( A, w );     // Computing only the eigenvalues of A (two arguments)
   eigen( A, w, V );  // Computing both the eigenvalues and eigenvectors of A (three arguments)
   \endcode

   \code
   using blaze::SymmetricMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > A( 5UL );  // The symmetric matrix A
   // ... Initialization

   DynamicVector<double,columnVector> w( 5UL );       // The vector for the real eigenvalues
   DynamicMatrix<double,rowMajor>     V( 5UL, 5UL );  // The matrix for the left eigenvectors

   w = eigen( A );    // Computing only the eigenvalues of A (one argument)
   eigen( A, w );     // Computing only the eigenvalues of A (two arguments)
   eigen( A, w, V );  // Computing both the eigenvalues and eigenvectors of A (three arguments)
   \endcode

   \code
   using blaze::HermitianMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> > A( 5UL );  // The Hermitian matrix A
   // ... Initialization

   DynamicVector<double,columnVector>      w( 5UL );       // The vector for the real eigenvalues
   DynamicMatrix<complex<double>,rowMajor> V( 5UL, 5UL );  // The matrix for the left eigenvectors

   w = eigen( A );    // Computing only the eigenvalues of A (one argument)
   eigen( A, w );     // Computing only the eigenvalues of A (two arguments)
   eigen( A, w, V );  // Computing both the eigenvalues and eigenvectors of A (three arguments)
   \endcode

// The one- and two-argument functions compute only the eigenvalues of the given \a n-by-\a n
// matrix, the three-argument function additionally computes the eigenvectors. The eigenvalues
// are returned in the given vector \a w and the eigenvectors are returned in the given matrix
// \a V, which are both resized to the correct dimensions (if possible and necessary).
//
// Depending on the given matrix type, the resulting eigenvalues are either of floating point
// or complex type: In case the given matrix is either a compile time symmetric matrix with
// floating point elements or an Hermitian matrix with complex elements, the resulting eigenvalues
// will be of floating point type and therefore the elements of the given eigenvalue vector are
// expected to be of floating point type. In all other cases they are expected to be of complex
// type. Please note that for complex eigenvalues no order of eigenvalues can be assumed, except
// that complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having
// the positive imaginary part first.
//
// In case \a A is a row-major matrix, \a V will contain the left eigenvectors, otherwise \a V
// will contain the right eigenvectors. In case \a V is a row-major matrix the eigenvectors are
// returned in the rows of \a V, in case \a V is a column-major matrix the eigenvectors are
// returned in the columns of \a V. In case the given matrix is a compile time symmetric matrix
// with floating point elements, the resulting eigenvectors will be of floating point type and
// therefore the elements of the given eigenvector matrix are expected to be of floating point
// type. In all other cases they are expected to be of complex type.
//
// The functions fail if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
//  - ... the eigenvalue computation fails.
//
// In all failure cases an exception is thrown.
//
// \note All \c eigen() functions can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The functions compute the eigenvalues and/or eigenvectors of a dense matrix by means of
// LAPACK kernels. Thus the functions can only be used if a fitting LAPACK library is available
// and linked to the executable. Otherwise a linker error will be created.
//
//
// \n \section matrix_operations_singularvalues Singular Values/Singular Vectors
// <hr>
//
// The singular value decomposition (SVD) of a dense matrix can be computed via the \c svd()
// functions. The following two examples give an impression of the computation of singular values
// and singular vectors for a general dense matrix with \c double and \c complex<double> element
// type, respectively:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor>  A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U;  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s;  // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V;  // The matrix for the right singular vectors

   s = svd( A );       // (1) Computing only the singular values of A
   svd( A, s );        // (2) Computing only the singular values of A
   svd( A, U, s, V );  // (3) Computing the singular values and vectors of A

   svd( A, s, 0.0, 1.0 );    // (4) Computing all singular values in the floating point range [0.0..1.0)
   svd( A, U, s, V, 0, 2 );  // (5) Computing the singular values and vectors in the index range [0..2]
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<complex<double>,rowMajor>  A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<complex<double>,rowMajor> U;  // The matrix for the left singular vectors
   DynamicVector<double,columnVector>      s;  // The vector for the singular values
   DynamicMatrix<complex<double>,rowMajor> V;  // The matrix for the right singular vectors

   s = svd( A );       // (1) Computing only the singular values of A
   svd( A, s );        // (2) Computing only the singular values of A
   svd( A, U, s, V );  // (3) Computing the singular values and vectors of A

   svd( A, s, 0.0, 1.0 );    // (4) Computing all singular values in the floating point range [0.0..1.0)
   svd( A, U, s, V, 0, 2 );  // (5) Computing the singular values and vectors in the index range [0..2]
   \endcode

// Functions (1), (2) and (4) compute only singular values of the given general \a m-by-\a n
// matrix, functions (3) and (5) additionally compute singular vectors. The resulting singular
// values are returned in the given vector \a s, the left singular vectors are returned in the
// given matrix \a U, and the right singular vectors are returned in the matrix \a V. \a s, \a U,
// and \a V are resized to the correct dimensions (if possible and necessary).
//
// Functions (4) and (5) allow for the specification of a subset of singular values and/or
// vectors. The number of singular values and vectors to be computed is specified by the lower
// bound \a low and the upper bound \a upp, which either form an integral or a floating point
// range.
//
// In case \a low and \a upp form are of integral type, the function computes all singular values
// in the index range \f$[low..upp]\f$. The \a num resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a \a num-dimensional vector. The resulting left singular vectors are stored
// in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-\a num matrix. The resulting right singular vectors are stored in the given matrix \a V,
// which is either resized (if possible) or expected to be a \a num-by-\a n matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all singular values
// in the half-open interval \f$(low..upp]\f$. The resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a min(\a m,\a n)-dimensional vector. The resulting left singular vectors are
// stored in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-min(\a m,\a n) matrix. The resulting right singular vectors are stored in the given
// matrix \a V, which is either resized (if possible) or expected to be a min(\a m,\a n)-by-\a n
// matrix.
//
// The functions fail if ...
//
//  - ... the given matrix \a U is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
//  - ... the given scalar values don't form a proper range;
//  - ... the singular value decomposition fails.
//
// In all failure cases an exception is thrown.
//
// \note All \c svd() functions can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type or with a sparse matrix results in a compile time error!
//
// \note The functions compute the singular values and/or singular vectors of a dense matrix by
// means of LAPACK kernels. Thus the functions can only be used if a fitting LAPACK library is
// available and linked to the executable. Otherwise a linker error will be created.
//
//
// \n Previous: \ref matrix_types &nbsp; &nbsp; Next: \ref adaptors
*/
//*************************************************************************************************


//**Adaptors***************************************************************************************
/*!\page adaptors Adaptors
//
// \tableofcontents
//
//
// \section adaptors_general General Concepts
// <hr>
//
// Adaptors act as wrappers around the general \ref matrix_types. They adapt the interface of the
// matrices such that certain invariants are preserved. Due to this adaptors can provide a compile
// time guarantee of certain properties, which can be exploited for optimized performance.
//
// The \b Blaze library provides a total of 9 different adaptors:
//
// <ul>
//    <li> \ref adaptors_symmetric_matrices </li>
//    <li> \ref adaptors_hermitian_matrices </li>
//    <li> \ref adaptors_triangular_matrices
//       <ul>
//          <li> \ref adaptors_triangular_matrices "Lower Triangular Matrices"
//             <ul>
//                <li> \ref adaptors_triangular_matrices_lowermatrix </li>
//                <li> \ref adaptors_triangular_matrices_unilowermatrix </li>
//                <li> \ref adaptors_triangular_matrices_strictlylowermatrix </li>
//             </ul>
//          </li>
//          <li> \ref adaptors_triangular_matrices "Upper Triangular Matrices"
//             <ul>
//                <li> \ref adaptors_triangular_matrices_uppermatrix </li>
//                <li> \ref adaptors_triangular_matrices_uniuppermatrix </li>
//                <li> \ref adaptors_triangular_matrices_strictlyuppermatrix </li>
//             </ul>
//          </li>
//          <li> \ref adaptors_triangular_matrices "Diagonal Matrices"
//             <ul>
//                <li> \ref adaptors_triangular_matrices_diagonalmatrix </li>
//             </ul>
//          </li>
//       </ul>
//    </li>
// </ul>
//
// In combination with the general matrix types, \b Blaze provides a total of 40 different matrix
// types that make it possible to exactly adapt the type of matrix to every specific problem.
//
//
// \n \section adaptors_examples Examples
// <hr>
//
// The following code examples give an impression on the use of adaptors. The first example shows
// the multiplication between two lower matrices:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   LowerMatrix< DynamicMatrix<double,rowMajor> > A;
   LowerMatrix< DynamicMatrix<double,columnMajor> > B;
   DynamicMatrix<double,columnMajor> C;

   // ... Resizing and initialization

   C = A * B;
   \endcode

// When multiplying two matrices, at least one of which is triangular, \b Blaze can exploit the
// fact that either the lower or upper part of the matrix contains only default elements and
// restrict the algorithm to the non-zero elements. Thus the adaptor provides a significant
// performance advantage in comparison to a general matrix multiplication, especially for large
// matrices.
//
// The second example shows the \c SymmetricMatrix adaptor in a row-major dense matrix/sparse
// vector multiplication:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > A;
   CompressedVector<double,columnVector> x;
   DynamicVector<double,columnVector> y;

   // ... Resizing and initialization

   y = A * x;
   \endcode

// In this example it is not intuitively apparent that using a row-major matrix is not the best
// possible choice in terms of performance since the computation cannot be vectorized. Choosing
// a column-major matrix instead, however, would enable a vectorized computation. Therefore
// \b Blaze exploits the fact that \c A is symmetric, selects the best suited storage order and
// evaluates the multiplication as

   \code
   y = trans( A ) * x;
   \endcode

// which significantly increases the performance.
//
// \n Previous: \ref matrix_operations &nbsp; &nbsp; Next: \ref adaptors_symmetric_matrices
*/
//*************************************************************************************************


//**Symmetric Matrices*****************************************************************************
/*!\page adaptors_symmetric_matrices Symmetric Matrices
//
// \tableofcontents
//
//
// \n \section adaptors_symmetric_matrices_general Symmetric Matrices
// <hr>
//
// In contrast to general matrices, which have no restriction in their number of rows and columns
// and whose elements can have any value, symmetric matrices provide the compile time guarantee
// to be square matrices with pair-wise identical values. Mathematically, this means that a
// symmetric matrix is always equal to its transpose (\f$ A = A^T \f$) and that all non-diagonal
// values have an identical counterpart (\f$ a_{ij} == a_{ji} \f$). This symmetry property can
// be exploited to provide higher efficiency and/or lower memory consumption. Within the \b Blaze
// library, symmetric matrices are realized by the \ref adaptors_symmetric_matrices_symmetricmatrix
// class template.
//
//
// \n \section adaptors_symmetric_matrices_symmetricmatrix SymmetricMatrix
// <hr>
//
// The SymmetricMatrix class template is an adapter for existing dense and sparse matrix types.
// It inherits the properties and the interface of the given matrix type \c MT and extends it
// by enforcing the additional invariant of symmetry (i.e. the matrix is always equal to its
// transpose \f$ A = A^T \f$). It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/SymmetricMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class SymmetricMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. SymmetricMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Note
// that the given matrix type must be either resizable (as for instance blaze::HybridMatrix or
// blaze::DynamicMatrix) or must be square at compile time (as for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible symmetric matrices:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of a 3x3 row-major dense symmetric matrix with static memory
   blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,rowMajor> > A;

   // Definition of a resizable column-major dense symmetric matrix based on HybridMatrix
   blaze::SymmetricMatrix< blaze::HybridMatrix<float,4UL,4UL,columnMajor> B;

   // Definition of a resizable row-major dense symmetric matrix based on DynamicMatrix
   blaze::SymmetricMatrix< blaze::DynamicMatrix<double,rowMajor> > C;

   // Definition of a fixed size row-major dense symmetric matrix based on CustomMatrix
   blaze::SymmetricMatrix< blaze::CustomMatrix<double,unaligned,unpadded,rowMajor> > D;

   // Definition of a compressed row-major single precision symmetric matrix
   blaze::SymmetricMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > E;
   \endcode

// The storage order of a symmetric matrix is depending on the storage order of the adapted matrix
// type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e. is specified as
// blaze::rowMajor), the symmetric matrix will also be a row-major matrix. Otherwise, if the
// adapted matrix is column-major (i.e. is specified as blaze::columnMajor), the symmetric matrix
// will also be a column-major matrix.
//
//
// \n \section adaptors_symmetric_matrices_special_properties Special Properties of Symmetric Matrices
// <hr>
//
// A symmetric matrix is used exactly like a matrix of the underlying, adapted matrix type \c MT.
// It also provides (nearly) the same interface as the underlying matrix type. However, there are
// some important exceptions resulting from the symmetry constraint:
//
//  -# <b>\ref adaptors_symmetric_matrices_square</b>
//  -# <b>\ref adaptors_symmetric_matrices_symmetry</b>
//  -# <b>\ref adaptors_symmetric_matrices_initialization</b>
//
// \n \subsection adaptors_symmetric_matrices_square Symmetric Matrices Must Always be Square!
//
// In case a resizable matrix is used (as for instance blaze::HybridMatrix, blaze::DynamicMatrix,
// or blaze::CompressedMatrix), this means that the according constructors, the \c resize() and
// the \c extend() functions only expect a single parameter, which specifies both the number of
// rows and columns, instead of two (one for the number of rows and one for the number of columns):

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;

   // Default constructed, default initialized, row-major 3x3 symmetric dynamic matrix
   SymmetricMatrix< DynamicMatrix<double,rowMajor> > A( 3 );

   // Resizing the matrix to 5x5
   A.resize( 5 );

   // Extending the number of rows and columns by 2, resulting in a 7x7 matrix
   A.extend( 2 );
   \endcode

// In case a matrix with a fixed size is used (as for instance blaze::StaticMatrix), the number
// of rows and number of columns must be specified equally:

   \code
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::columnMajor;

   // Correct setup of a fixed size column-major 3x3 symmetric static matrix
   SymmetricMatrix< StaticMatrix<int,3UL,3UL,columnMajor> > A;

   // Compilation error: the provided matrix type is not a square matrix type
   SymmetricMatrix< StaticMatrix<int,3UL,4UL,columnMajor> > B;
   \endcode

// \n \subsection adaptors_symmetric_matrices_symmetry The Symmetric Property is Always Enforced!
//
// This means that modifying the element \f$ a_{ij} \f$ of a symmetric matrix also modifies its
// counterpart element \f$ a_{ji} \f$. Also, it is only possible to assign matrices that are
// symmetric themselves:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;

   // Default constructed, row-major 3x3 symmetric compressed matrix
   SymmetricMatrix< CompressedMatrix<double,rowMajor> > A( 3 );

   // Initializing three elements via the function call operator
   A(0,0) = 1.0;  // Initialization of the diagonal element (0,0)
   A(0,2) = 2.0;  // Initialization of the elements (0,2) and (2,0)

   // Inserting three more elements via the insert() function
   A.insert( 1, 1, 3.0 );  // Inserting the diagonal element (1,1)
   A.insert( 1, 2, 4.0 );  // Inserting the elements (1,2) and (2,1)

   // Access via a non-const iterator
   *A.begin(1UL) = 10.0;  // Modifies both elements (1,0) and (0,1)

   // Erasing elements via the erase() function
   A.erase( 0, 0 );  // Erasing the diagonal element (0,0)
   A.erase( 0, 2 );  // Erasing the elements (0,2) and (2,0)

   // Construction from a symmetric dense matrix
   StaticMatrix<double,3UL,3UL> B{ {  3.0,  8.0, -2.0 },
                                   {  8.0,  0.0, -1.0 },
                                   { -2.0, -1.0,  4.0 } };

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > C( B );  // OK

   // Assignment of a non-symmetric dense matrix
   StaticMatrix<double,3UL,3UL> D{ {  3.0,  7.0, -2.0 },
                                   {  8.0,  0.0, -1.0 },
                                   { -2.0, -1.0,  4.0 } };

   C = D;  // Throws an exception; symmetric invariant would be violated!
   \endcode

// The same restriction also applies to the \c append() function for sparse matrices: Appending
// the element \f$ a_{ij} \f$ additionally inserts the element \f$ a_{ji} \f$ into the matrix.
// Despite the additional insertion, the \c append() function still provides the most efficient
// way to set up a symmetric sparse matrix. In order to achieve the maximum efficiency, the
// capacity of the individual rows/columns of the matrix should to be specifically prepared with
// \c reserve() calls:

   \code
   using blaze::CompressedMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;

   // Setup of the symmetric matrix
   //
   //       ( 0 1 3 )
   //   A = ( 1 2 0 )
   //       ( 3 0 0 )
   //
   SymmetricMatrix< CompressedMatrix<double,rowMajor> > A( 3 );

   A.reserve( 5 );         // Reserving enough space for 5 non-zero elements
   A.reserve( 0, 2 );      // Reserving two non-zero elements in the first row
   A.reserve( 1, 2 );      // Reserving two non-zero elements in the second row
   A.reserve( 2, 1 );      // Reserving a single non-zero element in the third row
   A.append( 0, 1, 1.0 );  // Appending the value 1 at position (0,1) and (1,0)
   A.append( 1, 1, 2.0 );  // Appending the value 2 at position (1,1)
   A.append( 2, 0, 3.0 );  // Appending the value 3 at position (2,0) and (0,2)
   \endcode

// The symmetry property is also enforced for symmetric custom matrices: In case the given array
// of elements does not represent a symmetric matrix, a \c std::invalid_argument exception is
// thrown:

   \code
   using blaze::CustomMatrix;
   using blaze::SymmetricMatrix;
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   using CustomSymmetric = SymmetricMatrix< CustomMatrix<double,unaligned,unpadded,rowMajor> >;

   // Creating a 3x3 symmetric custom matrix from a properly initialized array
   double array[9] = { 1.0, 2.0, 4.0,
                       2.0, 3.0, 5.0,
                       4.0, 5.0, 6.0 };
   CustomSymmetric A( array, 3UL );  // OK

   // Attempt to create a second 3x3 symmetric custom matrix from an uninitialized array
   std::unique_ptr<double[]> memory( new double[9UL] );
   CustomSymmetric B( memory.get(), 3UL );  // Throws an exception
   \endcode

// Finally, the symmetry property is enforced for views (rows, columns, submatrices, ...) on the
// symmetric matrix. The following example demonstrates that modifying the elements of an entire
// row of the symmetric matrix also affects the counterpart elements in the according column of
// the matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Setup of the symmetric matrix
   //
   //       ( 0 1 0 2 )
   //   A = ( 1 3 4 0 )
   //       ( 0 4 0 5 )
   //       ( 2 0 5 0 )
   //
   SymmetricMatrix< DynamicMatrix<int> > A( 4 );
   A(0,1) = 1;
   A(0,3) = 2;
   A(1,1) = 3;
   A(1,2) = 4;
   A(2,3) = 5;

   // Setting all elements in the 1st row to 0 results in the matrix
   //
   //       ( 0 0 0 2 )
   //   A = ( 0 0 0 0 )
   //       ( 0 0 0 5 )
   //       ( 2 0 5 0 )
   //
   row( A, 1 ) = 0;
   \endcode

// The next example demonstrates the (compound) assignment to submatrices of symmetric matrices.
// Since the modification of element \f$ a_{ij} \f$ of a symmetric matrix also modifies the
// element \f$ a_{ji} \f$, the matrix to be assigned must be structured such that the symmetry
// of the symmetric matrix is preserved. Otherwise a \c std::invalid_argument exception is
// thrown:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Setup of two default 4x4 symmetric matrices
   SymmetricMatrix< DynamicMatrix<int> > A1( 4 ), A2( 4 );

   // Setup of the 3x2 dynamic matrix
   //
   //       ( 1 2 )
   //   B = ( 3 4 )
   //       ( 5 6 )
   //
   DynamicMatrix<int> B{ { 1, 2 }, { 3, 4 }, { 5, 6 } };

   // OK: Assigning B to a submatrix of A1 such that the symmetry can be preserved
   //
   //        ( 0 0 1 2 )
   //   A1 = ( 0 0 3 4 )
   //        ( 1 3 5 6 )
   //        ( 2 4 6 0 )
   //
   submatrix( A1, 0UL, 2UL, 3UL, 2UL ) = B;  // OK

   // Error: Assigning B to a submatrix of A2 such that the symmetry cannot be preserved!
   //   The elements marked with X cannot be assigned unambiguously!
   //
   //        ( 0 1 2 0 )
   //   A2 = ( 1 3 X 0 )
   //        ( 2 X 6 0 )
   //        ( 0 0 0 0 )
   //
   submatrix( A2, 0UL, 1UL, 3UL, 2UL ) = B;  // Assignment throws an exception!
   \endcode

// \n \subsection adaptors_symmetric_matrices_initialization The Elements of a Dense Symmetric Matrix are Always Default Initialized!
//
// Although this results in a small loss of efficiency (especially in case all default values are
// overridden afterwards), this property is important since otherwise the symmetric property of
// dense symmetric matrices could not be guaranteed:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Uninitialized, 5x5 row-major dynamic matrix
   DynamicMatrix<int,rowMajor> A( 5, 5 );

   // Default initialized, 5x5 row-major symmetric dynamic matrix
   SymmetricMatrix< DynamicMatrix<int,rowMajor> > B( 5 );
   \endcode

// \n \section adaptors_symmetric_matrices_arithmetic_operations Arithmetic Operations
// <hr>
//
// A SymmetricMatrix matrix can participate in numerical operations in any way any other dense
// or sparse matrix can participate. It can also be combined with any other dense or sparse vector
// or matrix. The following code example gives an impression of the use of SymmetricMatrix within
// arithmetic operations:

   \code
   using blaze::SymmetricMatrix;
   using blaze::DynamicMatrix;
   using blaze::HybridMatrix;
   using blaze::StaticMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   DynamicMatrix<double,rowMajor> A( 3, 3 );
   CompressedMatrix<double,rowMajor> B( 3, 3 );

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > C( 3 );
   SymmetricMatrix< CompressedMatrix<double,rowMajor> > D( 3 );

   SymmetricMatrix< HybridMatrix<float,3UL,3UL,rowMajor> > E;
   SymmetricMatrix< StaticMatrix<float,3UL,3UL,columnMajor> > F;

   E = A + B;     // Matrix addition and assignment to a row-major symmetric matrix (includes runtime check)
   F = C - D;     // Matrix subtraction and assignment to a column-major symmetric matrix (only compile time check)
   F = A * D;     // Matrix multiplication between a dense and a sparse matrix (includes runtime check)

   C *= 2.0;      // In-place scaling of matrix C
   E  = 2.0 * B;  // Scaling of matrix B (includes runtime check)
   F  = C * 2.0;  // Scaling of matrix C (only compile time check)

   E += A - B;    // Addition assignment (includes runtime check)
   F -= C + D;    // Subtraction assignment (only compile time check)
   F *= A * D;    // Multiplication assignment (includes runtime check)
   \endcode

// Note that it is possible to assign any kind of matrix to a symmetric matrix. In case the matrix
// to be assigned is not symmetric at compile time, a runtime check is performed.
//
//
// \n \section adaptors_symmetric_matrices_block_matrices Symmetric Block Matrices
// <hr>
//
// It is also possible to use symmetric block matrices:

   \code
   using blaze::CompressedMatrix;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;

   // Definition of a 3x3 symmetric block matrix based on CompressedMatrix
   SymmetricMatrix< CompressedMatrix< StaticMatrix<int,3UL,3UL> > > A( 3 );
   \endcode

// Also in this case, the SymmetricMatrix class template enforces the invariant of symmetry and
// guarantees that a modifications of element \f$ a_{ij} \f$ of the adapted matrix is also
// applied to element \f$ a_{ji} \f$:

   \code
   // Inserting the elements (2,4) and (4,2)
   A.insert( 2, 4, StaticMatrix<int,3UL,3UL>{ { 1, -4,  5 },
                                              { 6,  8, -3 },
                                              { 2, -1,  2 } } );

   // Manipulating the elements (2,4) and (4,2)
   A(2,4)(1,1) = -5;
   \endcode

// For more information on block matrices, see the tutorial on \ref block_vectors_and_matrices.
//
//
// \n \section adaptors_symmetric_matrices_performance Performance Considerations
// <hr>
//
// When the symmetric property of a matrix is known beforehands using the SymmetricMatrix adaptor
// instead of a general matrix can be a considerable performance advantage. The \b Blaze library
// tries to exploit the properties of symmetric matrices whenever possible. However, there are
// also situations when using a symmetric matrix introduces some overhead. The following examples
// demonstrate several situations where symmetric matrices can positively or negatively impact
// performance.
//
// \n \subsection adaptors_symmetric_matrices_matrix_matrix_multiplication Positive Impact: Matrix/Matrix Multiplication
//
// When multiplying two matrices, at least one of which is symmetric, \b Blaze can exploit the fact
// that \f$ A = A^T \f$ and choose the fastest and most suited combination of storage orders for the
// multiplication. The following example demonstrates this by means of a dense matrix/sparse matrix
// multiplication:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > A;
   SymmetricMatrix< CompressedMatrix<double,columnMajor> > B;
   DynamicMatrix<double,columnMajor> C;

   // ... Resizing and initialization

   C = A * B;
   \endcode

// Intuitively, the chosen combination of a row-major and a column-major matrix is the most suited
// for maximum performance. However, \b Blaze evaluates the multiplication as

   \code
   C = A * trans( B );
   \endcode

// which significantly increases the performance since in contrast to the original formulation the
// optimized form can be vectorized. Therefore, in the context of matrix multiplications, using the
// SymmetricMatrix adapter is obviously an advantage.
//
// \n \subsection adaptors_symmetric_matrices_matrix_vector_multiplication Positive Impact: Matrix/Vector Multiplication
//
// A similar optimization is possible in case of matrix/vector multiplications:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > A;
   CompressedVector<double,columnVector> x;
   DynamicVector<double,columnVector> y;

   // ... Resizing and initialization

   y = A * x;
   \endcode

// In this example it is not intuitively apparent that using a row-major matrix is not the best
// possible choice in terms of performance since the computation cannot be vectorized. Choosing
// a column-major matrix instead, however, would enable a vectorized computation. Therefore
// \b Blaze exploits the fact that \c A is symmetric, selects the best suited storage order and
// evaluates the multiplication as

   \code
   y = trans( A ) * x;
   \endcode

// which also significantly increases the performance.
//
// \n \subsection adaptors_symmetric_matrices_views Positive Impact: Row/Column Views on Column/Row-Major Matrices
//
// Another example is the optimization of a row view on a column-major symmetric matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::columnMajor;

   SymmetricMatrix< DynamicMatrix<double,columnMajor> > A( 10UL );
   auto row5 = row( A, 5UL );
   \endcode

// Usually, a row view on a column-major matrix results in a considerable performance decrease in
// comparison to a row view on a row-major matrix due to the non-contiguous storage of the matrix
// elements. However, in case of symmetric matrices, \b Blaze instead uses the according column of
// the matrix, which provides the same performance as if the matrix would be row-major. Note that
// this also works for column views on row-major matrices, where \b Blaze can use the according
// row instead of a column in order to provide maximum performance.
//
// \n \subsection adaptors_symmetric_matrices_assignment Negative Impact: Assignment of a General Matrix
//
// In contrast to using a symmetric matrix on the right-hand side of an assignment (i.e. for read
// access), which introduces absolutely no performance penalty, using a symmetric matrix on the
// left-hand side of an assignment (i.e. for write access) may introduce additional overhead when
// it is assigned a general matrix, which is not symmetric at compile time:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   SymmetricMatrix< DynamicMatrix<double> > A, C;
   DynamicMatrix<double> B;

   B = A;  // Only read-access to the symmetric matrix; no performance penalty
   C = A;  // Assignment of a symmetric matrix to another symmetric matrix; no runtime overhead
   C = B;  // Assignment of a general matrix to a symmetric matrix; some runtime overhead
   \endcode

// When assigning a general, potentially not symmetric matrix to a symmetric matrix it is necessary
// to check whether the matrix is symmetric at runtime in order to guarantee the symmetry property
// of the symmetric matrix. In case it turns out to be symmetric, it is assigned as efficiently as
// possible, if it is not, an exception is thrown. In order to prevent this runtime overhead it is
// therefore generally advisable to assign symmetric matrices to other symmetric matrices.\n
// In this context it is especially noteworthy that in contrast to additions and subtractions the
// multiplication of two symmetric matrices does not necessarily result in another symmetric matrix:

   \code
   SymmetricMatrix< DynamicMatrix<double> > A, B, C;

   C = A + B;  // Results in a symmetric matrix; no runtime overhead
   C = A - B;  // Results in a symmetric matrix; no runtime overhead
   C = A * B;  // Is not guaranteed to result in a symmetric matrix; some runtime overhead
   \endcode

// \n Previous: \ref adaptors &nbsp; &nbsp; Next: \ref adaptors_hermitian_matrices
*/
//*************************************************************************************************


//**Hermitian Matrices*****************************************************************************
/*!\page adaptors_hermitian_matrices Hermitian Matrices
//
// \tableofcontents
//
//
// \n \section adaptors_hermitian_matrices_general Hermitian Matrices
// <hr>
//
// In addition to symmetric matrices, \b Blaze also provides an adaptor for Hermitian matrices.
// Hermitian matrices provide the compile time guarantee to be square matrices with pair-wise
// conjugate complex values. Mathematically, this means that an Hermitian matrix is always equal
// to its conjugate transpose (\f$ A = \overline{A^T} \f$) and that all non-diagonal values have
// a complex conjugate counterpart (\f$ a_{ij} == \overline{a_{ji}} \f$). Within the \b Blaze
// library, Hermitian matrices are realized by the \ref adaptors_hermitian_matrices_hermitianmatrix
// class template.
//
//
// \n \section adaptors_hermitian_matrices_hermitianmatrix HermitianMatrix
// <hr>
//
// The HermitianMatrix class template is an adapter for existing dense and sparse matrix types.
// It inherits the properties and the interface of the given matrix type \c MT and extends it by
// enforcing the additional invariant of Hermitian symmetry (i.e. the matrix is always equal to
// its conjugate transpose \f$ A = \overline{A^T} \f$). It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/HermitianMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class HermitianMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. HermitianMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Also,
// the given matrix type must have numeric element types (i.e. all integral types except \c bool,
// floating point and complex types). Note that the given matrix type must be either resizable (as
// for instance blaze::HybridMatrix or blaze::DynamicMatrix) or must be square at compile time (as
// for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible Hermitian matrices:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of a 3x3 row-major dense Hermitian matrix with static memory
   blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,rowMajor> > A;

   // Definition of a resizable column-major dense Hermitian matrix based on HybridMatrix
   blaze::HermitianMatrix< blaze::HybridMatrix<float,4UL,4UL,columnMajor> B;

   // Definition of a resizable row-major dense Hermitian matrix based on DynamicMatrix
   blaze::HermitianMatrix< blaze::DynamicMatrix<std::complex<double>,rowMajor> > C;

   // Definition of a fixed size row-major dense Hermitian matrix based on CustomMatrix
   blaze::HermitianMatrix< blaze::CustomMatrix<double,unaligned,unpadded,rowMajor> > D;

   // Definition of a compressed row-major single precision complex Hermitian matrix
   blaze::HermitianMatrix< blaze::CompressedMatrix<std::complex<float>,rowMajor> > E;
   \endcode

// The storage order of an Hermitian matrix is depending on the storage order of the adapted matrix
// type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e. is specified as
// blaze::rowMajor), the Hermitian matrix will also be a row-major matrix. Otherwise, if the
// adapted matrix is column-major (i.e. is specified as blaze::columnMajor), the Hermitian matrix
// will also be a column-major matrix.
//
//
// \n \section adaptors_hermitian_matrices_vs_symmetric_matrices Hermitian Matrices vs. Symmetric Matrices
//
// The blaze::HermitianMatrix adaptor and the blaze::SymmetricMatrix adaptor share several traits.
// However, there are a couple of differences, both from a mathematical point of view as well as
// from an implementation point of view.
//
// From a mathematical point of view, a matrix is called symmetric when it is equal to its
// transpose (\f$ A = A^T \f$) and it is called Hermitian when it is equal to its conjugate
// transpose (\f$ A = \overline{A^T} \f$). For matrices of real values, however, these two
// conditions coincide, which means that symmetric matrices of real values are also Hermitian
// and Hermitian matrices of real values are also symmetric.
//
// From an implementation point of view, \b Blaze restricts Hermitian matrices to numeric data
// types (i.e. all integral types except \c bool, floating point and complex types), whereas
// symmetric matrices can also be block matrices (i.e. can have vector or matrix elements).
// For built-in element types, the HermitianMatrix adaptor behaves exactly like the according
// SymmetricMatrix implementation. For complex element types, however, the Hermitian property
// is enforced (see also \ref adaptors_hermitian_matrices_hermitian).

	\code
	using blaze::DynamicMatrix;
	using blaze::DynamicVector;
	using blaze::HermitianMatrix;
	using blaze::SymmetricMatrix;

	// The following two matrices provide an identical experience (including performance)
	HermitianMatrix< DynamicMatrix<double> > A;  // Both Hermitian and symmetric
	SymmetricMatrix< DynamicMatrix<double> > B;  // Both Hermitian and symmetric

	// The following two matrices will behave differently
	HermitianMatrix< DynamicMatrix< complex<double> > > C;  // Only Hermitian
	SymmetricMatrix< DynamicMatrix< complex<double> > > D;  // Only symmetric

	// Hermitian block matrices are not allowed
	HermitianMatrix< DynamicMatrix< DynamicVector<double> > > E;  // Compilation error!
	SymmetricMatrix< DynamicMatrix< DynamicVector<double> > > F;  // Symmetric block matrix
	\endcode

// \n \section adaptors_hermitian_matrices_special_properties Special Properties of Hermitian Matrices
// <hr>
//
// An Hermitian matrix is used exactly like a matrix of the underlying, adapted matrix type \c MT.
// It also provides (nearly) the same interface as the underlying matrix type. However, there are
// some important exceptions resulting from the Hermitian symmetry constraint:
//
//  -# <b>\ref adaptors_hermitian_matrices_square</b>
//  -# <b>\ref adaptors_hermitian_matrices_hermitian</b>
//  -# <b>\ref adaptors_hermitian_matrices_initialization</b>
//
// \n \subsection adaptors_hermitian_matrices_square Hermitian Matrices Must Always be Square!
//
// In case a resizable matrix is used (as for instance blaze::HybridMatrix, blaze::DynamicMatrix,
// or blaze::CompressedMatrix), this means that the according constructors, the \c resize() and
// the \c extend() functions only expect a single parameter, which specifies both the number of
// rows and columns, instead of two (one for the number of rows and one for the number of columns):

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;
   using blaze::rowMajor;

   // Default constructed, default initialized, row-major 3x3 Hermitian dynamic matrix
   HermitianMatrix< DynamicMatrix<std::complex<double>,rowMajor> > A( 3 );

   // Resizing the matrix to 5x5
   A.resize( 5 );

   // Extending the number of rows and columns by 2, resulting in a 7x7 matrix
   A.extend( 2 );
   \endcode

// In case a matrix with a fixed size is used (as for instance blaze::StaticMatrix), the number
// of rows and number of columns must be specified equally:

   \code
   using blaze::StaticMatrix;
   using blaze::HermitianMatrix;
   using blaze::columnMajor;

   // Correct setup of a fixed size column-major 3x3 Hermitian static matrix
   HermitianMatrix< StaticMatrix<std::complex<float>,3UL,3UL,columnMajor> > A;

   // Compilation error: the provided matrix type is not a square matrix type
   HermitianMatrix< StaticMatrix<std::complex<float>,3UL,4UL,columnMajor> > B;
   \endcode

// \n \subsection adaptors_hermitian_matrices_hermitian The Hermitian Property is Always Enforced!
//
// This means that the following properties of an Hermitian matrix are always guaranteed:
//
//  - The diagonal elements are real numbers, i.e. the imaginary part is zero
//  - Element \f$ a_{ij} \f$ is always the complex conjugate of element \f$ a_{ji} \f$
//
// Thus modifying the element \f$ a_{ij} \f$ of an Hermitian matrix also modifies its
// counterpart element \f$ a_{ji} \f$. Also, it is only possible to assign matrices that
// are Hermitian themselves:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::HermitianMatrix;
   using blaze::rowMajor;

	using cplx = std::complex<double>;

   // Default constructed, row-major 3x3 Hermitian compressed matrix
   HermitianMatrix< CompressedMatrix<cplx,rowMajor> > A( 3 );

   // Initializing the matrix via the function call operator
	//
	//  ( (1, 0) (0,0) (2,1) )
	//  ( (0, 0) (0,0) (0,0) )
	//  ( (2,-1) (0,0) (0,0) )
   //
   A(0,0) = cplx( 1.0, 0.0 );  // Initialization of the diagonal element (0,0)
   A(0,2) = cplx( 2.0, 1.0 );  // Initialization of the elements (0,2) and (2,0)

   // Inserting three more elements via the insert() function
	//
	//  ( (1,-3) (0,0) (2, 1) )
	//  ( (0, 0) (2,0) (4,-2) )
	//  ( (2,-1) (4,2) (0, 0) )
   //
   A.insert( 1, 1, cplx( 2.0,  0.0 ) );  // Inserting the diagonal element (1,1)
   A.insert( 1, 2, cplx( 4.0, -2.0 ) );  // Inserting the elements (1,2) and (2,1)

   // Access via a non-const iterator
	//
	//  ( (1,-3) (8,1) (2, 1) )
	//  ( (8,-1) (2,0) (4,-2) )
	//  ( (2,-1) (4,2) (0, 0) )
   //
   *A.begin(1UL) = cplx( 8.0, -1.0 );  // Modifies both elements (1,0) and (0,1)

   // Erasing elements via the erase() function
	//
	//  ( (0, 0) (8,1) (0, 0) )
	//  ( (8,-1) (2,0) (4,-2) )
	//  ( (0, 0) (4,2) (0, 0) )
   //
   A.erase( 0, 0 );  // Erasing the diagonal element (0,0)
   A.erase( 0, 2 );  // Erasing the elements (0,2) and (2,0)

   // Construction from an Hermitian dense matrix
   StaticMatrix<cplx,3UL,3UL> B{ { cplx(  3.0,  0.0 ), cplx(  8.0, 2.0 ), cplx( -2.0,  2.0 ) },
                                 { cplx(  8.0,  1.0 ), cplx(  0.0, 0.0 ), cplx( -1.0, -1.0 ) },
                                 { cplx( -2.0, -2.0 ), cplx( -1.0, 1.0 ), cplx(  4.0,  0.0 ) } };

   HermitianMatrix< DynamicMatrix<double,rowMajor> > C( B );  // OK

   // Assignment of a non-Hermitian dense matrix
	StaticMatrix<cplx,3UL,3UL> D{ { cplx(  3.0, 0.0 ), cplx(  7.0, 2.0 ), cplx( 3.0, 2.0 ) },
                                 { cplx(  8.0, 1.0 ), cplx(  0.0, 0.0 ), cplx( 6.0, 4.0 ) },
                                 { cplx( -2.0, 2.0 ), cplx( -1.0, 1.0 ), cplx( 4.0, 0.0 ) } };

   C = D;  // Throws an exception; Hermitian invariant would be violated!
   \endcode

// The same restriction also applies to the \c append() function for sparse matrices: Appending
// the element \f$ a_{ij} \f$ additionally inserts the element \f$ a_{ji} \f$ into the matrix.
// Despite the additional insertion, the \c append() function still provides the most efficient
// way to set up an Hermitian sparse matrix. In order to achieve the maximum efficiency, the
// capacity of the individual rows/columns of the matrix should to be specifically prepared with
// \c reserve() calls:

   \code
   using blaze::CompressedMatrix;
   using blaze::HermitianMatrix;
   using blaze::rowMajor;

	using cplx = std::complex<double>;

   // Setup of the Hermitian matrix
   //
   //       ( (0, 0) (1,2) (3,-4) )
   //   A = ( (1,-2) (2,0) (0, 0) )
   //       ( (3, 4) (0,0) (0, 0) )
   //
   HermitianMatrix< CompressedMatrix<cplx,rowMajor> > A( 3 );

   A.reserve( 5 );         // Reserving enough space for 5 non-zero elements
   A.reserve( 0, 2 );      // Reserving two non-zero elements in the first row
   A.reserve( 1, 2 );      // Reserving two non-zero elements in the second row
   A.reserve( 2, 1 );      // Reserving a single non-zero element in the third row

   A.append( 0, 1, cplx( 1.0, 2.0 ) );  // Appending an element at position (0,1) and (1,0)
   A.append( 1, 1, cplx( 2.0, 0.0 ) );  // Appending an element at position (1,1)
   A.append( 2, 0, cplx( 3.0, 4.0 ) );  // Appending an element at position (2,0) and (0,2)
   \endcode

// The Hermitian property is also enforced for Hermitian custom matrices: In case the given array
// of elements does not represent an Hermitian matrix, a \c std::invalid_argument exception is
// thrown:

   \code
   using blaze::CustomMatrix;
   using blaze::HermitianMatrix;
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   using CustomHermitian = HermitianMatrix< CustomMatrix<double,unaligned,unpadded,rowMajor> >;

   // Creating a 3x3 Hermitian custom matrix from a properly initialized array
   double array[9] = { 1.0, 2.0, 4.0,
                       2.0, 3.0, 5.0,
                       4.0, 5.0, 6.0 };
   CustomHermitian A( array, 3UL );  // OK

   // Attempt to create a second 3x3 Hermitian custom matrix from an uninitialized array
   std::unique_ptr<double[]> memory( new double[9UL] );
   CustomHermitian B( memory.get(), 3UL );  // Throws an exception
   \endcode

// Finally, the Hermitian property is enforced for views (rows, columns, submatrices, ...) on the
// Hermitian matrix. The following example demonstrates that modifying the elements of an entire
// row of the Hermitian matrix also affects the counterpart elements in the according column of
// the matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermtianMatrix;

	using cplx = std::complex<double>;

   // Setup of the Hermitian matrix
   //
   //       ( (0, 0) (1,-1) (0,0) (2, 1) )
   //   A = ( (1, 1) (3, 0) (4,2) (0, 0) )
   //       ( (0, 0) (4,-2) (0,0) (5,-3) )
   //       ( (2,-1) (0, 0) (5,3) (0, 0) )
   //
   HermitianMatrix< DynamicMatrix<int> > A( 4 );
   A(0,1) = cplx( 1.0, -1.0 );
   A(0,3) = cplx( 2.0,  1.0 );
   A(1,1) = cplx( 3.0,  0.0 );
   A(1,2) = cplx( 4.0,  2.0 );
   A(2,3) = cplx( 5.0,  3.0 );

   // Setting all elements in the 1st row to 0 results in the matrix
   //
   //       ( (0, 0) (0,0) (0,0) (2, 1) )
   //   A = ( (0, 0) (0,0) (0,0) (0, 0) )
   //       ( (0, 0) (0,0) (0,0) (5,-3) )
   //       ( (2,-1) (0,0) (5,3) (0, 0) )
   //
   row( A, 1 ) = cplx( 0.0, 0.0 );
   \endcode

// The next example demonstrates the (compound) assignment to submatrices of Hermitian matrices.
// Since the modification of element \f$ a_{ij} \f$ of an Hermitian matrix also modifies the
// element \f$ a_{ji} \f$, the matrix to be assigned must be structured such that the Hermitian
// symmetry of the matrix is preserved. Otherwise a \c std::invalid_argument exception is thrown:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;

	std::complex<double>  cplx;

   // Setup of two default 4x4 Hermitian matrices
   HermitianMatrix< DynamicMatrix<cplx> > A1( 4 ), A2( 4 );

   // Setup of the 3x2 dynamic matrix
   //
   //       ( (1,-1) (2, 5) )
   //   B = ( (3, 0) (4,-6) )
   //       ( (5, 0) (6, 0) )
   //
   DynamicMatrix<int> B( 3UL, 2UL );
   B(0,0) = cplx( 1.0, -1.0 );
   B(0,1) = cplx( 2.0,  5.0 );
   B(1,0) = cplx( 3.0,  0.0 );
   B(1,1) = cplx( 4.0, -6.0 );
   B(2,1) = cplx( 5.0,  0.0 );
   B(2,2) = cplx( 6.0,  7.0 );

   // OK: Assigning B to a submatrix of A1 such that the Hermitian property is preserved
   //
   //        ( (0, 0) (0, 0) (1,-1) (2, 5) )
   //   A1 = ( (0, 0) (0, 0) (3, 0) (4,-6) )
   //        ( (1, 1) (3, 0) (5, 0) (6, 0) )
   //        ( (2,-5) (4, 6) (6, 0) (0, 0) )
   //
   submatrix( A1, 0UL, 2UL, 3UL, 2UL ) = B;  // OK

   // Error: Assigning B to a submatrix of A2 such that the Hermitian property isn't preserved!
   //   The elements marked with X cannot be assigned unambiguously!
   //
   //        ( (0, 0) (1,-1) (2,5) (0,0) )
   //   A2 = ( (1, 1) (3, 0) (X,X) (0,0) )
   //        ( (2,-5) (X, X) (6,0) (0,0) )
   //        ( (0, 0) (0, 0) (0,0) (0,0) )
   //
   submatrix( A2, 0UL, 1UL, 3UL, 2UL ) = B;  // Assignment throws an exception!
   \endcode

// \n \subsection adaptors_hermitian_matrices_initialization The Elements of a Dense Hermitian Matrix are Always Default Initialized!
//
// Although this results in a small loss of efficiency (especially in case all default values are
// overridden afterwards), this property is important since otherwise the Hermitian property of
// dense Hermitian matrices could not be guaranteed:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;

   // Uninitialized, 5x5 row-major dynamic matrix
   DynamicMatrix<int,rowMajor> A( 5, 5 );

   // Default initialized, 5x5 row-major Hermitian dynamic matrix
   HermitianMatrix< DynamicMatrix<int,rowMajor> > B( 5 );
   \endcode

// \n \section adaptors_hermitian_matrices_arithmetic_operations Arithmetic Operations
// <hr>
//
// An HermitianMatrix can be used within all numerical operations in any way any other dense or
// sparse matrix can be used. It can also be combined with any other dense or sparse vector or
// matrix. The following code example gives an impression of the use of HermitianMatrix within
// arithmetic operations:

   \code
   using blaze::HermitianMatrix;
   using blaze::DynamicMatrix;
   using blaze::HybridMatrix;
   using blaze::StaticMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

	using cplx = complex<float>;

   DynamicMatrix<cplx,rowMajor> A( 3, 3 );
   CompressedMatrix<cplx,rowMajor> B( 3, 3 );

   HermitianMatrix< DynamicMatrix<cplx,rowMajor> > C( 3 );
   HermitianMatrix< CompressedMatrix<cplx,rowMajor> > D( 3 );

   HermitianMatrix< HybridMatrix<cplx,3UL,3UL,rowMajor> > E;
   HermitianMatrix< StaticMatrix<cplx,3UL,3UL,columnMajor> > F;

   E = A + B;     // Matrix addition and assignment to a row-major Hermitian matrix (includes runtime check)
   F = C - D;     // Matrix subtraction and assignment to a column-major Hermitian matrix (only compile time check)
   F = A * D;     // Matrix multiplication between a dense and a sparse matrix (includes runtime check)

   C *= 2.0;      // In-place scaling of matrix C
   E  = 2.0 * B;  // Scaling of matrix B (includes runtime check)
   F  = C * 2.0;  // Scaling of matrix C (only compile time check)

   E += A - B;    // Addition assignment (includes runtime check)
   F -= C + D;    // Subtraction assignment (only compile time check)
   F *= A * D;    // Multiplication assignment (includes runtime check)
   \endcode

// Note that it is possible to assign any kind of matrix to an Hermitian matrix. In case the matrix
// to be assigned is not Hermitian at compile time, a runtime check is performed.
//
//
// \n \section adaptors_hermitian_matrices_performance Performance Considerations
// <hr>
//
// When the Hermitian property of a matrix is known beforehands using the HermitianMatrix adaptor
// instead of a general matrix can be a considerable performance advantage. This is particularly
// true in case the Hermitian matrix is also symmetric (i.e. has built-in element types). The
// \b Blaze library tries to exploit the properties of Hermitian (symmetric) matrices whenever
// possible. However, there are also situations when using an Hermitian matrix introduces some
// overhead. The following examples demonstrate several situations where Hermitian matrices can
// positively or negatively impact performance.
//
// \n \subsection adaptors_hermitian_matrices_matrix_matrix_multiplication Positive Impact: Matrix/Matrix Multiplication
//
// When multiplying two matrices, at least one of which is symmetric, \b Blaze can exploit the fact
// that \f$ A = A^T \f$ and choose the fastest and most suited combination of storage orders for the
// multiplication. The following example demonstrates this by means of a dense matrix/sparse matrix
// multiplication:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   HermitianMatrix< DynamicMatrix<double,rowMajor> > A;        // Both Hermitian and symmetric
   HermitianMatrix< CompressedMatrix<double,columnMajor> > B;  // Both Hermitian and symmetric
   DynamicMatrix<double,columnMajor> C;

   // ... Resizing and initialization

   C = A * B;
   \endcode

// Intuitively, the chosen combination of a row-major and a column-major matrix is the most suited
// for maximum performance. However, \b Blaze evaluates the multiplication as

   \code
   C = A * trans( B );
   \endcode

// which significantly increases the performance since in contrast to the original formulation the
// optimized form can be vectorized. Therefore, in the context of matrix multiplications, using a
// symmetric matrix is obviously an advantage.
//
// \n \subsection adaptors_hermitian_matrices_matrix_vector_multiplication Positive Impact: Matrix/Vector Multiplication
//
// A similar optimization is possible in case of matrix/vector multiplications:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
	using blaze::HermitianMatrix;
   using blaze::rowMajor;
   using blaze::columnVector;

   HermitianMatrix< DynamicMatrix<double,rowMajor> > A;  // Hermitian and symmetric
   CompressedVector<double,columnVector> x;
   DynamicVector<double,columnVector> y;

   // ... Resizing and initialization

   y = A * x;
   \endcode

// In this example it is not intuitively apparent that using a row-major matrix is not the best
// possible choice in terms of performance since the computation cannot be vectorized. Choosing
// a column-major matrix instead, however, would enable a vectorized computation. Therefore
// \b Blaze exploits the fact that \c A is symmetric, selects the best suited storage order and
// evaluates the multiplication as

   \code
   y = trans( A ) * x;
   \endcode

// which also significantly increases the performance.
//
// \n \subsection adaptors_hermitian_matrices_views Positive Impact: Row/Column Views on Column/Row-Major Matrices
//
// Another example is the optimization of a row view on a column-major symmetric matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;
   using blaze::columnMajor;

   HermitianMatrix< DynamicMatrix<double,columnMajor> > A( 10UL );  // Both Hermitian and symmetric
   auto row5 = row( A, 5UL );
   \endcode

// Usually, a row view on a column-major matrix results in a considerable performance decrease in
// comparison to a row view on a row-major matrix due to the non-contiguous storage of the matrix
// elements. However, in case of symmetric matrices, \b Blaze instead uses the according column of
// the matrix, which provides the same performance as if the matrix would be row-major. Note that
// this also works for column views on row-major matrices, where \b Blaze can use the according
// row instead of a column in order to provide maximum performance.
//
// \n \subsection adaptors_hermitian_matrices_assignment Negative Impact: Assignment of a General Matrix
//
// In contrast to using an Hermitian matrix on the right-hand side of an assignment (i.e. for read
// access), which introduces absolutely no performance penalty, using an Hermitian matrix on the
// left-hand side of an assignment (i.e. for write access) may introduce additional overhead when
// it is assigned a general matrix, which is not Hermitian at compile time:

   \code
   using blaze::DynamicMatrix;
   using blaze::HermitianMatrix;

   HermitianMatrix< DynamicMatrix< complex<double> > > A, C;
   DynamicMatrix<double> B;

   B = A;  // Only read-access to the Hermitian matrix; no performance penalty
   C = A;  // Assignment of an Hermitian matrix to another Hermitian matrix; no runtime overhead
   C = B;  // Assignment of a general matrix to an Hermitian matrix; some runtime overhead
   \endcode

// When assigning a general, potentially not Hermitian matrix to an Hermitian matrix it is necessary
// to check whether the matrix is Hermitian at runtime in order to guarantee the Hermitian property
// of the Hermitian matrix. In case it turns out to be Hermitian, it is assigned as efficiently as
// possible, if it is not, an exception is thrown. In order to prevent this runtime overhead it is
// therefore generally advisable to assign Hermitian matrices to other Hermitian matrices.\n
// In this context it is especially noteworthy that in contrast to additions and subtractions the
// multiplication of two Hermitian matrices does not necessarily result in another Hermitian matrix:

   \code
   HermitianMatrix< DynamicMatrix<double> > A, B, C;

   C = A + B;  // Results in an Hermitian matrix; no runtime overhead
   C = A - B;  // Results in an Hermitian matrix; no runtime overhead
   C = A * B;  // Is not guaranteed to result in an Hermitian matrix; some runtime overhead
   \endcode

// \n Previous: \ref adaptors_symmetric_matrices &nbsp; &nbsp; Next: \ref adaptors_triangular_matrices
*/
//*************************************************************************************************


//**Triangular Matrices****************************************************************************
/*!\page adaptors_triangular_matrices Triangular Matrices
//
// \tableofcontents
//
//
// \n \section adaptors_triangular_matrices_general Triangular Matrices
// <hr>
//
// Triangular matrices come in three flavors: Lower triangular matrices provide the compile time
// guarantee to be square matrices and that the upper part of the matrix contains only default
// elements that cannot be modified. Upper triangular matrices on the other hand provide the
// compile time guarantee to be square and that the lower part of the matrix contains only fixed
// default elements. Finally, diagonal matrices provide the compile time guarantee to be square
// and that both the lower and upper part of the matrix contain only immutable default elements.
// These properties can be exploited to gain higher performance and/or to save memory. Within the
// \b Blaze library, several kinds of lower and upper triangular and diagonal matrices are realized
// by the following class templates:
//
// Lower triangular matrices:
//  - <b>\ref adaptors_triangular_matrices_lowermatrix</b>
//  - <b>\ref adaptors_triangular_matrices_unilowermatrix</b>
//  - <b>\ref adaptors_triangular_matrices_strictlylowermatrix</b>
//
// Upper triangular matrices:
//  - <b>\ref adaptors_triangular_matrices_uppermatrix</b>
//  - <b>\ref adaptors_triangular_matrices_uniuppermatrix</b>
//  - <b>\ref adaptors_triangular_matrices_strictlyuppermatrix</b>
//
// Diagonal matrices
//  - <b>\ref adaptors_triangular_matrices_diagonalmatrix</b>
//
//
// \n \section adaptors_triangular_matrices_lowermatrix LowerMatrix
// <hr>
//
// The blaze::LowerMatrix class template is an adapter for existing dense and sparse matrix types.
// It inherits the properties and the interface of the given matrix type \c MT and extends it by
// enforcing the additional invariant that all matrix elements above the diagonal are 0 (lower
// triangular matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        l_{1,0} & l_{1,1} & 0       & \cdots & 0       \\
                        l_{2,0} & l_{2,1} & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/LowerMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class LowerMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::LowerMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Note
// that the given matrix type must be either resizable (as for instance blaze::HybridMatrix or
// blaze::DynamicMatrix) or must be square at compile time (as for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible lower matrices:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of a 3x3 row-major dense lower matrix with static memory
   blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,rowMajor> > A;

   // Definition of a resizable column-major dense lower matrix based on HybridMatrix
   blaze::LowerMatrix< blaze::HybridMatrix<float,4UL,4UL,columnMajor> B;

   // Definition of a resizable row-major dense lower matrix based on DynamicMatrix
   blaze::LowerMatrix< blaze::DynamicMatrix<double,rowMajor> > C;

   // Definition of a fixed size row-major dense lower matrix based on CustomMatrix
   blaze::LowerMatrix< blaze::CustomMatrix<double,unaligned,unpadded,rowMajor> > D;

   // Definition of a compressed row-major single precision lower matrix
   blaze::LowerMatrix< blaze::CompressedMatrix<float,rowMajor> > E;
   \endcode

// The storage order of a lower matrix is depending on the storage order of the adapted matrix
// type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e. is specified
// as blaze::rowMajor), the lower matrix will also be a row-major matrix. Otherwise, if the
// adapted matrix is column-major (i.e. is specified as blaze::columnMajor), the lower matrix
// will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_unilowermatrix UniLowerMatrix
// <hr>
//
// The blaze::UniLowerMatrix class template is an adapter for existing dense and sparse matrix
// types. It inherits the properties and the interface of the given matrix type \c MT and extends
// it by enforcing the additional invariant that all diagonal matrix elements are 1 and all matrix
// elements above the diagonal are 0 (lower unitriangular matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 1       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 1       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 1      \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/UniLowerMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class UniLowerMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::UniLowerMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Also,
// the given matrix type must have numeric element types (i.e. all integral types except \c bool,
// floating point and complex types). Note that the given matrix type must be either resizable (as
// for instance blaze::HybridMatrix or blaze::DynamicMatrix) or must be square at compile time (as
// for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible lower unitriangular matrices:

   \code
   // Definition of a 3x3 row-major dense unilower matrix with static memory
   blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense unilower matrix based on HybridMatrix
   blaze::UniLowerMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense unilower matrix based on DynamicMatrix
   blaze::UniLowerMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision unilower matrix
   blaze::UniLowerMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of a lower unitriangular matrix is depending on the storage order of the
// adapted matrix type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e.
// is specified as blaze::rowMajor), the unilower matrix will also be a row-major matrix.
// Otherwise if the adapted matrix is column-major (i.e. is specified as blaze::columnMajor),
// the unilower matrix will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_strictlylowermatrix StrictlyLowerMatrix
// <hr>
//
// The blaze::StrictlyLowerMatrix class template is an adapter for existing dense and sparse matrix
// types. It inherits the properties and the interface of the given matrix type \c MT and extends
// it by enforcing the additional invariant that all diagonal matrix elements and all matrix
// elements above the diagonal are 0 (strictly lower triangular matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 0       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 0       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 0      \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/StrictlyLowerMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class StrictlyLowerMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::StrictlyLowerMatrix can be used
// with any non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix
// type. Note that the given matrix type must be either resizable (as for instance
// blaze::HybridMatrix or blaze::DynamicMatrix) or must be square at compile time (as for instance
// blaze::StaticMatrix).
//
// The following examples give an impression of several possible strictly lower triangular matrices:

   \code
   // Definition of a 3x3 row-major dense strictly lower matrix with static memory
   blaze::StrictlyLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense strictly lower matrix based on HybridMatrix
   blaze::StrictlyLowerMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense strictly lower matrix based on DynamicMatrix
   blaze::StrictlyLowerMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision strictly lower matrix
   blaze::StrictlyLowerMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of a strictly lower triangular matrix is depending on the storage order of
// the adapted matrix type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e.
// is specified as blaze::rowMajor), the strictly lower matrix will also be a row-major matrix.
// Otherwise if the adapted matrix is column-major (i.e. is specified as blaze::columnMajor),
// the strictly lower matrix will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_uppermatrix UpperMatrix
// <hr>
//
// The blaze::UpperMatrix class template is an adapter for existing dense and sparse matrix types.
// It inherits the properties and the interface of the given matrix type \c MT and extends it by
// enforcing the additional invariant that all matrix elements below the diagonal are 0 (upper
// triangular matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        u_{0,0} & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & u_{1,1} & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & u_{2,2} & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & u_{N,N} \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/UpperMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class UpperMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::UpperMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Note
// that the given matrix type must be either resizable (as for instance blaze::HybridMatrix or
// blaze::DynamicMatrix) or must be square at compile time (as for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible upper matrices:

   \code
   // Definition of a 3x3 row-major dense upper matrix with static memory
   blaze::UpperMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense upper matrix based on HybridMatrix
   blaze::UpperMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense upper matrix based on DynamicMatrix
   blaze::UpperMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision upper matrix
   blaze::UpperMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of an upper matrix is depending on the storage order of the adapted matrix
// type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e. is specified
// as blaze::rowMajor), the upper matrix will also be a row-major matrix. Otherwise, if the
// adapted matrix is column-major (i.e. is specified as blaze::columnMajor), the upper matrix
// will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_uniuppermatrix UniUpperMatrix
// <hr>
//
// The blaze::UniUpperMatrix class template is an adapter for existing dense and sparse matrix
// types. It inherits the properties and the interface of the given matrix type \c MT and extends
// it by enforcing the additional invariant that all diagonal matrix elements are 1 and all matrix
// elements below the diagonal are 0 (upper unitriangular matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & 1       & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & 1       & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & 1       \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/UniUpperMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class UniUpperMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::UniUpperMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Also,
// the given matrix type must have numeric element types (i.e. all integral types except \c bool,
// floating point and complex types). Note that the given matrix type must be either resizable (as
// for instance blaze::HybridMatrix or blaze::DynamicMatrix) or must be square at compile time (as
// for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible upper unitriangular matrices:

   \code
   // Definition of a 3x3 row-major dense uniupper matrix with static memory
   blaze::UniUpperMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense uniupper matrix based on HybridMatrix
   blaze::UniUpperMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense uniupper matrix based on DynamicMatrix
   blaze::UniUpperMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision uniupper matrix
   blaze::UniUpperMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of an upper unitriangular matrix is depending on the storage order of the
// adapted matrix type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e.
// is specified as blaze::rowMajor), the uniupper matrix will also be a row-major matrix.
// Otherwise, if the adapted matrix is column-major (i.e. is specified as blaze::columnMajor),
// the uniupper matrix will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_strictlyuppermatrix StrictlyUpperMatrix
// <hr>
//
// The blaze::StrictlyUpperMatrix class template is an adapter for existing dense and sparse matrix
// types. It inherits the properties and the interface of the given matrix type \c MT and extends
// it by enforcing the additional invariant that all diagonal matrix elements and all matrix
// elements below the diagonal are 0 (strictly upper triangular matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & 0       & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & 0       & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & 0       \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/StrictlyUpperMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class StrictlyUpperMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::StrictlyUpperMatrix can be used
// with any non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix
// type. Note that the given matrix type must be either resizable (as for instance
// blaze::HybridMatrix or blaze::DynamicMatrix) or must be square at compile time (as for instance
// blaze::StaticMatrix).
//
// The following examples give an impression of several possible strictly upper triangular matrices:

   \code
   // Definition of a 3x3 row-major dense strictly upper matrix with static memory
   blaze::StrictlyUpperMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense strictly upper matrix based on HybridMatrix
   blaze::StrictlyUpperMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense strictly upper matrix based on DynamicMatrix
   blaze::StrictlyUpperMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision strictly upper matrix
   blaze::StrictlyUpperMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of a strictly upper triangular matrix is depending on the storage order of
// the adapted matrix type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e.
// is specified as blaze::rowMajor), the strictly upper matrix will also be a row-major matrix.
// Otherwise, if the adapted matrix is column-major (i.e. is specified as blaze::columnMajor),
// the strictly upper matrix will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_diagonalmatrix DiagonalMatrix
// <hr>
//
// The blaze::DiagonalMatrix class template is an adapter for existing dense and sparse matrix
// types. It inherits the properties and the interface of the given matrix type \c MT and extends
// it by enforcing the additional invariant that all matrix elements above and below the diagonal
// are 0 (diagonal matrix):

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        0       & l_{1,1} & 0       & \cdots & 0       \\
                        0       & 0       & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/DiagonalMatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The type of the adapted matrix can be specified via the first template parameter:

   \code
   namespace blaze {

   template< typename MT >
   class DiagonalMatrix;

   } // namespace blaze
   \endcode

// \c MT specifies the type of the matrix to be adapted. blaze::DiagonalMatrix can be used with any
// non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix type. Note
// that the given matrix type must be either resizable (as for instance blaze::HybridMatrix or
// blaze::DynamicMatrix) or must be square at compile time (as for instance blaze::StaticMatrix).
//
// The following examples give an impression of several possible diagonal matrices:

   \code
   // Definition of a 3x3 row-major dense diagonal matrix with static memory
   blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense diagonal matrix based on HybridMatrix
   blaze::DiagonalMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense diagonal matrix based on DynamicMatrix
   blaze::DiagonalMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision diagonal matrix
   blaze::DiagonalMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of a diagonal matrix is depending on the storage order of the adapted matrix
// type \c MT. In case the adapted matrix is stored in a row-wise fashion (i.e. is specified
// as blaze::rowMajor), the diagonal matrix will also be a row-major matrix. Otherwise, if the
// adapted matrix is column-major (i.e. is specified as blaze::columnMajor), the diagonal matrix
// will also be a column-major matrix.
//
//
// \n \section adaptors_triangular_matrices_special_properties Special Properties of Triangular Matrices
// <hr>
//
// A triangular matrix is used exactly like a matrix of the underlying, adapted matrix type \c MT.
// It also provides (nearly) the same interface as the underlying matrix type. However, there are
// some important exceptions resulting from the triangular matrix constraint:
//
//  -# <b>\ref adaptors_triangular_matrices_square</b>
//  -# <b>\ref adaptors_triangular_matrices_triangular</b>
//  -# <b>\ref adaptors_triangular_matrices_initialization</b>
//  -# <b>\ref adaptors_triangular_matrices_storage</b>
//  -# <b>\ref adaptors_triangular_matrices_scaling</b>
//
// \n \subsection adaptors_triangular_matrices_square Triangular Matrices Must Always be Square!
//
// In case a resizable matrix is used (as for instance blaze::HybridMatrix, blaze::DynamicMatrix,
// or blaze::CompressedMatrix), this means that the according constructors, the \c resize() and
// the \c extend() functions only expect a single parameter, which specifies both the number of
// rows and columns, instead of two (one for the number of rows and one for the number of columns):

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::rowMajor;

   // Default constructed, default initialized, row-major 3x3 lower dynamic matrix
   LowerMatrix< DynamicMatrix<double,rowMajor> > A( 3 );

   // Resizing the matrix to 5x5
   A.resize( 5 );

   // Extending the number of rows and columns by 2, resulting in a 7x7 matrix
   A.extend( 2 );
   \endcode

// In case a matrix with a fixed size is used (as for instance blaze::StaticMatrix), the number
// of rows and number of columns must be specified equally:

   \code
   using blaze::StaticMatrix;
   using blaze::LowerMatrix;
   using blaze::columnMajor;

   // Correct setup of a fixed size column-major 3x3 lower static matrix
   LowerMatrix< StaticMatrix<int,3UL,3UL,columnMajor> > A;

   // Compilation error: the provided matrix type is not a square matrix type
   LowerMatrix< StaticMatrix<int,3UL,4UL,columnMajor> > B;
   \endcode

// \n \subsection adaptors_triangular_matrices_triangular The Triangular Property is Always Enforced!
//
// This means that it is only allowed to modify elements in the lower part or the diagonal of
// a lower triangular matrix and in the upper part or the diagonal of an upper triangular matrix.
// Unitriangular and strictly triangular matrices are even more restrictive and don't allow the
// modification of diagonal elements. Also, triangular matrices can only be assigned matrices that
// don't violate their triangular property. The following example demonstrates this restriction
// by means of the blaze::LowerMatrix adaptor. For examples with other triangular matrix types
// see the according class documentations.

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::LowerMatrix;
   using blaze::rowMajor;

   using CompressedLower = LowerMatrix< CompressedMatrix<double,rowMajor> >;

   // Default constructed, row-major 3x3 lower compressed matrix
   CompressedLower A( 3 );

   // Initializing elements via the function call operator
   A(0,0) = 1.0;  // Initialization of the diagonal element (0,0)
   A(2,0) = 2.0;  // Initialization of the lower element (2,0)
   A(1,2) = 9.0;  // Throws an exception; invalid modification of upper element

   // Inserting two more elements via the insert() function
   A.insert( 1, 0, 3.0 );  // Inserting the lower element (1,0)
   A.insert( 2, 1, 4.0 );  // Inserting the lower element (2,1)
   A.insert( 0, 2, 9.0 );  // Throws an exception; invalid insertion of upper element

   // Appending an element via the append() function
   A.reserve( 1, 3 );      // Reserving enough capacity in row 1
   A.append( 1, 1, 5.0 );  // Appending the diagonal element (1,1)
   A.append( 1, 2, 9.0 );  // Throws an exception; appending an element in the upper part

   // Access via a non-const iterator
   CompressedLower::Iterator it = A.begin(1);
   *it = 6.0;  // Modifies the lower element (1,0)
   ++it;
   *it = 9.0;  // Modifies the diagonal element (1,1)

   // Erasing elements via the erase() function
   A.erase( 0, 0 );  // Erasing the diagonal element (0,0)
   A.erase( 2, 0 );  // Erasing the lower element (2,0)

   // Construction from a lower dense matrix
   StaticMatrix<double,3UL,3UL> B{ {  3.0,  0.0,  0.0 },
                                   {  8.0,  0.0,  0.0 },
                                   { -2.0, -1.0,  4.0 } };

   LowerMatrix< DynamicMatrix<double,rowMajor> > C( B );  // OK

   // Assignment of a non-lower dense matrix
   StaticMatrix<double,3UL,3UL> D{ {  3.0,  0.0, -2.0 },
                                   {  8.0,  0.0,  0.0 },
                                   { -2.0, -1.0,  4.0 } };

   C = D;  // Throws an exception; lower matrix invariant would be violated!
   \endcode

// The triangular property is also enforced during the construction of triangular custom matrices:
// In case the given array of elements does not represent the according triangular matrix type, a
// \c std::invalid_argument exception is thrown:

   \code
   using blaze::CustomMatrix;
   using blaze::LowerMatrix;
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   using CustomLower = LowerMatrix< CustomMatrix<double,unaligned,unpadded,rowMajor> >;

   // Creating a 3x3 lower custom matrix from a properly initialized array
   double array[9] = { 1.0, 0.0, 0.0,
                       2.0, 3.0, 0.0,
                       4.0, 5.0, 6.0 };
   CustomLower A( array, 3UL );  // OK

   // Attempt to create a second 3x3 lower custom matrix from an uninitialized array
   std::unique_ptr<double[]> memory( new double[9UL] );
   CustomLower B( memory.get(), 3UL );  // Throws an exception
   \endcode

// Finally, the triangular matrix property is enforced for views (rows, columns, submatrices, ...)
// on the triangular matrix. The following example demonstrates that modifying the elements of an
// entire row and submatrix of a lower matrix only affects the lower and diagonal matrix elements.
// Again, this example uses blaze::LowerMatrix, for examples with other triangular matrix types
// see the according class documentations.

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;

   // Setup of the lower matrix
   //
   //       ( 0 0 0 0 )
   //   A = ( 1 2 0 0 )
   //       ( 0 3 0 0 )
   //       ( 4 0 5 0 )
   //
   LowerMatrix< DynamicMatrix<int> > A( 4 );
   A(1,0) = 1;
   A(1,1) = 2;
   A(2,1) = 3;
   A(3,0) = 4;
   A(3,2) = 5;

   // Setting the lower and diagonal elements in the 2nd row to 9 results in the matrix
   //
   //       ( 0 0 0 0 )
   //   A = ( 1 2 0 0 )
   //       ( 9 9 9 0 )
   //       ( 4 0 5 0 )
   //
   row( A, 2 ) = 9;

   // Setting the lower and diagonal elements in the 1st and 2nd column to 7 results in
   //
   //       ( 0 0 0 0 )
   //   A = ( 1 7 0 0 )
   //       ( 9 7 7 0 )
   //       ( 4 7 7 0 )
   //
   submatrix( A, 0, 1, 4, 2 ) = 7;
   \endcode

// The next example demonstrates the (compound) assignment to rows/columns and submatrices of
// triangular matrices. Since only lower/upper and potentially diagonal elements may be modified
// the matrix to be assigned must be structured such that the triangular matrix invariant of the
// matrix is preserved. Otherwise a \c std::invalid_argument exception is thrown:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::LowerMatrix;
   using blaze::rowVector;

   // Setup of two default 4x4 lower matrices
   LowerMatrix< DynamicMatrix<int> > A1( 4 ), A2( 4 );

   // Setup of a 4-dimensional vector
   //
   //   v = ( 1 2 3 0 )
   //
   DynamicVector<int,rowVector> v{ 1, 2, 3, 0 };

   // OK: Assigning v to the 2nd row of A1 preserves the lower matrix invariant
   //
   //        ( 0 0 0 0 )
   //   A1 = ( 0 0 0 0 )
   //        ( 1 2 3 0 )
   //        ( 0 0 0 0 )
   //
   row( A1, 2 ) = v;  // OK

   // Error: Assigning v to the 1st row of A1 violates the lower matrix invariant! The element
   //   marked with X cannot be assigned and triggers an exception.
   //
   //        ( 0 0 0 0 )
   //   A1 = ( 1 2 X 0 )
   //        ( 1 2 3 0 )
   //        ( 0 0 0 0 )
   //
   row( A1, 1 ) = v;  // Assignment throws an exception!

   // Setup of the 3x2 dynamic matrix
   //
   //       ( 0 0 )
   //   B = ( 7 0 )
   //       ( 8 9 )
   //
   DynamicMatrix<int> B( 3UL, 2UL, 0 );
   B(1,0) = 7;
   B(2,0) = 8;
   B(2,1) = 9;

   // OK: Assigning B to a submatrix of A2 such that the lower matrix invariant can be preserved
   //
   //        ( 0 0 0 0 )
   //   A2 = ( 0 7 0 0 )
   //        ( 0 8 9 0 )
   //        ( 0 0 0 0 )
   //
   submatrix( A2, 0UL, 1UL, 3UL, 2UL ) = B;  // OK

   // Error: Assigning B to a submatrix of A2 such that the lower matrix invariant cannot be
   //   preserved! The elements marked with X cannot be assigned without violating the invariant!
   //
   //        ( 0 0 0 0 )
   //   A2 = ( 0 7 X 0 )
   //        ( 0 8 8 X )
   //        ( 0 0 0 0 )
   //
   submatrix( A2, 0UL, 2UL, 3UL, 2UL ) = B;  // Assignment throws an exception!
   \endcode

// \n \subsection adaptors_triangular_matrices_initialization The Elements of a Dense Triangular Matrix are Always Default Initialized!
//
// Although this results in a small loss of efficiency during the creation of a dense lower or
// upper matrix this initialization is important since otherwise the lower/upper matrix property
// of dense lower matrices would not be guaranteed:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;

   // Uninitialized, 5x5 row-major dynamic matrix
   DynamicMatrix<int,rowMajor> A( 5, 5 );

   // 5x5 row-major lower dynamic matrix with default initialized upper matrix
   LowerMatrix< DynamicMatrix<int,rowMajor> > B( 5 );

   // 7x7 column-major upper dynamic matrix with default initialized lower matrix
   UpperMatrix< DynamicMatrix<int,columnMajor> > C( 7 );

   // 3x3 row-major diagonal dynamic matrix with default initialized lower and upper matrix
   DiagonalMatrix< DynamicMatrix<int,rowMajor> > D( 3 );
   \endcode

// \n \subsection adaptors_triangular_matrices_storage Dense Triangular Matrices Store All Elements!
//
// All dense triangular matrices store all \f$ N \times N \f$ elements, including the immutable
// elements in the lower or upper part, respectively. Therefore dense triangular matrices don't
// provide any kind of memory reduction! There are two main reasons for this: First, storing also
// the zero elements guarantees maximum performance for many algorithms that perform vectorized
// operations on the triangular matrices, which is especially true for small dense matrices.
// Second, conceptually all triangular adaptors merely restrict the interface to the matrix type
// \c MT and do not change the data layout or the underlying matrix type.
//
// This property matters most for diagonal matrices. In order to achieve the perfect combination
// of performance and memory consumption for a diagonal matrix it is recommended to use dense
// matrices for small diagonal matrices and sparse matrices for large diagonal matrices:

   \code
   // Recommendation 1: use dense matrices for small diagonal matrices
   using SmallDiagonalMatrix = blaze::DiagonalMatrix< blaze::StaticMatrix<float,3UL,3UL> >;

   // Recommendation 2: use sparse matrices for large diagonal matrices
   using LargeDiagonalMatrix = blaze::DiagonalMatrix< blaze::CompressedMatrix<float> >;
   \endcode

// \n \subsection adaptors_triangular_matrices_scaling Unitriangular Matrices Cannot Be Scaled!
//
// Since the diagonal elements of a unitriangular matrix have a fixed value of 1 it is not possible
// to self-scale such a matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::UniLowerMatrix;

   UniLowerMatrix< DynamicMatrix<int> > A( 4 );

   A *= 2;        // Compilation error; Scale operation is not available on an unilower matrix
   A /= 2;        // Compilation error; Scale operation is not available on an unilower matrix
   A.scale( 2 );  // Compilation error; Scale function is not available on an unilower matrix

   A = A * 2;  // Throws an exception; Invalid assignment of non-unilower matrix
   A = A / 2;  // Throws an exception; Invalid assignment of non-unilower matrix
   \endcode

// \n \section adaptors_triangular_matrices_arithmetic_operations Arithmetic Operations
// <hr>
//
// A lower and upper triangular matrix can participate in numerical operations in any way any other
// dense or sparse matrix can participate. It can also be combined with any other dense or sparse
// vector or matrix. The following code example gives an impression of the use of blaze::LowerMatrix
// within arithmetic operations:

   \code
   using blaze::LowerMatrix;
   using blaze::DynamicMatrix;
   using blaze::HybridMatrix;
   using blaze::StaticMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   DynamicMatrix<double,rowMajor> A( 3, 3 );
   CompressedMatrix<double,rowMajor> B( 3, 3 );

   LowerMatrix< DynamicMatrix<double,rowMajor> > C( 3 );
   LowerMatrix< CompressedMatrix<double,rowMajor> > D( 3 );

   LowerMatrix< HybridMatrix<float,3UL,3UL,rowMajor> > E;
   LowerMatrix< StaticMatrix<float,3UL,3UL,columnMajor> > F;

   E = A + B;     // Matrix addition and assignment to a row-major lower matrix (includes runtime check)
   F = C - D;     // Matrix subtraction and assignment to a column-major lower matrix (only compile time check)
   F = A * D;     // Matrix multiplication between a dense and a sparse matrix (includes runtime check)

   C *= 2.0;      // In-place scaling of matrix C
   E  = 2.0 * B;  // Scaling of matrix B (includes runtime check)
   F  = C * 2.0;  // Scaling of matrix C (only compile time check)

   E += A - B;    // Addition assignment (includes runtime check)
   F -= C + D;    // Subtraction assignment (only compile time check)
   F *= A * D;    // Multiplication assignment (includes runtime check)
   \endcode

// Note that it is possible to assign any kind of matrix to a triangular matrix. In case the
// matrix to be assigned does not satisfy the invariants of the triangular matrix at compile
// time, a runtime check is performed. Also note that upper triangular, diagonal, unitriangular
// and strictly triangular matrix types can be used in the same way, but may pose some additional
// restrictions (see the according class documentations).
//
//
// \n \section adaptors_triangular_matrices_block_matrices Triangular Block Matrices
// <hr>
//
// It is also possible to use triangular block matrices:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;

   // Definition of a 5x5 lower block matrix based on DynamicMatrix
   LowerMatrix< DynamicMatrix< StaticMatrix<int,3UL,3UL> > > A( 5 );

   // Definition of a 7x7 upper block matrix based on CompressedMatrix
   UpperMatrix< CompressedMatrix< StaticMatrix<int,3UL,3UL> > > B( 7 );
   \endcode

// Also in this case the triangular matrix invariant is enforced, i.e. it is not possible to
// manipulate elements in the upper part (lower triangular matrix) or the lower part (upper
// triangular matrix) of the matrix:

   \code
   const StaticMatrix<int,3UL,3UL> C{ { 1, -4,  5 },
                                      { 6,  8, -3 },
                                      { 2, -1,  2 } };

   A(2,4)(1,1) = -5;     // Invalid manipulation of upper matrix element; Results in an exception
   B.insert( 4, 2, C );  // Invalid insertion of the elements (4,2); Results in an exception
   \endcode

// Note that unitriangular matrices are restricted to numeric element types and therefore cannot
// be used for block matrices:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::UniLowerMatrix;
   using blaze::UniUpperMatrix;

   // Compilation error: lower unitriangular matrices are restricted to numeric element types
   UniLowerMatrix< DynamicMatrix< StaticMatrix<int,3UL,3UL> > > A( 5 );

   // Compilation error: upper unitriangular matrices are restricted to numeric element types
   UniUpperMatrix< CompressedMatrix< StaticMatrix<int,3UL,3UL> > > B( 7 );
   \endcode

// For more information on block matrices, see the tutorial on \ref block_vectors_and_matrices.
//
//
// \n \section adaptors_triangular_matrices_performance Performance Considerations
// <hr>
//
// The \b Blaze library tries to exploit the properties of lower and upper triangular matrices
// whenever and wherever possible. Therefore using triangular matrices instead of a general
// matrices can result in a considerable performance improvement. However, there are also
// situations when using a triangular matrix introduces some overhead. The following examples
// demonstrate several common situations where triangular matrices can positively or negatively
// impact performance.
//
// \n \subsection adaptors_triangular_matrices_matrix_matrix_multiplication Positive Impact: Matrix/Matrix Multiplication
//
// When multiplying two matrices, at least one of which is triangular, \b Blaze can exploit the
// fact that either the lower or upper part of the matrix contains only default elements and
// restrict the algorithm to the non-zero elements. The following example demonstrates this by
// means of a dense matrix/dense matrix multiplication with lower triangular matrices:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   LowerMatrix< DynamicMatrix<double,rowMajor> > A;
   LowerMatrix< DynamicMatrix<double,columnMajor> > B;
   DynamicMatrix<double,columnMajor> C;

   // ... Resizing and initialization

   C = A * B;
   \endcode

// In comparison to a general matrix multiplication, the performance advantage is significant,
// especially for large matrices. Therefore is it highly recommended to use the blaze::LowerMatrix
// and blaze::UpperMatrix adaptors when a matrix is known to be lower or upper triangular,
// respectively. Note however that the performance advantage is most pronounced for dense matrices
// and much less so for sparse matrices.
//
// \n \subsection adaptors_triangular_matrices_matrix_vector_multiplication Positive Impact: Matrix/Vector Multiplication
//
// A similar performance improvement can be gained when using a triangular matrix in a matrix/vector
// multiplication:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   LowerMatrix< DynamicMatrix<double,rowMajor> > A;
   DynamicVector<double,columnVector> x, y;

   // ... Resizing and initialization

   y = A * x;
   \endcode

// In this example, \b Blaze also exploits the structure of the matrix and approx. halves the
// runtime of the multiplication. Also in case of matrix/vector multiplications the performance
// improvement is most pronounced for dense matrices and much less so for sparse matrices.
//
// \n \subsection adaptors_triangular_matrices_assignment Negative Impact: Assignment of a General Matrix
//
// In contrast to using a triangular matrix on the right-hand side of an assignment (i.e. for
// read access), which introduces absolutely no performance penalty, using a triangular matrix
// on the left-hand side of an assignment (i.e. for write access) may introduce additional
// overhead when it is assigned a general matrix, which is not triangular at compile time:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;

   LowerMatrix< DynamicMatrix<double> > A, C;
   DynamicMatrix<double> B;

   B = A;  // Only read-access to the lower matrix; no performance penalty
   C = A;  // Assignment of a lower matrix to another lower matrix; no runtime overhead
   C = B;  // Assignment of a general matrix to a lower matrix; some runtime overhead
   \endcode

// When assigning a general (potentially not lower triangular) matrix to a lower matrix or a
// general (potentially not upper triangular) matrix to an upper matrix it is necessary to check
// whether the matrix is lower or upper at runtime in order to guarantee the triangular property
// of the matrix. In case it turns out to be lower or upper, respectively, it is assigned as
// efficiently as possible, if it is not, an exception is thrown. In order to prevent this runtime
// overhead it is therefore generally advisable to assign lower or upper triangular matrices to
// other lower or upper triangular matrices.\n
// In this context it is especially noteworthy that the addition, subtraction, and multiplication
// of two triangular matrices of the same structure always results in another triangular matrix:

   \code
   LowerMatrix< DynamicMatrix<double> > A, B, C;

   C = A + B;  // Results in a lower matrix; no runtime overhead
   C = A - B;  // Results in a lower matrix; no runtime overhead
   C = A * B;  // Results in a lower matrix; no runtime overhead
   \endcode

   \code
   UpperMatrix< DynamicMatrix<double> > A, B, C;

   C = A + B;  // Results in an upper matrix; no runtime overhead
   C = A - B;  // Results in an upper matrix; no runtime overhead
   C = A * B;  // Results in an upper matrix; no runtime overhead
   \endcode

// \n Previous: \ref adaptors_hermitian_matrices &nbsp; &nbsp; Next: \ref views
*/
//*************************************************************************************************


//**Views******************************************************************************************
/*!\page views Views
//
// \tableofcontents
//
//
// \section views_general General Concepts
// <hr>
//
// Views represents parts of a vector or matrix, such as a subvector, a submatrix, or a specific
// row, column, or band of a matrix. As such, views act as a reference to specific elements of
// a vector or matrix. This reference is valid and can be used in every way as any other vector
// or matrix can be used as long as the referenced vector or matrix is not resized or entirely
// destroyed. Views also act as alias to the elements of the vector or matrix: Changes made to the
// elements (e.g. modifying values, inserting or erasing elements) via the view are immediately
// visible in the vector or matrix and changes made via the vector or matrix are immediately
// visible in the view.
//
// It is also possible to create nested views (compound views), such as for instance bands of
// submatrices or row selections on column selections. A compound view also acts as reference
// to specific elements of the underlying vector or matrix and is valid as long as the underlying,
// referenced vector or matrix is not resized or entirely destroyed.
//
// The \b Blaze library provides the following views on vectors and matrices:
//
// Vector views:
//  - \ref views_subvectors
//  - \ref views_element_selections
//
// Matrix views:
//  - \ref views_submatrices
//  - \ref views_rows
//  - \ref views_row_selections
//  - \ref views_columns
//  - \ref views_column_selections
//  - \ref views_bands
//
//
// \n \section views_examples Examples

   \code
   using blaze::DynamicMatrix;
   using blaze::StaticVector;

   // Setup of the 3x5 row-major matrix
   DynamicMatrix<int> A{ { 1,  0, -2,  3,  0 },
                         { 0,  2,  5, -1, -1 },
                         { 1,  0,  0,  2,  1 } };

   // Setup of the 2-dimensional row vector
   StaticVector<int,2UL,rowVector> vec{ 18, 19 };

   // Assigning to the elements (1,2) and (1,3) via a subvector of a row
   //
   //  ( 1  0 -2  3  0 )
   //  ( 0  2 18 19 -1 )
   //  ( 1  0  0  2  1 )
   //
   subvector( row( A, 1UL ), 2UL, 2UL ) = vec;

   // Switching rows 0 and 2 of A
   //
   //  ( 1  0  0  2  1 )
   //  ( 0  2 18 19 -1 )
   //  ( 1  0 -2  3  0 )
   //
   rows<0,2>( A ) = rows<2,0>( A );

   // Warning: It is the programmer's responsibility to ensure the view does not outlive
   //          the viewed vector or matrix (dangling reference)!
   auto row1 = row<1UL>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 } } );
   \endcode

// \n Previous: \ref adaptors_triangular_matrices &nbsp; &nbsp; Next: \ref views_subvectors
*/
//*************************************************************************************************


//**Subvectors*************************************************************************************
/*!\page views_subvectors Subvectors
//
// \tableofcontents
//
//
// Subvectors provide views on a specific part of a dense or sparse vector. As such, subvectors
// act as a reference to a specific range within a vector. This reference is valid and can be
// used in every way any other dense or sparse vector can be used as long as the vector containing
// the subvector is not resized or entirely destroyed. The subvector also acts as an alias to the
// vector elements in the specified range: Changes made to the elements (e.g. modifying values,
// inserting or erasing elements) are immediately visible in the vector and changes made via the
// vector are immediately visible in the subvector.
//
//
// \n \section views_subvectors_setup Setup of Subvectors
// <hr>
//
// A view on a dense or sparse subvector can be created very conveniently via the \c subvector()
// function. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Subvector.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The first parameter specifies the offset of the subvector within the underlying dense or sparse
// vector, the second parameter specifies the size of the subvector. The two parameters can be
// specified either at compile time or at runtime:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   // ... Resizing and initialization

   // Create a subvector from index 4 with a size of 12 (i.e. in the range [4..15]) (compile time arguments)
   auto sv1 = subvector<4UL,12UL>( x );

   // Create a subvector from index 8 with a size of 16 (i.e. in the range [8..23]) (runtime arguments)
   auto sv2 = subvector( x, 8UL, 16UL );
   \endcode

// The \c subvector() function returns an expression representing the subvector view. The type of
// this expression depends on the given subvector arguments, primarily the type of the vector and
// the compile time arguments. If the type is required, it can be determined via the \c decltype
// specifier:

   \code
   using VectorType = blaze::DynamicVector<int>;
   using SubvectorType = decltype( blaze::subvector<4UL,12UL>( std::declval<VectorType>() ) );
   \endcode

// The resulting view can be treated as any other dense or sparse vector, i.e. it can be assigned
// to, it can be copied from, and it can be used in arithmetic operations. A subvector created
// from a row vector can be used as any other row vector, a subvector created from a column vector
// can be used as any other column vector. The view can also be used on both sides of an assignment:
// The subvector can either be used as an alias to grant write access to a specific subvector of a
// vector primitive on the left-hand side of an assignment or to grant read-access to a specific
// subvector of a vector primitive or expression on the right-hand side of an assignment. The
// following example demonstrates this in detail:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   blaze::CompressedVector<double,blaze::rowVector> y;
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Create a subvector from index 0 with a size of 10 (i.e. in the range [0..9])
   auto sv = subvector( x, 0UL, 10UL );

   // Setting the first ten elements of x to the 2nd row of matrix A
   sv = row( A, 2UL );

   // Setting the second ten elements of x to y
   subvector( x, 10UL, 10UL ) = y;

   // Setting the 3rd row of A to a subvector of x
   row( A, 3UL ) = subvector( x, 3UL, 10UL );

   // Setting x to a subvector of the result of the addition between y and the 1st row of A
   x = subvector( y + row( A, 1UL ), 2UL, 5UL );
   \endcode

// \warning It is the programmer's responsibility to ensure the subvector does not outlive the
// viewed vector:

   \code
   // Creating a subvector on a temporary vector; results in a dangling reference!
   auto sv = subvector<1UL,3UL>( DynamicVector<int>{ 1, 2, 3, 4, 5 } );
   \endcode

// \n \section views_subvectors_element_access Element Access
// <hr>
//
// The elements of a subvector can be directly accessed via the subscript operator:

   \code
   blaze::DynamicVector<double,blaze::rowVector> v;
   // ... Resizing and initialization

   // Creating an 8-dimensional subvector, starting from index 4
   auto sv = subvector( v, 4UL, 8UL );

   // Setting the 1st element of the subvector, which corresponds to
   // the element at index 5 in vector v
   sv[1] = 2.0;
   \endcode

// The numbering of the subvector elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the specified size of the subvector. Alternatively, the elements of a subvector can
// be traversed via iterators. Just as with vectors, in case of non-const subvectors, \c begin()
// and \c end() return an iterator, which allows to manipulate the elements, in case of constant
// subvectors an iterator to immutable elements is returned:

   \code
   blaze::DynamicVector<int,blaze::rowVector> v( 256UL );
   // ... Resizing and initialization

   // Creating a reference to a specific subvector of vector v
   auto sv = subvector( v, 16UL, 64UL );

   // Traversing the elements via iterators to non-const elements
   for( auto it=sv.begin(); it!=sv.end(); ++it ) {
      *it = ...;  // OK: Write access to the dense subvector value.
      ... = *it;  // OK: Read access to the dense subvector value.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=sv.cbegin(); it!=sv.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense subvector value.
   }
   \endcode

   \code
   blaze::CompressedVector<int,blaze::rowVector> v( 256UL );
   // ... Resizing and initialization

   // Creating a reference to a specific subvector of vector v
   auto sv = subvector( v, 16UL, 64UL );

   // Traversing the elements via iterators to non-const elements
   for( auto it=sv.begin(); it!=sv.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=sv.cbegin(); it!=sv.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_subvectors_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse subvector can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   blaze::CompressedVector<double,blaze::rowVector> v( 256UL );  // Non-initialized vector of size 256

   auto sv = subvector( v, 10UL, 60UL );  // View on the range [10..69] of v

   // The subscript operator provides access to all possible elements of the sparse subvector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse subvector, the element is inserted into the
   // subvector.
   sv[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the subvector it is inserted into the subvector, if it is already contained
   // in the subvector its value is modified.
   sv.set( 45UL, -1.2 );

   // An alternative for inserting elements into the subvector is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the subvector.
   sv.insert( 50UL, 3.7 );

   // Just as in case of vectors, elements can also be inserted via the append() function. In
   // case of subvectors, append() also requires that the appended element's index is strictly
   // larger than the currently largest non-zero index of the subvector and that the subvector's
   // capacity is large enough to hold the new element. Note however that due to the nature of
   // a subvector, which may be an alias to the middle of a sparse vector, the append() function
   // does not work as efficiently for a subvector as it does for a vector.
   sv.reserve( 10UL );
   sv.append( 51UL, -2.1 );
   \endcode

// \n \section views_subvectors_common_operations Common Operations
// <hr>
//
// A subvector view can be used like any other dense or sparse vector. This means that with
// only a few exceptions all \ref vector_operations and \ref arithmetic_operations can be used.
// For instance, the current number of elements can be obtained via the \c size() function, the
// current capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since subvectors are references to a specific range of a
// vector, several operations are not possible, such as resizing and swapping. The following
// example shows this by means of a dense subvector view:

   \code
   blaze::DynamicVector<int,blaze::rowVector> v( 42UL );
   // ... Resizing and initialization

   // Creating a view on the range [5..15] of vector v
   auto sv = subvector( v, 5UL, 10UL );

   sv.size();          // Returns the number of elements in the subvector
   sv.capacity();      // Returns the capacity of the subvector
   sv.nonZeros();      // Returns the number of non-zero elements contained in the subvector

   sv.resize( 84UL );  // Compilation error: Cannot resize a subvector of a vector

   auto sv2 = subvector( v, 15UL, 10UL );
   swap( sv, sv2 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_subvectors_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse subvectors can be used in all arithmetic operations that any other dense
// or sparse vector can be used in. The following example gives an impression of the use of dense
// subvectors within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse subvectors with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::rowVector> d1, d2, d3;
   blaze::CompressedVector<double,blaze::rowVector> s1, s2;

   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> A;

   auto sv( subvector( d1, 0UL, 10UL ) );  // View on the range [0..9] of vector d1

   sv = d2;                           // Dense vector initialization of the range [0..9]
   subvector( d1, 10UL, 10UL ) = s1;  // Sparse vector initialization of the range [10..19]

   d3 = sv + d2;                           // Dense vector/dense vector addition
   s2 = s1 + subvector( d1, 10UL, 10UL );  // Sparse vector/dense vector addition
   d2 = sv * subvector( d1, 20UL, 10UL );  // Component-wise vector multiplication

   subvector( d1, 3UL, 4UL ) *= 2.0;      // In-place scaling of the range [3..6]
   d2 = subvector( d1, 7UL, 3UL ) * 2.0;  // Scaling of the range [7..9]
   d2 = 2.0 * subvector( d1, 7UL, 3UL );  // Scaling of the range [7..9]

   subvector( d1, 0UL , 10UL ) += d2;  // Addition assignment
   subvector( d1, 10UL, 10UL ) -= s2;  // Subtraction assignment
   subvector( d1, 20UL, 10UL ) *= sv;  // Multiplication assignment

   double scalar = subvector( d1, 5UL, 10UL ) * trans( s1 );  // Scalar/dot/inner product between two vectors

   A = trans( s1 ) * subvector( d1, 4UL, 16UL );  // Outer product between two vectors
   \endcode

// \n \section views_aligned_subvectors Aligned Subvectors
// <hr>
//
// Usually subvectors can be defined anywhere within a vector. They may start at any position and
// may have an arbitrary size (only restricted by the size of the underlying vector). However, in
// contrast to vectors themselves, which are always properly aligned in memory and therefore can
// provide maximum performance, this means that subvectors in general have to be considered to be
// unaligned. This can be made explicit by the \c blaze::unaligned flag:

   \code
   using blaze::unaligned;

   blaze::DynamicVector<double,blaze::rowVector> x;
   // ... Resizing and initialization

   // Identical creations of an unaligned subvector in the range [8..23]
   auto sv1 = subvector           ( x, 8UL, 16UL );
   auto sv2 = subvector<unaligned>( x, 8UL, 16UL );
   auto sv3 = subvector<8UL,16UL>          ( x );
   auto sv4 = subvector<unaligned,8UL,16UL>( x );
   \endcode

// All of these calls to the \c subvector() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned subvector. Whereas this may provide
// full flexibility in the creation of subvectors, this might result in performance disadvantages
// in comparison to vector primitives (even in case the specified subvector could be aligned).
// Whereas vector primitives are guaranteed to be properly aligned and therefore provide maximum
// performance in all operations, a general view on a vector might not be properly aligned. This
// may cause a performance penalty on some platforms and/or for some operations.
//
// However, it is also possible to create aligned subvectors. Aligned subvectors are identical to
// unaligned subvectors in all aspects, except that they may pose additional alignment restrictions
// and therefore have less flexibility during creation, but don't suffer from performance penalties
// and provide the same performance as the underlying vector. Aligned subvectors are created by
// explicitly specifying the \c blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned subvector in the range [8..23]
   auto sv1 = subvector<aligned>( x, 8UL, 16UL );
   auto sv2 = subvector<aligned,8UL,16UL>( x );
   \endcode

// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   using blaze::aligned;

   blaze::DynamicVector<double,blaze::columnVector> d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning, i.e. the first element is aligned
   auto dsv1 = subvector<aligned>( d, 0UL, 13UL );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   auto dsv2 = subvector<aligned>( d, 4UL, 7UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   auto dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   auto dsv4 = subvector<aligned>( d, 5UL, 8UL );
   \endcode

// Note that the discussed alignment restrictions are only valid for aligned dense subvectors.
// In contrast, aligned sparse subvectors at this time don't pose any additional restrictions.
// Therefore aligned and unaligned sparse subvectors are truly fully identical. Still, in case
// the \c blaze::aligned flag is specified during setup, an aligned subvector is created:

   \code
   using blaze::aligned;

   blaze::CompressedVector<double,blaze::rowVector> x;
   // ... Resizing and initialization

   // Creating an aligned subvector in the range [8..23]
   auto sv1 = subvector<aligned>( x, 8UL, 16UL );
   auto sv2 = subvector<aligned,8UL,16UL>( x );
   \endcode

// \n Previous: \ref views &nbsp; &nbsp; Next: \ref views_element_selections
*/
//*************************************************************************************************


//**Element Selections*****************************************************************************
/*!\page views_element_selections Element Selections
//
// \tableofcontents
//
//
// Element selections provide views on arbitrary compositions of elements of dense and sparse
// vectors. These views act as a reference to the selected elements and represent them as another
// dense or sparse vector. This reference is valid and can be used in every way any other dense
// or sparse vector can be used as long as the vector containing the elements is not resized or
// entirely destroyed. The element selection also acts as an alias to the vector elements in the
// specified range: Changes made to the elements (e.g. modifying values, inserting or erasing
// elements) are immediately visible in the vector and changes made via the vector are immediately
// visible in the elements.
//
//
// \n \section views_element_selections_setup Setup of Element Selections
//
// An element selection can be created very conveniently via the \c elements() function. It can
// be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Elements.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The indices of the elements to be selected can be specified either at compile time or at runtime
// (by means of an initializer list, array or vector):

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   // ... Resizing and initialization

   // Selecting the elements 4, 6, 8, and 10 (compile time arguments)
   auto e1 = elements<4UL,6UL,8UL,10UL>( x );

   // Selecting the elements 3, 2, and 1 (runtime arguments via an initializer list)
   const std::initializer_list<size_t> list{ 3UL, 2UL, 1UL };
   auto e2 = elements( x, { 3UL, 2UL, 1UL } );
   auto e3 = elements( x, list );

   // Selecting the elements 1, 2, 3, 3, 2, and 1 (runtime arguments via a std::array)
   const std::array<size_t> array{ 1UL, 2UL, 3UL, 3UL, 2UL, 1UL };
   auto e4 = elements( x, array );
   auto e5 = elements( x, array.data(), array.size() );

   // Selecting the element 4 fives times (runtime arguments via a std::vector)
   const std::vector<size_t> vector{ 4UL, 4UL, 4UL, 4UL, 4UL };
   auto e6 = elements( x, vector );
   auto e7 = elements( x, vector.data(), vector.size() );
   \endcode

// Note that it is possible to alias the elements of the underlying vector in any order. Also note
// that it is possible to use the same index multiple times.
//
// Alternatively it is possible to pass a callable such as a lambda or functor that produces the
// indices:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };

   // Selecting all even elements of the vector, i.e. selecting (0,2,4,6,8)
   auto e1 = elements( x, []( size_t i ){ return i*2UL; }, 5UL );

   // Selecting all odd elements of the vector, i.e. selecting (1,3,5,7)
   auto e2 = elements( x, []( size_t i ){ return i*2UL+1UL; }, 4UL );

   // Reversing the elements of the vector, i.e. selecting (8,7,6,5,4,3,2,1,0)
   auto e3 = elements( x, [max=v.size()-1UL]( size_t i ){ return max-i; }, 9UL );
   \endcode

// The \c elements() function returns an expression representing the view on the selected elements.
// The type of this expression depends on the given arguments, primarily the type of the vector and
// the compile time arguments. If the type is required, it can be determined via the \c decltype
// specifier:

   \code
   using VectorType = blaze::DynamicVector<int>;
   using ElementsType = decltype( blaze::elements<4UL,12UL>( std::declval<VectorType>() ) );
   \endcode

// The resulting view can be treated as any other dense or sparse vector, i.e. it can be assigned
// to, it can be copied from, and it can be used in arithmetic operations. An element selection
// created from a row vector can be used as any other row vector, an element selection created
// from a column vector can be used as any other column vector. The view can also be used on both
// sides of an assignment: It can either be used as an alias to grant write access to specific
// elements of a vector primitive on the left-hand side of an assignment or to grant read-access
// to specific elements of a vector primitive or expression on the right-hand side of an assignment.
// The following example demonstrates this in detail:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   blaze::CompressedVector<double,blaze::rowVector> y;
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Selecting the elements 1, 3, 5, and 7
   auto e = elements( x, { 1UL, 3UL, 5UL, 7UL } );

   // Setting the elements 1, 3, 5, and 7 of x to the 2nd row of matrix A
   e = row( A, 2UL );

   // Setting the elements 2, 4, 6, and 8 of x to y
   elements( x, { 2UL, 4UL, 6UL, 8UL } ) = y;

   // Setting the 3rd row of A to the elements 5, 4, 3, and 2 of x
   row( A, 3UL ) = elements( x, { 5UL, 4UL, 3UL, 2UL } );

   // Rotating the result of the addition between y and the 1st row of A
   x = elements( y + row( A, 1UL ), { 2UL, 3UL, 0UL, 1UL } )
   \endcode

// Please note that using an element selection, which refers to an index multiple times, on the
// left-hand side of an assignment leads to undefined behavior:

   \code
   blaze::DynamicVector<int,blaze::rowVector> a{ 1, 2, 3 };
   blaze::DynamicVector<int,blaze::rowVector> b{ 1, 2, 3, 4 };

   auto e = elements( a, { 1, 1, 1, 1 } );  // Selecting the element 1 four times
   e = b;  // Undefined behavior
   \endcode

// In this example both vectors have the same size, which results in a correct vector assignment,
// but the final value of the element at index 1 is unspecified.
//
// \warning It is the programmer's responsibility to ensure the element selection does not outlive
// the viewed vector:

   \code
   // Creating an element selection on a temporary vector; results in a dangling reference!
   auto e = elements<1UL,3UL>( DynamicVector<int>{ 1, 2, 3, 4, 5 } );
   \endcode

// \n \section views_element_selections_element_access Element Access
//
// The elements of an element selection can be directly accessed via the subscript operator:

   \code
   blaze::DynamicVector<double,blaze::rowVector> v;
   // ... Resizing and initialization

   // Selecting the elements 2, 4, 6, and 8
   auto e = elements( v, { 2UL, 4UL, 6UL, 8UL } );

   // Setting the 1st element of the element selection, which corresponds to
   // the element at index 4 in vector v
   e[1] = 2.0;
   \endcode

// The numbering of the selected elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of selected elements. Alternatively, the elements of an element selection
// can be traversed via iterators. Just as with vectors, in case of non-const element selections,
// \c begin() and \c end() return an iterator, which allows to manipulate the elements, in case of
// constant element selections an iterator to immutable elements is returned:

   \code
   blaze::DynamicVector<int,blaze::rowVector> v( 256UL );
   // ... Resizing and initialization

   // Creating an element selection including specific elements of dense vector v
   auto e = elements( v, { 0UL, 3UL, 6UL, 9UL, 12UL } );

   // Traversing the elements via iterators to non-const elements
   for( auto it=e.begin(); it!=e.end(); ++it ) {
      *it = ...;  // OK: Write access to the dense vector value.
      ... = *it;  // OK: Read access to the dense vector value.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=e.cbegin(); it!=e.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense vector value.
   }
   \endcode

   \code
   blaze::CompressedVector<int,blaze::rowVector> v( 256UL );
   // ... Resizing and initialization

   // Creating an element selection including specific elements of sparse vector v
   auto e = elements( v, { 0UL, 3UL, 6UL, 9UL, 12UL } );

   // Traversing the elements via iterators to non-const elements
   for( auto it=e.begin(); it!=e.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=e.cbegin(); it!=e.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_element_selections_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse element selection can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   blaze::CompressedVector<double,blaze::rowVector> v( 256UL );  // Non-initialized vector of size 256

   std::vector<size_t> indices;
   // ... Selecting indices of the sparse vector

   auto e = elements( v, indices );

   // The subscript operator provides access to the selected elements of the sparse vector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse vector, the element is inserted.
   e[42] = 2.0;

   // The second operation for inserting elements via the element selection is the set() function.
   // In case the element is not contained in the vector it is inserted into the vector, if it is
   // already contained in the vector its value is modified.
   e.set( 45UL, -1.2 );

   // An alternative for inserting elements into the vector is the insert() function. However, it
   // inserts the element only in case the element is not already contained in the vector.
   e.insert( 50UL, 3.7 );

   // Just as in case of vectors, elements can also be inserted via the append() function. In case
   // of element selections, append() also requires that the appended element's index is strictly
   // larger than the currently largest non-zero index of the selection and that the selections's
   // capacity is large enough to hold the new element. Note however that due to the nature of an
   // element selection, which is an alias to arbitrary elements of a sparse vector, the append()
   // function does not work as efficiently for an element selection as it does for a vector.
   e.reserve( 10UL );
   e.append( 51UL, -2.1 );
   \endcode

// \n \section views_element_selections_common_operations Common Operations
//
// An element selection can be used like any other dense or sparse vector. For instance, the
// number of selected elements can be obtained via the \c size() function, the current capacity
// via the \c capacity() function, and the number of non-zero elements via the \c nonZeros()
// function. However, since element selections are references to a specific range of a vector,
// several operations are not possible, such as resizing and swapping. The following example
// shows this by means of an element selection on a dense vector:

   \code
   blaze::DynamicVector<int,blaze::rowVector> v( 42UL );
   // ... Resizing and initialization

   // Selecting the elements 5 and 10
   auto e = elements( v, { 5UL, 10UL } );

   e.size();          // Returns the number of elements in the element selection
   e.capacity();      // Returns the capacity of the element selection
   e.nonZeros();      // Returns the number of non-zero elements contained in the element selection

   e.resize( 84UL );  // Compilation error: Cannot resize an element selection

   auto e2 = elements( v, { 15UL, 10UL } );
   swap( e, e2 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_element_selections_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse element selections can be used in all arithmetic operations that any other
// dense or sparse vector can be used in. The following example gives an impression of the use of
// dense element selections within arithmetic operations. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// element selections with fitting element types:

   \code
   blaze::DynamicVector<double,blaze::rowVector> d1, d2, d3;
   blaze::CompressedVector<double,blaze::rowVector> s1, s2;

   // ... Resizing and initialization

   blaze::DynamicMatrix<double,blaze::rowMajor> A;

   std::initializer_list<size_t> indices1{ 0UL, 3UL, 6UL,  9UL, 12UL, 15UL, 18UL, 21UL };
   std::initializer_list<size_t> indices2{ 1UL, 4UL, 7UL, 10UL, 13UL, 16UL, 19UL, 22UL };
   std::initializer_list<size_t> indices3{ 2UL, 5UL, 8UL, 11UL, 14UL, 17UL, 20UL, 23UL };

   auto e( elements( d1, indices1 ) );  // Selecting the every third element of d1 in the range [0..21]

   e = d2;                         // Dense vector assignment to the selected elements
   elements( d1, indices2 ) = s1;  // Sparse vector assignment to the selected elements

   d3 = e + d2;                         // Dense vector/dense vector addition
   s2 = s1 + elements( d1, indices2 );  // Sparse vector/dense vector addition
   d2 = e * elements( d1, indices3 );   // Component-wise vector multiplication

   elements( d1, indices2 ) *= 2.0;      // In-place scaling of the second selection of elements
   d2 = elements( d1, indices3 ) * 2.0;  // Scaling of the elements in the third selection of elements
   d2 = 2.0 * elements( d1, indices3 );  // Scaling of the elements in the third selection of elements

   elements( d1, indices1 ) += d2;  // Addition assignment
   elements( d1, indices2 ) -= s2;  // Subtraction assignment
   elements( d1, indices3 ) *= e;   // Multiplication assignment

   double scalar = elements( d1, indices2 ) * trans( s1 );  // Scalar/dot/inner product between two vectors

   A = trans( s1 ) * elements( d1, { 3UL, 6UL } );  // Outer product between two vectors
   \endcode

// \n Previous: \ref views_subvectors &nbsp; &nbsp; Next: \ref views_submatrices
*/
//*************************************************************************************************


//**Submatrices************************************************************************************
/*!\page views_submatrices Submatrices
//
// \tableofcontents
//
//
// Submatrices provide views on a specific part of a dense or sparse matrix just as subvectors
// provide views on specific parts of vectors. As such, submatrices act as a reference to a
// specific block within a matrix. This reference is valid and can be used in evary way any
// other dense or sparse matrix can be used as long as the matrix containing the submatrix is
// not resized or entirely destroyed. The submatrix also acts as an alias to the matrix elements
// in the specified block: Changes made to the elements (e.g. modifying values, inserting or
// erasing elements) are immediately visible in the matrix and changes made via the matrix are
// immediately visible in the submatrix.
//
//
// \n \section views_submatrices_setup Setup of Submatrices
// <hr>
//
// A view on a dense or sparse submatrix can be created very conveniently via the \c submatrix()
// function. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Submatrix.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The first and second parameter specify the row and column of the first element of the submatrix.
// The third and fourth parameter specify the number of rows and columns, respectively. The four
// parameters can be specified either at compile time or at runtime:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 4x8, starting in row 3 and column 0 (compile time arguments)
   auto sm1 = submatrix<3UL,0UL,4UL,8UL>( A );

   // Creating a dense submatrix of size 8x16, starting in row 0 and column 4 (runtime arguments)
   auto sm2 = submatrix( A, 0UL, 4UL, 8UL, 16UL );
   \endcode

// The \c submatrix() function returns an expression representing the submatrix view. The type of
// this expression depends on the given submatrix arguments, primarily the type of the matrix and
// the compile time arguments. If the type is required, it can be determined via the \c decltype
// specifier:

   \code
   using MatrixType = blaze::DynamicMatrix<int>;
   using SubmatrixType = decltype( blaze::submatrix<3UL,0UL,4UL,8UL>( std::declval<MatrixType>() ) );
   \endcode

// The resulting view can be treated as any other dense or sparse matrix, i.e. it can be assigned
// to, it can be copied from, and it can be used in arithmetic operations. A submatrix created from
// a row-major matrix will itself be a row-major matrix, a submatrix created from a column-major
// matrix will be a column-major matrix. The view can also be used on both sides of an assignment:
// The submatrix can either be used as an alias to grant write access to a specific submatrix
// of a matrix primitive on the left-hand side of an assignment or to grant read-access to
// a specific submatrix of a matrix primitive or expression on the right-hand side of an
// assignment. The following example demonstrates this in detail:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A, B;
   blaze::CompressedMatrix<double,blaze::rowMajor> C;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 2
   auto sm = submatrix( A, 0UL, 2UL, 8UL, 4UL );

   // Setting the submatrix of A to a 8x4 submatrix of B
   sm = submatrix( B, 0UL, 0UL, 8UL, 4UL );

   // Copying the sparse matrix C into another 8x4 submatrix of A
   submatrix( A, 8UL, 2UL, 8UL, 4UL ) = C;

   // Assigning part of the result of a matrix addition to the first submatrix
   sm = submatrix( B + C, 0UL, 0UL, 8UL, 4UL );
   \endcode

// \warning It is the programmer's responsibility to ensure the submatrix does not outlive the
// viewed matrix:

   \code
   // Creating a submatrix on a temporary matrix; results in a dangling reference!
   auto sm = submatrix<1UL,0UL,2UL,3UL>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } } );
   \endcode

// \n \section views_submatrices_element_access Element Access
// <hr>
//
// The elements of a submatrix can be directly accessed with the function call operator:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a 8x8 submatrix, starting from position (4,4)
   auto sm = submatrix( A, 4UL, 4UL, 8UL, 8UL );

   // Setting the element (0,0) of the submatrix, which corresponds to
   // the element at position (4,4) in matrix A
   sm(0,0) = 2.0;
   \endcode

// Alternatively, the elements of a submatrix can be traversed via (const) iterators. Just as
// with matrices, in case of non-const submatrices, \c begin() and \c end() return an iterator,
// which allows to manipuate the elements, in case of constant submatrices an iterator to
// immutable elements is returned:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a specific submatrix of matrix A
   auto sm = submatrix( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( auto it=sm.begin(0); it!=sm.end(0); ++it ) {
      *it = ...;  // OK: Write access to the dense submatrix value.
      ... = *it;  // OK: Read access to the dense submatrix value.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( auto it=sm.cbegin(1); it!=sm.cend(1); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense submatrix value.
   }
   \endcode

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a specific submatrix of matrix A
   auto sm = submatrix( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( auto it=sm.begin(0); it!=sm.end(0); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( auto it=sm.cbegin(1); it!=sm.cend(1); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_submatrices_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse submatrix can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 256UL, 512UL );  // Non-initialized matrix of size 256x512

   auto sm = submatrix( A, 10UL, 10UL, 16UL, 16UL );  // View on a 16x16 submatrix of A

   // The function call operator provides access to all possible elements of the sparse submatrix,
   // including the zero elements. In case the function call operator is used to access an element
   // that is currently not stored in the sparse submatrix, the element is inserted into the
   // submatrix.
   sm(2,4) = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the submatrix it is inserted into the submatrix, if it is already contained
   // in the submatrix its value is modified.
   sm.set( 2UL, 5UL, -1.2 );

   // An alternative for inserting elements into the submatrix is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the submatrix.
   sm.insert( 2UL, 6UL, 3.7 );

   // Just as in the case of sparse matrices, elements can also be inserted via the append()
   // function. In case of submatrices, append() also requires that the appended element's
   // index is strictly larger than the currently largest non-zero index in the according row
   // or column of the submatrix and that the according row's or column's capacity is large
   // enough to hold the new element. Note however that due to the nature of a submatrix, which
   // may be an alias to the middle of a sparse matrix, the append() function does not work as
   // efficiently for a submatrix as it does for a matrix.
   sm.reserve( 2UL, 10UL );
   sm.append( 2UL, 10UL, -2.1 );
   \endcode

// \n \section views_submatrices_common_operations Common Operations
// <hr>
//
// A submatrix view can be used like any other dense or sparse matrix. This means that with only
// a few exceptions all \ref matrix_operations and \ref arithmetic_operations can be used. For
// instance, the current size of the matrix, i.e. the number of rows or columns can be obtained
// via the \c rows() and \c columns() functions, the current total capacity via the \c capacity()
// function, and the number of non-zero elements via the \c nonZeros() function. However, since
// submatrices are views on a specific submatrix of a matrix, several operations are not possible,
// such as resizing and swapping:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a view on the a 8x12 submatrix of matrix A
   auto sm = submatrix( A, 0UL, 0UL, 8UL, 12UL );

   sm.rows();      // Returns the number of rows of the submatrix
   sm.columns();   // Returns the number of columns of the submatrix
   sm.capacity();  // Returns the capacity of the submatrix
   sm.nonZeros();  // Returns the number of non-zero elements contained in the submatrix

   sm.resize( 10UL, 8UL );  // Compilation error: Cannot resize a submatrix of a matrix

   auto sm2 = submatrix( A, 8UL, 0UL, 12UL, 8UL );
   swap( sm, sm2 );  // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_submatrices_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse submatrices can be used in all arithmetic operations that any other dense
// or sparse matrix can be used in. The following example gives an impression of the use of dense
// submatrices within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse matrices with
// fitting element types:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D1, D2, D3;
   blaze::CompressedMatrix<double,blaze::rowMajor> S1, S2;

   blaze::CompressedVector<double,blaze::columnVector> a, b;

   // ... Resizing and initialization

   auto sm = submatrix( D1, 0UL, 0UL, 8UL, 8UL );  // View on the 8x8 submatrix of matrix D1
                                                   // starting from row 0 and column 0

   submatrix( D1, 0UL, 8UL, 8UL, 8UL ) = D2;  // Dense matrix initialization of the 8x8 submatrix
                                              // starting in row 0 and column 8
   sm = S1;                                   // Sparse matrix initialization of the second 8x8 submatrix

   D3 = sm + D2;                                   // Dense matrix/dense matrix addition
   S2 = S1 - submatrix( D1, 8UL, 0UL, 8UL, 8UL );  // Sparse matrix/dense matrix subtraction
   D2 = sm * submatrix( D1, 8UL, 8UL, 8UL, 8UL );  // Dense matrix/dense matrix multiplication

   submatrix( D1, 8UL, 0UL, 8UL, 8UL ) *= 2.0;      // In-place scaling of a submatrix of D1
   D2 = submatrix( D1, 8UL, 8UL, 8UL, 8UL ) * 2.0;  // Scaling of the a submatrix of D1
   D2 = 2.0 * sm;                                   // Scaling of the a submatrix of D1

   submatrix( D1, 0UL, 8UL, 8UL, 8UL ) += D2;  // Addition assignment
   submatrix( D1, 8UL, 0UL, 8UL, 8UL ) -= S1;  // Subtraction assignment
   submatrix( D1, 8UL, 8UL, 8UL, 8UL ) *= sm;  // Multiplication assignment

   a = submatrix( D1, 4UL, 4UL, 8UL, 8UL ) * b;  // Dense matrix/sparse vector multiplication
   \endcode

// \n \section views_aligned_submatrices Aligned Submatrices
// <hr>
//
// Usually submatrices can be defined anywhere within a matrix. They may start at any position and
// may have an arbitrary extension (only restricted by the extension of the underlying matrix).
// However, in contrast to matrices themselves, which are always properly aligned in memory and
// therefore can provide maximum performance, this means that submatrices in general have to be
// considered to be unaligned. This can be made explicit by the \c blaze::unaligned flag:

   \code
   using blaze::unaligned;

   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Identical creations of an unaligned submatrix of size 8x8, starting in row 0 and column 0
   auto sm1 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   auto sm2 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   auto sm3 = submatrix<0UL,0UL,8UL,8UL>          ( A );
   auto sm4 = submatrix<unaligned,0UL,0UL,8UL,8UL>( A );
   \endcode

// All of these calls to the \c submatrix() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned submatrix. Whereas this may provide
// full flexibility in the creation of submatrices, this might result in performance disadvantages
// in comparison to matrix primitives (even in case the specified submatrix could be aligned).
// Whereas matrix primitives are guaranteed to be properly aligned and therefore provide maximum
// performance in all operations, a general view on a matrix might not be properly aligned. This
// may cause a performance penalty on some platforms and/or for some operations.
//
// However, it is also possible to create aligned submatrices. Aligned submatrices are identical to
// unaligned submatrices in all aspects, except that they may pose additional alignment restrictions
// and therefore have less flexibility during creation, but don't suffer from performance penalties
// and provide the same performance as the underlying matrix. Aligned submatrices are created by
// explicitly specifying the \c blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned submatrix of size 8x8, starting in row 0 and column 0
   auto sv1 = submatrix<aligned>( A, 0UL, 0UL, 8UL, 8UL );
   auto sv2 = submatrix<aligned,0UL,0UL,8UL,8UL>( A );
   \endcode

// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of each row/column of the submatrix must be aligned. The following source code
// gives some examples for a double precision row-major dynamic matrix, assuming that padding is
// enabled and that AVX is available, which packs 4 \c double values into a SIMD vector:

   \code
   using blaze::aligned;

   blaze::DynamicMatrix<double,blaze::rowMajor> D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0), i.e. the first element of each row is aligned (due to padding)
   auto dsm1 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 11UL );

   // OK: First column is a multiple of 4, i.e. the first element of each row is aligned (due to padding)
   auto dsm2 = submatrix<aligned>( D, 3UL, 12UL, 8UL, 16UL );

   // OK: First column is a multiple of 4 and the submatrix includes the last row and column
   auto dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First column is not a multiple of 4, i.e. the first element is not aligned
   auto dsm4 = submatrix<aligned>( D, 2UL, 3UL, 12UL, 12UL );
   \endcode

// Note that the discussed alignment restrictions are only valid for aligned dense submatrices.
// In contrast, aligned sparse submatrices at this time don't pose any additional restrictions.
// Therefore aligned and unaligned sparse submatrices are truly fully identical. Still, in case
// the \c blaze::aligned flag is specified during setup, an aligned submatrix is created:

   \code
   using blaze::aligned;

   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating an aligned submatrix of size 8x8, starting in row 0 and column 0
   auto sv = submatrix<aligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// \n \section views_submatrices_on_symmetric_matrices Submatrices on Symmetric Matrices
//
// Submatrices can also be created on symmetric matrices (see the \c SymmetricMatrix class template):

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Setup of a 16x16 symmetric matrix
   SymmetricMatrix< DynamicMatrix<int> > A( 16UL );

   // Creating a dense submatrix of size 8x12, starting in row 2 and column 4
   auto sm = submatrix( A, 2UL, 4UL, 8UL, 12UL );
   \endcode

// It is important to note, however, that (compound) assignments to such submatrices have a
// special restriction: The symmetry of the underlying symmetric matrix must not be broken!
// Since the modification of element \f$ a_{ij} \f$ of a symmetric matrix also modifies the
// element \f$ a_{ji} \f$, the matrix to be assigned must be structured such that the symmetry
// of the symmetric matrix is preserved. Otherwise a \c std::invalid_argument exception is
// thrown:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Setup of two default 4x4 symmetric matrices
   SymmetricMatrix< DynamicMatrix<int> > A1( 4 ), A2( 4 );

   // Setup of the 3x2 dynamic matrix
   //
   //       ( 1 2 )
   //   B = ( 3 4 )
   //       ( 5 6 )
   //
   DynamicMatrix<int> B{ { 1, 2 }, { 3, 4 }, { 5, 6 } };

   // OK: Assigning B to a submatrix of A1 such that the symmetry can be preserved
   //
   //        ( 0 0 1 2 )
   //   A1 = ( 0 0 3 4 )
   //        ( 1 3 5 6 )
   //        ( 2 4 6 0 )
   //
   submatrix( A1, 0UL, 2UL, 3UL, 2UL ) = B;  // OK

   // Error: Assigning B to a submatrix of A2 such that the symmetry cannot be preserved!
   //   The elements marked with X cannot be assigned unambiguously!
   //
   //        ( 0 1 2 0 )
   //   A2 = ( 1 3 X 0 )
   //        ( 2 X 6 0 )
   //        ( 0 0 0 0 )
   //
   submatrix( A2, 0UL, 1UL, 3UL, 2UL ) = B;  // Assignment throws an exception!
   \endcode

// \n Previous: \ref views_element_selections &nbsp; &nbsp; Next: \ref views_rows
*/
//*************************************************************************************************


//**Rows*******************************************************************************************
/*!\page views_rows Rows
//
// \tableofcontents
//
//
// Rows provide views on a specific row of a dense or sparse matrix. As such, rows act as a
// reference to a specific row. This reference is valid and can be used in every way any other
// row vector can be used as long as the matrix containing the row is not resized or entirely
// destroyed. The row also acts as an alias to the row elements: Changes made to the elements
// (e.g. modifying values, inserting or erasing elements) are immediately visible in the matrix
// and changes made via the matrix are immediately visible in the row.
//
//
// \n \section views_rows_setup Setup of Rows
// <hr>
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
// time arguments. If the type is required, it can be determined via the \c decltype specifier:

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

// \warning It is the programmer's responsibility to ensure the row does not outlive the viewed
// matrix:

   \code
   // Creating a row on a temporary matrix; results in a dangling reference!
   auto row1 = row<1UL>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } } );
   \endcode

// \n \section views_rows_element_access Element Access
// <hr>
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

// \n \section views_rows_element_insertion Element Insertion
// <hr>
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

// \n \section views_rows_common_operations Common Operations
// <hr>
//
// A row view can be used like any other row vector. This means that with only a few exceptions
// all \ref vector_operations and \ref arithmetic_operations can be used. For instance, the
// current number of elements can be obtained via the \c size() function, the current capacity
// via the \c capacity() function, and the number of non-zero elements via the \c nonZeros()
// function. However, since rows are references to specific rows of a matrix, several operations
// are not possible on views, such as resizing and swapping. The following example shows this by
// means of a dense row view:

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

// \n \section views_rows_arithmetic_operations Arithmetic Operations
// <hr>
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

// \n \section views_rows_non_fitting_storage_order Views on Matrices with Non-Fitting Storage Order
// <hr>
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

// Although \b Blaze performs the resulting vector/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix \c A would result in a more efficient evaluation.
//
// \n Previous: \ref views_submatrices &nbsp; &nbsp; Next: \ref views_row_selections
*/
//*************************************************************************************************


//**Row Selections*********************************************************************************
/*!\page views_row_selections Row Selections
//
// \tableofcontents
//
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
// \n \section views_row_selections_setup Setup of Row Selections
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
   auto rs2 = rows( A, []( size_t i ){ return i*2UL+1UL; }, 4UL );

   // Reversing the rows of the matrix, i.e. selecting the rows 8, 7, 6, 5, 4, 3, 2, 1, and 0
   auto rs3 = rows( A, [max=A.rows()-1UL]( size_t i ){ return max-i; }, 9UL );
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

// \warning It is the programmer's responsibility to ensure the row selection does not outlive the
// viewed matrix:

   \code
   // Creating a row selection on a temporary matrix; results in a dangling reference!
   auto rs = rows<2UL,0UL>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } } );
   \endcode

// \n \section views_row_selections_element_access Element Access
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

// \n \section views_row_selections_element_insertion Element Insertion
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

   // An alternative for inserting elements into the row selection is the insert() function.
   // However, it inserts the element only in case the element is not already contained in the
   // row selection.
   rs.insert( 2UL, 6UL, 3.7 );

   // Just as in the case of sparse matrices, elements can also be inserted via the append()
   // function. In case of row selections, append() also requires that the appended element's
   // index is strictly larger than the currently largest non-zero index in the according row
   // of the row selection and that the according row's capacity is large enough to hold the new
   // element. Note however that due to the nature of a row selection, which may be an alias to
   // an arbitrary collection of rows, the append() function does not work as efficiently for
   // a row selection as it does for a matrix.
   rs.reserve( 2UL, 10UL );
   rs.append( 2UL, 10UL, -2.1 );
   \endcode

// \n \section views_row_selections_common_operations Common Operations
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

// \n \section views_row_selections_arithmetic_operations Arithmetic Operations
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

// \n \section views_row_selections_on_column_major_matrix Row Selections on Column-Major Matrices
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

// Although \b Blaze performs the resulting matrix/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix \c A would result in a more efficient evaluation.
//
// \n Previous: \ref views_rows &nbsp; &nbsp; Next: \ref views_columns
*/
//*************************************************************************************************


//**Columns****************************************************************************************
/*!\page views_columns Columns
//
// \tableofcontents
//
//
// Just as rows provide a view on a specific row of a matrix, columns provide views on a specific
// column of a dense or sparse matrix. As such, columns act as a reference to a specific column.
// This reference is valid an can be used in every way any other column vector can be used as long
// as the matrix containing the column is not resized or entirely destroyed. Changes made to the
// elements (e.g. modifying values, inserting or erasing elements) are immediately visible in the
// matrix and changes made via the matrix are immediately visible in the column.
//
//
// \n \section views_colums_setup Setup of Columns
// <hr>
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
// compile time arguments. If the type is required, it can be determined via the \c decltype
// specifier:

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

// \warning It is the programmer's responsibility to ensure the column does not outlive the
// viewed matrix:

   \code
   // Creating a column on a temporary matrix; results in a dangling reference!
   auto col1 = column<1UL>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } } );
   \endcode

// \n \section views_columns_element_access Element Access
// <hr>
//
// The elements of a column can be directly accessed with the subscript operator.

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

// \n \section views_columns_element_insertion Element Insertion
// <hr>
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

// \n \section views_columns_common_operations Common Operations
// <hr>
//
// A column view can be used like any other column vector. This means that with only a few
// exceptions all \ref vector_operations and \ref arithmetic_operations can be used. For instance,
// the current number of elements can be obtained via the \c size() function, the current capacity
// via the \c capacity() function, and the number of non-zero elements via the \c nonZeros()
// function. However, since columns are references to specific columns of a matrix, several
// operations are not possible on views, such as resizing and swapping. The following example
// shows this by means of a dense column view:

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

// \n \section views_columns_arithmetic_operations Arithmetic Operations
// <hr>
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

// \n \section views_columns_non_fitting_storage_order Views on Matrices with Non-Fitting Storage Order
// <hr>
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

// Although \b Blaze performs the resulting matrix/vector multiplication as efficiently as possible
// using a column-major storage order for matrix \c B would result in a more efficient evaluation.
//
// \n Previous: \ref views_row_selections &nbsp; &nbsp; Next: \ref views_column_selections
*/
//*************************************************************************************************


//**Column Selections******************************************************************************
/*!\page views_column_selections Column Selections
//
// \tableofcontents
//
//
// Column selections provide views on arbitrary compositions of columns of dense and sparse
// matrices. These views act as a reference to the selected columns and represent them as another
// dense or sparse matrix. This reference is valid and can be used in every way any other dense
// or sparse matrix can be used as long as the matrix containing the columns is not resized or
// entirely destroyed. The column selection also acts as an alias to the matrix elements in the
// specified range: Changes made to the columns (e.g. modifying values, inserting or erasing
// elements) are immediately visible in the matrix and changes made via the matrix are immediately
// visible in the columns.
//
//
// \n \section views_column_selections_setup Setup of Column Selections
//
// A column selection can be created very conveniently via the \c columns() function. It can be
// included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Columns.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The indices of the columns to be selected can be specified either at compile time or at runtime
// (by means of an initializer list, array or vector):

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   // Selecting the columns 4, 6, 8, and 10 (compile time arguments)
   auto cs1 = columns<4UL,6UL,8UL,10UL>( A );

   // Selecting the columns 3, 2, and 1 (runtime arguments via an initializer list)
   const std::initializer_list<size_t> list{ 3UL, 2UL, 1UL };
   auto cs2 = columns( A, { 3UL, 2UL, 1UL } );
   auto cs3 = columns( A, list );

   // Selecting the columns 1, 2, 3, 3, 2, and 1 (runtime arguments via a std::array)
   const std::array<size_t> array{ 1UL, 2UL, 3UL, 3UL, 2UL, 1UL };
   auto cs4 = columns( A, array );
   auto cs5 = columns( A, array.data(), array.size() );

   // Selecting the column 4 fives times (runtime arguments via a std::vector)
   const std::vector<size_t> vector{ 4UL, 4UL, 4UL, 4UL, 4UL };
   auto cs6 = columns( A, vector );
   auto cs7 = columns( A, vector.data(), vector.size() );
   \endcode

// Note that it is possible to alias the columns of the underlying matrix in any order. Also note
// that it is possible to use the same index multiple times.
//
// Alternatively it is possible to pass a callable such as a lambda or functor that produces the
// indices:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A( 18UL, 9UL );

   // Selecting all even columns of the matrix, i.e. selecting the columns 0, 2, 4, 6, and 8
   auto cs1 = columns( A, []( size_t i ){ return i*2UL; }, 5UL );

   // Selecting all odd columns of the matrix, i.e. selecting the columns 1, 3, 5, and 7
   auto cs2 = columns( A, []( size_t i ){ return i*2UL+1UL; }, 4UL );

   // Reversing the columns of the matrix, i.e. selecting the columns 8, 7, 6, 5, 4, 3, 2, 1, and 0
   auto cs3 = columns( A, [max=A.columns()-1UL]( size_t i ){ return max-i; }, 9UL );
   \endcode

// The \c columns() function returns an expression representing the view on the selected columns.
// The type of this expression depends on the given arguments, primarily the type of the matrix
// and the compile time arguments. If the type is required, it can be determined via the \c decltype
// specifier:

   \code
   using MatrixType = blaze::DynamicMatrix<int>;
   using ColumnsType = decltype( blaze::columns<3UL,0UL,4UL,8UL>( std::declval<MatrixType>() ) );
   \endcode

// The resulting view can be treated as any other dense or sparse matrix, i.e. it can be assigned
// to, it can be copied from, and it can be used in arithmetic operations. Note, however, that a
// column selection will always be treated as a column-major matrix, regardless of the storage
// order of the matrix containing the columns. The view can also be used on both sides of an
// assignment: It can either be used as an alias to grant write access to specific columns of a
// matrix primitive on the left-hand side of an assignment or to grant read-access to specific
// columns of a matrix primitive or expression on the right-hand side of an assignment. The
// following example demonstrates this in detail:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A;
   blaze::DynamicMatrix<double,blaze::rowMajor> B;
   blaze::CompressedMatrix<double,blaze::columnMajor> C;
   // ... Resizing and initialization

   // Selecting the columns 1, 3, 5, and 7 of A
   auto cs = columns( A, { 1UL, 3UL, 5UL, 7UL } );

   // Setting columns 1, 3, 5, and 7 of A to column 4 of B
   cs = columns( B, { 4UL, 4UL, 4UL, 4UL } );

   // Setting the columns 2, 4, 6, and 8 of A to C
   columns( A, { 2UL, 4UL, 6UL, 8UL } ) = C;

   // Setting the first 4 columns of A to the columns 5, 4, 3, and 2 of C
   submatrix( A, 0UL, 0UL, A.rows(), 4UL ) = columns( C, { 5UL, 4UL, 3UL, 2UL } );

   // Rotating the result of the addition between columns 1, 3, 5, and 7 of A and C
   B = columns( cs + C, { 2UL, 3UL, 0UL, 1UL } );
   \endcode

// \warning It is the programmer's responsibility to ensure the column selection does not outlive
// the viewed matrix:

   \code
   // Creating a column selection on a temporary matrix; results in a dangling reference!
   auto cs = columns<2UL,0UL>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } } );
   \endcode

// \n \section views_column_selections_element_access Element Access
//
// The elements of a column selection can be directly accessed via the function call operator:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   // Creating a view on the first four columns of A in reverse order
   auto cs = columns( A, { 3UL, 2UL, 1UL, 0UL } );

   // Setting the element (0,0) of the column selection, which corresponds
   // to the element at position (0,3) in matrix A
   cs(0,0) = 2.0;
   \endcode

// Alternatively, the elements of a column selection can be traversed via (const) iterators.
// Just as with matrices, in case of non-const column selection, \c begin() and \c end() return
// an iterator, which allows to manipuate the elements, in case of constant column selection an
// iterator to immutable elements is returned:

   \code
   blaze::DynamicMatrix<int,blaze::columnMajor> A( 512UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to a selection of columns of matrix A
   auto cs = columns( A, { 16UL, 32UL, 64UL, 128UL } );

   // Traversing the elements of the 0th column via iterators to non-const elements
   for( auto it=cs.begin(0); it!=cs.end(0); ++it ) {
      *it = ...;  // OK: Write access to the dense value.
      ... = *it;  // OK: Read access to the dense value.
   }

   // Traversing the elements of the 1st column via iterators to const elements
   for( auto it=cs.cbegin(1); it!=cs.cend(1); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense value.
   }
   \endcode

   \code
   blaze::CompressedMatrix<int,blaze::columnMajor> A( 512UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to a selection of columns of matrix A
   auto cs = columns( A, { 16UL, 32UL, 64UL, 128UL } );

   // Traversing the elements of the 0th column via iterators to non-const elements
   for( auto it=cs.begin(0); it!=cs.end(0); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements of the 1st column via iterators to const elements
   for( auto it=cs.cbegin(1); it!=cs.cend(1); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_column_selections_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse column selection can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A( 512UL, 256UL );  // Non-initialized matrix of size 512x256

   auto cs = columns( A, { 10UL, 20UL, 30UL, 40UL } );  // View on the columns 10, 20, 30, and 40 of A

   // The function call operator provides access to all possible elements of the sparse column
   // selection, including the zero elements. In case the function call operator is used to
   // access an element that is currently not stored in the sparse column selection, the element
   // is inserted into the column selection.
   cs(2,4) = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the column selection it is inserted into the column selection, if it is
   // already contained in the column selection its value is modified.
   cs.set( 2UL, 5UL, -1.2 );

   // An alternative for inserting elements into the column selection is the insert() function.
   // However, it inserts the element only in case the element is not already contained in the
   // column selection.
   cs.insert( 2UL, 6UL, 3.7 );

   // Just as in the case of sparse matrices, elements can also be inserted via the append()
   // function. In case of column selections, append() also requires that the appended element's
   // index is strictly larger than the currently largest non-zero index in the according column
   // of the column selection and that the according column's capacity is large enough to hold the
   // new element. Note however that due to the nature of a column selection, which may be an alias
   // to an arbitrary collection of columns, the append() function does not work as efficiently
   // for a column selection as it does for a matrix.
   cs.reserve( 2UL, 10UL );
   cs.append( 2UL, 10UL, -2.1 );
   \endcode

// \n \section views_column_selections_common_operations Common Operations
//
// A view on specific columns of a matrix can be used like any other dense or sparse matrix. For
// instance, the current size of the matrix, i.e. the number of rows or columns can be obtained
// via the \c rows() and \c columns() functions, the current total capacity via the \c capacity()
// function, and the number of non-zero elements via the \c nonZeros() function. However, since
// column selections are views on specific columns of a matrix, several operations are not possible,
// such as resizing and swapping:

   \code
   blaze::DynamicMatrix<int,blaze::columnMajor> A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a view on the columns 8, 16, 24, and 32 of matrix A
   auto cs = columns( A, { 8UL, 16UL, 24UL, 32UL } );

   cs.rows();      // Returns the number of rows of the column selection
   cs.columns();   // Returns the number of columns of the column selection
   cs.capacity();  // Returns the capacity of the column selection
   cs.nonZeros();  // Returns the number of non-zero elements contained in the column selection

   cs.resize( 10UL, 8UL );  // Compilation error: Cannot resize a column selection

   auto cs2 = columns( A, 9UL, 17UL, 25UL, 33UL );
   swap( cs, cs2 );  // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_column_selections_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse column selections can be used in all arithmetic operations that any other
// dense or sparse matrix can be used in. The following example gives an impression of the use of
// dense column selctions within arithmetic operations. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and
// sparse matrices with fitting element types:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> D1, D2, D3;
   blaze::CompressedMatrix<double,blaze::columnMajor> S1, S2;

   blaze::CompressedVector<double,blaze::columnVector> a, b;

   // ... Resizing and initialization

   std::initializer_list<size_t> indices1{ 0UL, 3UL, 6UL,  9UL, 12UL, 15UL, 18UL, 21UL };
   std::initializer_list<size_t> indices2{ 1UL, 4UL, 7UL, 10UL, 13UL, 16UL, 19UL, 22UL };
   std::initializer_list<size_t> indices3{ 2UL, 5UL, 8UL, 11UL, 14UL, 17UL, 20UL, 23UL };

   auto cs = columns( D1, indices1 );  // Selecting the every third column of D1 in the range [0..21]

   cs = D2;                       // Dense matrix assignment to the selected columns
   columns( D1, indices2 ) = S1;  // Sparse matrix assignment to the selected columns

   D3 = cs + D2;                       // Dense matrix/dense matrix addition
   S2 = S1 - columns( D1, indices2 );  // Sparse matrix/dense matrix subtraction
   D2 = cs % columns( D1, indices3 );  // Dense matrix/dense matrix Schur product
   D2 = columns( D1, indices2 ) * D1;  // Dense matrix/dense matrix multiplication

   columns( D1, indices2 ) *= 2.0;      // In-place scaling of the second selection of columns
   D2 = columns( D1, indices3 ) * 2.0;  // Scaling of the elements in the third selection of columns
   D2 = 2.0 * columns( D1, indices3 );  // Scaling of the elements in the third selection of columns

   columns( D1, indices1 ) += D2;  // Addition assignment
   columns( D1, indices2 ) -= S1;  // Subtraction assignment
   columns( D1, indices3 ) %= cs;  // Schur product assignment

   a = columns( D1, indices1 ) * b;  // Dense matrix/sparse vector multiplication
   \endcode

// \n \section views_column_selections_on_row_major_matrix Column Selections on a Row-Major Matrix
//
// Especially noteworthy is that column selections can be created for both row-major and
// column-major matrices. Whereas the interface of a row-major matrix only allows to traverse a
// row directly and the interface of a column-major matrix only allows to traverse a column, via
// views it is possible to traverse a row of a column-major matrix or a column of a row-major
// matrix. For instance:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 1st and 3rd column of a column-major matrix A
   auto cs = columns( A, { 1UL, 3UL } );

   // Traversing column 0 of the selection, which corresponds to the 1st column of matrix A
   for( auto it=cs.begin( 0UL ); it!=cs.end( 0UL ); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a column selection on a matrix stored in a row-major fashion
// can result in a considerable performance decrease in comparison to a column selection on a
// matrix with column-major storage format. This is due to the non-contiguous storage of the
// matrix elements. Therefore care has to be taken in the choice of the most suitable storage
// order:

   \code
   // Setup of two row-major matrices
   blaze::DynamicMatrix<double,blaze::rowMajor> A( 128UL, 128UL );
   blaze::DynamicMatrix<double,blaze::rowMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th, 30th, and 45th column of the multiplication between A and B ...
   blaze::DynamicMatrix<double,blaze::columnMajor> x = columns( A * B, { 15UL, 30UL, 45UL } );

   // ... is essentially the same as the following computation, which multiplies
   // A with the 15th, 30th, and 45th column of the row-major matrix B.
   blaze::DynamicMatrix<double,blaze::columnMajor> x = A * column( B, { 15UL, 30UL, 45UL } );
   \endcode

// Although \b Blaze performs the resulting matrix/matrix multiplication as efficiently as possible
// using a column-major storage order for matrix \c A would result in a more efficient evaluation.
//
// \n Previous: \ref views_columns &nbsp; &nbsp; Next: \ref views_bands
*/
//*************************************************************************************************


//**Bands******************************************************************************************
/*!\page views_bands Bands
//
// \tableofcontents
//
//
// Bands provide views on a specific band of a dense or sparse matrix (e.g. the diagonal, the
// subdiagonal, ...). As such, bands act as a reference to a specific band. This reference
// is valid and can be used in every way any other vector can be used as long as the matrix
// containing the band is not resized or entirely destroyed. The band also acts as an alias to
// the band elements: Changes made to the elements (e.g. modifying values, inserting or erasing
// elements) are immediately visible in the matrix and changes made via the matrix are immediately
// visible in the band.
//
//
// \n \section views_bands_setup Setup of Bands
// <hr>
//
// \image html band.png
// \image latex band.eps "Band view" width=250pt
//
// A reference to a dense or sparse band can be created very conveniently via the \c band()
// function. It can be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Band.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The band index must be in the range from \f$[min(0,1-M)..max(0,N-1)]\f$, where \c M is the
// total number of rows and \c N is the total number of columns, and can be specified both at
// compile time or at runtime:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a reference to the 1st lower band of matrix A (compile time index)
   auto band1 = band<-1L>( A );

   // Creating a reference to the 2nd upper band of matrix A (runtime index)
   auto band2 = band( A, 2L );
   \endcode

// In addition, the \c diagonal() function provides a convenient shortcut for the setup of a view
// on the diagonal of a dense or sparse matrix. It has the same effect as calling the \c band()
// function with a compile time index of 0:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a reference to the diagonal of matrix A via the band() and diagonal() functions
   auto diag1 = band<0L>( A );
   auto diag2 = diagonal( A );

   static_assert( blaze::IsSame< decltype(diag1), decltype(diag2) >::value, "Non-identical types detected" );
   \endcode

// Both the \c band() and the \c diagonal() function return an expression representing the band
// view. The type of this expression depends on the given arguments, primarily the type of the
// matrix and the compile time arguments. If the type is required, it can be determined via
// \c decltype specifier:

   \code
   using MatrixType = blaze::DynamicMatrix<int>;
   using BandType = decltype( blaze::band<1L>( std::declval<MatrixType>() ) );
   using DiagonalType = decltype( blaze::diagonal( std::declval<MatrixType>() ) );
   \endcode

// This resulting view can be treated as any other vector, i.e. it can be assigned to, it can
// be copied from, and it can be used in arithmetic operations. By default, bands are considered
// column vectors, but this setting can be changed via the \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch
// (see \ref transpose_flag). The reference can also be used on both sides of an assignment: The
// band can either be used as an alias to grant write access to a specific band of a matrix
// primitive on the left-hand side of an assignment or to grant read-access to a specific band of
// a matrix primitive or expression on the right-hand side of an assignment. The following example
// demonstrates this in detail:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   blaze::CompressedVector<double,blaze::rowVector> y;
   blaze::DynamicMatrix<double,blaze::rowMajor> A, B;
   blaze::CompressedMatrix<double,blaze::rowMajor> C, D;
   // ... Resizing and initialization

   // Setting the 2nd upper band of matrix A to x
   auto band2 = band( A, 2L );
   band2 = x;

   // Setting the 3rd upper band of matrix B to y
   band( B, 3L ) = y;

   // Setting x to the 2nd lower band of the result of the matrix multiplication
   x = band( A * B, -2L );

   // Setting y to the 2nd upper band of the result of the sparse matrix multiplication
   y = band( C * D, 2L );
   \endcode

// \warning It is the programmer's responsibility to ensure the band does not outlive the viewed
// matrix:

   \code
   // Creating a band on a temporary matrix; results in a dangling reference!
   auto band1 = band<1L>( DynamicMatrix<int>{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } } );
   \endcode

// \n \section views_bands_element_access Element Access
// <hr>
//
// The elements of a band can be directly accessed with the subscript operator:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Creating a view on the 4th upper band of matrix A
   auto band4 = band( A, 4L );

   // Setting the 1st element of the dense band, which corresponds
   // to the 1st element in the 4th upper band of matrix A
   band4[1] = 2.0;
   \endcode

// The numbering of the band elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of elements of the referenced band. Alternatively, the elements of a band
// can be traversed via iterators. Just as with vectors, in case of non-const band, \c begin() and
// \c end() return an iterator, which allows to manipulate the elements, in case of constant bands
// an iterator to immutable elements is returned:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 5th upper band of matrix A
   auto band5 = band( A, 5L );

   // Traversing the elements via iterators to non-const elements
   for( auto it=band5.begin(); it!=band5.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense band value
      ... = *it;  // OK: Read access to the dense band value.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=band5.cbegin(); it!=band5.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense band value.
   }
   \endcode

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 5th band of matrix A
   auto band5 = band( A, 5L );

   // Traversing the elements via iterators to non-const elements
   for( auto it=band5.begin(); it!=band5.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=band5.cbegin(); it!=band5.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_bands_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse band can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 10UL, 100UL );  // Non-initialized 10x100 matrix

   auto diag( band( A, 0L ) );  // Reference to the diagonal of A

   // The subscript operator provides access to all possible elements of the sparse band,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse band, the element is inserted into the band.
   diag[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the band it is inserted into the band, if it is already contained in
   // the band its value is modified.
   diag.set( 45UL, -1.2 );

   // An alternative for inserting elements into the band is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the band.
   diag.insert( 50UL, 3.7 );
   \endcode

// \n \section views_bands_common_operations Common Operations
// <hr>
//
// A band view can be used like any other column vector. This means that with only a few
// exceptions all \ref vector_operations and \ref arithmetic_operations can be used. For instance,
// the current number of band elements can be obtained via the \c size() function, the current
// capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since bands are references to specific bands of a matrix,
// several operations are not possible, such as resizing and swapping. The following example
// shows this by means of a dense band view:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd upper band of matrix A
   auto band2 = band( A, 2L );

   band2.size();          // Returns the number of elements in the band
   band2.capacity();      // Returns the capacity of the band
   band2.nonZeros();      // Returns the number of non-zero elements contained in the band

   band2.resize( 84UL );  // Compilation error: Cannot resize a single band of a matrix

   auto band3 = band( A, 3L );
   swap( band2, band3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_bands_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse bands can be used in all arithmetic operations that any other dense or
// sparse vector can be used in. The following example gives an impression of the use of dense
// bands within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse bands with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::columnVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::columnVector> c( 2UL );
   c[1] = 3.0;

   blaze::DynamicMatrix<double,blaze::rowMajor> A( 4UL, 2UL );  // Non-initialized 4x2 matrix

   auto band1( band( A, 1L ) );  // Reference to the 1st upper band of A
   auto diag ( band( A, 0L ) );  // Reference to the diagonal of A

   band1[0] = 0.0;      // Manual initialization of the 1st upper band of A
   diag = 1.0;          // Homogeneous initialization of the diagonal of A
   band( A, -1L ) = a;  // Dense vector initialization of the 1st lower band of A
   band( A, -2L ) = c;  // Sparse vector initialization of the 2nd lower band of A

   b = diag + a;               // Dense vector/dense vector addition
   b = c + band( A, -1L );     // Sparse vector/dense vector addition
   b = diag * band( A, -2L );  // Component-wise vector multiplication

   band( A, -1L ) *= 2.0;     // In-place scaling of the 1st upper band
   b = band( A, -1L ) * 2.0;  // Scaling of the 1st upper band
   b = 2.0 * band( A, -1L );  // Scaling of the 1st upper band

   band( A, -2L ) += a;              // Addition assignment
   band( A, -2L ) -= c;              // Subtraction assignment
   band( A, -2L ) *= band( A, 0L );  // Multiplication assignment

   double scalar = trans( c ) * band( A, -1L );  // Scalar/dot/inner product between two vectors

   A = band( A, -1L ) * trans( c );  // Outer product between two vectors
   \endcode

// \n Previous: \ref views_column_selections &nbsp; &nbsp; Next: \ref arithmetic_operations
*/
//*************************************************************************************************


//**Arithmetic Operations**************************************************************************
/*!\page arithmetic_operations Arithmetic Operations
//
// \tableofcontents
//
//
// \b Blaze provides the following arithmetic operations for vectors and matrices:
//
// <ul>
//    <li> \ref addition
//       <ul>
//          <li> \ref vector_vector_addition </li>
//          <li> \ref matrix_matrix_addition </li>
//          <li> \ref scalar_addition </li>
//       </ul>
//    </li>
//    <li> \ref subtraction
//       <ul>
//          <li> \ref vector_vector_subtraction </li>
//          <li> \ref matrix_matrix_subtraction </li>
//          <li> \ref scalar_subtraction </li>
//       </ul>
//    </li>
//    <li> \ref scalar_multiplication </li>
//    <li> \ref vector_vector_multiplication
//       <ul>
//          <li> \ref componentwise_multiplication </li>
//          <li> \ref inner_product </li>
//          <li> \ref outer_product </li>
//          <li> \ref cross_product </li>
//          <li> \ref vector_kronecker_product </li>
//       </ul>
//    </li>
//    <li> \ref vector_vector_division </li>
//    <li> \ref matrix_vector_multiplication </li>
//    <li> \ref matrix_matrix_multiplication
//       <ul>
//          <li> \ref schur_product </li>
//          <li> \ref matrix_product </li>
//          <li> \ref matrix_kronecker_product </li>
//       </ul>
//    </li>
// </ul>
//
// \n Previous: \ref views_bands &nbsp; &nbsp; Next: \ref addition
*/
//*************************************************************************************************


//**Addition***************************************************************************************
/*!\page addition Addition
//
// \n \section vector_vector_addition Vector/Vector Addition
// <hr>
//
// The addition of vectors is as intuitive as the addition of scalar values. For the addition of
// any two vectors the addition operator (i.e. \c operator+()) can be used. It even enables the
// addition of dense and sparse vectors:

   \code
   blaze::DynamicVector<int>      v1( 5UL ), v3;
   blaze::CompressedVector<float> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 + v2;  // Addition of a dense and a sparse column vector of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to add vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<int,columnVector>   v1( 5UL );
   blaze::CompressedVector<float,rowVector> v2( 5UL );

   v1 + v2;           // Compilation error: Cannot add a column vector and a row vector
   v1 + trans( v2 );  // OK: Addition of two column vectors
   \endcode

// Also note that the addition of two vectors with the same element type is favorable due to
// possible vectorization of the operation:

   \code
   blaze::DynamicVector<double> v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 + v2;  // Vectorized addition of two double precision vectors
   \endcode

// \n \section outer_sum Outer Sum
// <hr>
//
// The addition between a column vector and a row vector results in the outer sum of the two
// vectors:

   \code
   blaze::StaticVector<int,3UL,columnVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   // Results in the matrix
   //
   //       (  1  5  0  6 )
   //   A = (  4  8  3  9 )
   //       ( -2  2 -3  3 )
   //
   blaze::StaticMatrix<int,3UL,4UL> M1 = v1 + v2;
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   blaze::StaticMatrix<int,3UL,4UL> M1 = trans( v1 ) + v2;
   \endcode

// \n \section matrix_matrix_addition Matrix/Matrix Addition
// <hr>
//
// For the addition of any two matrices the addition operator (i.e. \c operator+()) can be used.
// It even enables the addition of dense and sparse matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedMatrix<size_t,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<float,rowMajor>        M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 + M2;  // Addition of a sparse column-major and a dense row-major matrix of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to add row-major and column-major matrices.
// Note however that in favor of performance the addition of two matrices with the same storage
// order is favorable. The same argument holds for the element type: In case two matrices with
// the same element type are added, the performance can be much higher due to vectorization of
// the operation.

   \code
   blaze::DynamicMatrix<float> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 + M2;  // Vectorized addition of two row-major, single precision dense matrices
   \endcode

// \n \section scalar_addition Scalar Addition
// <hr>
//
// For convenience it is also possible to add a scalar value to a dense vector or dense matrix,
// which has the same effect as adding a uniform vector or matrix. In \b Blaze it is possible to
// use all built-in/fundamental data types except bool as scalar values. Additionally, it is
// possible to use \c std::complex values with the same built-in data types as element type.
// Examples:

   \code
   blaze::StaticVector<int,3UL> v1{ 3, 2, 5, -4, 1, 6 };

   blaze::DynamicVector<int>    v2 = v1 + 2;  // Results in { 5, 4, 7, -2, 3, 8 }
   blaze::CompressedVector<int> v3 = 3 + v1;  // Results in { 6, 5, 8, -1, 4, 9 }
   \endcode

   \code
   blaze::StaticMatrix<int,2UL,3UL> M1{ {  3, 2, 5 },
                                        { -4, 1, 6 } };

   blaze::DynamicMatrix<int>    M2 = M1 + 2;  // Results in { { 5, 4, 7 }, { -2, 3, 8 } }
   blaze::CompressedMatrix<int> M3 = 3 + M1;  // Results in { { 6, 5, 8 }, { -1, 4, 9 } }
   \endcode

// \n Previous: \ref arithmetic_operations &nbsp; &nbsp; Next: \ref subtraction
*/
//*************************************************************************************************


//**Subtraction************************************************************************************
/*!\page subtraction Subtraction
//
// \n \section vector_vector_subtraction Vector/Vector Subtraction
// <hr>
//
// The subtraction of vectors works exactly as intuitive as the addition, but with the subtraction
// operator (i.e. \c operator-()). It also enables the subtraction of dense and sparse vectors:

   \code
   blaze::DynamicVector<int>      v1( 5UL ), v3;
   blaze::CompressedVector<float> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 - v2;  // Subtraction of a dense and a sparse column vector of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that in case of vectors it is only possible to
// subtract vectors with the same transpose flag:

   \code
   blaze::DynamicVector<int,columnVector>   v1( 5UL );
   blaze::CompressedVector<float,rowVector> v2( 5UL );

   v1 - v2;           // Compilation error: Cannot subtract a row vector from a column vector
   v1 - trans( v2 );  // OK: Subtraction of two column vectors
   \endcode

// Also note that the subtraction of two vectors with the same element type is favorable due to
// possible vectorization of the operation:

   \code
   blaze::DynamicVector<double> v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 - v2;  // Vectorized subtraction of two double precision vectors
   \endcode

// \n \section outer_difference Outer Difference
// <hr>
//
// The subtraction between a column vector and a row vector results in the outer difference of
// the two vectors:

   \code
   blaze::StaticVector<int,3UL,columnVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   // Results in the matrix
   //
   //       ( 3 -1  4 -2 )
   //   A = ( 6  2  7  1 )
   //       ( 0 -4  1 -5 )
   //
   StaticMatrix<int,3UL,3UL> M1 = v1 - v2;
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   blaze::StaticMatrix<int,3UL,4UL> M1 = trans( v1 ) - v2;
   \endcode

// \n \section matrix_matrix_subtraction Matrix/Matrix Subtraction
// <hr>
//
// For the subtraction of any two matrices the subtraction operator (i.e. \c operator-()) can be
// used. It even enables the subtraction of dense and sparse matrices:

   \code
   blaze::DynamicMatrix<float,rowMajor>        M1( 7UL, 3UL );
   blaze::CompressedMatrix<size_t,columnMajor> M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 - M2;  // Subtraction of a row-major and a column-major matrix of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to subtract row-major and column-major
// matrices. Note however that in favor of performance the subtraction of two matrices with the
// same storage order is favorable. The same argument holds for the element type: In case two
// matrices with the same element type are subtracted, the performance can be much higher due
// to vectorization of the operation.

   \code
   blaze::DynamicMatrix<float> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 - M2;  // Vectorized subtraction of two row-major, single precision dense matrices
   \endcode

// \n \section scalar_subtraction Scalar Subtraction
// <hr>
//
// For convenience it is also possible to subtract a scalar value from a dense vector or dense
// matrix, which has the same effect as subtracting a uniform vector or matrix. In \b Blaze it is
// possible to use all built-in/fundamental data types except bool as scalar values. Additionally,
// it is possible to use \c std::complex values with the same built-in data types as element type.
// Examples:

   \code
   blaze::StaticVector<int,3UL> v1{ 3, 2, 5, -4, 1, 6 };

   blaze::DynamicVector<int>    v2 = v1 - 2;  // Results in { 1, 0, 3, -6, -1, 4 }
   blaze::CompressedVector<int> v3 = 3 - v1;  // Results in { 0, 1, -2, 7, 2, -3 }
   \endcode

   \code
   blaze::StaticMatrix<int,2UL,3UL> M1{ {  3, 2, 5 },
                                        { -4, 1, 6 } };

   blaze::DynamicMatrix<int>    M2 = M1 - 2;  // Results in { { 1, 0, 3 }, { -6, -1, 4 } }
   blaze::CompressedMatrix<int> M3 = 3 - M1;  // Results in { { 0, 1, -2 }, { 7, 2, -3 } }
   \endcode

// \n Previous: \ref addition &nbsp; &nbsp; Next: \ref scalar_multiplication
*/
//*************************************************************************************************


//**Scalar Multiplication**************************************************************************
/*!\page scalar_multiplication Scalar Multiplication
//
// The scalar multiplication is the multiplication of vector or a matrix with a scalar value.
// Alternatively it is also possible to divide a vector or a matrix by a scalar value. In \b Blaze
// it is possible to use all built-in/fundamental data types except bool as scalar values.
// Additionally, it is possible to use \c std::complex values with the same built-in data types
// as element type.

   \code
   blaze::StaticVector<int,3UL> v1{ 1, 2, 3 };

   blaze::DynamicVector<double>   v2 = v1 * 1.2;    // Scalar multiplication
   blaze::CompressedVector<float> v3 = -0.3F * v1;  // Scalar multiplication
   blaze::DynamicVector<double>   v4 = v1 / 1.2;    // Scalar division
   blaze::CompressedVector<float> v5 = 12.0F / v1;  // Scalar division (only dense vectors)
   \endcode

   \code
   blaze::StaticMatrix<int,3UL,2UL> M1{ { 1, 2 }, { 3, 4 }, { 5, 6 } };

   blaze::DynamicMatrix<double>   M2 = M1 * 1.2;    // Scalar multiplication
   blaze::CompressedMatrix<float> M3 = -0.3F * M1;  // Scalar multiplication
   blaze::DynamicMatrix<double>   M4 = M1 / 1.2;    // Scalar division
   blaze::CompressedMatrix<float> M5 = 12.0F / M1;  // Scalar division (only dense matrices)
   \endcode

// Vectors and matrices cannot be used for as scalar value for scalar multiplications or divisions
// (see the following example). However, each vector and matrix provides the \c scale() function,
// which can be used to scale a vector or matrix element-wise with arbitrary scalar data types:

   \code
   blaze::CompressedMatrix< blaze::StaticMatrix<int,3UL,3UL> > M1;
   blaze::StaticMatrix<int,3UL,3UL> scalar;

   M1 * scalar;  // No scalar multiplication, but matrix/matrix multiplication

   M1.scale( scalar );  // Scalar multiplication
   \endcode

// \n Previous: \ref subtraction &nbsp; &nbsp; Next: \ref componentwise_multiplication
*/
//*************************************************************************************************


//**Vector/Vector Multiplication*******************************************************************
/*!\page vector_vector_multiplication Vector/Vector Multiplication
//
// \n \section componentwise_multiplication Componentwise Multiplication
// <hr>
//
// Multiplying two vectors with the same transpose flag (i.e. either blaze::columnVector or
// blaze::rowVector) via the multiplication operator results in a componentwise multiplication
// of the two vectors:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   CompressedVector<int,columnVector> v1( 17UL );
   DynamicVector<int,columnVector>    v2( 17UL );

   StaticVector<double,10UL,rowVector> v3;
   DynamicVector<double,rowVector>     v4( 10UL );

   // ... Initialization of the vectors

   CompressedVector<int,columnVector> v5( v1 * v2 );  // Componentwise multiplication of a sparse and
                                                      // a dense column vector. The result is a sparse
                                                      // column vector.
   DynamicVector<double,rowVector>    v6( v3 * v4 );  // Componentwise multiplication of two dense row
                                                      // vectors. The result is a dense row vector.
   \endcode

// \n \section inner_product Inner Product / Scalar Product / Dot Product
// <hr>
//
// The multiplication between a row vector and a column vector results in an inner product between
// the two vectors:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{  2, 5, -1 };
   blaze::DynamicVector<int,columnVector> v2{ -1, 3, -2 };

   int result = v1 * v2;  // Results in the value 15
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{  2, 5, -1 };
   blaze::StaticVector<int,3UL,rowVector> v2{ -1, 3, -2 };

   int result = v1 * trans( v2 );  // Also results in the value 15
   \endcode

// Alternatively, either the \c inner() function, the \c dot() function or the comma operator can
// be used for any combination of vectors (row or column vectors) to perform an inner product:

   \code
   blaze::StaticVector<int,3UL,columnVector> v1{  2, 5, -1 };
   blaze::StaticVector<int,3UL,rowVector>    v2{ -1, 3, -2 };

   // All alternatives for the inner product between a column vector and a row vector
   int result1 = trans( v1 ) * trans( v2 );
   int result2 = inner( v1, v2 );
   int result3 = dot( v1, v2 );
   int result4 = (v1,v2);
   \endcode

// When using the comma operator, please note the brackets embracing the inner product expression.
// Due to the low precedence of the comma operator (lower even than the assignment operator) these
// brackets are strictly required for a correct evaluation of the inner product.
//
//
// \n \section outer_product Outer Product
// <hr>
//
// The multiplication between a column vector and a row vector results in the outer product of
// the two vectors:

   \code
   blaze::StaticVector<int,3UL,columnVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   // Results in the matrix
   //
   //       ( -2  6  -4  8 )
   //   A = ( -5 15 -10 20 )
   //       (  1 -3   2 -4 )
   //
   StaticMatrix<int,3UL,3UL> M1 = v1 * v2;
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{  2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   blaze::StaticMatrix<int,3UL,4UL> M1 = trans( v1 ) * v2;
   \endcode

// Alternatively, the \c outer() function can be used for any combination of vectors (row or column
// vectors) to perform an outer product:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{  2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 3, -2, 4 };

   blaze::StaticMatrix<int,3UL,4UL> M1 = outer( v1, v2 );  // Outer product between two row vectors
   \endcode

// \n \section cross_product Cross Product
// <hr>
//
// Two vectors with the same transpose flag can be multiplied via the cross product. The cross
// product between two vectors \f$ a \f$ and \f$ b \f$ is defined as

   \f[
   \left(\begin{array}{*{1}{c}}
   c_0 \\
   c_1 \\
   c_2 \\
   \end{array}\right)
   =
   \left(\begin{array}{*{1}{c}}
   a_1 b_2 - a_2 b_1 \\
   a_2 b_0 - a_0 b_2 \\
   a_0 b_1 - a_1 b_0 \\
   \end{array}\right).
   \f]

// Due to the absence of a \f$ \times \f$ operator in the C++ language, the cross product is
// realized via the \c cross() function. Alternatively, the modulo operator (i.e. \c operator%)
// can be used in case infix notation is required:

   \code
   blaze::StaticVector<int,3UL,columnVector> v1{  2, 5, -1 };
   blaze::DynamicVector<int,columnVector>    v2{ -1, 3, -2 };

   blaze::StaticVector<int,3UL,columnVector> v3( cross( v1, v2 ) );
   blaze::StaticVector<int,3UL,columnVector> v4( v1 % v2 );
   \endcode

// Please note that the cross product is restricted to three dimensional (dense and sparse)
// column vectors.
//
//
// \n \section vector_kronecker_product Kronecker Product
// <hr>
//
// The Kronecker product of two vectors with the same transpose flag can be computed via the
// \a kron() function:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   DynamicVector<double>   v1( 28UL );
   CompressedVector<float> v2( 17UL );

   // ... Initialization of the vectors

   CompressedVector<double> v3 = kron( v1, v2 );
   \endcode

// Both dense and sparse vectors can be used for a Kronecker product. It is possible to multiply
// two vectors with different element type, as long as the element types themselves can be
// multiplied.
//
// \n Previous: \ref scalar_multiplication &nbsp; &nbsp; Next: \ref vector_vector_division
*/
//*************************************************************************************************


//**Vector/Vector Division*************************************************************************
/*!\page vector_vector_division Vector/Vector Division
//
// \n \section componentwise_division Componentwise Division
// <hr>
//
// Dividing a vector by a dense vector with the same transpose flag (i.e. either blaze::columnVector
// or blaze::rowVector) via the division operator results in a componentwise division:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   CompressedVector<int,columnVector> v1( 17UL );
   DynamicVector<int,columnVector>    v2( 17UL );

   StaticVector<double,10UL,rowVector> v3;
   DynamicVector<double,rowVector>     v4( 10UL );

   // ... Initialization of the vectors

   CompressedVector<int,columnVector> v5( v1 / v2 );  // Componentwise division of a sparse and a
                                                      // dense column vector. The result is a sparse
                                                      // column vector.
   DynamicVector<double,rowVector>    v6( v3 / v4 );  // Componentwise division of two dense row
                                                      // vectors. The result is a dense row vector.
   \endcode

// Note that all values of the divisor must be non-zero and that no checks are performed to assert
// this precondition!
//
//
// \n \section outer_quotient Outer Quotient
// <hr>
//
// The division between a column vector and a row vector results in the outer quotient of the
// two vectors:

   \code
   blaze::StaticVector<double,3UL,columnVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<double,rowVector> v2{ -1, 5, -2, 4 };

   // Results in the matrix
   //
   //       ( -2  0.4   -1   0.5 )
   //   A = ( -5    1 -2.5  1.25 )
   //       (  1 -0.2  0.5 -0.25 )
   //
   blaze::StaticMatrix<int,3UL,4UL> M1 = v1 / v2;
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1{ 2, 5, -1 };
   blaze::DynamicVector<int,rowVector> v2{ -1, 5, -2, 4 };

   blaze::StaticMatrix<int,3UL,4UL> M1 = trans( v1 ) / v2;
   \endcode

// Note that all values of the divisor must be non-zero and that no checks are performed to assert
// this precondition!
//
// \n Previous: \ref vector_vector_multiplication &nbsp; &nbsp; Next: \ref matrix_vector_multiplication
*/
//*************************************************************************************************


//**Matrix/Vector Multiplication*******************************************************************
/*!\page matrix_vector_multiplication Matrix/Vector Multiplication
//
// In \b Blaze matrix/vector multiplications can be as intuitively formulated as in mathematical
// textbooks. Just as in textbooks there are two different multiplications between a matrix and
// a vector: a matrix/column vector multiplication and a row vector/matrix multiplication:

   \code
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;

   DynamicMatrix<int>                  M1( 39UL, 12UL );
   StaticVector<int,12UL,columnVector> v1;

   // ... Initialization of the matrix and the vector

   DynamicVector<int,columnVector> v2 = M1 * v1;           // Matrix/column vector multiplication
   DynamicVector<int,rowVector>    v3 = trans( v1 ) * M1;  // Row vector/matrix multiplication
   \endcode

// Note that the storage order of the matrix poses no restrictions on the operation. Also note,
// that the highest performance for a multiplication between a dense matrix and a dense vector can
// be achieved if both the matrix and the vector have the same scalar element type.
//
// \n Previous: \ref vector_vector_division &nbsp; &nbsp; Next: \ref matrix_matrix_multiplication
*/
//*************************************************************************************************


//**Matrix/Matrix Multiplication*******************************************************************
/*!\page matrix_matrix_multiplication Matrix/Matrix Multiplication
//
// \n \section schur_product Componentwise Multiplication / Schur Product
// <hr>
//
// Multiplying two matrices with the same dimensions (i.e. the same number of rows and columns)
// via the modulo operator results in a componentwise multiplication (Schur product) of the two
// matrices:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   DynamicMatrix<double>   M1( 28UL, 35UL );
   CompressedMatrix<float> M2( 28UL, 35UL );

   // ... Initialization of the matrices

   DynamicMatrix<double> M3 = M1 % M2;
   \endcode

// Both dense and sparse matrices can be used for a Schur product. The storage order of the two
// matrices poses no restrictions on the operation, all variations are possible. It is also
// possible to multiply two matrices with different element type, as long as the element types
// themselves can be multiplied.
//
//
// \n \section matrix_product Matrix Product
// <hr>
//
// The matrix/matrix product can be formulated exactly as in mathematical textbooks:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   DynamicMatrix<double>   M1( 45UL, 85UL );
   CompressedMatrix<float> M2( 85UL, 37UL );

   // ... Initialization of the matrices

   DynamicMatrix<double> M3 = M1 * M2;
   \endcode

// The storage order of the two matrices poses no restrictions on the operation, all variations
// are possible. It is also possible to multiply two matrices with different element type, as
// long as the element types themselves can be multiplied and added. Note however that the
// highest performance for a multiplication between two matrices can be expected for two
// matrices with the same scalar element type.
//
// In case the resulting matrix is known to be symmetric, Hermitian, lower triangular, upper
// triangular, or diagonal, the computation can be optimized by explicitly declaring the
// multiplication as symmetric, Hermitian, lower triangular, upper triangular, or diagonal by
// means of the \ref matrix_operations_declaration_operations :

   \code
   using blaze::DynamicMatrix;

   DynamicMatrix<double> M1, M2, M3;

   // ... Initialization of the square matrices

   M3 = declsym ( M1 * M2 );  // Declare the result of the matrix multiplication as symmetric
   M3 = declherm( M1 * M2 );  // Declare the result of the matrix multiplication as Hermitian
   M3 = decllow ( M1 * M2 );  // Declare the result of the matrix multiplication as lower triangular
   M3 = declupp ( M1 * M2 );  // Declare the result of the matrix multiplication as upper triangular
   M3 = decldiag( M1 * M2 );  // Declare the result of the matrix multiplication as diagonal
   \endcode

// Using a declaration operation on the a multiplication expression can speed up the computation
// by a factor of 2. Note however that the caller of the according declaration operation takes
// full responsibility for the correctness of the declaration. Falsely declaring a multiplication
// as symmetric, Hermitian, lower triangular, upper triangular, or diagonal leads to undefined
// behavior!
//
//
// \n \section matrix_kronecker_product Kronecker Product
// <hr>
//
// The Kronecker product of two matrices can be computed via the \a kron() function:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   DynamicMatrix<double>   M1( 28UL, 35UL );
   CompressedMatrix<float> M2( 17UL, 11UL );

   // ... Initialization of the matrices

   CompressedMatrix<double> M3 = kron( M1, M2 );
   \endcode

// Both dense and sparse matrices can be used for a Kronecker product. The storage order of the
// two matrices poses no restrictions on the operation, all variations are possible. It is also
// possible to multiply two matrices with different element type, as long as the element types
// themselves can be multiplied.
//
// \n Previous: \ref matrix_vector_multiplication &nbsp; &nbsp; Next: \ref bitwise_operations
*/
//*************************************************************************************************


//**Bitwise Operations*****************************************************************************
/*!\page bitwise_operations Bitwise Operations
//
// \tableofcontents
//
//
// \b Blaze provides the following bitwise operations for vectors and matrices:
//
// <ul>
//    <li> \ref bitwise_shift
//       <ul>
//          <li> \ref vector_vector_shift </li>
//          <li> \ref matrix_matrix_shift </li>
//          <li> \ref scalar_shift </li>
//       </ul>
//    </li>
//    <li> \ref bitwise_and
//       <ul>
//          <li> \ref vector_vector_bitand </li>
//          <li> \ref matrix_matrix_bitand </li>
//          <li> \ref scalar_bitand </li>
//       </ul>
//    </li>
//    <li> \ref bitwise_or
//       <ul>
//          <li> \ref vector_vector_bitor </li>
//          <li> \ref matrix_matrix_bitor </li>
//          <li> \ref scalar_bitor </li>
//       </ul>
//    </li>
//    <li> \ref bitwise_xor
//       <ul>
//          <li> \ref vector_vector_bitxor </li>
//          <li> \ref matrix_matrix_bitxor </li>
//          <li> \ref scalar_bitxor </li>
//       </ul>
//    </li>
// </ul>
//
// \n Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref bitwise_shift
*/
//*************************************************************************************************


//**Bitwise Shift**********************************************************************************
/*!\page bitwise_shift Bitwise Shift
//
// \n \section vector_vector_shift Vector/Vector Shift
// <hr>
//
// Via the left-shift operator (i.e. operator<<()) and the right-shift operator (i.e. operator>>())
// it is possible to perform an elementwise shift of a dense vector:

   \code
   blaze::DynamicVector<unsigned int> v1( 5UL ), v3;
   blaze::DynamicVector<unsigned short> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 << v2;  // Elementwise left-shift of a dense column vector
   v3 = v1 >> v2;  // Elementwise right-shift of a dense column vector
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to shift vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<unsigned int,columnVector> v1( 5UL );
   blaze::DynamicVector<unsigned int,rowVector>    v2( 5UL );

   v1 << v2;           // Compilation error: Cannot shift a column vector by a row vector
   v1 << trans( v2 );  // OK: Shifting a column vector by another column vector
   \endcode

// Furthermore, it is possible to use different element types in the two vector operands, but
// shifting two vectors with the same element type is favorable due to possible vectorization
// of the operation:

   \code
   blaze::DynamicVector<unsigned int> v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 << v2;  // Vectorized left-shift of an unsigned int vector
   \endcode

// \n \section matrix_matrix_shift Matrix/Matrix Shift
// <hr>
//
// The left-shift operator (i.e. operator<<()) and the right-shift operator (i.e. operator>>())
// can also be used to perform an elementwise shift of a dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<unsigned int,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<unsigned short,rowMajor>  M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 << M2;  // Elementwise left-shift of a dense column-major matrix
   M3 = M1 >> M2;  // Elementwise right-shift of a dense column-major matrix
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to use any combination of row-major and
// column-major matrices. Note however that in favor of performance using two matrices with the
// same storage order is favorable. The same argument holds for the element type: While it is
// possible to use matrices with different element type, using two matrices with the same element
// type potentially leads to better performance due to vectorization of the operation.

   \code
   blaze::DynamicMatrix<unsigned int> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 << M2;  // Vectorized left-shift of an unsigned int matrix
   \endcode

// \n \section scalar_shift Scalar Shift
// <hr>
//
// It is also possible to uniformly shift all elements of a dense vector or dense matrix by means
// of a scalar, which has the same effect as shifting by means of a uniform vector or matrix (see
// \ref vector_types_uniform_vector and \ref matrix_types_uniform_matrix). In \b Blaze it is
// possible to use all built-in/fundamental data types except bool as scalar values. Examples:

   \code
   blaze::DynamicVector<unsigned int> v1{ 3, 2, 5, 4, 1, 6 };

   // Uniform left-shift by one bit of all elements of v1; Results in
   //
   //    ( 6, 4, 10, 8, 2, 12 )
   //
   blaze::DynamicVector<int> v2( v1 << 1U );
   \endcode

   \code
   blaze::DynamicMatrix<unsigned int> M1{ { 3, 2, 5 },
                                          { 4, 1, 6 } };

   // Uniform left-shift by one bit of all elements of M1; Results in
   //
   //    ( 6, 4, 10 )
   //    ( 8, 2, 12 )
   //
   blaze::DynamicMatrix<unsigned int> M2( M1 << 1U );
   \endcode

// \n Previous: \ref bitwise_operations &nbsp; &nbsp; Next: \ref bitwise_and
*/
//*************************************************************************************************


//**Bitwise AND************************************************************************************
/*!\page bitwise_and Bitwise AND
//
// \n \section vector_vector_bitand Vector/Vector Bitwise AND
// <hr>
//
// Via the bitwise AND operator (i.e. operator&()) it is possible to perform an elementwise
// bitwise AND with dense vectors:

   \code
   blaze::DynamicVector<unsigned int> v1( 5UL ), v3;
   blaze::DynamicVector<unsigned short> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 & v2;  // Elementwise bitwise AND of two dense column vectors of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to use vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<unsigned int,columnVector> v1( 5UL );
   blaze::DynamicVector<unsigned int,rowVector>    v2( 5UL );

   v1 & v2;           // Compilation error: Cannot AND a column vector and a row vector
   v1 & trans( v2 );  // OK: Bitwise AND of two column vectors
   \endcode

// Furthermore, it is possible to use different element types in the two vector operands, but a
// bitwise AND of two vectors with the same element type is favorable due to possible vectorization
// of the operation:

   \code
   blaze::DynamicVector<unsigned int> v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 & v2;  // Vectorized bitwise AND of an unsigned int vector
   \endcode

// \n \section matrix_matrix_bitand Matrix/Matrix Bitwise AND
// <hr>
//
// The bitwise AND operator (i.e. operator&()) can also be used to perform an elementwise bitwise
// AND with dense matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<unsigned int,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<unsigned short,rowMajor>  M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 & M2;  // Elementwise bitwise AND of two dense matrices of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to use any combination of row-major and
// column-major matrices. Note however that in favor of performance using two matrices with the
// same storage order is favorable. The same argument holds for the element type: While it is
// possible to use matrices with different element type, using two matrices with the same element
// type potentially leads to better performance due to vectorization of the operation.

   \code
   blaze::DynamicMatrix<unsigned int> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 & M2;  // Vectorized bitwise AND of two row-major, unsigned int dense matrices
   \endcode

// \n \section scalar_bitand Scalar Bitwise AND
// <hr>
//
// Is is also possible to perform a bitwise AND between a dense vector or dense matrix and a
// scalar value, which has the same effect as performing a bitwise AND by means of a uniform
// vector or matrix (see \ref vector_types_uniform_vector and \ref matrix_types_uniform_matrix).
// In \b Blaze it is possible to use all built-in/fundamental data types except bool as scalar
// values. Examples:

   \code
   blaze::DynamicVector<unsigned int> v1{ 3U, 2U, 5U, 4U, 1U, 6U };

   // Perform a bitwise AND with all elements of v1; Results in
   //
   //    ( 3, 2, 1, 0, 1, 2 )
   //
   blaze::DynamicVector<int> v2( v1 & 3U );
   \endcode

   \code
   blaze::DynamicMatrix<unsigned int> M1{ { 3U, 2U, 5U },
                                          { 4U, 1U, 6U } };

   // Perform a bitwise AND with all elements of M1; Results in
   //
   //    ( 3, 2, 1 )
   //    ( 0, 1, 2 )
   //
   blaze::DynamicMatrix<unsigned int> M2( M1 & 3U );
   \endcode

// \n Previous: \ref bitwise_shift &nbsp; &nbsp; Next: \ref bitwise_or
*/
//*************************************************************************************************


//**Bitwise OR*************************************************************************************
/*!\page bitwise_or Bitwise OR
//
// \n \section vector_vector_bitor Vector/Vector Bitwise OR
// <hr>
//
// Via the bitwise OR operator (i.e. operator|()) it is possible to perform an elementwise
// bitwise OR with dense vectors:

   \code
   blaze::DynamicVector<unsigned int> v1( 5UL ), v3;
   blaze::DynamicVector<unsigned short> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 | v2;  // Elementwise bitwise OR of two dense column vectors of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to use vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<unsigned int,columnVector> v1( 5UL );
   blaze::DynamicVector<unsigned int,rowVector>    v2( 5UL );

   v1 | v2;           // Compilation error: Cannot OR a column vector and a row vector
   v1 | trans( v2 );  // OK: Bitwise OR of two column vectors
   \endcode

// Furthermore, it is possible to use different element types in the two vector operands, but a
// bitwise OR of two vectors with the same element type is favorable due to possible vectorization
// of the operation:

   \code
   blaze::DynamicVector<unsigned int> v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 | v2;  // Vectorized bitwise OR of an unsigned int vector
   \endcode

// \n \section matrix_matrix_bitor Matrix/Matrix Bitwise OR
// <hr>
//
// The bitwise OR operator (i.e. operator|()) can also be used to perform an elementwise bitwise
// OR with dense matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<unsigned int,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<unsigned short,rowMajor>  M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 | M2;  // Elementwise bitwise OR of two dense matrices of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to use any combination of row-major and
// column-major matrices. Note however that in favor of performance using two matrices with the
// same storage order is favorable. The same argument holds for the element type: While it is
// possible to use matrices with different element type, using two matrices with the same element
// type potentially leads to better performance due to vectorization of the operation.

   \code
   blaze::DynamicMatrix<unsigned int> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 | M2;  // Vectorized bitwise OR of two row-major, unsigned int dense matrices
   \endcode

// \n \section scalar_bitor Scalar Bitwise OR
// <hr>
//
// Is is also possible to perform a bitwise OR between a dense vector or dense matrix and a
// scalar value, which has the same effect as performing a bitwise OR by means of a uniform
// vector or matrix (see \ref vector_types_uniform_vector and \ref matrix_types_uniform_matrix).
// In \b Blaze it is possible to use all built-in/fundamental data types except bool as scalar
// values. Examples:

   \code
   blaze::DynamicVector<unsigned int> v1{ 3U, 2U, 5U, 4U, 1U, 6U };

   // Perform a bitwise OR with all elements of v1; Results in
   //
   //    ( 3, 3, 7, 7, 3, 3 )
   //
   blaze::DynamicVector<int> v2( v1 | 3U );
   \endcode

   \code
   blaze::DynamicMatrix<unsigned int> M1{ { 3U, 2U, 5U },
                                          { 4U, 1U, 6U } };

   // Perform a bitwise OR with all elements of M1; Results in
   //
   //    ( 3, 3, 7 )
   //    ( 7, 3, 3 )
   //
   blaze::DynamicMatrix<unsigned int> M2( M1 | 3U );
   \endcode

// \n Previous: \ref bitwise_and &nbsp; &nbsp; Next: \ref bitwise_xor
*/
//*************************************************************************************************


//**Bitwise XOR************************************************************************************
/*!\page bitwise_xor Bitwise XOR
//
// \n \section vector_vector_bitxor Vector/Vector Bitwise XOR
// <hr>
//
// Via the bitwise XOR operator (i.e. operator^()) it is possible to perform an elementwise
// bitwise XOR with dense vectors:

   \code
   blaze::DynamicVector<unsigned int> v1( 5UL ), v3;
   blaze::DynamicVector<unsigned short> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 ^ v2;  // Elementwise bitwise XOR of two dense column vectors of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to use vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<unsigned int,columnVector> v1( 5UL );
   blaze::DynamicVector<unsigned int,rowVector>    v2( 5UL );

   v1 ^ v2;           // Compilation error: Cannot XOR a column vector and a row vector
   v1 ^ trans( v2 );  // OK: Bitwise XOR of two column vectors
   \endcode

// Furthermore, it is possible to use different element types in the two vector operands, but a
// bitwise XOR of two vectors with the same element type is favorable due to possible vectorization
// of the operation:

   \code
   blaze::DynamicVector<unsigned int> v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 ^ v2;  // Vectorized bitwise XOR of an unsigned int vector
   \endcode

// \n \section matrix_matrix_bitxor Matrix/Matrix Bitwise XOR
// <hr>
//
// The bitwise XOR operator (i.e. operator^()) can also be used to perform an elementwise bitwise
// XOR with dense matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<unsigned int,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<unsigned short,rowMajor>  M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 ^ M2;  // Elementwise bitwise XOR of two dense matrices of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to use any combination of row-major and
// column-major matrices. Note however that in favor of performance using two matrices with the
// same storage order is favorable. The same argument holds for the element type: While it is
// possible to use matrices with different element type, using two matrices with the same element
// type potentially leads to better performance due to vectorization of the operation.

   \code
   blaze::DynamicMatrix<unsigned int> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 ^ M2;  // Vectorized bitwise XOR of two row-major, unsigned int dense matrices
   \endcode

// \n \section scalar_bitxor Scalar Bitwise XOR
// <hr>
//
// Is is also possible to perform a bitwise XOR between a dense vector or dense matrix and a
// scalar value, which has the same effect as performing a bitwise XOR by means of a uniform
// vector or matrix (see \ref vector_types_uniform_vector and \ref matrix_types_uniform_matrix).
// In \b Blaze it is possible to use all built-in/fundamental data types except bool as scalar
// values. Examples:

   \code
   blaze::DynamicVector<unsigned int> v1{ 3U, 2U, 5U, 4U, 1U, 6U };

   // Perform a bitwise XOR with all elements of v1; Results in
   //
   //    ( 0, 1, 6, 7, 2, 5 )
   //
   blaze::DynamicVector<int> v2( v1 ^ 3U );
   \endcode

   \code
   blaze::DynamicMatrix<unsigned int> M1{ { 3U, 2U, 5U },
                                          { 4U, 1U, 6U } };

   // Perform a bitwise XOR with all elements of M1; Results in
   //
   //    ( 0, 1, 6 )
   //    ( 7, 2, 5 )
   //
   blaze::DynamicMatrix<unsigned int> M2( M1 ^ 3U );
   \endcode

// \n Previous: \ref bitwise_or &nbsp; &nbsp; Next: \ref logical_operations
*/
//*************************************************************************************************


//**Logical Operations*****************************************************************************
/*!\page logical_operations Logical Operations
//
// \tableofcontents
//
//
// \b Blaze provides the following logical operations for vectors and matrices:
//
// <ul>
//    <li> \ref logical_not
//       <ul>
//          <li> \ref vector_vector_not </li>
//          <li> \ref matrix_matrix_not </li>
//       </ul>
//    </li>
//    <li> \ref logical_and
//       <ul>
//          <li> \ref vector_vector_and </li>
//          <li> \ref matrix_matrix_and </li>
//       </ul>
//    </li>
//    <li> \ref logical_or
//       <ul>
//          <li> \ref vector_vector_or </li>
//          <li> \ref matrix_matrix_or </li>
//       </ul>
//    </li>
// </ul>
//
// \n Previous: \ref bitwise_xor &nbsp; &nbsp; Next: \ref logical_not
*/
//*************************************************************************************************


//**Logical NOT************************************************************************************
/*!\page logical_not Logical NOT
//
// \n \section vector_vector_not Vector/Vector Logical NOT
// <hr>
//
// Via the logical NOT operator (i.e. operator!()) it is possible to compute an elementwise
// logical NOT of a dense vector:

   \code
   blaze::DynamicVector<bool> v1( 5UL ), v2;

   // ... Initializing the vectors

   v2 = !v1;  // Elementwise logical NOT of a dense column vector
   \endcode

// \n \section matrix_matrix_not Matrix/Matrix Logical NOT
// <hr>
//
// The logical NOT operator (i.e. operator!()) can also be used to compute an elementwise logical
// NOT with dense matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<bool,rowMajor> M1( 7UL, 3UL ), M2;

   // ... Initializing the matrices

   M2 = !M1;  // Elementwise logical NOT of a dense row-major matrix
   \endcode

// \n Previous: \ref logical_operations &nbsp; &nbsp; Next: \ref logical_and
*/
//*************************************************************************************************


//**Logical AND************************************************************************************
/*!\page logical_and Logical AND
//
// \n \section vector_vector_and Vector/Vector Logical AND
// <hr>
//
// Via the logical AND operator (i.e. operator&&()) it is possible to compute an elementwise
// logical AND with dense vectors:

   \code
   blaze::DynamicVector<bool> v1( 5UL ), v3;
   blaze::DynamicVector<bool> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 && v2;  // Elementwise logical AND of two dense column vectors
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to use vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<bool,columnVector> v1( 5UL );
   blaze::DynamicVector<bool,rowVector>    v2( 5UL );

   v1 && v2;           // Compilation error: Cannot AND a column vector and a row vector
   v1 && trans( v2 );  // OK: Logical AND of two column vectors
   \endcode

// \n \section matrix_matrix_and Matrix/Matrix Logical AND
// <hr>
//
// The logical AND operator (i.e. operator&&()) can also be used to compute an elementwise logical
// AND with dense matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<bool,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<bool,rowMajor>    M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 && M2;  // Elementwise logical AND of two dense matrices
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to use any combination of row-major and
// column-major matrices. Note however that in favor of performance using two matrices with the
// same storage order is favorable.
//
// \n Previous: \ref logical_not &nbsp; &nbsp; Next: \ref logical_or
*/
//*************************************************************************************************


//**Logical OR*************************************************************************************
/*!\page logical_or Logical OR
//
// \n \section vector_vector_or Vector/Vector Logical OR
// <hr>
//
// Via the logical OR operator (i.e. operator||()) it is possible to perform an elementwise
// logical OR with dense vectors:

   \code
   blaze::DynamicVector<bool> v1( 5UL ), v3;
   blaze::DynamicVector<bool> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 || v2;  // Elementwise logical OR of two dense column vectors
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that it is only possible to use vectors with
// the same transpose flag:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<unsigned int,columnVector> v1( 5UL );
   blaze::DynamicVector<unsigned int,rowVector>    v2( 5UL );

   v1 || v2;           // Compilation error: Cannot OR a column vector and a row vector
   v1 || trans( v2 );  // OK: Logical OR of two column vectors
   \endcode

// \n \section matrix_matrix_or Matrix/Matrix Logical OR
// <hr>
//
// The logical OR operator (i.e. operator||()) can also be used to perform an elementwise logical
// OR with dense matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<bool,columnMajor> M1( 7UL, 3UL );
   blaze::DynamicMatrix<bool,rowMajor>    M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 || M2;  // Elementwise logical OR of two dense matrices
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. It is possible to use any combination of row-major and
// column-major matrices. Note however that in favor of performance using two matrices with the
// same storage order is favorable.
//
// \n Previous: \ref logical_and &nbsp; &nbsp; Next: \ref shared_memory_parallelization
*/
//*************************************************************************************************


//**Shared Memory Parallelization******************************************************************
/*!\page shared_memory_parallelization Shared Memory Parallelization
//
// For all possible operations \b Blaze tries to achieve maximum performance on a single CPU
// core. However, today's CPUs are not single core anymore, but provide several (homogeneous
// or heterogeneous) compute cores. In order to fully exploit the performance potential of a
// multicore CPU, computations have to be parallelized across all available cores of a CPU.
// For this purpose, \b Blaze provides four different shared memory parallelization techniques:
//
//  - \ref hpx_parallelization
//  - \ref cpp_threads_parallelization
//  - \ref boost_threads_parallelization
//  - \ref openmp_parallelization
//
// When any of the shared memory parallelization techniques is activated, all arithmetic
// operations on dense vectors and matrices (including additions, subtractions, multiplications,
// divisions, and all componentwise arithmetic operations) and most operations on sparse vectors
// and matrices are automatically run in parallel. However, in addition, \b Blaze provides means
// to enforce the serial execution of specific operations:
//
//  - \ref serial_execution
//
// \n Previous: \ref logical_or &nbsp; &nbsp; Next: \ref hpx_parallelization
*/
//*************************************************************************************************


//**HPX Parallelization****************************************************************************
/*!\page hpx_parallelization HPX Parallelization
//
// \tableofcontents
//
//
// The first shared memory parallelization provided with \b Blaze is based on
// <a href="http://stellar.cct.lsu.edu/projects/hpx/">HPX</a>.
//
//
// \n \section hpx_setup HPX Setup
// <hr>
//
// In order to enable the HPX-based parallelization, the following steps have to be taken: First,
// the \c BLAZE_USE_HPX_THREADS command line argument has to be explicitly specified during
// compilation:

   \code
   ... -DBLAZE_USE_HPX_THREADS ...
   \endcode

// Second, the HPX library and depending libraries such as Boost, hwloc, etc. have to be linked.
// And third, the HPX threads have to be initialized by a call to the \c hpx::init() function (see
// the <a href="http://stellar.cct.lsu.edu/files/hpx_0.9.0/docs/hpx/tutorial.html">HPX tutorial</a>
// for further details). These three actions will cause the \b Blaze library to automatically try
// to run all operations in parallel with the specified number of HPX threads.
//
// Note that the HPX-based parallelization has priority over the OpenMP-based, C++11 thread-based,
// and Boost thread-based parallelizations, i.e. is preferred in case multiple parallelizations
// are enabled in combination with the HPX thread parallelization.
//
// The number of threads used by the HPX backend has to be specified via the command line:

   \code
   ... --hpx:threads 4 ...
   \endcode

// Please note that the \b Blaze library does not limit the available number of threads. Therefore
// it is in YOUR responsibility to choose an appropriate number of threads. The best performance,
// though, can be expected if the specified number of threads matches the available number of
// cores.
//
// In order to query the number of threads used for the parallelization of operations, the
// \c getNumThreads() function can be used:

   \code
   const size_t threads = blaze::getNumThreads();
   \endcode

// In the context of HPX threads, the function will return the actual number of threads used by
// the HPX subsystem.
//
//
// \n \section hpx_configuration HPX Configuration
// <hr>
//
// As in case of the other shared memory parallelizations \b Blaze is not unconditionally running
// an operation in parallel (see for instance \ref openmp_parallelization). Only in case a given
// operation is large enough and exceeds a certain threshold the operation is executed in parallel.
// All thresholds related to the HPX-based parallelization are contained within the configuration
// file <tt><blaze/config/Thresholds.h></tt>.
//
// Please note that these thresholds are highly sensitiv to the used system architecture and
// the shared memory parallelization technique. Therefore the default values cannot guarantee
// maximum performance for all possible situations and configurations. They merely provide a
// reasonable standard for the current CPU generation. Also note that the provided defaults
// have been determined using the OpenMP parallelization and require individual adaption for
// the HPX-based parallelization.
//
// \n Previous: \ref shared_memory_parallelization &nbsp; &nbsp; Next: \ref cpp_threads_parallelization
*/
//*************************************************************************************************


//**C++11 Thread Parallelization*******************************************************************
/*!\page cpp_threads_parallelization C++11 Thread Parallelization
//
// \tableofcontents
//
//
// In addition to the HPX-based shared memory parallelization, starting with \b Blaze 2.1,
// \b Blaze also provides a shared memory parallelization based on C++11 threads.
//
//
// \n \section cpp_threads_setup C++11 Thread Setup
// <hr>
//
// In order to enable the C++11 thread-based parallelization, first the according C++11-specific
// compiler flags have to be used and second the \c BLAZE_USE_CPP_THREADS command line argument
// has to be explicitly specified. For instance, in case of the GNU C++ and Clang compilers the
// compiler flags have to be extended by

   \code
   ... -std=c++11 -DBLAZE_USE_CPP_THREADS ...
   \endcode

// This simple action will cause the \b Blaze library to automatically try to run all operations
// in parallel with the specified number of C++11 threads. Note that in case both HPX and C++11
// threads are enabled on the command line, the HPX-based parallelization has priority and is
// preferred.
//
// The number of threads can be either specified via the environment variable \c BLAZE_NUM_THREADS

   \code
   export BLAZE_NUM_THREADS=4  // Unix systems
   set BLAZE_NUM_THREADS=4     // Windows systems
   \endcode

// or alternatively via the \c setNumThreads() function provided by the \b Blaze library:

   \code
   blaze::setNumThreads( 4 );
   \endcode

// Please note that the \b Blaze library does not limit the available number of threads. Therefore
// it is in YOUR responsibility to choose an appropriate number of threads. The best performance,
// though, can be expected if the specified number of threads matches the available number of
// cores.
//
// In order to query the number of threads used for the parallelization of operations, the
// \c getNumThreads() function can be used:

   \code
   const size_t threads = blaze::getNumThreads();
   \endcode

// In the context of C++11 threads, the function will return the previously specified number of
// threads.
//
//
// \n \section cpp_threads_configuration C++11 Thread Configuration
// <hr>
//
// As in case of the OpenMP-based parallelization \b Blaze is not unconditionally running an
// operation in parallel. In case \b Blaze deems the parallel execution as counterproductive for
// the overall performance, the operation is executed serially. One of the main reasons for not
// executing an operation in parallel is the size of the operands. For instance, a vector addition
// is only executed in parallel if the size of both vector operands exceeds a certain threshold.
// Otherwise, the performance could seriously decrease due to the overhead caused by the thread
// setup. However, in order to be able to adjust the \b Blaze library to a specific system, it
// is possible to configure these thresholds manually. All thresholds are contained within the
// configuration file <tt><blaze/config/Thresholds.h></tt>.
//
// Please note that these thresholds are highly sensitiv to the used system architecture and
// the shared memory parallelization technique. Therefore the default values cannot guarantee
// maximum performance for all possible situations and configurations. They merely provide a
// reasonable standard for the current CPU generation. Also note that the provided defaults
// have been determined using the OpenMP parallelization and require individual adaption for
// the C++11 thread parallelization.
//
//
// \n \section cpp_threads_known_issues Known Issues
// <hr>
//
// There is a known issue in Visual Studio 2012 and 2013 that may cause C++11 threads to hang
// if their destructor is executed after the \c main() function:
//
//    http://connect.microsoft.com/VisualStudio/feedback/details/747145
//
// Unfortunately, the C++11 parallelization of the \b Blaze library is affected from this bug.
// In order to circumvent this problem, \b Blaze provides the \c shutDownThreads() function,
// which can be used to manually destroy all threads at the end of the \c main() function:

   \code
   int main()
   {
      // ... Using the C++11 thread parallelization of Blaze

      shutDownThreads();
   }
   \endcode

// Please note that this function may only be used at the end of the \c main() function. After
// this function no further computation may be executed! Also note that this function has an
// effect for Visual Studio compilers only and doesn't need to be used with any other compiler.
//
// \n Previous: \ref hpx_parallelization &nbsp; &nbsp; Next: \ref boost_threads_parallelization
*/
//*************************************************************************************************


//**Boost Thread Parallelization*******************************************************************
/*!\page boost_threads_parallelization Boost Thread Parallelization
//
// \tableofcontents
//
//
// The third available shared memory parallelization provided with \b Blaze is based
// on <a href="https://www.boost.org/doc/libs/1_68_0/doc/html/thread.html">Boost threads</a>.
//
//
// \n \section boost_threads_setup Boost Thread Setup
// <hr>
//
// In order to enable the Boost thread-based parallelization, two steps have to be taken: First,
// the \c BLAZE_USE_BOOST_THREADS command line argument has to be explicitly specified during
// compilation:

   \code
   ... -DBLAZE_USE_BOOST_THREADS ...
   \endcode

// Second, the according Boost libraries have to be linked. These two simple actions will cause
// the \b Blaze library to automatically try to run all operations in parallel with the specified
// number of Boost threads. Note that the HPX-based and C++11 thread-based parallelizations have
// priority, i.e. are preferred in case either is enabled in combination with the Boost thread
// parallelization.
//
// The number of threads can be either specified via the environment variable \c BLAZE_NUM_THREADS

   \code
   export BLAZE_NUM_THREADS=4  // Unix systems
   set BLAZE_NUM_THREADS=4     // Windows systems
   \endcode

// or alternatively via the \c setNumThreads() function provided by the \b Blaze library:

   \code
   blaze::setNumThreads( 4 );
   \endcode

// Please note that the \b Blaze library does not limit the available number of threads. Therefore
// it is in YOUR responsibility to choose an appropriate number of threads. The best performance,
// though, can be expected if the specified number of threads matches the available number of
// cores.
//
// In order to query the number of threads used for the parallelization of operations, the
// \c getNumThreads() function can be used:

   \code
   const size_t threads = blaze::getNumThreads();
   \endcode

// In the context of Boost threads, the function will return the previously specified number of
// threads.
//
//
// \n \section boost_threads_configuration Boost Thread Configuration
// <hr>
//
// As in case of the other shared memory parallelizations \b Blaze is not unconditionally running
// an operation in parallel (see \ref openmp_parallelization or \ref cpp_threads_parallelization).
// All thresholds related to the Boost thread parallelization are also contained within the
// configuration file <tt><blaze/config/Thresholds.h></tt>.
//
// Please note that these thresholds are highly sensitiv to the used system architecture and
// the shared memory parallelization technique. Therefore the default values cannot guarantee
// maximum performance for all possible situations and configurations. They merely provide a
// reasonable standard for the current CPU generation. Also note that the provided defaults
// have been determined using the OpenMP parallelization and require individual adaption for
// the Boost thread parallelization.
//
// \n Previous: \ref cpp_threads_parallelization &nbsp; &nbsp; Next: \ref openmp_parallelization
*/
//*************************************************************************************************


//**OpenMP Parallelization*************************************************************************
/*!\page openmp_parallelization OpenMP Parallelization
//
// \tableofcontents
//
//
// The fourth and final shared memory parallelization provided with \b Blaze is based on
// <a href="https://www.openmp.org">OpenMP</a>.
//
//
// \n \section openmp_setup OpenMP Setup
// <hr>
//
// To enable the OpenMP-based parallelization, all that needs to be done is to explicitly specify
// the use of OpenMP on the command line:

   \code
   -fopenmp   // GNU/Clang C++ compiler
   -openmp    // Intel C++ compiler
   /openmp    // Visual Studio
   \endcode

// This simple action will cause the \b Blaze library to automatically try to run all operations
// in parallel with the specified number of threads. Note however that the HPX-based, the C++11
// thread-based, and the Boost thread-based parallelizations have priority, i.e. are preferred in
// case either is enabled in combination with the OpenMP thread parallelization.
//
// As common for OpenMP, the number of threads can be specified either via an environment variable

   \code
   export OMP_NUM_THREADS=4  // Unix systems
   set OMP_NUM_THREADS=4     // Windows systems
   \endcode

// or via an explicit call to the \c omp_set_num_threads() function:

   \code
   omp_set_num_threads( 4 );
   \endcode

// Alternatively, the number of threads can also be specified via the \c setNumThreads() function
// provided by the \b Blaze library:

   \code
   blaze::setNumThreads( 4 );
   \endcode

// Please note that the \b Blaze library does not limit the available number of threads. Therefore
// it is in YOUR responsibility to choose an appropriate number of threads. The best performance,
// though, can be expected if the specified number of threads matches the available number of
// cores.
//
// In order to query the number of threads used for the parallelization of operations, the
// \c getNumThreads() function can be used:

   \code
   const size_t threads = blaze::getNumThreads();
   \endcode

// In the context of OpenMP, the function returns the maximum number of threads OpenMP will use
// within a parallel region and is therefore equivalent to the \c omp_get_max_threads() function.
//
//
// \n \section openmp_configuration OpenMP Configuration
// <hr>
//
// Note that \b Blaze is not unconditionally running an operation in parallel. In case \b Blaze
// deems the parallel execution as counterproductive for the overall performance, the operation
// is executed serially. One of the main reasons for not executing an operation in parallel is
// the size of the operands. For instance, a vector addition is only executed in parallel if the
// size of both vector operands exceeds a certain threshold. Otherwise, the performance could
// seriously decrease due to the overhead caused by the thread setup. However, in order to be
// able to adjust the \b Blaze library to a specific system, it is possible to configure these
// thresholds manually. All shared memory thresholds are contained within the configuration file
// <tt><blaze/config/Thresholds.h></tt>.
//
// Please note that these thresholds are highly sensitiv to the used system architecture and
// the shared memory parallelization technique (see also \ref cpp_threads_parallelization and
// \ref boost_threads_parallelization). Therefore the default values cannot guarantee maximum
// performance for all possible situations and configurations. They merely provide a reasonable
// standard for the current CPU generation.
//
//
// \n \section openmp_first_touch First Touch Policy
// <hr>
//
// So far the \b Blaze library does not (yet) automatically initialize dynamic memory according
// to the first touch principle. Consider for instance the following vector triad example:

   \code
   using blaze::columnVector;

   const size_t N( 1000000UL );

   blaze::DynamicVector<double,columnVector> a( N ), b( N ), c( N ), d( N );

   // Initialization of the vectors b, c, and d
   for( size_t i=0UL; i<N; ++i ) {
      b[i] = rand<double>();
      c[i] = rand<double>();
      d[i] = rand<double>();
   }

   // Performing a vector triad
   a = b + c * d;
   \endcode

// If this code, which is prototypical for many OpenMP applications that have not been optimized
// for ccNUMA architectures, is run across several locality domains (LD), it will not scale
// beyond the maximum performance achievable on a single LD if the working set does not fit into
// the cache. This is because the initialization loop is executed by a single thread, writing to
// \c b, \c c, and \c d for the first time. Hence, all memory pages belonging to those arrays will
// be mapped into a single LD.
//
// As mentioned above, this problem can be solved by performing vector initialization in parallel:

   \code
   // ...

   // Initialization of the vectors b, c, and d
   #pragma omp parallel for
   for( size_t i=0UL; i<N; ++i ) {
      b[i] = rand<double>();
      c[i] = rand<double>();
      d[i] = rand<double>();
   }

   // ...
   \endcode

// This simple modification makes a huge difference on ccNUMA in memory-bound situations (as for
// instance in all BLAS level 1 operations and partially BLAS level 2 operations). Therefore, in
// order to achieve the maximum possible performance, it is imperative to initialize the memory
// according to the later use of the data structures.
//
//
// \n \section openmp_limitations Limitations of the OpenMP Parallelization
// <hr>
//
// There are a few important limitations to the current \b Blaze OpenMP parallelization. The first
// one involves the explicit use of an OpenMP parallel region (see \ref openmp_parallel), the
// other one the OpenMP \c sections directive (see \ref openmp_sections).
//
//
// \n \subsection openmp_parallel The Parallel Directive
//
// In OpenMP threads are explicitly spawned via the an OpenMP parallel directive:

   \code
   // Serial region, executed by a single thread

   #pragma omp parallel
   {
      // Parallel region, executed by the specified number of threads
   }

   // Serial region, executed by a single thread
   \endcode

// Conceptually, the specified number of threads (see \ref openmp_setup) is created every time a
// parallel directive is encountered. Therefore, from a performance point of view, it seems to be
// beneficial to use a single OpenMP parallel directive for several operations:

   \code
   blaze::DynamicVector<double> x, y1, y2;
   blaze::DynamicMatrix<double> A, B;

   #pragma omp parallel
   {
      y1 = A * x;
      y2 = B * x;
   }
   \endcode

// Unfortunately, this optimization approach is not allowed within the \b Blaze library. More
// explicitly, it is not allowed to put an operation into a parallel region. The reason is that
// the entire code contained within a parallel region is executed by all threads. Although this
// appears to just comprise the contained computations, a computation (or more specifically the
// assignment of an expression to a vector or matrix) can contain additional logic that must not
// be handled by multiple threads (as for instance memory allocations, setup of temporaries, etc.).
// Therefore it is not possible to manually start a parallel region for several operations, but
// \b Blaze will spawn threads automatically, depending on the specifics of the operation at hand
// and the given operands.
//
// \n \subsection openmp_sections The Sections Directive
//
// OpenMP provides several work-sharing construct to distribute work among threads. One of these
// constructs is the \c sections directive:

   \code
   blaze::DynamicVector<double> x, y1, y2;
   blaze::DynamicMatrix<double> A, B;

   // ... Resizing and initialization

   #pragma omp sections
   {
   #pragma omp section

      y1 = A * x;

   #pragma omp section

      y2 = B * x;

   }
   \endcode

// In this example, two threads are used to compute two distinct matrix/vector multiplications
// concurrently. Thereby each of the \c sections is executed by exactly one thread.
//
// Unfortunately \b Blaze does not support concurrent parallel computations and therefore this
// approach does not work with any of the \b Blaze parallelization techniques. All techniques
// (including the C++11 and Boost thread parallelizations; see \ref cpp_threads_parallelization
// and \ref boost_threads_parallelization) are optimized for the parallel computation of an
// operation within a single thread of execution. This means that \b Blaze tries to use all
// available threads to compute the result of a single operation as efficiently as possible.
// Therefore, for this special case, it is advisable to disable all \b Blaze parallelizations
// and to let \b Blaze compute all operations within a \c sections directive in serial. This can
// be done by either completely disabling the \b Blaze parallelization (see \ref serial_execution)
// or by selectively serializing all operations within a \c sections directive via the \c serial()
// function:

   \code
   blaze::DynamicVector<double> x, y1, y2;
   blaze::DynamicMatrix<double> A, B;

   // ... Resizing and initialization

   #pragma omp sections
   {
   #pragma omp section

      y1 = serial( A * x );

   #pragma omp section

      y2 = serial( B * x );

   }
   \endcode

// Please note that the use of the \c BLAZE_SERIAL_SECTION (see also \ref serial_execution) does
// NOT work in this context!
//
// \n Previous: \ref boost_threads_parallelization &nbsp; &nbsp; Next: \ref serial_execution
*/
//*************************************************************************************************


//**Serial Execution*******************************************************************************
/*!\page serial_execution Serial Execution
//
// Sometimes it may be necessary to enforce the serial execution of specific operations. For this
// purpose, the \b Blaze library offers three possible options: the serialization of a single
// expression via the \c serial() function, the serialization of a block of expressions via the
// \c BLAZE_SERIAL_SECTION, and the general deactivation of the parallel execution.
//
//
// \n \section serial_execution_serial_expression Option 1: Serialization of a Single Expression
// <hr>
//
// The first option is the serialization of a specific operation via the \c serial() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = serial( A + B );
   \endcode

// \c serial() enforces the serial evaluation of the enclosed expression. It can be used on any
// kind of dense or sparse vector or matrix expression.
//
//
// \n \section serial_execution_serial_section Option 2: Serialization of Multiple Expressions
// <hr>
//
// The second option is the temporary and local enforcement of a serial execution via the
// \c BLAZE_SERIAL_SECTION:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<double,rowMajor> A;
   blaze::DynamicVector<double,columnVector> b, c, d, x, y, z;

   // ... Resizing and initialization

   // Parallel execution
   // If possible and beneficial for performance the following operation is executed in parallel.
   x = A * b;

   // Serial execution
   // All operations executed within the serial section are guaranteed to be executed in
   // serial (even if a parallel execution would be possible and/or beneficial).
   BLAZE_SERIAL_SECTION
   {
      y = A * c;
      z = A * d;
   }

   // Parallel execution continued
   // ...
   \endcode

// Within the scope of the \c BLAZE_SERIAL_SECTION, all operations are guaranteed to run in serial.
// Outside the scope of the serial section, all operations are run in parallel (if beneficial for
// the performance).
//
// Note that the \c BLAZE_SERIAL_SECTION must only be used within a single thread of execution.
// The use of the serial section within several concurrent threads will result undefined behavior!
//
//
// \n \section serial_execution_deactivate_parallelism Option 3: Deactivation of Parallel Execution
// <hr>
//
// The third option is the general deactivation of the parallel execution (even in case OpenMP is
// enabled on the command line). This can be achieved via the \c BLAZE_USE_SHARED_MEMORY_PARALLELIZATION
// switch in the <tt>./blaze/config/SMP.h</tt> configuration file:

   \code
   #define BLAZE_USE_SHARED_MEMORY_PARALLELIZATION 1
   \endcode

// In case the \c BLAZE_USE_SHARED_MEMORY_PARALLELIZATION switch is set to 0, the shared memory
// parallelization is deactivated altogether.
//
// \n Previous: \ref openmp_parallelization &nbsp; &nbsp; Next: \ref serialization
*/
//*************************************************************************************************


//**Serialization**********************************************************************************
/*!\page serialization Serialization
//
// Sometimes it is necessary to store vector and/or matrices on disk, for instance for storing
// results or for sharing specific setups with other people. The \b Blaze math serialization
// module provides the according functionality to create platform independent, portable, binary
// representations of vectors and matrices that can be used to store the \b Blaze data structures
// without loss of precision and to reliably transfer them from one machine to another.
//
// The following two pages explain how to serialize vectors and matrices:
//
//  - \ref vector_serialization
//  - \ref matrix_serialization
//
// \n Previous: \ref serial_execution &nbsp; &nbsp; Next: \ref vector_serialization
*/
//*************************************************************************************************


//**Vector Serialization***************************************************************************
/*!\page vector_serialization Vector Serialization
//
// The following example demonstrates the (de-)serialization of dense and sparse vectors:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Serialization of both vectors
   {
      blaze::StaticVector<double,5UL,rowVector> d;
      blaze::CompressedVector<int,columnVector> s;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vectors.blaze"
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

      // Serialization of both vectors into the same archive. Note that d lies before s!
      archive << d << s;
   }

   // Reconstitution of both vectors
   {
      blaze::DynamicVector<double,rowVector> d1;
      blaze::DynamicVector<int,rowVector> d2;

      // Creating an archive that reads from the file "vectors.blaze"
      blaze::Archive<std::ifstream> archive( "vectors.blaze" );

      // Reconstituting the former d vector into d1. Note that it is possible to reconstitute
      // the vector into a differrent kind of vector (StaticVector -> DynamicVector), but that
      // the type of elements has to be the same.
      archive >> d1;

      // Reconstituting the former s vector into d2. Note that is is even possible to reconstitute
      // a sparse vector as a dense vector (also the reverse is possible) and that a column vector
      // can be reconstituted as row vector (and vice versa). Note however that also in this case
      // the type of elements is the same!
      archive >> d2
   }
   \endcode

// The (de-)serialization of vectors is not restricted to vectors of built-in data type, but can
// also be used for vectors with vector or matrix element type:

   \code
   // Serialization
   {
      blaze::CompressedVector< blaze::DynamicVector< blaze::complex<double> > > vec;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vector.blaze"
      blaze::Archive<std::ofstream> archive( "vector.blaze" );

      // Serialization of the vector into the archive
      archive << vec;
   }

   // Deserialization
   {
      blaze::CompressedVector< blaze::DynamicVector< blaze::complex<double> > > vec;

      // Creating an archive that reads from the file "vector.blaze"
      blaze::Archive<std::ifstream> archive( "vector.blaze" );

      // Reconstitution of the vector from the archive
      archive >> vec;
   }
   \endcode

// As the examples demonstrates, the vector serialization offers an enormous flexibility. However,
// several actions result in errors:
//
//  - vectors cannot be reconstituted as matrices (and vice versa)
//  - the element type of the serialized and reconstituted vector must match, which means
//    that on the source and destination platform the general type (signed/unsigned integral
//    or floating point) and the size of the type must be exactly the same
//  - when reconstituting a \c StaticVector, its size must match the size of the serialized vector
//
// In case an error is encountered during (de-)serialization, a \c std::runtime_exception is
// thrown.
//
// \n Previous: \ref serialization &nbsp; &nbsp; Next: \ref matrix_serialization
*/
//*************************************************************************************************


//**Matrix Serialization***************************************************************************
/*!\page matrix_serialization Matrix Serialization
//
// The serialization of matrices works in the same manner as the serialization of vectors. The
// following example demonstrates the (de-)serialization of dense and sparse matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Serialization of both matrices
   {
      blaze::StaticMatrix<double,3UL,5UL,rowMajor> D;
      blaze::CompressedMatrix<int,columnMajor> S;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "matrices.blaze"
      blaze::Archive<std::ofstream> archive( "matrices.blaze" );

      // Serialization of both matrices into the same archive. Note that D lies before S!
      archive << D << S;
   }

   // Reconstitution of both matrices
   {
      blaze::DynamicMatrix<double,rowMajor> D1;
      blaze::DynamicMatrix<int,rowMajor> D2;

      // Creating an archive that reads from the file "matrices.blaze"
      blaze::Archive<std::ifstream> archive( "matrices.blaze" );

      // Reconstituting the former D matrix into D1. Note that it is possible to reconstitute
      // the matrix into a differrent kind of matrix (StaticMatrix -> DynamicMatrix), but that
      // the type of elements has to be the same.
      archive >> D1;

      // Reconstituting the former S matrix into D2. Note that is is even possible to reconstitute
      // a sparse matrix as a dense matrix (also the reverse is possible) and that a column-major
      // matrix can be reconstituted as row-major matrix (and vice versa). Note however that also
      // in this case the type of elements is the same!
      archive >> D2
   }
   \endcode

// Note that also in case of matrices it is possible to (de-)serialize matrices with vector or
// matrix elements:

   \code
   // Serialization
   {
      blaze::CompressedMatrix< blaze::DynamicMatrix< blaze::complex<double> > > mat;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "matrix.blaze"
      blaze::Archive<std::ofstream> archive( "matrix.blaze" );

      // Serialization of the matrix into the archive
      archive << mat;
   }

   // Deserialization
   {
      blaze::CompressedMatrix< blaze::DynamicMatrix< blaze::complex<double> > > mat;

      // Creating an archive that reads from the file "matrix.blaze"
      blaze::Archive<std::ifstream> archive( "matrix.blaze" );

      // Reconstitution of the matrix from the archive
      archive >> mat;
   }
   \endcode

// Note that just as the vector serialization, the matrix serialization is restricted by a
// few important rules:
//
//  - matrices cannot be reconstituted as vectors (and vice versa)
//  - the element type of the serialized and reconstituted matrix must match, which means
//    that on the source and destination platform the general type (signed/unsigned integral
//    or floating point) and the size of the type must be exactly the same
//  - when reconstituting a \c StaticMatrix, the number of rows and columns must match those
//    of the serialized matrix
//
// In case an error is encountered during (de-)serialization, a \c std::runtime_exception is
// thrown.
//
// \n Previous: \ref vector_serialization &nbsp; &nbsp; Next: \ref customization \n
*/
//*************************************************************************************************


//**Customization**********************************************************************************
/*!\page customization Customization
//
// Although \b Blaze tries to work out of the box for every possible setting, still it may be
// necessary to adapt the library to specific requirements. The following three pages explain
// how to customize the \b Blaze library to your own needs:
//
//  - \ref configuration_files
//  - \ref vector_and_matrix_customization
//  - \ref grouping_tagging
//  - \ref error_reporting_customization
//
// \n Previous: \ref matrix_serialization &nbsp; &nbsp; Next: \ref configuration_files
*/
//*************************************************************************************************


//**Configuration Files****************************************************************************
/*!\page configuration_files Configuration Files
//
// \tableofcontents
//
//
// Sometimes it is necessary to adapt \b Blaze to specific requirements. For this purpose
// \b Blaze provides several configuration files in the <tt>./blaze/config/</tt> subdirectory,
// which provide ample opportunity to customize internal settings, behavior, and thresholds.
// This chapter explains the most important of these configuration files. For a complete
// overview of all customization opportunities, please go to the configuration files in the
// <tt>./blaze/config/</tt> subdirectory or see the complete \b Blaze documentation.
//
//
// \n \section transpose_flag Default Vector Storage
// <hr>
//
// The \b Blaze default is that all vectors are created as column vectors (if not specified
// explicitly):

   \code
   blaze::StaticVector<double,3UL> x;  // Creates a 3-dimensional static column vector
   \endcode

// The header file <tt>./blaze/config/TransposeFlag.h</tt> allows the configuration of the default
// vector storage (i.e. the default transpose flag) of all vectors within the \b Blaze library.
// The default transpose flag is specified via the \c BLAZE_DEFAULT_TRANSPOSE_FLAG macro:

   \code
   #define BLAZE_DEFAULT_TRANSPOSE_FLAG blaze::columnVector
   \endcode

// Alternatively the default transpose flag can be specified via command line or by or defining
// this symbol manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_DEFAULT_TRANSPOSE_FLAG=blaze::columnVector ...
   \endcode

   \code
   #define BLAZE_DEFAULT_TRANSPOSE_FLAG blaze::columnVector
   #include <blaze/Blaze.h>
   \endcode

// Valid settings for \c BLAZE_DEFAULT_TRANSPOSE_FLAG are blaze::rowVector and blaze::columnVector.
//
//
// \n \section storage_order Default Matrix Storage
// <hr>
//
// Matrices are by default created as row-major matrices:

   \code
   blaze::StaticMatrix<double,3UL,3UL>  A;  // Creates a 3x3 row-major matrix
   \endcode

// The header file <tt>./blaze/config/StorageOrder.h</tt> allows the configuration of the default
// matrix storage order. Via the \c BLAZE_DEFAULT_STORAGE_ORDER macro the default storage order
// for all matrices of the \b Blaze library can be specified.

   \code
   #define BLAZE_DEFAULT_STORAGE_ORDER blaze::rowMajor
   \endcode

// Alternatively the default storage order can be specified via command line or by or defining
// this symbol manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_DEFAULT_STORAGE_ORDER=blaze::rowMajor ...
   \endcode

   \code
   #define BLAZE_DEFAULT_STORAGE_ORDER blaze::rowMajor
   #include <blaze/Blaze.h>
   \endcode

// Valid settings for \c BLAZE_DEFAULT_STORAGE_ORDER are blaze::rowMajor and blaze::columnMajor.
//
//
// \n \section blas_mode BLAS Mode
// <hr>
//
// In order to achieve maximum performance for multiplications with dense matrices, \b Blaze can
// be configured to use a BLAS library. Via the following compilation switch in the configuration
// file <tt>./blaze/config/BLAS.h</tt> BLAS can be enabled:

   \code
   #define BLAZE_BLAS_MODE 1
   \endcode

// By default, \b Blaze assumes a 32-bit BLAS library. Via the \c BLAZE_BLAS_IS_64BIT compilation
// switch, the 64-bit BLAS mode can be selected:

   \code
   #define BLAZE_BLAS_IS_64BIT 1
   \endcode

// Note that the \c BLAZE_BLAS_IS_64BIT switch also has an effect on the \ref lapack_functions.
// Please also note that it might additionally be necessary to use a compilation switch to put
// the BLAS/LAPACK library into 64-bit mode (e.g. \c MKL_ILP64 for the Intel MKL library).
//
// In case the selected BLAS library provides parallel execution, the \c BLAZE_BLAS_IS_PARALLEL
// switch should be activated to prevent \b Blaze from parallelizing on its own:

   \code
   #define BLAZE_BLAS_IS_PARALLEL 1
   \endcode

// Additionally, it is possible to specify the name of the BLAS include file via the
// \c BLAZE_BLAS_INCLUDE_FILE switch. The default setting is <tt><cblas.h></tt>:

   \code
   #define BLAZE_BLAS_INCLUDE_FILE <cblas.h>
   \endcode

// Alternatively, all settings can be specified via command line or by or defining the symbols
// manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_BLAS_MODE=1 -DBLAZE_BLAS_IS_64BIT=1 -DBLAZE_BLAS_IS_PARALLEL=1 -DBLAZE_BLAS_INCLUDE_FILE='<cblas.h>' ...
   \endcode

   \code
   #define BLAZE_BLAS_MODE 1
   #define BLAZE_BLAS_IS_64BIT 1
   #define BLAZE_BLAS_IS_PARALLEL 1
   #define BLAZE_BLAS_INCLUDE_FILE <cblas.h>
   #include <blaze/Blaze.h>
   \endcode

// In case no BLAS library is available, \b Blaze will still work and will not be reduced in
// functionality, but performance may be limited.
//
//
// \n \section cache_size Cache Size
// <hr>
//
// The optimization of several \b Blaze compute kernels depends on the cache size of the target
// architecture. By default, \b Blaze assumes a cache size of 3 MiByte. However, for optimal
// speed the exact cache size of the system should be provided via the \c cacheSize value in the
// <tt>./blaze/config/CacheSize.h</tt> configuration file:

   \code
   #define BLAZE_CACHE_SIZE 3145728UL;
   \endcode

// The cache size can also be specified via command line or by defining this symbol manually
// before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_CACHE_SIZE=3145728 ...
   \endcode

   \code
   #define BLAZE_CACHE_SIZE 3145728
   #include <blaze/Blaze.h>
   \endcode

// \n \section vectorization Vectorization
// <hr>
//
// In order to achieve maximum performance and to exploit the compute power of a target platform
// the \b Blaze library attempts to vectorize all linear algebra operations by SSE, AVX, and/or
// AVX-512 intrinsics, depending on which instruction set is available. However, it is possible
// to disable the vectorization entirely by the compile time switch in the configuration file
// <tt>./blaze/config/Vectorization.h</tt>:

   \code
   #define BLAZE_USE_VECTORIZATION 1
   \endcode

// It is also possible to (de-)activate vectorization via command line or by defining this symbol
// manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_USE_VECTORIZATION=1 ...
   \endcode

   \code
   #define BLAZE_USE_VECTORIZATION 1
   #include <blaze/Blaze.h>
   \endcode

// In case the switch is set to 1, vectorization is enabled and the \b Blaze library is allowed
// to use intrinsics to speed up computations. In case the switch is set to 0, vectorization is
// disabled entirely and the \b Blaze library chooses default, non-vectorized functionality for
// the operations. Note that deactivating the vectorization may pose a severe performance
// limitation for a large number of operations!
//
//
// \n \section sleef Sleef
// <hr>
//
// For several complex operations \b Blaze can make use of the Sleef library for vectorization
// (https://github.com/shibatch/sleef). This compilation switch enables/disables the vectorization
// by means of Sleef. In case the switch is set to 1, \b Blaze uses Sleef for instance for the
// vectorized computation of trigonometric functions (i.e. \c sin(), \c cos(), \c tan(), etc.)
// and exponential functions (i.e. \c exp(), \c log(), ...).

   \code
   #define BLAZE_USE_SLEEF 1
   \endcode

// It is also possible to enable/disable Sleef vectorization via command line or by defining this
// symbol manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_USE_SLEEF=1 ...
   \endcode

   \code
   #define BLAZE_USE_SLEEF 1
   #include <blaze/Blaze.h>
   \endcode

// \n \section thresholds Thresholds
// <hr>
//
// For many computations \b Blaze distinguishes between small and large vectors and matrices.
// This separation is especially important for the parallel execution of computations, since
// the use of several threads only pays off for sufficiently large vectors and matrices.
// Additionally, it also enables \b Blaze to select kernels that are optimized for a specific
// size.
//
// In order to distinguish between small and large data structures \b Blaze provides several
// thresholds that can be adapted to the characteristics of the target platform. For instance,
// the \c DMATDVECMULT_THRESHOLD specifies the threshold between the application of the custom
// \b Blaze kernels for small dense matrix/dense vector multiplications and the BLAS kernels
// for large multiplications. All thresholds, including the thresholds for the OpenMP- and
// thread-based parallelization, are contained within the configuration file
// <tt><blaze/config/Thresholds.h></tt>.
//
//
// \n \section alignment Alignment
// <hr>
//
// For performance reasons, the vector types \ref vector_types_static_vector and
// \ref vector_types_hybrid_vector and the matrix types \ref matrix_types_static_matrix and
// \ref matrix_types_hybrid_matrix by default make use of aligned memory. Via the configuration
// file <tt>./blaze/config/Alignment.h</tt> it is possible to define the default alignment flag:

   \code
   #define BLAZE_DEFAULT_ALIGNMENT_FLAG blaze::aligned
   \endcode

// Alternatively it is possible set the default alignment flag via command line or by defining
// this symbol manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_DEFAULT_ALIGNMENT_FLAG=blaze::aligned ...
   \endcode

   \code
   #define BLAZE_DEFAULT_ALIGNMENT_FLAG blaze::aligned
   #include <blaze/Blaze.h>
   \endcode

// If \c BLAZE_DEFAULT_ALIGNMENT_FLAG is set to \c blaze::aligned then \ref vector_types_static_vector,
// \ref vector_types_hybrid_vector, \ref matrix_types_static_matrix, and \ref matrix_types_hybrid_matrix
// use aligned memory by default. If it is set to \c blaze::unaligned they don't enforce aligned
// memory. Note however that disabling alignment can considerably reduce the performance of all
// operations with these vector and matrix types!
//
//
// \n \section padding Padding
// <hr>
//
// By default the \b Blaze library uses padding for the vector types \ref vector_types_static_vector
// and \ref vector_types_hybrid_vector and the matrix types \ref matrix_types_static_matrix and
// \ref matrix_types_hybrid_matrix in order to achieve maximum performance in all operations. Due
// to padding, the proper alignment of data elements can be guaranteed and the need for remainder
// loops is minimized. However, on the downside padding introduces an additional memory overhead,
// which can be large depending on the used data type.
//
// The configuration file <tt>./blaze/config/Padding.h</tt> provides a compile time switch that
// can be used to define the default padding flag:

   \code
   #define BLAZE_DEFAULT_PADDING_FLAG blaze::padded
   \endcode

// Alternatively it is possible to define the default padding flag via command line or by defining
// this symbol manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_DEFAULT_PADDING_FLAG=blaze::padded ...
   \endcode

   \code
   #define BLAZE_DEFAULT_PADDING_FLAG blaze::padded
   #include <blaze/Blaze.h>
   \endcode

// If \c BLAZE_DEFAULT_ALIGNMENT_FLAG is set to \c blaze::padded, by default padding is enabled
// for \ref vector_types_static_vector, \ref vector_types_hybrid_vector, \ref matrix_types_static_matrix
// and \ref matrix_types_hybrid_matrix. If it is set to \c blaze::unpadded, then padding is by
// default disabled. Note however that disabling padding can considerably reduce the performance
// of all dense vector and matrix operations!
//
//
// \n \section streaming Streaming (Non-Temporal Stores)
// <hr>
//
// For vectors and matrices that don't fit into the cache anymore non-temporal stores can provide
// a significant performance advantage of about 20%. However, this advantage is only in effect in
// case the memory bandwidth of the target architecture is maxed out. If the target architecture's
// memory bandwidth cannot be exhausted the use of non-temporal stores can decrease performance
// instead of increasing it.
//
// The configuration file <tt>./blaze/config/Optimizations.h</tt> provides a compile time switch
// that can be used to (de-)activate streaming:

   \code
   #define BLAZE_USE_STREAMING 1
   \endcode

// Alternatively streaming can be (de-)activated via command line or by defining this symbol
// manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_USE_STREAMING=1 ...
   \endcode

   \code
   #define BLAZE_USE_STREAMING 1
   #include <blaze/Blaze.h>
   \endcode

// If \c BLAZE_USE_STREAMING is set to 1 streaming is enabled, if it is set to 0 streaming is
// disabled. It is recommended to consult the target architecture's white papers to decide whether
// streaming is beneficial or hurtful for performance.
//
//
// \n Previous: \ref customization &nbsp; &nbsp; Next: \ref vector_and_matrix_customization \n
*/
//*************************************************************************************************


//**Customization of Vectors and Matrices**********************************************************
/*!\page vector_and_matrix_customization Customization of Vectors and Matrices
//
// \tableofcontents
//
//
// \n \section custom_data_members Custom Data Members
// <hr>
//
// So far the \b Blaze library does not provide a lot of flexibility to customize the data
// members of existing \ref vector_types and \ref matrix_types. However, to some extend it is
// possible to customize vectors and matrices by inheritance. The following example gives an
// impression on how to create a simple variation of \ref matrix_types_custom_matrix, which
// automatically takes care of acquiring and releasing custom memory.

   \code
   template< typename Type                    // Data type of the matrix
           , bool SO = defaultStorageOrder >  // Storage order
   class MyCustomMatrix
      : public CustomMatrix< Type, unaligned, unpadded, SO >
   {
    public:
      explicit inline MyCustomMatrix( size_t m, size_t n )
         : CustomMatrix<Type,unaligned,unpadded,SO>()
         , array_( new Type[m*n] )
      {
         this->reset( array_.get(), m, n );
      }

    private:
      std::unique_ptr<Type[]> array_;
   };
   \endcode

// Please note that this is a simplified example with the intent to show the general approach.
// The number of constructors, the memory acquisition, and the kind of memory management can of
// course be adapted to specific requirements. Also, please note that since none of the \b Blaze
// vectors and matrices have virtual destructors polymorphic destruction cannot be used.
//
//
// \n \section custom_operations Custom Operations
// <hr>
//
// There are two approaches to extend \b Blaze with custom operations. First, the \c map()
// functions provide the possibility to execute componentwise custom operations on vectors and
// matrices. Second, it is possible to add customized free functions.
//
// \n \subsection custom_operations_map The map() Functions
//
// Via the unary and binary \c map() functions it is possible to execute componentwise custom
// operations on vectors and matrices. The unary \c map() function can be used to apply a custom
// operation on each single element of a dense vector or matrix or each non-zero element of a
// sparse vector or matrix. For instance, the following example demonstrates a custom square
// root computation on a dense matrix:

   \code
   blaze::DynamicMatrix<double> A, B;

   B = map( A, []( double d ) { return std::sqrt( d ); } );
   \endcode

// The binary \c map() function can be used to apply an operation pairwise to the elements of
// two dense vectors or two dense matrices. The following example demonstrates the merging of
// two matrices of double precision values into a matrix of double precision complex numbers:

   \code
   blaze::DynamicMatrix<double> real{ { 2.1, -4.2 }, { 1.0,  0.6 } };
   blaze::DynamicMatrix<double> imag{ { 0.3,  1.4 }, { 2.9, -3.4 } };

   blaze::DynamicMatrix< complex<double> > cplx;

   // Creating the matrix
   //    ( ( 2.1,  0.3) (-4.2,  1.4) )
   //    ( ( 1.0,  2.9) ( 0.6, -3.4) )
   cplx = map( real, imag, []( double r, double i ){ return complex<double>( r, i ); } );
   \endcode

// These examples demonstrate the most convenient way of defining a unary custom operation by
// passing a lambda to the \c map() function. Alternatively, it is possible to pass a custom
// functor:

   \code
   struct Sqrt
   {
      double operator()( double a ) const
      {
         return std::sqrt( a );
      }
   };

   B = map( A, Sqrt() );
   \endcode

// In order for the functor to work in a call to \c map() it must define a function call operator,
// which accepts arguments of the type of the according vector or matrix elements.
//
// Although the operation is automatically parallelized depending on the size of the vector or
// matrix, no automatic vectorization is possible. In order to enable vectorization, a \c load()
// function can be added to the functor, which handles the vectorized computation. Depending on
// the data type this function is passed one of the following \b Blaze SIMD data types:
//
// <ul>
//    <li>SIMD data types for fundamental data types
//       <ul>
//          <li>\c blaze::SIMDint8: Packed SIMD type for 8-bit signed integral data types</li>
//          <li>\c blaze::SIMDuint8: Packed SIMD type for 8-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDint16: Packed SIMD type for 16-bit signed integral data types</li>
//          <li>\c blaze::SIMDuint16: Packed SIMD type for 16-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDint32: Packed SIMD type for 32-bit signed integral data types</li>
//          <li>\c blaze::SIMDuint32: Packed SIMD type for 32-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDint64: Packed SIMD type for 64-bit signed integral data types</li>
//          <li>\c blaze::SIMDuint64: Packed SIMD type for 64-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDfloat: Packed SIMD type for single precision floating point data</li>
//          <li>\c blaze::SIMDdouble: Packed SIMD type for double precision floating point data</li>
//       </ul>
//    </li>
//    <li>SIMD data types for complex data types
//       <ul>
//          <li>\c blaze::SIMDcint8: Packed SIMD type for complex 8-bit signed integral data types</li>
//          <li>\c blaze::SIMDcuint8: Packed SIMD type for complex 8-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDcint16: Packed SIMD type for complex 16-bit signed integral data types</li>
//          <li>\c blaze::SIMDcuint16: Packed SIMD type for complex 16-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDcint32: Packed SIMD type for complex 32-bit signed integral data types</li>
//          <li>\c blaze::SIMDcuint32: Packed SIMD type for complex 32-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDcint64: Packed SIMD type for complex 64-bit signed integral data types</li>
//          <li>\c blaze::SIMDcuint64: Packed SIMD type for complex 64-bit unsigned integral data types</li>
//          <li>\c blaze::SIMDcfloat: Packed SIMD type for complex single precision floating point data</li>
//          <li>\c blaze::SIMDcdouble: Packed SIMD type for complex double precision floating point data</li>
//       </ul>
//    </li>
// </ul>
//
// All SIMD types provide the \c value data member for a direct access to the underlying intrinsic
// data element. In the following example, this intrinsic element is passed to the AVX function
// \c _mm256_sqrt_pd():

   \code
   struct Sqrt
   {
      double operator()( double a ) const
      {
         return std::sqrt( a );
      }

      SIMDdouble load( const SIMDdouble& a ) const
      {
         return _mm256_sqrt_pd( a.value );
      }
   };
   \endcode

// In this example, whenever vectorization is generally applicable, the \c load() function is
// called instead of the function call operator for as long as the number of remaining elements
// is larger-or-equal to the width of the packed SIMD type. In all other cases (which also
// includes peel-off and remainder loops) the scalar operation is used.
//
// Please note that this example has two drawbacks: First, it will only compile in case the
// intrinsic \c _mm256_sqrt_pd() function is available (i.e. when AVX is active). Second, the
// availability of AVX is not taken into account. The first drawback can be alleviated by making
// the \c load() function a function template. The second drawback can be dealt with by adding a
// \c simdEnabled() function template to the functor:

   \code
   struct Sqrt
   {
      double operator()( double a ) const
      {
         return std::sqrt( a );
      }

      template< typename T >
      T load( const T& a ) const
      {
         return _mm256_sqrt_pd( a.value );
      }

      template< typename T >
      static constexpr bool simdEnabled() {
#if defined(__AVX__)
         return true;
#else
         return false;
#endif
      }
   };
   \endcode

// The \c simdEnabled() function must be a \c static, \c constexpr function and must return whether
// or not vectorization is available for the given data type \c T. In case the function returns
// \c true, the \c load() function is used for a vectorized evaluation, in case the function
// returns \c false, \c load() is neither called nor instantiated.
//
// By default the \c map() function uses peel-off and remainder loops if the number of elements is
// not a multiple of the width of the packed SIMD type. However, all dense vector and matrix types
// in \b Blaze provide padding as an optimization. In case the custom operation preserves the
// value zero of the padding elements, it is possible to omit the peel-off and remainder loops,
// include the padding elements in the computation and by that increase performance. For that
// purpose the \c paddingEnabled() function can be added to the functor:

   \code
   struct Sqrt
   {
      // ...

      static constexpr bool paddingEnabled() { return true; }
   };
   \endcode

// Also the \c paddingEnabled() function must be a \c static, \c constexpr function and must
// return whether padding elements can be used in the custom operation. In case the function
// returns \c true, the padding elements are used during a vectorized operation, in case the
// function returns \c false, the padding elements are not used.
//
// Note that this is a simplified example that is only working when used for dense vectors and
// matrices with double precision floating point elements. The following code shows the complete
// implementation of the according functor that is used within the \b Blaze library. The \b Blaze
// \c Sqrt functor is working for all data types that are providing a square root operation:

   \code
   namespace blaze {

   struct Sqrt
   {
      template< typename T >
      BLAZE_ALWAYS_INLINE auto operator()( const T& a ) const
      {
         return sqrt( a );
      }

      template< typename T >
      static constexpr bool simdEnabled() { return HasSIMDSqrt<T>::value; }

      static constexpr bool paddingEnabled() { return true; }

      template< typename T >
      BLAZE_ALWAYS_INLINE auto load( const T& a ) const
      {
         BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T );
         return sqrt( a );
      }
   };

   } // namespace blaze
   \endcode

// The same approach can be taken for binary custom operations. The following code demonstrates
// the \c Min functor of the \b Blaze library, which is working for all data types that provide
// a \c min() operation:

   \code
   struct Min
   {
      explicit inline Min()
      {}

      template< typename T1, typename T2 >
      BLAZE_ALWAYS_INLINE decltype(auto) operator()( const T1& a, const T2& b ) const
      {
         return min( a, b );
      }

      template< typename T1, typename T2 >
      static constexpr bool simdEnabled() { return HasSIMDMin<T1,T2>::value; }

      static constexpr bool paddingEnabled() { return true; }

      template< typename T1, typename T2 >
      BLAZE_ALWAYS_INLINE decltype(auto) load( const T1& a, const T2& b ) const
      {
         BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T1 );
         BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T2 );
         return min( a, b );
      }
   };
   \endcode

// For more information on the available \b Blaze SIMD data types and functions, please see the
// SIMD module in the complete \b Blaze documentation.
//
// \n \subsection custom_operations_free_functions Free Functions
//
// In order to extend \b Blaze with new functionality it is possible to add free functions. Free
// functions can be used either as wrappers around calls to the map() function or to implement
// general, non-componentwise operations. The following two examples will demonstrate both ideas.
//
// The first example shows the \c setToZero() function, which resets a sparse matrix to zero
// without affecting the sparsity pattern. It is implemented as a convenience wrapper around
// the map() function:

   \code
   template< typename MT  // Type of the sparse matrix
           , bool SO >    // Storage order
   void setToZero( blaze::SparseMatrix<MT,SO>& mat )
   {
      (~mat) = blaze::map( ~mat, []( const auto& value ){ return decltype(value){}; } );
   }
   \endcode

// The blaze::SparseMatrix class template is the base class for all kinds of sparse matrices and
// provides an abstraction from the actual type \c MT of the sparse matrix. However, due to the
// <a href="https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern">Curiously Recurring Template Pattern (CRTP)</a>
// it also enables a conversion back to the actual type. This downcast is performed via the tilde
// operator (i.e. \c operator~()). The template parameter \c SO represents the storage order
// (blaze::rowMajor or blaze::columnMajor) of the matrix.
//
// The second example shows the \c countZeros() function, which counts the number of values, which
// are exactly zero, in a dense, row-major matrix:

   \code
   template< typename MT >
   size_t countZeros( blaze::DenseMatrix<MT,rowMajor>& mat )
   {
      const size_t M( (~mat).rows() );
      const size_t N( (~mat).columns() );
      size_t count( 0UL );

      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=0UL; j<N; ++j ) {
            if( blaze::isDefault<strict>( (~mat)(i,j) ) )
               ++count;
         }
      }

      return count;
   }
   \endcode

// The blaze::DenseMatrix class template is the base class for all kinds of dense matrices. Again,
// it is possible to perform the conversion to the actual type via the tilde operator.
//
// The following two listings show the declarations of all vector and matrix base classes, which
// can be used for custom free functions:

   \code
   template< typename VT  // Concrete type of the dense or sparse vector
           , bool TF >    // Transpose flag (blaze::columnVector or blaze::rowVector)
   class Vector;

   template< typename VT  // Concrete type of the dense vector
           , bool TF >    // Transpose flag (blaze::columnVector or blaze::rowVector)
   class DenseVector;

   template< typename VT  // Concrete type of the sparse vector
           , bool TF >    // Transpose flag (blaze::columnVector or blaze::rowVector)
   class SparseVector;
   \endcode

   \code
   template< typename MT  // Concrete type of the dense or sparse matrix
           , bool SO >    // Storage order (blaze::rowMajor or blaze::columnMajor)
   class Matrix;

   template< typename MT  // Concrete type of the dense matrix
           , bool SO >    // Storage order (blaze::rowMajor or blaze::columnMajor)
   class DenseMatrix;

   template< typename MT  // Concrete type of the sparse matrix
           , bool SO >    // Storage order (blaze::rowMajor or blaze::columnMajor)
   class SparseMatrix;
   \endcode

// \n \section custom_data_types Custom Data Types
// <hr>
//
// \subsection custom_data_types_introduction Introduction
//
// The \b Blaze library is not restricted to integral, floating point and complex data types
// (called numeric types in \b Blaze), but it supports custom data types. For instance, the
// following example demonstrates that it is possible to use \c std::string as data type:

   \code
   blaze::DynamicVector<std::string> a{ "Hello, ", "Blaze " , "Expression" };
   blaze::DynamicVector<std::string> b{ "World"  , "Library", " Templates" };

   const auto c( evaluate( a + b ) );
   std::cout <<  "c =\n" << c << "\n\n";

   const std::string maxString( max( c ) );
   std::cout << "maxString = " << std::quoted(maxString) << "\n";
   \endcode

// Output:

   \code
   c =
   ( Hello, World )
   ( Blaze Library )
   ( Expression Templates )

   maxString = "Hello, World"
   \endcode

// \b Blaze tries hard to make the use of custom data types as convenient, easy and intuitive as
// possible. In order to work flawlessly with \b Blaze, custom data types are required to provide
// a certain interface (depending on the operations that the type is used for). The following
// sections give an overview of the necessary steps to enable the use of the hypothetical custom
// data type \c custom::double_t for vector and matrix operations.

   \code
   namespace custom {

   struct double_t
   {
      constexpr double_t() = default;
      constexpr double_t( double i ) : value( i ) {}
      double value{};
   };

   } // namespace custom
   \endcode

// \subsection custom_data_types_arithmetic_operations Arithmetic Operations
//
// The \b Blaze library assumes that a custom data type provides \c operator<<() for streaming,
// \c operator+=() and \c operator+() for additions (which for instance includes additions inside
// matrix/vector multiplications, matrix/matrix multiplications, reduction or norm operations),
// \c operator-=() and \c operator-() for subtractions, \c operator*=() and \c operator*() for
// multiplications and \c operator/=() and \c operator/() for divisions:

   \code
   namespace custom {

   constexpr double_t& operator+=( double_t& lhs, double_t rhs ) noexcept { lhs.value += rhs.value; return lhs; }
   constexpr double_t& operator-=( double_t& lhs, double_t rhs ) noexcept { lhs.value -= rhs.value; return lhs; }
   constexpr double_t& operator*=( double_t& lhs, double_t rhs ) noexcept { lhs.value *= rhs.value; return lhs; }
   constexpr double_t& operator/=( double_t& lhs, double_t rhs ) noexcept { lhs.value /= rhs.value; return lhs; }

   constexpr double_t operator+( double_t lhs, double_t rhs ) noexcept { return double_t{ lhs.value + rhs.value }; }
   constexpr double_t operator-( double_t lhs, double_t rhs ) noexcept { return double_t{ lhs.value - rhs.value }; }
   constexpr double_t operator*( double_t lhs, double_t rhs ) noexcept { return double_t{ lhs.value * rhs.value }; }
   constexpr double_t operator/( double_t lhs, double_t rhs ) noexcept { return double_t{ lhs.value / rhs.value }; }

   inline std::ostream& operator<<( std::ostream& os, double_t d )
   {
      return os << d.value;
   }

   } // namespace custom
   \endcode

// Example:

   \code
   int main()
   {
      blaze::DynamicVector<custom::double_t> a{ 1.0, 2.0, 3.0, 4.0 };
      blaze::DynamicVector<custom::double_t> b{ 0.1, 0.2, 0.3, 0.4 };

      std::cout << "a + b =\n" << ( a + b ) << "\n";
      std::cout << "a * b =\n" << ( a * b ) << "\n";

      std::cout << "sum(a) = " << sum(a) << "\n"
                << "prod(a) = " << prod(a) << "\n";
   }
   \endcode

// Output:

   \code
   a + b =
   (         1.1 )
   (         2.2 )
   (         3.3 )
   (         4.4 )

   a * b =
   (         0.1 )
   (         0.4 )
   (         0.9 )
   (         1.6 )

   sum(a) = 10
   prod(a) = 24
   \endcode

// Note that similar steps are necessary if several custom data types are combined (as for instance
// \c custom::double_t and \c custom::float_t). Note that in this case both permutations need to
// be taken into account:

   \code
   custom::double_t operator+( const custom::double_t& a, const custom::float_t& b );
   custom::double_t operator+( const custom::float_t& a, const custom::double_t& b );
   // ...
   \endcode

// Please note that only built-in data types apply for vectorization and thus custom data types
// cannot achieve maximum performance!
//
// \subsection custom_data_types_relational_operations Relational Operations
//
// In order to compare the element type, \b Blaze expects the equality operator (i.e. \c operator==())
// and the inequality operator (i.e. \c operator!=()). Alternatively it is possible to provide an
// \c equal() function, which distinguishes between strict and relaxed comparison:

   \code
   namespace custom {

   constexpr bool operator==( double_t lhs, double_t rhs ) noexcept { return lhs.value == rhs.value; }
   constexpr bool operator!=( double_t lhs, double_t rhs ) noexcept { return !( lhs == rhs ); }

   template< blaze::RelaxationFlag RF >
   constexpr bool equal( double_t lhs, double_t rhs ) noexcept { return blaze::equal<RF>( lhs.value, rhs.value ); }

   } // namespace custom
   \endcode

// Example:

   \code
   int main()
   {
      blaze::DynamicVector<custom::double_t> a{ 1.0, 2.0, 3.0, 4.0 };
      blaze::DynamicVector<custom::double_t> b{ 0.1, 0.2, 0.3, 0.4 };

      std::cout << "a == b: " << ( a == b ) << "\n"
                << "a != b: " << ( a != b ) << "\n";
   }
   \endcode

// Output:

   \code
   a == b: 0
   a != b: 1
   \endcode

// \subsection custom_data_types_elementwise_operations Elementwise Operations
//
// For the different kinds of elementwise operations on vectors and matrices (\c abs(), \c sin(),
// \c cos(), \c sqrt(), \c log(), \c exp(), \c min(), \c max(), ...), the custom type is required
// to provide the according function overload. Note that the \c sqrt() operation may also be
// required for several norm computations. Also, for any inversion operation, the type is required
// to suport the \c inv() function:

   \code
   namespace custom {

   inline    double_t abs ( double_t d ) noexcept { return double_t{ std::abs ( d.value ) }; }
   inline    double_t sin ( double_t d ) noexcept { return double_t{ std::sin ( d.value ) }; }
   inline    double_t cos ( double_t d ) noexcept { return double_t{ std::cos ( d.value ) }; }
   inline    double_t sqrt( double_t d ) noexcept { return double_t{ std::sqrt( d.value ) }; }
   inline    double_t log ( double_t d ) noexcept { return double_t{ std::log ( d.value ) }; }
   inline    double_t exp ( double_t d ) noexcept { return double_t{ std::exp ( d.value ) }; }
   constexpr double_t inv ( double_t d ) noexcept { return double_t{ 1.0/d.value }; }

   constexpr double_t min( double_t lhs, double_t rhs ) noexcept { return double_t{ blaze::min( lhs.value, rhs.value ) }; }
   constexpr double_t max( double_t lhs, double_t rhs ) noexcept { return double_t{ blaze::max( lhs.value, rhs.value ) }; }

   } // namespace custom
   \endcode

// Example:

   \code
   int main()
   {
      blaze::DynamicVector<custom::double_t> a{ 1.0, 2.0, 3.0, 4.0 };
      blaze::DynamicVector<custom::double_t> b{ 0.1, 0.2, 0.3, 0.4 };

      std::cout << "abs(a) =\n" << abs(a) << "\n";
      std::cout << "sin(a) =\n" << sin(a) << "\n";
      std::cout << "cos(a) =\n" << cos(a) << "\n";
      std::cout << "sqrt(a) =\n" << sqrt(a) << "\n";
      std::cout << "log(a) =\n" << log(a) << "\n";
      std::cout << "exp(a) =\n" << exp(a) << "\n\n";
      std::cout << "min(a) =\n" << min(a) <<  "\n";
      std::cout << "max(a) =\n" << max(a) << "\n\n";
      std::cout << "min(a,b) =\n" << min(a,b) << "\n";
      std::cout << "max(a,b) =\n" << max(a,b) << "\n";
      std::cout << "norm(a) = " << norm(a) << "\n";
   }
   \endcode

// Output:

   \code
   abs(a) =
   (           1 )
   (           2 )
   (           3 )
   (           4 )

   sin(a) =
   (    0.841471 )
   (    0.909297 )
   (     0.14112 )
   (   -0.756802 )

   cos(a) =
   (    0.540302 )
   (   -0.416147 )
   (   -0.989992 )
   (   -0.653644 )

   sqrt(a) =
   (           1 )
   (     1.41421 )
   (     1.73205 )
   (           2 )

   log(a) =
   (           0 )
   (    0.693147 )
   (     1.09861 )
   (     1.38629 )

   exp(a) =
   (     2.71828 )
   (     7.38906 )
   (     20.0855 )
   (     54.5982 )

   min(a) = 1
   max(a) = 4

   min(a,b) =
   (         0.1 )
   (         0.2 )
   (         0.3 )
   (         0.4 )

   max(a,b) =
   (           1 )
   (           2 )
   (           3 )
   (           4 )

   norm(a) = 5.47723
   \endcode

// \subsection custom_data_types_adaptors Adaptors
//
// If the custom data type is used in the context of the HermitianMatrix, UniLowerMatrix, or
// UniUpperMatrix adaptors, it will be necessary to provide overloads of the \c isZero(),
// \c isOne(), and \c isReal() functions:

   \code
   namespace custom {

   template< blaze::RelaxationFlag RF >
   constexpr bool isZero( double_t d ) { return blaze::isZero<RF>( d.value ); }

   template< blaze::RelaxationFlag RF >
   constexpr bool isOne ( double_t d ) { return blaze::isOne<RF> ( d.value ); }

   template< blaze::RelaxationFlag RF >
   constexpr bool isReal( double_t d ) { MAYBE_UNUSED( d ); return true; }

   } // namespace custom
   \endcode

// Example:

   \code
   int main()
   {
      blaze::UniLowerMatrix< blaze::DynamicMatrix<custom::double_t> > L
         { { 1.0, 0.0, 0.0 },
           { 2.0, 1.0, 0.0 },
           { 3.0, 4.0, 1.0 } };

      blaze::UniUpperMatrix< blaze::DynamicMatrix<custom::double_t> > U
         { { 1.0, 2.0, 3.0 },
           { 0.0, 1.0, 4.0 },
           { 0.0, 0.0, 1.0 } };

      const auto A( evaluate( L * U ) );
      std::cout << "A =\n" << A << "\n";
   }
   \endcode

// Output:

   \code
   A =
   (            1            2            3 )
   (            2            5           10 )
   (            3           10           26 )
   \endcode

// \n Previous: \ref configuration_files &nbsp; &nbsp; Next: \ref grouping_tagging \n
*/
//*************************************************************************************************


//**Grouping/Tagging*******************************************************************************
/*!\page grouping_tagging Grouping/Tagging
//
// \tableofcontents
//
//
// \n \section grouping_tagging_tagging_and_groups Tagging and Groups
// <hr>
//
// Sometimes it may be desirable to separate two or more distinct groups of vectors and matrices,
// for instance in order to allow operations only within a group and to prevent operations across
// groups. This goal can be achieved by means of tags. All vector and matrix classes provide a
// template parameter to specify a tag (for instance, the fourth template parameter for
// blaze::DynamicVector and the sixth template parameter for blaze::StaticVector):

   \code
   template< typename Type, bool TF, typename Alloc, typename Tag >
   class DynamicVector;

   template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
   class StaticVector;
   \endcode

// By default, all vectors and matrices are associated with blaze::Group0 (i.e. the tag is set
// to blaze::Group0). However, it is possible to explicitly associate vectors and matrices with
// different groups:

   \code
   using blaze::DynamicVector;
   using blaze::AlignedAllocator;
   using blaze::Group0;
   using blaze::Group1;
   using blaze::columnVector;

   DynamicVector<int,columnVector,AlignedAllocator<int>,Group0> a0, b0;
   DynamicVector<int,columnVector,AlignedAllocator<int>,Group1> a1, b1;

   a0 + b0;  // Compiles, a0 and b0 are in the same group (Group0)
   a1 + b1;  // Compiles, a1 and b1 are in the same group (Group1)
   a0 + b1;  // Compilation error: a0 and b1 are not in the same group
   \endcode

// All vectors or matrices that are associated with the same group can be freely combined with any
// other vector or matrix from the same group. The attempt to combine vectors and matrices from
// different groups results in a compilation error.
//
//
// \n \section grouping_tagging_creating_new_groups Creating New Groups
// <hr>
//
// \b Blaze provides the tags for the ten predefined groups blaze::Group0 through blaze::Group9.
// In order to create further groups, all that needs to be done is to create new instances of the
// blaze::GroupTag class template:

   \code
   using Group10 = blaze::GroupTag<10>;
   using Group11 = blaze::GroupTag<11>;
   // ... further groups
   \endcode

// All groups based on the blaze::GroupTag class template will be treated as separate groups just
// as the ten predefined groups.
//
//
// \n \section grouping_tagging_custom_tags Custom Tags
// <hr>
//
// Sometimes it is not enough to separate vectors and matrices into different groups, but it is
// required to define the interaction between different groups. This situation for instance occurs
// if a vector or matrix is associated with a physical quantity. This problem can be solved by
// using custom tags. The following example gives an impression on how to define the physics on
// meters (represented by the \c Meter tag) and seconds (represented by the \c Second tag):

   \code
   struct Meter {};           // Definition of the 'Meter' tag
   struct Second {};          // Definition of the 'Second' tag
   struct SquareMeter {};     // Definition of the 'SquareMeter' tag
   struct MeterPerSecond {};  // Definition of the 'MeterPerSecond' tag
   \endcode

// The \c Meter and \c Second tags are not associated with the blaze::GroupTag class template. For
// that reason, by default, it is not possible to perform any operation on an accordingly tagged
// vector or matrix. All required operations need to be declared explicitly in order to specify
// the resulting tag of an operation. In the following code example, this happens by declaring
// both the addition for the \c Meter tag and the \c Second tag, the multiplication between two
// \c Meter tags and the division between \c Meter and \c Second. Note that it is enough to
// declare the operations, it is not necessary to define them!

   \code
   Meter          operator+( Meter , Meter  );  // Enabling addition between 'Meter'
   Second         operator+( Second, Second );  // Enabling addition between 'Second'
   SquareMeter    operator*( Meter , Meter  );  // Enabling multiplication between 'Meter'
   MeterPerSecond operator/( Meter , Second );  // Enabling division between 'Meter' and 'Second'
   \endcode

// With these declarations it is now possible to add meters and seconds, but not to subtract them
// (no subtraction operator was declared). Also, it is possible to multiply meters and to divide
// meters and seconds:

   \code
   const DynamicVector<int,rowVector,AlignedAllocator<int>,Meter> m1{ 1, 2, 3 };
   const DynamicVector<int,rowVector,AlignedAllocator<int>,Meter> m2{ 4, 5, 6 };

   const DynamicVector<int,rowVector,AlignedAllocator<int>,Second> s1{ 1, 2, 3 };
   const DynamicVector<int,rowVector,AlignedAllocator<int>,Second> s2{ 4, 5, 6 };

   m1 + m2;  // Compiles and results in vector tagged with 'Meter'
   s1 + s2;  // Compiles and results in vector tagged with 'Second'

   m1 - m2;  // Compilation error: No subtraction defined for 'Meter'!
   m1 + s2;  // Compilation error: No addition between 'Meter' and 'Second' defined!

   m1 * m2;  // Compiles and results in vector tagged with 'SquareMeter'
   m1 / s1;  // Compiles and results in vector tagged with 'MeterPerSecond'
   \endcode

// At this point it is possible to use the \c pow2() function for vectors and matrices tagged with
// \c Meter since \c pow2() is based on multiplication, which has already been declared. However,
// it is not possible to use the \c abs() function:

   \code
   pow2( m1 );  // Compiles and results in vector tagged with 'SquareMeter'
   abs ( m1 );  // Compilation error: No 'abs()' declared for the 'Meter' tag
   \endcode

// In order to enable the \c abs() function it also needs to be explicitly declared for the
// \c Meter tag:

   \code
   Meter abs( Meter );  // Enabling the 'abs()' function on 'Meter'

   abs ( m1 );  // Compiles and results in vector tagged with 'Meter'
   \endcode

// \n Previous: \ref vector_and_matrix_customization &nbsp; &nbsp; Next: \ref error_reporting_customization \n
*/
//*************************************************************************************************


//**Customization of the Error Reporting Mechanism*************************************************
/*!\page error_reporting_customization Customization of the Error Reporting Mechanism
//
// \tableofcontents
//
//
// \n \section error_reporting_background Background
// <hr>
//
// The default way of \b Blaze to report errors of any kind is to throw a standard exception.
// However, although in general this approach works well, in certain environments and under
// special circumstances exceptions may not be the mechanism of choice and a different error
// reporting mechanism may be desirable. For this reason, \b Blaze provides several macros,
// which enable the customization of the error reporting mechanism. Via these macros it is
// possible to replace the standard exceptions by some other exception type or a completely
// different approach to report errors.
//
//
// \n \section error_reporting_general_customization Customization of the Reporting Mechanism
// <hr>
//
// In some cases it might be necessary to adapt the entire error reporting mechanism and to
// replace it by some other means to signal failure. The primary macro for this purpose is the
// \c BLAZE_THROW macro:

   \code
   #define BLAZE_THROW( EXCEPTION ) \
      throw EXCEPTION
   \endcode

// This macro represents the default mechanism of the \b Blaze library to report errors of any
// kind. In order to customize the error reporing mechanism all that needs to be done is to
// define the macro prior to including any \b Blaze header file. This will cause the \b Blaze
// specific mechanism to be overridden. The following example demonstrates this by replacing
// exceptions by a call to a \c log() function and a direct call to abort:

   \code
   #define BLAZE_THROW( EXCEPTION ) \
      log( "..." ); \
      abort()

   #include <blaze/Blaze.h>
   \endcode

// Doing this will trigger a call to \c log() and an abort instead of throwing an exception
// whenever an error (such as an invalid argument) is detected.
//
// \note It is possible to execute several statements instead of executing a single statement to
// throw an exception. Also note that it is recommended to define the macro such that a subsequent
// semicolon is required!
//
// \warning This macro is provided with the intention to assist in adapting \b Blaze to special
// conditions and environments. However, the customization of the error reporting mechanism via
// this macro can have a significant effect on the library. Thus be advised to use the macro
// with due care!
//
//
// \n \section error_reporting_exception_customization Customization of the Type of Exceptions
// <hr>
//
// In addition to the customization of the entire error reporting mechanism it is also possible
// to customize the type of exceptions being thrown. This can be achieved by customizing any
// number of the following macros:

   \code
   #define BLAZE_THROW_BAD_ALLOC \
      BLAZE_THROW( std::bad_alloc() )

   #define BLAZE_THROW_LOGIC_ERROR( MESSAGE ) \
      BLAZE_THROW( std::logic_error( MESSAGE ) )

   #define BLAZE_THROW_INVALID_ARGUMENT( MESSAGE ) \
      BLAZE_THROW( std::invalid_argument( MESSAGE ) )

   #define BLAZE_THROW_LENGTH_ERROR( MESSAGE ) \
      BLAZE_THROW( std::length_error( MESSAGE ) )

   #define BLAZE_THROW_OUT_OF_RANGE( MESSAGE ) \
      BLAZE_THROW( std::out_of_range( MESSAGE ) )

   #define BLAZE_THROW_RUNTIME_ERROR( MESSAGE ) \
      BLAZE_THROW( std::runtime_error( MESSAGE ) )
   \endcode

// In order to customize the type of exception the according macro has to be defined prior to
// including any \b Blaze header file. This will override the \b Blaze default behavior. The
// following example demonstrates this by replacing \c std::invalid_argument by a custom
// exception type:

   \code
   class InvalidArgument
   {
    public:
      InvalidArgument();
      explicit InvalidArgument( const std::string& message );
      // ...
   };

   #define BLAZE_THROW_INVALID_ARGUMENT( MESSAGE ) \
      BLAZE_THROW( InvalidArgument( MESSAGE ) )

   #include <blaze/Blaze.h>
   \endcode

// By manually defining the macro, an \c InvalidArgument exception is thrown instead of a
// \c std::invalid_argument exception. Note that it is recommended to define the macro such
// that a subsequent semicolon is required!
//
// \warning These macros are provided with the intention to assist in adapting \b Blaze to
// special conditions and environments. However, the customization of the type of an exception
// via this macro may have an effect on the library. Thus be advised to use the macro with due
// care!
//
//
// \n \section error_reporting_special_errors Customization of Special Errors
// <hr>
//
// Last but not least it is possible to customize the error reporting for special kinds of errors.
// This can be achieved by customizing any number of the following macros:

   \code
   #define BLAZE_THROW_DIVISION_BY_ZERO( MESSAGE ) \
      BLAZE_THROW_RUNTIME_ERROR( MESSAGE )

   #define BLAZE_THROW_LAPACK_ERROR( MESSAGE ) \
      BLAZE_THROW_RUNTIME_ERROR( MESSAGE )
   \endcode

// As explained in the previous sections, in order to customize the handling of special errors
// the according macro has to be defined prior to including any \b Blaze header file. This will
// override the \b Blaze default behavior.
//
//
// \n Previous: \ref grouping_tagging &nbsp; &nbsp; Next: \ref blas_functions \n
*/
//*************************************************************************************************


//**BLAS Functions*********************************************************************************
/*!\page blas_functions BLAS Functions
//
// \tableofcontents
//
//
// For vector/vector, matrix/vector and matrix/matrix multiplications with large dense matrices
// \b Blaze relies on the efficiency of BLAS libraries. For this purpose, \b Blaze implements
// several convenient C++ wrapper functions for several BLAS functions. The following sections
// give a complete overview of all available BLAS level 1, 2 and 3 functions.
//
//
// \n \section blas_level_1 BLAS Level 1
// <hr>
//
// \subsection blas_level_1_dotu Dot Product (dotu)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// dot product of two dense vectors (\c cblas_sdot(), \c cblas_ddot(), \c cblas_cdotu_sub(), and
// \c cblas_zdotu_sub()):

   \code
   namespace blaze {

   float dotu( blas_int_t n, const float* x, blas_int_t incX, const float* y, blas_int_t incY );

   double dotu( blas_int_t n, const double* x, blas_int_t incX, const double* y, blas_int_t incY );

   complex<float> dotu( blas_int_t n, const complex<float>* x, blas_int_t incX,
                        const complex<float>* y, blas_int_t incY );

   complex<double> dotu( blas_int_t n, const complex<double>* x, blas_int_t incX,
                         const complex<double>* y, blas_int_t incY );

   template< typename VT1, bool TF1, typename VT2, bool TF2 >
   ElementType_<VT1> dotu( const DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& y );

   } // namespace blaze
   \endcode

// \subsection blas_level_1_dotc Complex Conjugate Dot Product (dotc)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// complex conjugate dot product of two dense vectors (\c cblas_sdot(), \c cblas_ddot(),
// \c cblas_cdotc_sub(), and \c cblas_zdotc_sub()):

   \code
   namespace blaze {

   float dotc( blas_int_t n, const float* x, blas_int_t incX, const float* y, blas_int_t incY );

   double dotc( blas_int_t n, const double* x, blas_int_t incX, const double* y, blas_int_t incY );

   complex<float> dotc( blas_int_t n, const complex<float>* x, blas_int_t incX,
                        const complex<float>* y, blas_int_t incY );

   complex<double> dotc( blas_int_t n, const complex<double>* x, blas_int_t incX,
                         const complex<double>* y, blas_int_t incY );

   template< typename VT1, bool TF1, typename VT2, bool TF2 >
   ElementType_<VT1> dotc( const DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& y );

   } // namespace blaze
   \endcode

// \subsection blas_level_1_axpy Axpy Product (axpy)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// axpy product of two dense vectors (\c cblas_saxpy(), \c cblas_daxpy(), \c cblas_caxpy(), and
// \c cblas_zaxpy()):

   \code
   namespace blaze {

   void axpy( blas_int_t n, float alpha, const float* x, blas_int_t incX, float* y, blas_int_t incY );

   void axpy( blas_int_t n, double alpha, const double* x, blas_int_t incX, double* y, blas_int_t incY );

   void axpy( blas_int_t n, complex<float> alpha, const complex<float>* x,
              blas_int_t incX, complex<float>* y, blas_int_t incY );

   void axpy( blas_int_t n, complex<double> alpha, const complex<double>* x,
              blas_int_t incX, complex<double>* y, blas_int_t incY );

   template< typename VT1, bool TF1, typename VT2, bool TF2, typename ST >
   void axpy( const DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& y, ST alpha );

   } // namespace blaze
   \endcode

// \n \section blas_level_2 BLAS Level 2
// <hr>
//
// \subsection blas_level_2_gemv General Matrix/Vector Multiplication (gemv)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// general matrix/vector multiplication (\c cblas_sgemv(), \c cblas_dgemv(), \c cblas_cgemv(),
// and \c cblas_zgemv()):

   \code
   namespace blaze {

   void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
              float alpha, const float* A, blas_int_t lda, const float* x, blas_int_t incX,
              float beta, float* y, blas_int_t incY );

   void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
              double alpha, const double* A, blas_int_t lda, const double* x, blas_int_t incX,
              double beta, double* y, blas_int_t incY );

   void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
              complex<float> alpha, const complex<float>* A, blas_int_t lda,
              const complex<float>* x, blas_int_t incX, complex<float> beta,
              complex<float>* y, blas_int_t incY );

   void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
              complex<double> alpha, const complex<double>* A, blas_int_t lda,
              const complex<double>* x, blas_int_t incX, complex<double> beta,
              complex<double>* y, blas_int_t incY );

   } // namespace blaze
   \endcode

// \n \subsection blas_level_2_trmv Triangular Matrix/Vector Multiplication (trmv)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// matrix/vector multiplication with a triangular matrix (\c cblas_strmv(), \c cblas_dtrmv(),
// \c cblas_ctrmv(), and \c cblas_ztrmv()):

   \code
   namespace blaze {

   void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA, CBLAS_DIAG diag,
              blas_int_t n, const float* A, blas_int_t lda, float* x, blas_int_t incX );

   void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA, CBLAS_DIAG diag,
              blas_int_t n, const double* A, blas_int_t lda, double* x, blas_int_t incX );

   void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA, CBLAS_DIAG diag,
              blas_int_t n, const complex<float>* A, blas_int_t lda, complex<float>* x, blas_int_t incX );

   void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA, CBLAS_DIAG diag,
              blas_int_t n, const complex<double>* A, blas_int_t lda, complex<double>* x, blas_int_t incX );

   template< typename VT, typename MT, bool SO >
   void trmv( DenseVector<VT,false>& x, const DenseMatrix<MT,SO>& A, CBLAS_UPLO uplo );

   template< typename VT, typename MT, bool SO >
   void trmv( DenseVector<VT,true>& x, const DenseMatrix<MT,SO>& A, CBLAS_UPLO uplo );

   } // namespace blaze
   \endcode

// \n \section blas_level_3 BLAS Level 3
// <hr>
//
// \subsection blas_level_3_gemm General Matrix/Matrix Multiplication (gemm)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// general matrix/matrix multiplication (\c cblas_sgemm(), \c cblas_dgemm(), \c cblas_cgemm(),
// and \c cblas_zgemm()):

   \code
   namespace blaze {

   void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
              blas_int_t m, blas_int_t n, blas_int_t k, float alpha, const float* A,
              blas_int_t lda, const float* B, blas_int_t ldb, float beta, float* C,
              blas_int_t ldc );

   void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
              blas_int_t m, blas_int_t n, blas_int_t k, double alpha, const double* A,
              blas_int_t lda, const double* B, blas_int_t ldb, double beta, float* C,
              blas_int_t ldc );

   void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
              blas_int_t m, blas_int_t n, blas_int_t k, complex<float> alpha,
              const complex<float>* A, blas_int_t lda, const complex<float>* B,
              blas_int_t ldb, complex<float> beta, float* C, blas_int_t ldc );

   void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
              blas_int_t m, blas_int_t n, blas_int_t k, complex<double> alpha,
              const complex<double>* A, blas_int_t lda, const complex<double>* B,
              blas_int_t ldb, complex<double> beta, float* C, blas_int_t ldc );x

   } // namespace blaze
   \endcode

// \n \subsection blas_level_3_trmm Triangular Matrix/Matrix Multiplication (trmm)
//
// The following wrapper functions provide a generic interface for the BLAS functions for the
// matrix/matrix multiplication with a triangular matrix (\c cblas_strmm(), \c cblas_dtrmm(),
// \c cblas_ctrmm(), and \c cblas_ztrmm()):

   \code
   namespace blaze {

   void trmm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, float alpha, const float* A,
              blas_int_t lda, float* B, blas_int_t ldb );

   void trmm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, double alpha, const double* A,
              blas_int_t lda, double* B, blas_int_t ldb );

   void trmm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<float> alpha,
              const complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb );

   void trmm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<double> alpha,
              const complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb );

   template< typename MT1, bool SO1, typename MT2, bool SO2, typename ST >
   void trmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
              CBLAS_SIDE side, CBLAS_UPLO uplo, ST alpha );

   } // namespace blaze
   \endcode

// \n \subsection blas_level_3_trsm Triangular System Solver (trsm)
//
// The following wrapper functions provide a generic interface for the BLAS functions for solving
// a triangular system of equations (\c cblas_strsm(), \c cblas_dtrsm(), \c cblas_ctrsm(), and
// \c cblas_ztrsm()):

   \code
   namespace blaze {

   void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, float alpha, const float* A,
              blas_int_t lda, float* B, blas_int_t ldb );

   void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, double alpha, const double* A,
              blas_int_t lda, double* B, blas_int_t ldb );

   void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<float> alpha,
              const complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb );

   void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
              CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<double> alpha,
              const complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb );

   template< typename MT, bool SO, typename VT, bool TF, typename ST >
   void trsm( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b,
              CBLAS_SIDE side, CBLAS_UPLO uplo, ST alpha );

   template< typename MT1, bool SO1, typename MT2, bool SO2, typename ST >
   void trsm( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
              CBLAS_SIDE side, CBLAS_UPLO uplo, ST alpha );

   } // namespace blaze
   \endcode

// \n Previous: \ref error_reporting_customization &nbsp; &nbsp; Next: \ref lapack_functions \n
*/
//*************************************************************************************************


//**LAPACK Functions*******************************************************************************
/*!\page lapack_functions LAPACK Functions
//
// \tableofcontents
//
//
// \n \section lapack_introction Introduction
// <hr>
//
// The \b Blaze library makes extensive use of the LAPACK functionality for various compute tasks
// (including the decomposition, inversion and the computation of the determinant of dense matrices).
// For this purpose, \b Blaze implements several convenient C++ wrapper functions for all required
// LAPACK functions. The following sections give a complete overview of all available LAPACK wrapper
// functions. For more details on the individual LAPACK functions see the \b Blaze function
// documentation or the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// Most of the wrapper functions are implemented as thin wrappers around LAPACK functions. They
// provide the parameters of the original LAPACK functions and thus provide maximum flexibility:

   \code
   using blaze::blas_int_t;

   constexpr size_t N( 100UL );

   blaze::DynamicMatrix<double,blaze::columnMajor> A( N, N );
   // ... Initializing the matrix

   const blas_int_t m    ( numeric_cast<blas_int_t>( A.rows()    ) );  // == N
   const blas_int_t n    ( numeric_cast<blas_int_t>( A.columns() ) );  // == N
   const blas_int_t lda  ( numeric_cast<blas_int_t>( A.spacing() ) );  // >= N
   const blas_int_t lwork( n*lda );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );  // No initialization required
   const std::unique_ptr<double[]> work( new double[N] );          // No initialization required

   blas_int_t info( 0 );

   getrf( m, n, A.data(), lda, ipiv.get(), &info );                  // Reports failure via 'info'
   getri( n, A.data(), lda, ipiv.get(), work.get(), lwork, &info );  // Reports failure via 'info'
   \endcode

// In this context, \c blas_int_t is either a 32-bit or 64-bit signed integral type, depending
// on the setting of the \c BLAZE_BLAS_IS_64BIT compilation switch (see \ref blas_mode).
//
// Additionally, \b Blaze provides wrappers that provide a higher level of abstraction. These
// wrappers provide a maximum of convenience:

   \code
   using blaze::blas_int_t;

   constexpr size_t N( 100UL );

   blaze::DynamicMatrix<double,blaze::columnMajor> A( N, N );
   // ... Initializing the matrix

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );  // No initialization required

   getrf( A, ipiv.get() );  // Cannot fail
   getri( A, ipiv.get() );  // Reports failure via exception
   \endcode

// \note All functions only work for general, non-adapted matrices with \c float, \c double,
// \c complex<float>, or \c complex<double> element type. The attempt to call the function with
// adaptors or matrices of any other element type results in a compile time error!
//
// \note All functions can only be used if a fitting LAPACK library is available and linked to
// the final executable. Otherwise a call to this function will result in a linker error.
//
// \note For performance reasons all functions do only provide the basic exception safety guarantee,
// i.e. in case an exception is thrown the given matrix may already have been modified.
//
//
// \n \section lapack_decomposition Matrix Decomposition
// <hr>
//
// The following functions decompose/factorize the given dense matrix. Based on this decomposition
// the matrix can be inverted or used to solve a linear system of equations.
//
//
// \n \subsection lapack_lu_decomposition LU Decomposition
//
// The following functions provide an interface for the LAPACK functions \c sgetrf(), \c dgetrf(),
// \c cgetrf(), and \c zgetrf(), which compute the LU decomposition for the given general matrix:

   \code
   namespace blaze {

   void getrf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda, blas_int_t* ipiv, blas_int_t* info );

   void getrf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda, blas_int_t* ipiv, blas_int_t* info );

   void getrf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, blas_int_t* ipiv, blas_int_t* info );

   void getrf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, blas_int_t* ipiv, blas_int_t* info );

   template< typename MT, bool SO >
   void getrf( DenseMatrix<MT,SO>& A, blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// The decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed.
//
// \note The LU decomposition will never fail, even for singular matrices. However, in case of a
// singular matrix the resulting decomposition cannot be used for a matrix inversion or solving
// a linear system of equations.
//
//
// \n \subsection lapack_ldlt_decomposition LDLT Decomposition
//
// The following functions provide an interface for the LAPACK functions \c ssytrf(), \c dsytrf(),
// \c csytrf(), and \c zsytrf(), which compute the LDLT (Bunch-Kaufman) decomposition for the given
// symmetric indefinite matrix:

   \code
   namespace blaze {

   void sytrf( char uplo, blas_int_t n, float* A, blas_int_t lda, blas_int_t* ipiv, float* work, blas_int_t lwork, blas_int_t* info );

   void sytrf( char uplo, blas_int_t n, double* A, blas_int_t lda, blas_int_t* ipiv, double* work, blas_int_t lwork, blas_int_t* info );

   void sytrf( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, blas_int_t* ipiv, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void sytrf( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, blas_int_t* ipiv, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void sytrf( DenseMatrix<MT,SO>& A, char uplo, blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// The decomposition has the form

                      \f[ A = U D U^{T} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// \note The Bunch-Kaufman decomposition will never fail, even for singular matrices. However, in
// case of a singular matrix the resulting decomposition cannot be used for a matrix inversion or
// solving a linear system of equations.
//
//
// \n \subsection lapack_ldlh_decomposition LDLH Decomposition
//
// The following functions provide an interface for the LAPACK functions \c chetrf() and \c zsytrf(),
// which compute the LDLH (Bunch-Kaufman) decomposition for the given Hermitian indefinite matrix:

   \code
   namespace blaze {

   void hetrf( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, blas_int_t* ipiv, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void hetrf( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, blas_int_t* ipiv, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void hetrf( DenseMatrix<MT,SO>& A, char uplo, blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// The decomposition has the form

                      \f[ A = U D U^{H} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{H} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// \note The Bunch-Kaufman decomposition will never fail, even for singular matrices. However, in
// case of a singular matrix the resulting decomposition cannot be used for a matrix inversion or
// solving a linear system of equations.
//
//
// \n \subsection lapack_llh_decomposition Cholesky Decomposition
//
// The following functions provide an interface for the LAPACK functions \c spotrf(), \c dpotrf(),
// \c cpotrf(), and \c zpotrf(), which compute the Cholesky (LLH) decomposition for the given
// positive definite matrix:

   \code
   namespace blaze {

   void potrf( char uplo, blas_int_t n, float* A, blas_int_t lda, blas_int_t* info );

   void potrf( char uplo, blas_int_t n, double* A, blas_int_t lda, blas_int_t* info );

   void potrf( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, blas_int_t* info );

   void potrf( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, blas_int_t* info );

   template< typename MT, bool SO >
   void potrf( DenseMatrix<MT,SO>& A, char uplo );

   } // namespace blaze
   \endcode

// The decomposition has the form

                      \f[ A = U^{T} U \texttt{ (if uplo = 'U'), or }
                          A = L L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U is an upper triangular matrix and \c L is a lower triangular matrix. The Cholesky
// decomposition fails if the given matrix \a A is not a positive definite matrix. In this case
// a \c std::invalid_argument exception is thrown.
//
//
// \n \subsection lapack_qr_decomposition QR Decomposition
//
// The following functions provide an interface for the LAPACK functions \c sgeqrf(), \c dgeqrf(),
// \c cgeqrf(), and \c zgeqrf(), which compute the QR decomposition of the given general matrix:

   \code
   namespace blaze {

   void geqrf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void geqrf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   void geqrf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void geqrf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void geqrf( DenseMatrix<MT,SO>& A, typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The decomposition has the form

                              \f[ A = Q \cdot R, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(1) H(2) . . . H(k) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:m)</tt> is stored on exit in <tt>A(i+1:m,i)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and above the diagonal of the matrix contain the
// min(\a m,\a n)-by-\a n upper trapezoidal matrix \c R (\c R is upper triangular if \a m >= \a n);
// the elements below the diagonal, with the array \c tau, represent the orthogonal matrix \c Q as
// a product of min(\a m,\a n) elementary reflectors.
//
// The following functions provide an interface for the LAPACK functions \c sorgqr(), \c dorgqr(),
// \c sorg2r(), \c dorg2r(), \c cungqr(), \c zunqqr(), \c cung2r(), and \c zung2r(), which
// reconstruct the \c Q matrix from a QR decomposition:

   \code
   namespace blaze {

   void orgqr( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void orgqr( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void orgqr( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void org2r( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t* info );

   void org2r( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t* info );

   template< typename MT, bool SO >
   void org2r( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ungqr( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void ungqr( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void ungqr( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ung2r( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t* info );

   void ung2r( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t* info );

   template< typename MT, bool SO >
   void ung2r( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The following functions provide an interface for the LAPACK functions \c sormqr(), \c dormqr(),
// \c cunmqr(), and \c zunmqr(), which can be used to multiply a matrix with the \c Q matrix from
// a QR decomposition:

   \code
   namespace blaze {

   void ormqr( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const float* A, blas_int_t lda, const float* tau, float* C, blas_int_t ldc, float* work, blas_int_t lwork, blas_int_t* info );

   void ormqr( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const double* A, blas_int_t lda, const double* tau, double* C, blas_int_t ldc, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void ormqr( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A, char side, char trans, const ElementType_<MT2>* tau );


   void unmqr( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* C, blas_int_t ldc, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void unmqr( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* C, blas_int_t ldc, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO, typename MT2 >
   void unmqr( DenseMatrix<MT1,SO>& C, DenseMatrix<MT2,SO>& A, char side, char trans, ElementType_<MT2>* tau );

   } // namespace blaze
   \endcode

// \n \subsection lapack_rq_decomposition RQ Decomposition
//
// The following functions provide an interface for the LAPACK functions \c sgerqf(), \c dgerqf(),
// \c cgerqf(), and \c zgerqf(), which compute the RQ decomposition of the given general matrix:

   \code
   namespace blaze {

   void gerqf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void gerqf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   void gerqf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void gerqf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void gerqf( DenseMatrix<MT,SO>& A, typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The decomposition has the form

                              \f[ A = R \cdot Q, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(1) H(2) . . . H(k) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(n-k+i+1:n) = 0</tt> and
// <tt>v(n-k+i) = 1</tt>. <tt>v(1:n-k+i-1)</tt> is stored on exit in <tt>A(m-k+i,1:n-k+i-1)</tt>,
// and \c tau in \c tau(i). Thus in case \a m <= \a n, the upper triangle of the subarray
// <tt>A(1:m,n-m+1:n)</tt> contains the \a m-by-\a m upper triangular matrix \c R and in case
// \a m >= \a n, the elements on and above the (\a m-\a n)-th subdiagonal contain the \a m-by-\a n
// upper trapezoidal matrix \c R; the remaining elements in combination with the array \c tau
// represent the orthogonal matrix \c Q as a product of min(\a m,\a n) elementary reflectors.
//
// The following functions provide an interface for the LAPACK functions \c sorgrq(), \c dorgrq(),
// \c sorgr2(), \c dorgr2(), \c cungrq(), \c zunqrq(), \c cungr2(), and \c zunqr2(), which
// reconstruct the \c Q matrix from a RQ decomposition:

   \code
   namespace blaze {

   void orgrq( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void orgrq( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void orgrq( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void orgr2( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t* info );

   void orgr2( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t* info );

   template< typename MT, bool SO >
   void orgr2( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ungrq( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void ungrq( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void ungrq( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ungr2( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t* info );

   void ungr2( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t* info );

   template< typename MT, bool SO >
   void ungr2( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The following functions provide an interface for the LAPACK functions \c sormrq(), \c dormrq(),
// \c cunmrq(), and \c zunmrq(), which can be used to multiply a matrix with the \c Q matrix from
// a RQ decomposition:

   \code
   namespace blaze {

   void ormrq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const float* A, blas_int_t lda, const float* tau, float* C, blas_int_t ldc, float* work, blas_int_t lwork, blas_int_t* info );

   void ormrq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const double* A, blas_int_t lda, const double* tau, double* C, blas_int_t ldc, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void ormrq( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A, char side, char trans, const ElementType_<MT2>* tau );


   void unmrq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* C, blas_int_t ldc, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void unmrq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* C, blas_int_t ldc, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO, typename MT2 >
   void unmrq( DenseMatrix<MT1,SO>& C, DenseMatrix<MT2,SO>& A, char side, char trans, ElementType_<MT2>* tau );

   } // namespace blaze
   \endcode

// \n \subsection lapack_ql_decomposition QL Decomposition
//
// The following functions provide an interface for the LAPACK functions \c sgeqlf(), \c dgeqlf(),
// \c cgeqlf(), and \c zgeqlf(), which compute the QL decomposition of the given general matrix:

   \code
   namespace blaze {

   void geqlf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void geqlf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   void geqlf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void geqlf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void geqlf( DenseMatrix<MT,SO>& A, typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The decomposition has the form

                              \f[ A = Q \cdot L, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(k) . . . H(2) H(1) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(m-k+i+1:m) = 0</tt> and
// <tt>v(m-k+i) = 1</tt>. <tt>v(1:m-k+i-1)</tt> is stored on exit in <tt>A(1:m-k+i-1,n-k+i)</tt>,
// and \c tau in \c tau(i). Thus in case \a m >= \a n, the lower triangle of the subarray
// A(m-n+1:m,1:n) contains the \a n-by-\a n lower triangular matrix \c L and in case \a m <= \a n,
// the elements on and below the (\a n-\a m)-th subdiagonal contain the \a m-by-\a n lower
// trapezoidal matrix \c L; the remaining elements in combination with the array \c tau represent
// the orthogonal matrix \c Q as a product of min(\a m,\a n) elementary reflectors.
//
// The following functions provide an interface for the LAPACK functions \c sorgql(), \c dorgql(),
// \c sorg2l(), \c dorg2l(), \c cungql(), \c zungql(), \c cung2l(), and \c zung2l(), which
// reconstruct the \c Q matrix from an QL decomposition:

   \code
   namespace blaze {

   void orgql( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void orgql( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void orgql( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void org2l( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t* info );

   void org2l( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t* info );

   template< typename MT, bool SO >
   void org2l( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ungql( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void ungql( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void ungql( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ung2l( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t* info );

   void ung2l( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t* info );

   template< typename MT, bool SO >
   void ung2l( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The following functions provide an interface for the LAPACK functions \c sormql(), \c dormql(),
// \c cunmql(), and \c zunmql(), which can be used to multiply a matrix with the \c Q matrix from
// a QL decomposition:

   \code
   namespace blaze {

   void ormql( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const float* A, blas_int_t lda, const float* tau, float* C, blas_int_t ldc, float* work, blas_int_t lwork, blas_int_t* info );

   void ormql( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const double* A, blas_int_t lda, const double* tau, double* C, blas_int_t ldc, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void ormql( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A, char side, char trans, const ElementType_<MT2>* tau );


   void unmql( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* C, blas_int_t ldc, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void unmql( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* C, blas_int_t ldc, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO, typename MT2 >
   void unmql( DenseMatrix<MT1,SO>& C, DenseMatrix<MT2,SO>& A, char side, char trans, ElementType_<MT2>* tau );

   } // namespace blaze
   \endcode

// \n \subsection lapack_lq_decomposition LQ Decomposition
//
// The following functions provide an interface for the LAPACK functions \c sgelqf(), \c dgelqf(),
// \c cgelqf(), and \c zgelqf(), which compute the LQ decomposition of the given general matrix:

   \code
   namespace blaze {

   void gelqf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void gelqf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   void gelqf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void gelqf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void gelqf( DenseMatrix<MT,SO>& A, typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The decomposition has the form

                              \f[ A = L \cdot Q, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(k) . . . H(2) H(1) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:n)</tt> is stored on exit in <tt>A(i,i+1:n)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and below the diagonal of the matrix contain the
// \a m-by-min(\a m,\a n) lower trapezoidal matrix \c L (\c L is lower triangular if \a m <= \a n);
// the elements above the diagonal, with the array \c tau, represent the orthogonal matrix \c Q
// as a product of min(\a m,\a n) elementary reflectors.
//
// The following functions provide an interface for the LAPACK functions \c sorglq(), \c dorglq(),
// \c sorgl2(), \c dorgl2(), \c cunglq(), \c zunqlq(), \c cungl2(), and \c zunql2(), which
// reconstruct the \c Q matrix from an LQ decomposition:

   \code
   namespace blaze {

   void orglq( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t lwork, blas_int_t* info );

   void orglq( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void orglq( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void orgl2( blas_int_t m, blas_int_t n, blas_int_t k, float* A, blas_int_t lda, const float* tau, float* work, blas_int_t* info );

   void orgl2( blas_int_t m, blas_int_t n, blas_int_t k, double* A, blas_int_t lda, const double* tau, double* work, blas_int_t* info );

   template< typename MT, bool SO >
   void orgl2( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void unglq( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void unglq( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void unglq( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );


   void ungl2( blas_int_t m, blas_int_t n, blas_int_t k, complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* work, blas_int_t* info );

   void ungl2( blas_int_t m, blas_int_t n, blas_int_t k, complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* work, blas_int_t* info );

   template< typename MT, bool SO >
   void ungl2( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );

   } // namespace blaze
   \endcode

// The following functions provide an interface for the LAPACK functions \c sormlq(), \c dormlq(),
// \c cunmlq(), and \c zunmlq(), which can be used to multiply a matrix with the \c Q matrix from
// a LQ decomposition:

   \code
   namespace blaze {

   void ormlq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const float* A, blas_int_t lda, const float* tau, float* C, blas_int_t ldc, float* work, blas_int_t lwork, blas_int_t* info );

   void ormlq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const double* A, blas_int_t lda, const double* tau, double* C, blas_int_t ldc, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void ormlq( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A, char side, char trans, const ElementType_<MT2>* tau );


   void unmlq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<float>* A, blas_int_t lda, const complex<float>* tau, complex<float>* C, blas_int_t ldc, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void unmlq( char side, char trans, blas_int_t m, blas_int_t n, blas_int_t k, const complex<double>* A, blas_int_t lda, const complex<double>* tau, complex<double>* C, blas_int_t ldc, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT1, bool SO, typename MT2 >
   void unmlq( DenseMatrix<MT1,SO>& C, DenseMatrix<MT2,SO>& A, char side, char trans, ElementType_<MT2>* tau );

   } // namespace blaze
   \endcode

// \n \section lapack_inversion Matrix Inversion
// <hr>
//
// Given a matrix that has already been decomposed, the following functions can be used to invert
// the matrix in-place.
//
//
// \n \subsection lapack_lu_inversion LU-based Inversion
//
// The following functions provide an interface for the LAPACK functions \c sgetri(), \c dgetri(),
// \c cgetri(), and \c zgetri(), which invert a general matrix that has already been decomposed by
// an \ref lapack_lu_decomposition :

   \code
   namespace blaze {

   void getri( blas_int_t n, float* A, blas_int_t lda, const blas_int_t* ipiv, float* work, blas_int_t lwork, blas_int_t* info );

   void getri( blas_int_t n, double* A, blas_int_t lda, const blas_int_t* ipiv, double* work, blas_int_t lwork, blas_int_t* info );

   void getri( blas_int_t n, complex<float>* A, blas_int_t lda, const blas_int_t* ipiv, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void getri( blas_int_t n, complex<double>* A, blas_int_t lda, const blas_int_t* ipiv, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO >
   void getri( DenseMatrix<MT,SO>& A, const blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// The functions fail if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_ldlt_inversion LDLT-based Inversion
//
// The following functions provide an interface for the LAPACK functions \c ssytri(), \c dsytri(),
// \c csytri(), and \c zsytri(), which invert a symmetric indefinite matrix that has already been
// decomposed by an \ref lapack_ldlt_decomposition :

   \code
   namespace blaze {

   void sytri( char uplo, blas_int_t n, float* A, blas_int_t lda, const blas_int_t* ipiv, float* work, blas_int_t* info );

   void sytri( char uplo, blas_int_t n, double* A, blas_int_t lda, const blas_int_t* ipiv, double* work, blas_int_t* info );

   void sytri( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, const blas_int_t* ipiv, complex<float>* work, blas_int_t* info );

   void sytri( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, const blas_int_t* ipiv, complex<double>* work, blas_int_t* info );

   template< typename MT, bool SO >
   void sytri( DenseMatrix<MT,SO>& A, char uplo, const blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// The functions fail if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_ldlh_inversion LDLH-based Inversion
//
// The following functions provide an interface for the LAPACK functions \c chetri() and
// \c zhetri(), which invert an Hermitian indefinite matrix that has already been decomposed by
// an \ref lapack_ldlh_decomposition :

   \code
   namespace blaze {

   void hetri( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, const blas_int_t* ipiv, complex<float>* work, blas_int_t* info );

   void hetri( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, const blas_int_t* ipiv, complex<double>* work, blas_int_t* info );

   template< typename MT, bool SO >
   void hetri( DenseMatrix<MT,SO>& A, char uplo, const blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// The functions fail if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_llh_inversion Cholesky-based Inversion
//
// The following functions provide an interface for the LAPACK functions \c spotri(), \c dpotri(),
// \c cpotri(), and \c zpotri(), which invert a positive definite matrix that has already been
// decomposed by an \ref lapack_llh_decomposition :

   \code
   namespace blaze {

   void potri( char uplo, blas_int_t n, float* A, blas_int_t lda, blas_int_t* info );

   void potri( char uplo, blas_int_t n, double* A, blas_int_t lda, blas_int_t* info );

   void potri( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, blas_int_t* info );

   void potri( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, blas_int_t* info );

   template< typename MT, bool SO >
   void potri( DenseMatrix<MT,SO>& A, char uplo );

   } // namespace blaze
   \endcode

// The functions fail if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the given matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_triangular_inversion Inversion of Triangular Matrices
//
// The following functions provide an interface for the LAPACK functions \c strtri(), \c dtrtri(),
// \c ctrtri(), and \c ztrtri(), which invert the given triangular matrix in-place:

   \code
   namespace blaze {

   void trtri( char uplo, char diag, blas_int_t n, float* A, blas_int_t lda, blas_int_t* info );

   void trtri( char uplo, char diag, blas_int_t n, double* A, blas_int_t lda, blas_int_t* info );

   void trtri( char uplo, char diag, blas_int_t n, complex<float>* A, blas_int_t lda, blas_int_t* info );

   void trtri( char uplo, char diag, blas_int_t n, complex<double>* A, blas_int_t lda, blas_int_t* info );

   template< typename MT, bool SO >
   void trtri( DenseMatrix<MT,SO>& A, char uplo, char diag );

   } // namespace blaze
   \endcode

// The functions fail if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the given \a diag argument is neither 'U' nor 'N';
//  - ... the given matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \section lapack_substitution Substitution
// <hr>
//
// Given a matrix that has already been decomposed the following functions can be used to perform
// the forward/backward substitution step to compute the solution to a system of linear equations.
// Note that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems:
//
// Single right-hand side:
//  - \f$ A  *x=b \f$ if \a A is column-major
//  - \f$ A^T*x=b \f$ if \a A is row-major
//
// Multiple right-hand sides:
//  - \f$ A  *X  =B   \f$ if both \a A and \a B are column-major
//  - \f$ A^T*X  =B   \f$ if \a A is row-major and \a B is column-major
//  - \f$ A  *X^T=B^T \f$ if \a A is column-major and \a B is row-major
//  - \f$ A^T*X^T=B^T \f$ if both \a A and \a B are row-major
//
// In this context the general system matrix \a A is a n-by-n matrix that has already been
// factorized by the according decomposition function, \a x and \a b are n-dimensional vectors
// and \a X and \a B are either row-major m-by-n matrices or column-major n-by-m matrices.
//
//
// \n \subsection lapack_lu_substitution LU-based Substitution
//
// The following functions provide an interface for the LAPACK functions \c sgetrs(), \c dgetrs(),
// \c cgetrs(), and \c zgetrs(), which perform the substitution step for a general matrix that has
// already been decomposed by an \ref lapack_lu_decomposition :

   \code
   namespace blaze {

   void getrs( char trans, blas_int_t n, blas_int_t nrhs, const float* A, blas_int_t lda, const blas_int_t* ipiv, float* B, blas_int_t ldb, blas_int_t* info );

   void getrs( char trans, blas_int_t n, blas_int_t nrhs, const double* A, blas_int_t lda, const blas_int_t* ipiv, double* B, blas_int_t ldb, blas_int_t* info );

   void getrs( char trans, blas_int_t n, const complex<float>* A, blas_int_t lda, const blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void getrs( char trans, blas_int_t n, const complex<double>* A, blas_int_t lda, const blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void getrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char trans, const blas_int_t* ipiv );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void getrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char trans, const blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems (see \ref lapack_substitution). If the function exits
// successfully, the vector \a b or the matrix \a B contain the solution(s) of the linear system of
// equations. The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a trans argument is neither 'N' nor 'T' nor 'C';
//  - ... the sizes of the two given matrices do not match.
//
// The first four functions report failure via the \c info argument, the last two functions throw
// a \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_ldlt_substitution LDLT-based Substitution
//
// The following functions provide an interface for the LAPACK functions \c ssytrs(), \c dsytrs(),
// \c csytrs(), and \c zsytrs(), which perform the substitution step for a symmetric indefinite
// matrix that has already been decomposed by an \ref lapack_ldlt_decomposition :

   \code
   namespace blaze {

   void sytrs( char uplo, blas_int_t n, blas_int_t nrhs, const float* A, blas_int_t lda, const blas_int_t* ipiv, float* B, blas_int_t ldb, blas_int_t* info );

   void sytrs( char uplo, blas_int_t n, blas_int_t nrhs, const double* A, blas_int_t lda, const blas_int_t* ipiv, double* B, blas_int_t ldb, blas_int_t* info );

   void sytrs( char uplo, blas_int_t n, blas_int_t nrhs, const complex<float>* A, blas_int_t lda, const blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void sytrs( char uplo, blas_int_t n, blas_int_t nrhs, const complex<double>* A, blas_int_t lda, const blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void sytrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo, const blas_int_t* ipiv );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void sytrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo, const blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems (see \ref lapack_substitution). If the function exits
// successfully, the vector \a b or the matrix \a B contain the solution(s) of the linear system of
// equations. The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the sizes of the two given matrices do not match.
//
// The first four functions report failure via the \c info argument, the last two functions throw
// a \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_ldlh_substitution LDLH-based Substitution
//
// The following functions provide an interface for the LAPACK functions \c chetrs(), and \c zhetrs(),
// which perform the substitution step for an Hermitian indefinite matrix that has already been
// decomposed by an \ref lapack_ldlh_decomposition :

   \code
   namespace blaze {

   void hetrs( char uplo, blas_int_t n, blas_int_t nrhs, const complex<float>* A, blas_int_t lda, const blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void hetrs( char uplo, blas_int_t n, blas_int_t nrhs, const complex<double>* A, blas_int_t lda, const blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void hetrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo, const blas_int_t* ipiv );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void hetrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo, const blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems (see \ref lapack_substitution). If the function exits
// successfully, the vector \a b or the matrix \a B contain the solution(s) of the linear system of
// equations. The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the sizes of the two given matrices do not match.
//
// The first two functions report failure via the \c info argument, the last two functions throw
// a \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_llh_substitution Cholesky-based Substitution
//
// The following functions provide an interface for the LAPACK functions \c spotrs(), \c dpotrs(),
// \c cpotrs(), and \c zpotrs(), which perform the substitution step for a positive definite matrix
// that has already been decomposed by an \ref lapack_llh_decomposition :

   \code
   namespace blaze {

   void potrs( char uplo, blas_int_t n, blas_int_t nrhs, const float* A, blas_int_t lda, float* B, blas_int_t ldb, blas_int_t* info );

   void potrs( char uplo, blas_int_t n, blas_int_t nrhs, const double* A, blas_int_t lda, double* B, blas_int_t ldb, blas_int_t* info );

   void potrs( char uplo, blas_int_t n, blas_int_t nrhs, const complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void potrs( char uplo, blas_int_t n, blas_int_t nrhs, const complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void potrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void potrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems (see \ref lapack_substitution). If the function exits
// successfully, the vector \a b or the matrix \a B contain the solution(s) of the linear system of
// equations. The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the sizes of the two given matrices do not match.
//
// The first two functions report failure via the \c info argument, the last two functions throw
// a \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_triangular_substitution Substitution for Triangular Matrices
//
// The following functions provide an interface for the LAPACK functions \c strtrs(), \c dtrtrs(),
// \c ctrtrs(), and \c ztrtrs(), which perform the substitution step for a triangular matrix:

   \code
   namespace blaze {

   void trtrs( char uplo, char trans, char diag, blas_int_t n, blas_int_t nrhs, const float* A, blas_int_t lda, float* B, blas_int_t ldb, blas_int_t* info );

   void trtrs( char uplo, char trans, char diag, blas_int_t n, blas_int_t nrhs, const double* A, blas_int_t lda, double* B, blas_int_t ldb, blas_int_t* info );

   void trtrs( char uplo, char trans, char diag, blas_int_t n, blas_int_t nrhs, const complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void trtrs( char uplo, char trans, char diag, blas_int_t n, blas_int_t nrhs, const complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void trtrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo, char trans, char diag );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void trtrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo, char trans, char diag );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems (see \ref lapack_substitution). If the function exits
// successfully, the vector \a b or the matrix \a B contain the solution(s) of the linear system of
// equations. The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the given \a trans argument is neither 'N' nor 'T' nor 'C';
//  - ... the given \a diag argument is neither 'U' nor 'N';
//  - ... the sizes of the two given matrices do not match.
//
// The first four functions report failure via the \c info argument, the last two functions throw
// a \c std::invalid_argument exception in case of an error.
//
//
// \n \section lapack_linear_system_solver Linear System Solver
// <hr>
//
// The following functions represent compound functions that perform both the decomposition step
// as well as the substitution step to compute the solution to a system of linear equations. Note
// that depending on the storage order of the system matrix and the given right-hand side the
// functions solve different equation systems:
//
// Single right-hand side:
//  - \f$ A  *x=b \f$ if \a A is column-major
//  - \f$ A^T*x=b \f$ if \a A is row-major
//
// Multiple right-hand sides:
//  - \f$ A  *X  =B   \f$ if both \a A and \a B are column-major
//  - \f$ A^T*X  =B   \f$ if \a A is row-major and \a B is column-major
//  - \f$ A  *X^T=B^T \f$ if \a A is column-major and \a B is row-major
//  - \f$ A^T*X^T=B^T \f$ if both \a A and \a B are row-major
//
// In this context the general system matrix \a A is a n-by-n matrix that has already been
// factorized by the according decomposition function, \a x and \a b are n-dimensional vectors
// and \a X and \a B are either row-major m-by-n matrices or column-major n-by-m matrices.
//
//
// \subsection lapack_lu_linear_system_solver LU-based Linear System Solver
//
// The following functions provide an interface for the LAPACK functions \c sgesv(), \c dgesv(),
// \c cgesv(), and \c zgesv(), which combine an \ref lapack_lu_decomposition and the according
// \ref lapack_lu_substitution :

   \code
   namespace blaze {

   void gesv( blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda, blas_int_t* ipiv, float* B, blas_int_t ldb, blas_int_t* info );

   void gesv( blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda, blas_int_t* ipiv, double* B, blas_int_t ldb, blas_int_t* info );

   void gesv( blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda, blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void gesv( blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda, blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void gesv( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, blas_int_t* ipiv );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void gesv( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side
// the functions solve different equation systems (see \ref lapack_linear_system_solver). If
// the function exits successfully, the vector \a b or the matrix \a B contain the solution(s)
// of the linear system of equations and \a A has been decomposed by means of an
// \ref lapack_lu_decomposition.
//
// The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given system matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_ldlt_linear_system_solver LDLT-based Linear System Solver
//
// The following functions provide an interface for the LAPACK functions \c ssysv(), \c dsysv(),
// \c csysv(), and \c zsysv(), which combine an \ref lapack_ldlt_decomposition and the according
// \ref lapack_ldlt_substitution :

   \code
   namespace blaze {

   void sysv( char uplo, blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda, blas_int_t* ipiv, float* B, blas_int_t ldb, float* work, blas_int_t lwork, blas_int_t* info );

   void sysv( char uplo, blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda, blas_int_t* ipiv, double* B, blas_int_t ldb, double* work, blas_int_t lwork, blas_int_t* info );

   void sysv( char uplo, blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda, blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void sysv( char uplo, blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda, blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void sysv( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo, blas_int_t* ipiv );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void sysv( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo, blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side
// the functions solve different equation systems (see \ref lapack_linear_system_solver). If
// the function exits successfully, the vector \a b or the matrix \a B contain the solution(s)
// of the linear system of equations and \a A has been decomposed by means of an
// \ref lapack_ldlt_decomposition.
//
// The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the sizes of the two given matrices do not match;
//  - ... the given system matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_ldlh_linear_system_solver LDLH-based Linear System Solver
//
// The following functions provide an interface for the LAPACK functions \c shesv(), \c dhesv(),
// \c chesv(), and \c zhesv(), which combine an \ref lapack_ldlh_decomposition and the according
// \ref lapack_ldlh_substitution :

   \code
   namespace blaze {

   void hesv( char uplo, blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda, blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, complex<float>* work, blas_int_t lwork, blas_int_t* info );

   void hesv( char uplo, blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda, blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, complex<double>* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void hesv( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo, blas_int_t* ipiv );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void hesv( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo, blas_int_t* ipiv );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side
// the functions solve different equation systems (see \ref lapack_linear_system_solver). If
// the function exits successfully, the vector \a b or the matrix \a B contain the solution(s)
// of the linear system of equations and \a A has been decomposed by means of an
// \ref lapack_ldlh_decomposition.
//
// The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the sizes of the two given matrices do not match;
//  - ... the given system matrix is singular and not invertible.
//
// The first two functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_llh_linear_system_solver Cholesky-based Linear System Solver
//
// The following functions provide an interface for the LAPACK functions \c sposv(), \c dposv(),
// \c cposv(), and \c zposv(), which combine an \ref lapack_llh_decomposition and the according
// \ref lapack_llh_substitution :

   \code
   namespace blaze {

   void posv( char uplo, blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda, float* B, blas_int_t ldb, blas_int_t* info );

   void posv( char uplo, blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda, double* B, blas_int_t ldb, blas_int_t* info );

   void posv( char uplo, blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb, blas_int_t* info );

   void posv( char uplo, blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void posv( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo );

   template< typename MT1, bool SO1, typename MT2, bool SO2 >
   void posv( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char uplo );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side
// the functions solve different equation systems (see \ref lapack_linear_system_solver). If
// the function exits successfully, the vector \a b or the matrix \a B contain the solution(s)
// of the linear system of equations and \a A has been decomposed by means of an
// \ref lapack_llh_decomposition.
//
// The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the sizes of the two given matrices do not match;
//  - ... the given system matrix is singular and not invertible.
//
// The first four functions report failure via the \c info argument, the fifth function throws a
// \c std::invalid_argument exception in case of an error.
//
//
// \n \subsection lapack_triangular_linear_system_solver Linear System Solver for Triangular Matrices
//
// The following functions provide an interface for the LAPACK functions \c strsv(), \c dtrsv(),
// \c ctrsv(), and \c ztrsv():

   \code
   namespace blaze {

   void trsv( char uplo, char trans, char diag, blas_int_t n, const float* A, blas_int_t lda, float* x, blas_int_t incX );

   void trsv( char uplo, char trans, char diag, blas_int_t n, const double* A, blas_int_t lda, double* x, blas_int_t incX );

   void trsv( char uplo, char trans, char diag, blas_int_t n, const complex<float>* A, blas_int_t lda, complex<float>* x, blas_int_t incX );

   void trsv( char uplo, char trans, char diag, blas_int_t n, const complex<double>* A, blas_int_t lda, complex<double>* x, blas_int_t incX );

   template< typename MT, bool SO, typename VT, bool TF >
   void trsv( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char uplo, char trans, char diag );

   } // namespace blaze
   \endcode

// Note that depending on the storage order of the system matrix and the given right-hand side
// the functions solve different equation systems (see \ref lapack_linear_system_solver). If the
// function exits successfully, the vector \a b or the matrix \a B contain the solution(s) of the
// linear system of equations.
//
// The functions fail if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the given \a trans argument is neither 'N' nor 'T' nor 'C';
//  - ... the given \a diag argument is neither 'U' nor 'N'.
//
// The last function throws a \c std::invalid_argument exception in case of an error. Note that
// none of the functions does perform any test for singularity or near-singularity. Such tests
// must be performed prior to calling this function!
//
//
// \n \section lapack_eigenvalues Eigenvalues/Eigenvectors
//
// \subsection lapack_eigenvalues_general General Matrices
//
// The following functions provide an interface for the LAPACK functions \c sgeev(), \c dgeev(),
// \c cgeev(), and \c zgeev(), which compute the eigenvalues and optionally the eigenvectors of
// the given general matrix:

   \code
   namespace blaze {

   void geev( char jobvl, char jobvr, blas_int_t n, float* A, blas_int_t lda, float* wr, float* wi, float* VL, blas_int_t ldvl, float* VR, blas_int_t ldvr, float* work, blas_int_t lwork, blas_int_t* info );

   void geev( char jobvl, char jobvr, blas_int_t n, double* A, blas_int_t lda, double* wr, double* wi, double* VL, blas_int_t ldvl, double* VR, blas_int_t ldvr, double* work, blas_int_t lwork, blas_int_t* info );

   void geev( char jobvl, char jobvr, blas_int_t n, complex<float>* A, blas_int_t lda, complex<float>* w, complex<float>* VL, blas_int_t ldvl, complex<float>* VR, blas_int_t ldvr, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* info );

   void geev( char jobvl, char jobvr, blas_int_t n, complex<double>* A, blas_int_t lda, complex<double>* w, complex<double>* VL, blas_int_t ldvl, complex<double>* VR, blas_int_t ldvr, complex<double>* work, blas_int_t lwork, double* rwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void geev( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w );

   template< typename MT1, bool SO1, typename MT2, bool SO2, typename VT, bool TF >
   void geev( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL, DenseVector<VT,TF>& w );

   template< typename MT1, bool SO1, typename VT, bool TF, typename MT2, bool SO2 >
   void geev( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& VR );

   template< typename MT1, bool SO1, typename MT2, bool SO2, typename VT, bool TF, typename MT3, bool SO3 >
   void geev( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL, DenseVector<VT,TF>& w, DenseMatrix<MT3,SO3>& VR );

   } // namespace blaze
   \endcode

// The complex eigenvalues of the given matrix \a A are returned in the given vector \a w.
// Please note that no order of eigenvalues can be assumed, except that complex conjugate pairs
// of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part
// first.
//
// If \a VR is provided as an argument, the right eigenvectors are returned in the rows of \a VR
// in case \a VR is a row-major matrix and in the columns of \a VR in case \a VR is a column-major
// matrix. The right eigenvector \f$v[j]\f$ of \a A satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue.
//
// If \a VL is provided as an argument, the left eigenvectors are returned in the rows of \a VL
// in case \a VL is a row-major matrix and in the columns of \a VL in case \a VL is a column-major
// matrix. The left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// \a w, \a VL, and \a VR are resized to the correct dimensions (if possible and necessary). The
// functions fail if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given matrix \a VL is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VR is a fixed size matrix and the dimensions don't match;
//  - ... the eigenvalue computation fails.
//
// The first four functions report failure via the \c info argument, the last four functions throw
// an exception in case of an error.
//
//
// \n \subsection lapack_eigenvalues_symmetric Symmetric Matrices
//
// The following functions provide an interface for the LAPACK functions \c ssyev() and \c dsyev(),
// which compute the eigenvalues and eigenvectors of the given symmetric matrix:

   \code
   namespace blaze {

   void syev( char jobz, char uplo, blas_int_t n, float* A, blas_int_t lda, float* w, float* work, blas_int_t lwork, blas_int_t* info );

   void syev( char jobz, char uplo, blas_int_t n, double* A, blas_int_t lda, double* w, double* work, blas_int_t lwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void syev( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char jobz, char uplo );

   } // namespace blaze
   \endcode

// Alternatively, the following functions can be used, which provide an interface to the LAPACK
// functions \c ssyevd() and \c dsyevd(). In contrast to the \c syev() functions they use a
// divide-and-conquer strategy for the computation of the left and right eigenvectors:

   \code
   namespace blaze {

   void syevd( char jobz, char uplo, blas_int_t n, float* A, blas_int_t lda, float* w, float* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t liwork, blas_int_t* info );

   void syevd( char jobz, char uplo, blas_int_t n, double* A, blas_int_t lda, double* w, double* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t liwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void syevd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char jobz, char uplo );

   } // namespace blaze
   \endcode

// The real eigenvalues are returned in ascending order in the given vector \a w. \a w is resized
// to the correct size (if possible and necessary). In case \a A is a row-major matrix, the left
// eigenvectors are returned in the rows of \a A, in case \a A is a column-major matrix, the right
// eigenvectors are returned in the columns of \a A.
//
// The functions fail if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given \a jobz argument is neither \c 'V' nor \c 'N';
//  - ... the given \a uplo argument is neither \c 'L' nor \c 'U';
//  - ... the eigenvalue computation fails.
//
// The first two functions report failure via the \c info argument, the last function throws an
// exception in case of an error.
//
// Via the following functions, which wrap the LAPACK functions \c ssyevx() and \c dsyevx(), it
// is possible to compute a subset of eigenvalues and/or eigenvectors of a symmetric matrix:

   \code
   namespace blaze {

   void syevx( char jobz, char range, char uplo, blas_int_t n, float* A, blas_int_t lda, float vl, float vu, blas_int_t il, blas_int_t iu, float abstol, blas_int_t* m, float* w, float* Z, blas_int_t ldz, float* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info );

   void syevx( char jobz, char range, char uplo, blas_int_t n, double* A, blas_int_t lda, double vl, double vu, blas_int_t il, blas_int_t iu, double abstol, blas_int_t* m, double* w, double* Z, blas_int_t ldz, double* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   size_t syevx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char uplo );

   template< typename MT, bool SO, typename VT, bool TF, typename ST >
   size_t syevx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char uplo, ST low, ST upp );

   template< typename MT1, bool SO1, typename VT, bool TF, typename MT2, bool SO2 >
   size_t syevx( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& Z, char uplo );

   template< typename MT1, bool SO1, typename VT, bool TF, typename MT2, bool SO2, typename ST >
   size_t syevx( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& Z, char uplo, ST low, ST upp );

   } // namespace blaze
   \endcode

// The number of eigenvalues to be computed is specified by the lower bound \c low and the upper
// bound \c upp, which either form an integral or a floating point range.
//
// In case \a low and \a upp are of integral type, the function computes all eigenvalues in the
// index range \f$[low..upp]\f$. The \a num resulting real eigenvalues are stored in ascending
// order in the given vector \a w, which is either resized (if possible) or expected to be a
// \a num-dimensional vector. The eigenvectors are returned in the rows of \a Z in case \a Z is
// row-major matrix and in the columns of \a Z in case \a Z is a column-major matrix. \a Z is
// resized (if possible) or expected to be a \a num-by-\a n row-major matrix or a \a n-by-\a num
// column-major matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all eigenvalues
// in the half-open interval \f$(low..upp]\f$. The resulting real eigenvalues are stored in
// ascending order in the given vector \a w, which is either resized (if possible) or expected
// to be an \a n-dimensional vector. The eigenvectors are returned in the rows of \a Z in case
// \a Z is a row-major matrix and in the columns of \a Z in case \a Z is a column-major matrix.
// \a Z is resized (if possible) or expected to be a \a n-by-\a n matrix.
//
// The functions fail if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a Z is a fixed size matrix and the dimensions don't match;
//  - ... the given \a uplo argument is neither \c 'L' nor \c 'U';
//  - ... the eigenvalue computation fails.
//
// The first two functions report failure via the \c info argument, the last four functions throw
// an exception in case of an error.
//
//
// \n \subsection lapack_eigenvalues_hermitian Hermitian Matrices
//
// The following functions provide an interface for the LAPACK functions \c cheev() and \c zheev(),
// which compute the eigenvalues and eigenvectors of the given Hermitian matrix:

   \code
   namespace blaze {

   void heev( char jobz, char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, float* w, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* info );

   void heev( char jobz, char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, double* w, complex<double>* work, blas_int_t lwork, float* rwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void heev( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char jobz, char uplo );

   } // namespace blaze
   \endcode

// Alternatively, the following functions can be used, which provide an interface to the LAPACK
// functions \c cheevd() and \c zheevd(). In contrast to the \c heev() functions they use a
// divide-and-conquer strategy for the computation of the left and right eigenvectors:

   \code
   namespace blaze {

   void heevd( char jobz, char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, float* w, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* lrwork, blas_int_t* iwork, blas_int_t* liwork, blas_int_t* info );

   void heevd( char jobz, char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, double* w, complex<double>* work, blas_int_t lwork, double* rwork, blas_int_t lrwork, blas_int_t* iwork, blas_int_t* liwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void heevd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char jobz, char uplo );

   } // namespace blaze
   \endcode

// The real eigenvalues are returned in ascending order in the given vector \a w. \a w is resized
// to the correct size (if possible and necessary). In case \a A is a row-major matrix, the left
// eigenvectors are returned in the rows of \a A, in case \a A is a column-major matrix, the right
// eigenvectors are returned in the columns of \a A.
//
// The functions fail if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given \a jobz argument is neither \c 'V' nor \c 'N';
//  - ... the given \a uplo argument is neither \c 'L' nor \c 'U';
//  - ... the eigenvalue computation fails.
//
// The first two functions report failure via the \c info argument, the last function throws an
// exception in case of an error.
//
// Via the following functions, which wrap the LAPACK functions \c cheevx() and \c zheevx(), it
// is possible to compute a subset of eigenvalues and/or eigenvectors of an Hermitian matrix:

   \code
   namespace blaze {

   void heevx( char jobz, char range, char uplo, blas_int_t n, complex<float>* A, blas_int_t lda, float vl, float vu, blas_int_t il, blas_int_t iu, float abstol, blas_int_t* m, float* w, complex<float>* Z, blas_int_t ldz, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info );

   void heevx( char jobz, char range, char uplo, blas_int_t n, complex<double>* A, blas_int_t lda, double vl, double vu, blas_int_t il, blas_int_t iu, double abstol, blas_int_t* m, double* w, complex<double>* Z, blas_int_t ldz, complex<double>* work, blas_int_t lwork, double* rwork, blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   size_t heevx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char uplo );

   template< typename MT, bool SO, typename VT, bool TF, typename ST >
   size_t heevx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char uplo, ST low, ST upp );

   template< typename MT1, bool SO1, typename VT, bool TF, typename MT2, bool SO2 >
   size_t heevx( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& Z, char uplo );

   template< typename MT1, bool SO1, typename VT, bool TF, typename MT2, bool SO2, typename ST >
   size_t heevx( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& Z, char uplo, ST low, ST upp );

   } // namespace blaze
   \endcode

// The number of eigenvalues to be computed is specified by the lower bound \c low and the upper
// bound \c upp, which either form an integral or a floating point range.
//
// In case \a low and \a upp are of integral type, the function computes all eigenvalues in the
// index range \f$[low..upp]\f$. The \a num resulting real eigenvalues are stored in ascending
// order in the given vector \a w, which is either resized (if possible) or expected to be a
// \a num-dimensional vector. The eigenvectors are returned in the rows of \a Z in case \a Z is
// row-major matrix and in the columns of \a Z in case \a Z is a column-major matrix. \a Z is
// resized (if possible) or expected to be a \a num-by-\a n row-major matrix or a \a n-by-\a num
// column-major matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all eigenvalues
// in the half-open interval \f$(low..upp]\f$. The resulting real eigenvalues are stored in
// ascending order in the given vector \a w, which is either resized (if possible) or expected
// to be an \a n-dimensional vector. The eigenvectors are returned in the rows of \a Z in case
// \a Z is a row-major matrix and in the columns of \a Z in case \a Z is a column-major matrix.
// \a Z is resized (if possible) or expected to be a \a n-by-\a n matrix.
//
// The functions fail if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a Z is a fixed size matrix and the dimensions don't match;
//  - ... the given \a uplo argument is neither \c 'L' nor \c 'U';
//  - ... the eigenvalue computation fails.
//
// The first two functions report failure via the \c info argument, the last four functions throw
// an exception in case of an error.
//
//
// \n \section lapack_singular_values Singular Values/Singular Vectors
//
// The following functions provide an interface for the LAPACK functions \c sgesvd(), \c dgesvd(),
// \c cgesvd(), and \c zgesvd(), which perform a singular value decomposition (SVD) on the given
// general matrix:

   \code
   namespace blaze {

   void gesvd( char jobu, char jobv, blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float* s, float* U, blas_int_t ldu, float* V, blas_int_t ldv, float* work, blas_int_t lwork, blas_int_t* info );

   void gesvd( char jobu, char jobv, blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double* s, double* U, blas_int_t ldu, double* V, blas_int_t ldv, double* work, blas_int_t lwork, blas_int_t* info );

   void gesvd( char jobu, char jobv, blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, float* s, complex<float>* U, blas_int_t ldu, complex<float>* V, blas_int_t ldv, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* info );

   void gesvd( char jobu, char jobv, blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, double* s, complex<double>* U, blas_int_t ldu, complex<double>* V, blas_int_t ldv, complex<double>* work, blas_int_t lwork, double* rwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void gesvd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s, char jobu, char jobv );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
   void gesvd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, char jobu, char jobv );

   template< typename MT1, bool SO, typename VT, bool TF, typename MT2 >
   void gesvd( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V, char jobu, char jobv );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3 >
   void gesvd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, char jobu, char jobv );

   } // namespace blaze
   \endcode

// Alternatively, the following functions can be used, which provide an interface to the LAPACK
// functions \c sgesdd(), \c dgesdd(), \c cgesdd(), and \c zgesdd(). In contrast to the \c gesvd()
// functions they compute the singular value decomposition (SVD) of the given general matrix by
// applying a divide-and-conquer strategy for the computation of the left and right singular
// vectors:

   \code
   namespace blaze {

   void gesdd( char jobz, blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float* s, float* U, blas_int_t ldu, float* V, blas_int_t ldv, float* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t* info );

   void gesdd( char jobz, blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double* s, double* U, blas_int_t ldu, double* V, blas_int_t ldv, double* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t* info );

   void gesdd( char jobz, blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, float* s, complex<float>* U, blas_int_t ldu, complex<float>* V, blas_int_t ldv, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* iwork, blas_int_t* info );

   void gesdd( char jobz, blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, double* s, complex<double>* U, blas_int_t ldu, complex<double>* V, blas_int_t ldv, complex<double>* work, blas_int_t lwork, double* rwork, blas_int_t* iwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   void gesdd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
   void gesdd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, char jobz );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
   void gesdd( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V, char jobz );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3 >
   void gesdd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, char jobz );

   } // namespace blaze
   \endcode

// The resulting decomposition has the form

                          \f[ A = U \cdot S \cdot V, \f]

// where \a S is a \a m-by-\a n matrix, which is zero except for its min(\a m,\a n) diagonal
// elements, \a U is an \a m-by-\a m orthogonal matrix, and \a V is a \a n-by-\a n orthogonal
// matrix. The diagonal elements of \a S are the singular values of \a A, the first min(\a m,\a n)
// columns of \a U and rows of \a V are the left and right singular vectors of \a A, respectively.
//
// The resulting min(\a m,\a n) real and non-negative singular values are returned in descending
// order in the vector \a s, which is resized to the correct size (if possible and necessary).
//
// Via the following functions, which wrap the LAPACK functions \c sgesvdx(), \c dgesvdx(),
// \c cgesvdx(), and \c zgesvdx(), it is possible to compute a subset of singular values and/or
// vectors:

   \code
   namespace blaze {

   void gesvdx( char jobu, char jobv, char range, blas_int_t m, blas_int_t n, float* A, blas_int_t lda, float vl, float vu, blas_int_t il, blas_int_t iu, blas_int_t* ns, float* s, float* U, blas_int_t ldu, float* V, blas_int_t ldv, float* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t* info );

   void gesvdx( char jobu, char jobv, char range, blas_int_t m, blas_int_t n, double* A, blas_int_t lda, double vl, double vu, blas_int_t il, blas_int_t iu, blas_int_t* ns, double* s, double* U, blas_int_t ldu, double* V, blas_int_t ldv, double* work, blas_int_t lwork, blas_int_t* iwork, blas_int_t* info );

   void gesvdx( char jobu, char jobv, char range, blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda, float vl, float vu, blas_int_t il, blas_int_t iu, blas_int_t* ns, float* s, complex<float>* U, blas_int_t ldu, complex<float>* V, blas_int_t ldv, complex<float>* work, blas_int_t lwork, float* rwork, blas_int_t* iwork, blas_int_t* info );

   void gesvdx( char jobu, char jobv, char range, blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda, double vl, double vu, blas_int_t il, blas_int_t iu, blas_int_t* ns, double* s, complex<double>* U, blas_int_t ldu, complex<double>* V, blas_int_t ldv, complex<double>* work, blas_int_t lwork, double* rwork, blas_int_t* iwork, blas_int_t* info );

   template< typename MT, bool SO, typename VT, bool TF >
   size_t gesvdx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s );

   template< typename MT, bool SO, typename VT, bool TF, typename ST >
   size_t gesvdx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s, ST low, ST upp );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
   size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename ST >
   size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, ST low, ST upp );

   template< typename MT1, bool SO, typename VT, bool TF, typename MT2 >
   size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V );

   template< typename MT1, bool SO, typename VT, bool TF, typename MT2, typename ST >
   size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V, ST low, ST upp );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3 >
   size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V );

   template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3, typename ST >
   size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, ST low, ST upp );

   } // namespace blaze
   \endcode

// The number of singular values to be computed is specified by the lower bound \a low and the
// upper bound \a upp, which either form an integral or a floating point range.
//
// In case \a low and \a upp form are of integral type, the function computes all singular values
// in the index range \f$[low..upp]\f$. The \a num resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a \a num-dimensional vector. The resulting left singular vectors are stored
// in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-\a num matrix. The resulting right singular vectors are stored in the given matrix \a V,
// which is either resized (if possible) or expected to be a \a num-by-\a n matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all singular values
// in the half-open interval \f$(low..upp]\f$. The resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a min(\a m,\a n)-dimensional vector. The resulting left singular vectors are
// stored in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-min(\a m,\a n) matrix. The resulting right singular vectors are stored in the given
// matrix \a V, which is either resized (if possible) or expected to be a min(\a m,\a n)-by-\a n
// matrix.
//
// The functions fail if ...
//
//  - ... the given matrix \a U is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
//  - ... the given scalar values don't form a proper range;
//  - ... the singular value decomposition fails.
//
// The first four functions report failure via the \c info argument, the remaining functions throw
// an exception in case of an error.
//
//
// \n Previous: \ref blas_functions &nbsp; &nbsp; Next: \ref block_vectors_and_matrices \n
*/
//*************************************************************************************************


//**Block Vectors and Matrices*********************************************************************
/*!\page block_vectors_and_matrices Block Vectors and Matrices
//
// \tableofcontents
//
//
// \n \section block_vectors_and_matrices_general General Concepts
// <hr>
//
// In addition to fundamental element types, the \b Blaze library supports vectors and matrices
// with non-fundamental element type. For instance, it is possible to define block matrices by
// using a matrix type as the element type:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix< DynamicMatrix<double,rowMajor>, rowMajor > A;
   DynamicVector< DynamicVector<double,columnVector >, columnVector > x, y;

   // ... Resizing and initialization

   y = A * x;
   \endcode

// The matrix/vector multiplication in this example runs fully parallel and uses vectorization
// for every inner matrix/vector multiplication and vector addition.
//
//
// \n \section block_vectors_and_matrices_pitfalls Pitfalls
// <hr>
//
// The only thing to keep in mind when using non-fundamental element types is that all operations
// between the elements have to be well defined. More specifically, the size of vector and matrix
// elements has to match. The attempt to combine two non-matching elements results in either a
// compilation error (in case of statically sized elements) or an exception (for dynamically sized
// elements):

   \code
   DynamicVector< StaticVector<int,2UL> > a;
   DynamicVector< StaticVector<int,3UL> > b;

   DynamicVector< DynamicVector<int> > c( a + b );  // Compilation error: element size doesn't match
   \endcode

// Therefore please don't forget that dynamically sized elements (e.g. \c blaze::DynamicVector,
// \c blaze::HybridVector, \c blaze::DynamicMatrix, \c blaze::HybridMatrix, ...) need to be sized
// accordingly upfront.
//
//
// \n \section block_vectors_and_matrices_examples Examples
// <hr>
//
// The first example demonstrates the multiplication between a statically sized block matrix
// and a block vector:

   \code
   using namespace blaze;

   // ( ( 1 1 )  ( 2 2 ) )   ( ( 1 ) )   ( ( 10 ) )
   // ( ( 1 1 )  ( 2 2 ) )   ( ( 1 ) )   ( ( 10 ) )
   // (                  ) * (       ) = (        )
   // ( ( 3 3 )  ( 4 4 ) )   ( ( 2 ) )   ( ( 22 ) )
   // ( ( 3 3 )  ( 4 4 ) )   ( ( 2 ) )   ( ( 22 ) )

   using M2x2 = StaticMatrix<int,2UL,2UL,rowMajor>;
   using V2   = StaticVector<int,2UL,columnVector>;

   DynamicMatrix<M2x2,rowMajor> A{ { M2x2(1), M2x2(2) },
                                   { M2x2(3), M2x2(4) } };

   DynamicVector<V2,columnVector> x{ V2(1), V2(2) };

   DynamicVector<V2,columnVector> y( A * x );
   \endcode

// The second example shows the multiplication between a compressed block matrix with blocks of
// varying size and a compressed block vector:

   \code
   using namespace blaze;

   // ( ( 1 -2  3 )         ( 5 -1 ) )   ( (  1 ) )   ( ( -3 ) )
   // ( ( 4  1  0 )         ( 1  2 ) )   ( (  0 ) )   ( (  7 ) )
   // ( ( 0  2  4 )         ( 3  1 ) )   ( (  1 ) )   ( (  3 ) )
   // (                              )   (        )   (        )
   // (              ( 1 )           ) * ( (  2 ) ) = ( (  2 ) )
   // (                              )   (        )   (        )
   // ( ( 0 -1  1 )         ( 1  0 ) )   ( ( -1 ) )   ( (  0 ) )
   // ( ( 2 -1  2 )         ( 0  1 ) )   ( (  2 ) )   ( (  6 ) )

   using M3x3 = HybridMatrix<int,3UL,3UL,rowMajor>;
   using V3   = HybridVector<int,3UL,columnVector>;

   CompressedMatrix<M3x3,rowMajor> A( 3UL, 3UL, 5UL );
   A(0,0) = M3x3{ { 1, -2, 3 }, { 4, 1, 0 }, { 0, 2, 4 } };
   A(0,2) = M3x3{ { 5, -1 }, { 1, 2 }, { 3, 1 } };
   A(1,1) = M3x3{ { 1 } };
   A(2,0) = M3x3{ { 0, -1, 1 }, { 2, -1, 2 } };
   A(2,2) = M3x3{ { 1, 0 }, { 0, 1 } };

   CompressedVector<V3,columnVector> x( 3UL, 3UL );
   x[0] = V3{ 1, 0, 1 };
   x[1] = V3{ 2 };
   x[2] = V3{ -1, 2 };

   CompressedVector<V3,columnVector> y( A * x );
   \endcode

// \n Previous: \ref lapack_functions &nbsp; &nbsp; Next: \ref intra_statement_optimization \n
*/
//*************************************************************************************************


//**Intra-Statement Optimization*******************************************************************
/*!\page intra_statement_optimization Intra-Statement Optimization
//
// One of the prime features of the \b Blaze library is the automatic intra-statement optimization.
// In order to optimize the overall performance of every single statement \b Blaze attempts to
// rearrange the operands based on their types. For instance, the following addition of dense and
// sparse vectors

   \code
   blaze::DynamicVector<double> d1, d2, d3;
   blaze::CompressedVector<double> s1;

   // ... Resizing and initialization

   d3 = d1 + s1 + d2;
   \endcode

// is automatically rearranged and evaluated as

   \code
   // ...
   d3 = d1 + d2 + s1;  // <- Note that s1 and d2 have been rearranged
   \endcode

// This order of operands is highly favorable for the overall performance since the addition of
// the two dense vectors \c d1 and \c d2 can be handled much more efficiently in a vectorized
// fashion.
//
// This intra-statement optimization can have a tremendous effect on the performance of a statement.
// Consider for instance the following computation:

   \code
   blaze::DynamicMatrix<double> A, B;
   blaze::DynamicVector<double> x, y;

   // ... Resizing and initialization

   y = A * B * x;
   \endcode

// Since multiplications are evaluated from left to right, this statement would result in a
// matrix/matrix multiplication, followed by a matrix/vector multiplication. However, if the
// right subexpression is evaluated first, the performance can be dramatically improved since the
// matrix/matrix multiplication can be avoided in favor of a second matrix/vector multiplication.
// The \b Blaze library exploits this by automatically restructuring the expression such that the
// right multiplication is evaluated first:

   \code
   // ...
   y = A * ( B * x );
   \endcode

// Note however that although this intra-statement optimization may result in a measurable or
// even significant performance improvement, this behavior may be undesirable for several reasons,
// for instance because of numerical stability. Therefore, in case the order of evaluation matters,
// the best solution is to be explicit and to separate a statement into several statements:

   \code
   blaze::DynamicVector<double> d1, d2, d3;
   blaze::CompressedVector<double> s1;

   // ... Resizing and initialization

   d3  = d1 + s1;  // Compute the dense vector/sparse vector addition first ...
   d3 += d2;       // ... and afterwards add the second dense vector
   \endcode

   \code
   // ...
   blaze::DynamicMatrix<double> A, B, C;
   blaze::DynamicVector<double> x, y;

   // ... Resizing and initialization

   C = A * B;  // Compute the left-hand side matrix-matrix multiplication first ...
   y = C * x;  // ... before the right-hand side matrix-vector multiplication
   \endcode

// Alternatively, it is also possible to use the \c eval() function to fix the order of evaluation:

   \code
   blaze::DynamicVector<double> d1, d2, d3;
   blaze::CompressedVector<double> s1;

   // ... Resizing and initialization

   d3 = d1 + eval( s1 + d2 );
   \endcode

   \code
   blaze::DynamicMatrix<double> A, B;
   blaze::DynamicVector<double> x, y;

   // ... Resizing and initialization

   y = eval( A * B ) * x;
   \endcode

// \n Previous: \ref block_vectors_and_matrices &nbsp; &nbsp; Next: \ref faq \n
*/
//*************************************************************************************************


//**FAQ********************************************************************************************
/*!\page faq Frequently Asked Questions (FAQ)
//
// \tableofcontents
//
//
// <hr>
// \section faq_padding A StaticVector/StaticMatrix is larger than expected. Is this a bug?
//
// The size of a \ref vector_types_static_vector, \ref matrix_types_static_matrix,
// \ref vector_types_hybrid_vector, or \ref matrix_types_hybrid_matrix can indeed be larger
// than expected:

   \code
   StaticVector<int,3> a;
   StaticMatrix<int,3,3> A;

   sizeof( a );  // Evaluates to 16, 32, or even 64, but not 12
   sizeof( A );  // Evaluates to 48, 96, or even 144, but not 36
   \endcode

// In order to achieve the maximum possible performance the \b Blaze library tries to enable
// SIMD vectorization even for small vectors. For that reason \b Blaze by default uses padding
// elements for all dense vectors and matrices to guarantee that at least a single SIMD vector
// can be loaded. Depending on the used SIMD technology that can significantly increase the size
// of a \ref vector_types_static_vector, \ref matrix_types_static_matrix,
// \ref vector_types_hybrid_vector, or \ref matrix_types_hybrid_matrix :

   \code
   StaticVector<int,3> a;
   StaticMatrix<int,3,3> A;

   sizeof( a );  // Evaluates to 16 in case of SSE, 32 in case of AVX, and 64 in case of AVX-512
                 // (under the assumption that an integer occupies 4 bytes)
   sizeof( A );  // Evaluates to 48 in case of SSE, 96 in case of AVX, and 144 in case of AVX-512
                 // (under the assumption that an integer occupies 4 bytes)
   \endcode

// The configuration files <tt>./blaze/config/Padding.h</tt> provides a compile time switch
// that can be used to (de-)activate padding:

   \code
   #define BLAZE_DEFAULT_PADDING_FLAG blaze::padded
   \endcode

// Alternatively it is possible to (de-)activate padding via command line or by defining this
// symbol manually before including any \b Blaze header file:

   \code
   g++ ... -BLAZE_DEFAULT_PADDING_FLAG=blaze::padded ...
   \endcode

   \code
   #define BLAZE_DEFAULT_PADDING_FLAG blaze::padded
   #include <blaze/Blaze.h>
   \endcode

// If \c BLAZE_DEFAULT_PADDING_FLAG is set to \c blaze::padded, by default padding is enabled for
// \ref vector_types_static_vector, \ref vector_types_hybrid_vector, \ref matrix_types_static_matrix,
// and \ref matrix_types_hybrid_matrix. If it is set to \c blaze::unpadded, then padding is by
// default disabled. Note however that disabling padding can considerably reduce the performance
// of all dense vector and matrix operations!
//
//
// <hr>
// \section faq_alignment Despite disabling padding, a StaticVector/StaticMatrix is still larger than expected. Is this a bug?
//
// Despite disabling padding via the \c BLAZE_DEFAULT_PADDING_FLAG compile time switch (see
// \ref faq_padding), the size of a \ref vector_types_static_vector, \ref matrix_types_static_matrix,
// \ref vector_types_hybrid_vector, or \ref matrix_types_hybrid_matrix can still be larger than
// expected:

   \code
   #define BLAZE_DEFAULT_PADDING_FLAG blaze::unpadded
   #include <blaze/Blaze.h>

   StaticVector<int,3> a;
   StaticVector<int,5> b;

   sizeof( a );  // Always evaluates to 12
   sizeof( b );  // Evaluates to 32 with SSE (larger than expected) and to 20 with AVX or AVX-512 (expected)
   \endcode

// The reason for this behavior is the used SIMD technology. If SSE is used, which provides 128
// bit wide registers, a single SIMD pack can usually hold 4 integers (128 bit divided by 32 bit).
// Since the second vector contains enough elements is possible to benefit from vectorization.
// However, SSE requires an alignment of 16 bytes, which ultimately results in a total size of
// 32 bytes for the \c StaticVector (2 times 16 bytes due to 5 integer elements). If AVX or AVX-512
// is used, which provide 256 bit or 512 bit wide registers, a single SIMD vector can hold 8 or 16
// integers, respectively. Even the second vector does not hold enough elements to benefit from
// vectorization, which is why \b Blaze does not enforce a 32 byte (for AVX) or even 64 byte
// alignment (for AVX-512).
//
// It is possible to disable the SIMD-specific alignment for \ref vector_types_static_vector,
// \ref matrix_types_static_matrix, \ref vector_types_hybrid_vector, or \ref matrix_types_hybrid_matrix
// via the compile time switch in the <tt>./blaze/config/Alignment.h</tt> configuration file:

   \code
   #define BLAZE_DEFAULT_ALIGNMENT_FLAG blaze::aligned
   \endcode

// Alternatively it is possible set the default alignment flag via command line or by defining
// this symbol manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_DEFAULT_ALIGNMENT_FLAG=blaze::aligned ...
   \endcode

   \code
   #define BLAZE_DEFAULT_ALIGNMENT_FLAG blaze::aligned
   #include <blaze/Blaze.h>
   \endcode

// If \c BLAZE_DEFAULT_ALIGNMENT_FLAG is set to \c blaze::aligned then \ref vector_types_static_vector,
// \ref vector_types_hybrid_vector, \ref matrix_types_static_matrix, and \ref matrix_types_hybrid_matrix
// use aligned memory by default. If it is set to \c blaze::unaligned they don't enforce aligned
// memory. Note however that disabling alignment can considerably reduce the performance of all
// operations with these vector and matrix types!
//
// Alternatively it is possible to disable the vectorization entirely by the compile time switch
// in the <tt>./blaze/config/Vectorization.h</tt> configuration file:

   \code
   #define BLAZE_USE_VECTORIZATION 1
   \endcode

// It is also possible to (de-)activate vectorization via command line or by defining this symbol
// manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_USE_VECTORIZATION=1 ...
   \endcode

   \code
   #define BLAZE_USE_VECTORIZATION 1
   #include <blaze/Blaze.h>
   \endcode

// In case the switch is set to 1, vectorization is enabled and the \b Blaze library is allowed
// to use intrinsics and the necessary alignment to speed up computations. In case the switch is
// set to 0, vectorization is disabled entirely and the \b Blaze library chooses default,
// non-vectorized functionality for the operations. Note that deactivating the vectorization may
// pose a severe performance limitation for a large number of operations!
//
//
// <hr>
// \section faq_std_vector I experience crashes when using StaticVector/StaticMatrix in a std::vector. Is this a bug?
//
// With active vectorization the elements of a \ref vector_types_static_vector,
// \ref vector_types_hybrid_vector, \ref matrix_types_static_matrix, and \ref matrix_types_hybrid_matrix
// are possibly over-aligned to meet the alignment requirements of the available instruction set
// (SSE, AVX, AVX-512, ...). The alignment for fundamental types (\c short, \c int, \c float,
// \c double, ...) and complex types (\c complex<float>, \c complex<double>, ...) is 16 bytes
// for SSE, 32 bytes for AVX, and 64 bytes for AVX-512. All other types are aligned according to
// their intrinsic alignment:

   \code
   struct Int { int i; };

   using VT1 = blaze::StaticVector<double,3UL>;
   using VT2 = blaze::StaticVector<complex<float>,2UL>;
   using VT3 = blaze::StaticVector<Int,5UL>;

   alignof( VT1 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( VT2 );  // Evaluates to 16 for SSE, 32 for AVX, and 64 for AVX-512
   alignof( VT3 );  // Evaluates to 'alignof( Int )'
   \endcode

// For this reason \ref vector_types_static_vector, \ref vector_types_hybrid_vector,
// \ref matrix_types_static_matrix, and \ref matrix_types_hybrid_matrix cannot be used in
// containers using dynamic memory such as \c std::vector without additionally providing an
// allocator that can provide over-aligned memory:

   \code
   using Type = blaze::StaticVector<double,3UL>;
   using Allocator = blaze::AlignedAllocator<Type>;

   std::vector<Type> v1;  // Might be misaligned for AVX or AVX-512
   std::vector<Type,Allocator> v2;  // Properly aligned for AVX or AVX-512
   \endcode

// It is possible to disable the vectorization entirely by the compile time switch in the
// <tt>./blaze/config/Vectorization.h</tt> configuration file:

   \code
   #define BLAZE_USE_VECTORIZATION 1
   \endcode

// It is also possible to (de-)activate vectorization via command line or by defining this symbol
// manually before including any \b Blaze header file:

   \code
   g++ ... -DBLAZE_USE_VECTORIZATION=1 ...
   \endcode

   \code
   #define BLAZE_USE_VECTORIZATION 1
   #include <blaze/Blaze.h>
   \endcode

// In case the switch is set to 1, vectorization is enabled and the \b Blaze library is allowed
// to use intrinsics and the necessary alignment to speed up computations. In case the switch is
// set to 0, vectorization is disabled entirely and the \b Blaze library chooses default,
// non-vectorized functionality for the operations. Note that deactivating the vectorization may
// pose a severe performance limitation for a large number of operations!
//
//
// <hr>
// \section faq_blas To which extend does Blaze make use of BLAS functions under the hood?
//
// Currently the only BLAS functions that are utilized by \b Blaze are the \c gemm() functions
// for the multiplication of two dense matrices (i.e. \c sgemm(), \c dgemm(), \c cgemm(), and
// \c zgemm()). All other operations are always and unconditionally performed by native \b Blaze
// kernels.
//
// The \c BLAZE_BLAS_MODE config switch (see <tt>./blaze/config/BLAS.h</tt>) determines whether
// \b Blaze is allowed to use BLAS kernels. If \c BLAZE_BLAS_MODE is set to 0 then \b Blaze
// does not utilize the BLAS kernels and unconditionally uses its own custom kernels. If
// \c BLAZE_BLAS_MODE is set to 1 then \b Blaze is allowed to choose between using BLAS kernels
// or its own custom kernels. In case of the dense matrix multiplication this decision is based
// on the size of the dense matrices. For large matrices, \b Blaze uses the BLAS kernels, for
// small matrices it uses its own custom kernels. The threshold for this decision can be
// configured via the \c BLAZE_DMATDMATMULT_THRESHOLD, \c BLAZE_DMATTDMATMULT_THRESHOLD,
// \c BLAZE_TDMATDMATMULT_THRESHOLD and \c BLAZE_TDMATTDMATMULT_THRESHOLD config switches
// (see <tt>./blaze/config/Thresholds.h</tt>).
//
// Please note that the extend to which \b Blaze uses BLAS kernels can change in future releases
// of \b Blaze!
//
//
// <hr>
// \section faq_lapack To which extend does Blaze make use of LAPACK functions under the hood?
//
// \b Blaze uses LAPACK functions for matrix decomposition, matrix inversion, computing the
// determinants and eigenvalues, and the SVD. In contrast to the BLAS functionality (see
// \ref faq_blas), you cannot disable LAPACK or switch to custom kernels. In case you try to
// use any of these functionalities, but do not provide (i.e. link) a LAPACK library you will
// get link time errors.
//
// Please note that the extend to which \b Blaze uses LAPACK kernels can change in future releases
// of \b Blaze!
//
//
// <hr>
// \section faq_sparse_matrix_setup What is the fastest way to setup a very large sparse matrix?
//
// The following examples give an overview of different approaches to setup a sparse, row-major NxN
// matrix with the following pattern, where all values on the diagonal and the two sub-diagonals
// are filled:

   \f[\left(\begin{array}{*{9}{c}}
   1      & 1      & 0      & 0      & 0      & \cdots & 0      & 0      & 0      \\
   1      & 1      & 1      & 0      & 0      & \cdots & 0      & 0      & 0      \\
   0      & 1      & 1      & 1      & 0      & \cdots & 0      & 0      & 0      \\
   0      & 0      & 1      & 1      & 1      & \cdots & 0      & 0      & 0      \\
   0      & 0      & 0      & 1      & 1      & \cdots & 0      & 0      & 0      \\
   \vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
   0      & 0      & 0      & 0      & 0      & \cdots & 1      & 1      & 0      \\
   0      & 0      & 0      & 0      & 0      & \cdots & 1      & 1      & 1      \\
   0      & 0      & 0      & 0      & 0      & \cdots & 0      & 1      & 1      \\
   \end{array}\right)\f]

// Special emphasis is given to the runtime until the matrix setup is complete. In all cases the
// runtime is benchmarked with Clang-9.0 (compilation flags \c -O2 and \c -DNDEBUG) for \c N=200000.
//
//
// <b>Approach 1: Using the function call operator</b>
//
// In this approach the function call operator (i.e. \c operator()) is used to insert the according
// elements into the matrix:

   \code
   blaze::CompressedMatrix<int,rowMajor> A( N, N );

   A.reserve( N*3UL-2UL );  // Optional: Reserve capacity for all elements upfront

   for( size_t i=0; i<N; ++i ) {
      const size_t jbegin( i == 0UL ? 0UL : i-1UL );
      const size_t jend  ( i == N-1UL ? N-1UL : i+1UL );
      for( size_t j=jbegin; j<=jend; ++j ) {
         A(i,j) = 1;
      }
   }
   \endcode

// This approach is the most general and convenient, but also the slowest of all (approx. \b 64
// seconds). With every call to \c operator(), a new element is inserted at the specified position.
// This implies shifting all subsequent elements and adapting every subsequent row. Since all
// non-zero elements are stored in a single array inside a \c CompressedMatrix, this approach is
// similar to inserting elements at the front of a \c std::vector; all subsequent elements have
// to be shifted.
//
//
// <b>Approach 2: Rowwise reserve and insert</b>
//
// The next approach performs a rowwise reservation of capacity:

   \code
   blaze::CompressedMatrix<int,rowMajor> A( N, N );

   A.reserve( N*3UL );                // Allocate the total amount of memory
   A.reserve( 0UL, 2UL );             // Reserve a capacity of 2 for row 0
   for( size_t i=1; i<N-1UL; ++i ) {
      A.reserve( i, 3UL );            // Reserve a capacity of 3 for row i
   }
   A.reserve( N-1UL, 2UL );           // Reserve a capacity of 2 for the last row

   for( size_t i=0; i<N; ++i ) {
      const size_t jbegin( i == 0UL ? 0UL : i-1UL );
      const size_t jend  ( i == N-1UL ? N-1UL : i+1UL );
      for( size_t j=jbegin; j<=jend; ++j ) {
         A.insert( i, j, 1 );
      }
   }
   \endcode

// The first call to reserve() performs the memory allocation for the entire matrix. The complete
// matrix now holds the entire capacity, but each single row has a capacity of 0. Therefore the
// subsequent calls to \c reserve() divide the existing capacity to all rows.
//
// Unfortunately, also this approach is rather slow. The runtime is approx. \b 30 seconds. The
// downside of this approach is that changing the capacity of a single row causes a change in
// all following rows. Therefore this approach is similar to the first approach.
//
//
// <b>Approach 3: reserve/append/finalize</b>
//
// As the wiki explains, the most efficient way to fill a sparse matrix is a combination of
// \c reserve(), \c append() and \c finalize():

   \code
   CompressedMatrix<int,rowMajor> A( N, N );

   A.reserve( N*3UL );
   for( size_t i=0; i<N; ++i ) {
      const size_t jbegin( i == 0UL ? 0UL : i-1UL );
      const size_t jend  ( i == N-1UL ? N-1UL : i+1UL );
      for( size_t j=jbegin; j<=jend; ++j ) {
         A.append( i, j, 1 );
      }
      A.finalize( i );
   }
   \endcode

// The initial call to \c reserve() allocates enough memory for all non-zero elements of the
// entire matrix. \c append() and \c finalize() are then used to insert the elements and to mark
// the end of each single row. This is a very low-level approach and very similar to writing to
// an array manually, which results in a mere \b 0.026 seconds. The \c append() function writes
// the new element to the next memory location, and at the end of each row or column the
// \c finalize() function sets the internal pointers accordingly. It is very important to note
// that the \c finalize() function has to be explicitly called for each row, even for empty ones!
// Else the internal data structure will be corrupt! Also note that although \c append() does not
// allocate new memory, it still invalidates all iterators returned by the \c end() functions!
//
//
// <b>Approach 4: Reservation via the constructor</b>
//
// In case the number of non-zero elements is known upfront, it is also possible to perform the
// reservation via the constructor of \c CompressedMatrix. For that purpose \c CompressedMatrix
// provides a constructor taking a \c std::vector<size_t>:

   \code
   std::vector<size_t> nonzeros( N, 3UL );  // Create a vector of N elements with value 3
   nonzeros[  0] = 2UL;                     // We need only 2 elements in the first row ...
   nonzeros[N-1] = 2UL;                     //  ... and last row.

   CompressedMatrix<int,rowMajor> A( N, N, nonzeros );

   //std::cerr << " Inserting values...\n";
   for( size_t i=0; i<N; ++i ) {
      const size_t jbegin( i == 0UL ? 0UL : i-1UL );
      const size_t jend  ( i == N-1UL ? N-1UL : i+1UL );
      for( size_t j=jbegin; j<=jend; ++j ) {
         A.insert( i, j, 1 );
      }
   }
   \endcode

// The runtime for this approach is \b 0.027 seconds.
//
//
// <hr>
// \section faq_compile_times The compile time is too high if I include <blaze/Blaze.h>. Can I reduce it?
//
// The include file <tt><blaze/Blaze.h></tt> includes the entire functionality of the \b Blaze
// library, which by now is several hundred thousand lines of source code. That means that a lot
// of source code has to be parsed whenever <tt><blaze/Blaze.h></tt> is encountered. However, it
// is rare that everything is required within a single compilation unit. Therefore it is easily
// possible to reduce compile times by including only those \b Blaze features that are used within
// the compilation unit. For instance, instead of including <tt><blaze/Blaze.h></tt> it could be
// enough to include <tt><blaze/math/DynamicVector.h></tt>, which would reduce the compilation
// times by about 20%.
//
// Additionally we are taking care to implement new \b Blaze functionality such that compile times
// do not explode and try to reduce the compile times of existing features. Thus newer releases of
// \b Blaze can also improve compile times.
//
//
// <hr>
// \section faq_custom_operations Blaze does not provide feature XYZ. What can I do?
//
// In some cases you might be able to implement the required functionality very conveniently by
// building on the existing \c map() functions (see \ref custom_operations_map). For instance,
// the following code demonstrates the addition of a function that merges two vectors of floating
// point type into a vector of complex numbers:

   \code
   template< typename VT1, typename VT2, bool TF >
   decltype(auto) zip( const blaze::DenseVector<VT1,TF>& lhs, const blaze::DenseVector<VT2,TF>& rhs )
   {
      return blaze::map( ~lhs, ~rhs, []( const auto& r, const auto& i ) {
         using ET1 = ElementType_t<VT1>;
         using ET2 = ElementType_t<VT2>;
         return std::complex<std::common_type_t<ET1,ET2>>( r, i );
      } );
   }
   \endcode

// You will find a summary of the necessary steps to create custom features in \ref customization.
//
// Sometimes, however, the available customization points might not be sufficient. In this case
// you are cordially invited to create a pull request that provides the implementation of a
// feature or to create an issue according to our \ref issue_creation_guidelines. Please try
// to explain the feature as descriptive as possible, for instance by providing conceptual code
// examples.
//
// \n Previous: \ref intra_statement_optimization &nbsp; &nbsp; Next: \ref issue_creation_guidelines \n
*/
//*************************************************************************************************


//**FAQ********************************************************************************************
/*!\page issue_creation_guidelines Issue Creation Guidelines
//
// \tableofcontents
//
//
// One of the most important aspects of the \b Blaze project is the
// <a href="https://bitbucket.org/blaze-lib/blaze/issues">issue management</a> on the official
// \b Blaze Bitbucket page. We cordially invite all \b Blaze users to submit feature requests
// and bug reports, as we believe that this is a significant part of making \b Blaze a better
// library. However, we are asking to follow a small set of guidelines when creating an issue
// to facilitate the issue management on our side and also to make issues more useful for users
// of \b Blaze.
//
//
// <hr>
// \section issues_title Title
//
// The title is the most important detail of an issue. A well chosen title makes it easy to grasp
// the idea of an issue and improves the discoverability. Therefore, please choose a title that
// is ...
//
//  - ... as descriptive as possible;
//  - ... as concise as possible;
//  - ... as unambiguous as possible.
//
// Also, please create a separate issue for each idea/problem/etc. A very general title or an
// \"and\" in the title could be an indication that the issue is not specific enough and should
// be split into several issues.
//
// \subsection issues_title_good_examples Good Examples
//
//  - \"Provide support for AVX-512 SIMD operations\"
//  - \"Add support for the Boost Multiprecision Library\"
//  - \"Introduce reduction operations into Blaze\"
//  - \"Compilation error on KNL with -march=knl\"
//
// \subsection issues_title_bad_examples Bad Examples
//
//  - \"Several requests\" (instead create separate issues for each single request)
//  - \"Improve the performance\" (instead specify which operation should perform better)
//  - \"Blaze library compilation error\" (instead try to be more specific)
//
//
// <hr>
// \section issues_description Description
//
// The description should help us to understand your idea or problem in as much detail as possible.
// Also, it helps to clearly spell out your expectations (how a feature is supposed to work, how
// the behavior should be, etc.). Please spend a couple of minutes to try to make the description
// as comprehensive as possible.
//
//
// <hr>
// \section issues_assignee Assignee
//
// There is no need to assign the issue to a particular person. It is perfectly ok if you just
// ignore this setting.
//
//
// <hr>
// \section issues_kind Kind of Issue
//
// There are four kinds of issues available in the Bitbucket issue tracker: \ref issues_kind_bug,
// \ref issues_kind_enhancement, \ref issues_kind_proposal, and \ref issues_kind_task. In the
// following we try to give guidelines on which kind to choose for a particular issue:
//
// \subsection issues_kind_bug Bug
//
// Please choose the category \ref issues_kind_bug if ...
//
//  - ... you experience a compilation error despite your best efforts to get it right;
//  - ... you experience a crash/failure despite your best efforts to get it right;
//  - ... you experience problems when combining features;
//  - ... a feature does not work as specified/documented (i.e. can be considered broken).
//
// Please \b don't choose the category \ref issues_kind_bug if ...
//
//  - ... you feel a feature should work differently than it currently does (instead create a
//        \ref issues_kind_proposal with a convincing title and description);
//  - ... you are not sure how to use a feature (instead create an \ref issues_kind_enhancement
//        issue to extend the documentation);
//  - ... you are missing a feature (instead create a \ref issues_kind_proposal or
//        \ref issues_kind_enhancement issue).
//
// If you select the category \ref issues_kind_bug, please also try to provide a minimum example
// that fails. That helps us to minimize the time to resolve the bug.
//
// As we try to keep \b Blaze bug-free, we will always prioritize bug issues. However, we will
// also quickly close bug issues as \"wontfix\" if the described issue is not a bug (i.e. one of
// the problems mentioned above). We will \b not relabel a bug issue to \ref issues_kind_enhancement
// or \ref issues_kind_proposal, even if they would be reasonable extensions to \b Blaze.
//
// \subsection issues_kind_enhancement Enhancement
//
// Please choose the category \ref issues_kind_enhancement if ...
//
//  - ... you need an add-on to an existing feature;
//  - ... you need an extension of an existing feature;
//  - ... you need an extended documentation for an existing feature.
//
// \ref issues_kind_enhancement is very similar to \ref issues_kind_proposal, so we don't mind
// if an \ref issues_kind_enhancement is labeled as a \ref issues_kind_proposal or vice versa.
// Just make sure you don't request an extension or new feature as a \ref issues_kind_bug.
//
// \subsection issues_kind_proposal Proposal
//
// Please choose the category \ref issues_kind_proposal if ...
//
//  - ... you want to request a new feature;
//  - ... you want to change an existing feature.
//
// \ref issues_kind_proposal is very similar to \ref issues_kind_enhancement, so we don't mind if
// a \ref issues_kind_proposal is labeled as an \ref issues_kind_enhancement or vice versa. Just
// make sure you don't request an extension or new feature as a \ref issues_kind_bug.
//
// \subsection issues_kind_task Task
//
// Please choose the category \ref issues_kind_task if ...
//
//  - ... you want us to do something not feature related;
//  - ... you have something else in mind which does not fall in the other three categories.
//
//
// <hr>
// \section issues_priority Priority
//
// Via the priority of an issue you can tell us how important the issue is to you. Therefore the
// priority can have an influence on when we will deal with the issue. However, unfortunately we
// don't have an infinite amount of time and we can not deal with an arbitrary amount of issues
// at the same time. We will therefore take the priority into account, but mainly schedule the
// issues based on impact to all \b Blaze users and the estimated time to resolve it.
//
// You can choose between \ref issues_priority_blocker, \ref issues_priority_critical,
// \ref issues_priority_major, \ref issues_priority_minor, and \ref issues_priority_trivial.
//
// \subsection issues_priority_blocker Blocker
//
// Please choose a \ref issues_priority_blocker priority if ...
//
//  - ... you cannot work with \b Blaze due to the described \ref issues_kind_bug;
//  - ... the \ref issues_kind_bug likely has an influence on \b all \b Blaze users.
//
// Please note that the categories \ref issues_kind_enhancement or \ref issues_kind_proposal
// should never be a \ref issues_priority_blocker!
//
// \subsection issues_priority_critical Critical
//
// Please choose a \ref issues_priority_critical priority if ...
//
//  - ... you can work around a \ref issues_kind_bug, but the workaround is (much) slower or awful;
//  - ... you cannot use \b Blaze without the proposed feature;
//  - ... you consider it to be essential for \b all \b Blaze users.
//
// \subsection issues_priority_major Major
//
// Please choose a \ref issues_priority_major priority if ...
//
//  - ... a \ref issues_kind_bug or feature request is not \ref issues_priority_critical, but
//        still very important to you;
//  - ... you consider it to have a \ref issues_priority_major impact on most \b Blaze users.
//
// The \ref issues_priority_major category is the default setting in Bitbucket and we therefore
// consider it as the default priority for issues.
//
// \subsection issues_priority_minor Minor
//
// Please choose a \ref issues_priority_minor priority if ...
//
//  - ... a \ref issues_kind_bug does not affect many \b Blaze users;
//  - ... a feature request would only be useful for a small number of \b Blaze users;
//  - ... a feature would be nice to have, but is not particularly important.
//
// \subsection issues_priority_trivial Trivial
//
// Please choose a \ref issues_priority_trivial priority if ...
//
//  - ... a \ref issues_kind_bug hardly affects anyone;
//  - ... a feature request would only be useful for very few \b Blaze users;
//  - ... the expected time to resolve an issue is very small.
//
//
// <hr>
// \section issues_attachment Attachments
//
// You can always provide us with additional information in the form of attachments. Feel free
// to attach something to the issue if ...
//
//  - ... it can help us to analyze a \ref issues_kind_bug;
//  - ... you have some source code that demonstrates a problem;
//  - ... you already have a working prototype that sketches the idea;
//  - ... you have additional resources that could help us.
//
// We appreciate anything that simplifies our work and speeds up our progress.
//
// \n Previous: \ref faq &nbsp; &nbsp; Next: \ref blaze_references \n
*/
//*************************************************************************************************


//**Blaze References*******************************************************************************
/*!\page blaze_references Blaze References
//
// In case you need references to the \b Blaze library (for papers or other publications), please
// feel free to use one of the following references:

   \code
   @misc{blazelib,
      author       = "Klaus {Iglberger}",
      title        = "Blaze C++ Linear Algebra Library",
      howpublished = "https://bitbucket.org/blaze-lib",
      year         = 2012
   }
   \endcode

   \code
   @article{iglberger2012_1,
      author  = "Klaus {Iglberger} and Georg {Hager} and Jan {Treibig} and Ulrich {R{\"u}de}",
      title   = "Expression Templates Revisited: A Performance Analysis of Current Methodologies",
      journal = "SIAM Journal on Scientific Computing",
      year    = 2012,
      volume  = 34(2),
      pages   = C42--C69
   }
   \endcode

   \code
   @inproceedings{iglberger2012_2,
      author    = "Klaus {Iglberger} and Georg {Hager} and Jan {Treibig} and Ulrich {R{\"u}de}",
      title     = "High Performance Smart Expression Template Math Libraries",
      booktitle = "Proceedings of the 2nd International Workshop on New Algorithms and Programming Models for the Manycore Era (APMM 2012) at HPCS 2012",
      year      = 2012
   }
   \endcode

// \n Previous: \ref issue_creation_guidelines
*/
//*************************************************************************************************

#endif
