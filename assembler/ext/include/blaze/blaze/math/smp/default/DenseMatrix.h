//=================================================================================================
/*!
//  \file blaze/math/smp/default/DenseMatrix.h
//  \brief Header file for the default dense matrix SMP implementation
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

#ifndef _BLAZE_MATH_SMP_DEFAULT_DENSEMATRIX_H_
#define _BLAZE_MATH_SMP_DEFAULT_DENSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/system/SMP.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Dense matrix SMP functions */
//@{
template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpAddAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpSubAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpSchurAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the SMP assignment of a matrix to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be assigned.
// \return void
//
// This function implements the default SMP assignment of a matrix to a dense matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   assign( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the SMP addition assignment of a matrix to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function implements the default SMP addition assignment of a matrix to a dense matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpAddAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   addAssign( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the SMP subtraction assignment of a matrix to dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function implements the default SMP subtraction assignment of a matrix to a dense matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpSubAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   subAssign( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the SMP Schur product assignment of a matrix to dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix for the Schur product.
// \return void
//
// This function implements the default SMP Schur product assignment of a matrix to a dense
// matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpSchurAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   schurAssign( *lhs, *rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( !BLAZE_HPX_PARALLEL_MODE           );
BLAZE_STATIC_ASSERT( !BLAZE_CPP_THREADS_PARALLEL_MODE   );
BLAZE_STATIC_ASSERT( !BLAZE_BOOST_THREADS_PARALLEL_MODE );
BLAZE_STATIC_ASSERT( !BLAZE_OPENMP_PARALLEL_MODE        );

}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
