//=================================================================================================
/*!
//  \file blaze/math/InversionFlag.h
//  \brief Header file for the dense matrix inversion flags
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

#ifndef _BLAZE_MATH_INVERSIONFLAG_H_
#define _BLAZE_MATH_INVERSIONFLAG_H_


namespace blaze {

//=================================================================================================
//
//  INVERSION FLAG VALUES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Inversion flag.
// \ingroup math
//
// The InversionFlag type enumeration represents the different types of matrix inversion algorithms
// that are available within the Blaze library. The following flags are available:
//
//  - \c byLU: The default inversion algorithm for general square matrices. It uses the LU
//          algorithm to decompose a matrix into a lower unitriangular matrix \c L, an upper
//          triangular matrix \c U, and a permutation matrix \c P (\f$ A = P L U \f$). If no
//          permutations are required, \c P is the identity matrix.
//  - \c byLDLT: The Bunch-Kaufman inversion algorithm for symmetric indefinite matrices. It
//          decomposes the given matrix into either \f$ A = U D U^{T} \f$ or \f$ A = L D L^{T} \f$,
//          where \c U (or \c L) is a product of permutation and unit upper (lower) triangular
//          matrices, and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal
//          blocks.
//  - \c byLDLH: The Bunch-Kaufman inversion algorithm for Hermitian indefinite matrices. It
//          decomposes the given matrix into either \f$ A = U D U^{H} \f$ or \f$ A = L D L^{H} \f$,
//          where \c U (or \c L) is a product of permutation and unit upper (lower) triangular
//          matrices, and \c D is Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal
//          blocks.
//  - \c byLLH: The Cholesky inversion algorithm for Hermitian positive definite matrices. It
//          decomposes a given matrix into either \f$ A = L L^H \f$, where \c L is a lower
//          triangular matrix, or \f$ A = U^H U \f$, where \c U is an upper triangular matrix.
//
// Alternatively, the type of the matrix can be specified, leaving it to the Blaze library to
// select the appropriate matrix inversion algorithm. The following flags are available:
//
//  - \c asGeneral: This flag selects the best suited inversion algorithm for general square
//          matrices. In case no further compile time information is available, this will imply
//          the use of the LU decomposition algorithm (see the \c byLU flag).
//  - \c asSymmetric: This flag selects the most suited inversion algorithm for symmetric matrices.
//          In case no further compile time information is available, the Bunch-Kaufman matrix
//          decomposition algorithm will be used (see the \c byLDLT flag).
//  - \c asHermitian: This flag selects the most suited inversion algorithm for Hermitian matrices.
//          In case no further compile time information is available, the Bunch-Kaufman matrix
//          decomposition algorithm will be used (see the \c byLDLH flag).
//  - \c asLower: This flag selects the most suited inversion algorithm for lower triangular
//          matrices. In case no further compile time information is available, the inversion will
//          be performed by a forward substitution. No matrix decomposition will be performed.
//  - \c asUniLower: This flag selects the most suited inversion algorithm for lower unitriangular
//          matrices. In case no further compile time information is available, the inversion will
//          be performed by a forward substitution. No matrix decomposition will be performed.
//  - \c asUpper: This flag selects the most suited inversion algorithm for upper triangular
//          matrices. In case no further compile time information is available, the inversion will
//          be performed by a back substitution. No matrix decomposition will be performed.
//  - \c asUniUpper: This flag selects the most suited inversion algorithm for upper unitriangular
//          matrices. In case no further compile time information is available, the inversion will
//          be performed by a back substitution. No matrix decomposition will be performed.
//  - \c asDiagonal: This flag selects the most suited inversion algorithm for diagonal matrices.
//          In case no further compile time information is available, the inversion will be
//          performed by directly computing the reciprocal of each diagonal element. No matrix
//          decomposition will be performed.
*/
enum InversionFlag
{
   byLU        =  0,  //!< Flag for the LU-based matrix inversion.
   byLDLT      =  1,  //!< Flag for the Bunch-Kaufman-based inversion for symmetric matrices.
   byLDLH      =  2,  //!< Flag for the Bunch-Kaufman-based inversion for Hermitian matrices.
   byLLH       =  3,  //!< Flag for the Cholesky-based inversion for positive-definite matrices.

   asGeneral   =  4,  //!< Flag for the inversion of a general matrix (same as byLU).
   asSymmetric =  5,  //!< Flag for the inversion of a symmetric matrix (same as byLDLT).
   asHermitian =  6,  //!< Flag for the inversion of a Hermitian matrix (same as byLDLH).
   asLower     =  7,  //!< Flag for the inversion of a lower triangular matrix.
   asUniLower  =  8,  //!< Flag for the inversion of a lower unitriangular matrix.
   asUpper     =  9,  //!< Flag for the inversion of a upper triangular matrix.
   asUniUpper  = 10,  //!< Flag for the inversion of a upper unitriangular matrix.
   asDiagonal  = 11   //!< Flag for the inversion of a diagonal matrix.
};
//*************************************************************************************************

} // namespace blaze

#endif
