//=================================================================================================
/*!
//  \file blaze/math/functors/Forward.h
//  \brief Header file for all functor forward declarations
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

#ifndef _BLAZE_MATH_FUNCTORS_FORWARD_H_
#define _BLAZE_MATH_FUNCTORS_FORWARD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

struct Abs;
struct Acos;
struct Acosh;
struct Add;
struct AddAssign;
struct And;
struct AndAssign;
struct Arg;
struct Asin;
struct Asinh;
struct Assign;
struct Atan;
struct Atan2;
struct Atanh;
template< typename, typename > struct Bind1st;
template< typename, typename > struct Bind2nd;
template< typename, typename > struct Bind3rd;
struct Bitand;
struct Bitor;
struct Bitxor;
struct Cbrt;
struct Ceil;
struct Clamp;
struct Clear;
struct Conj;
struct Cos;
struct Cosh;
struct CTrans;
struct DeclDiag;
struct DeclHerm;
struct DeclId;
struct DeclLow;
struct DeclStrLow;
struct DeclStrUpp;
struct DeclSym;
struct DeclUniLow;
struct DeclUniUpp;
struct DeclUpp;
struct DeclZero;
struct Div;
struct DivAssign;
struct Erf;
struct Erfc;
struct Eval;
struct Exp;
struct Exp2;
struct Exp10;
struct Floor;
struct Greater;
struct Hypot;
struct Imag;
struct Inv;
struct InvAdd;
struct InvCbrt;
struct InvSqrt;
template< typename > struct Join;
struct Kron;
struct L1Norm;
struct L2Norm;
struct L3Norm;
struct L4Norm;
struct LeftShiftAssign;
struct Less;
struct LGamma;
struct Log;
struct Log1p;
struct Log2;
struct Log10;
template< size_t... > struct LpNorm;
struct MAC;
struct MakePair;
struct Max;
struct Min;
struct Minmax;
struct ModuloAssign;
struct Mult;
struct MultAssign;
struct NoAlias;
struct Noop;
struct NoSIMD;
struct Not;
struct Or;
struct OrAssign;
struct Pow;
struct Pow2;
struct Pow3;
struct Pow4;
struct Qdrt;
struct Real;
struct Reset;
struct RightShiftAssign;
struct Round;
struct Schur;
struct Serial;
struct ShiftLI;
struct ShiftLV;
struct ShiftRI;
struct ShiftRV;
struct Sign;
struct Sin;
struct Sinh;
struct SqrAbs;
struct Sqrt;
struct Sub;
struct SubAssign;
struct Tan;
struct Tanh;
struct Trans;
struct Trunc;
struct XorAssign;

} // namespace blaze

#endif
