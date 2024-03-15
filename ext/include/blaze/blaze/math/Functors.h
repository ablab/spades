//=================================================================================================
/*!
//  \file blaze/math/Functors.h
//  \brief Header file for all functors
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

#ifndef _BLAZE_MATH_FUNCTORS_H_
#define _BLAZE_MATH_FUNCTORS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/functors/Abs.h>
#include <blaze/math/functors/Acos.h>
#include <blaze/math/functors/Acosh.h>
#include <blaze/math/functors/Add.h>
#include <blaze/math/functors/AddAssign.h>
#include <blaze/math/functors/And.h>
#include <blaze/math/functors/AndAssign.h>
#include <blaze/math/functors/Arg.h>
#include <blaze/math/functors/Asin.h>
#include <blaze/math/functors/Asinh.h>
#include <blaze/math/functors/Assign.h>
#include <blaze/math/functors/Atan.h>
#include <blaze/math/functors/Atan2.h>
#include <blaze/math/functors/Atanh.h>
#include <blaze/math/functors/Bind1st.h>
#include <blaze/math/functors/Bind2nd.h>
#include <blaze/math/functors/Bind3rd.h>
#include <blaze/math/functors/Bitand.h>
#include <blaze/math/functors/Bitor.h>
#include <blaze/math/functors/Bitxor.h>
#include <blaze/math/functors/Cbrt.h>
#include <blaze/math/functors/Ceil.h>
#include <blaze/math/functors/Clamp.h>
#include <blaze/math/functors/Clear.h>
#include <blaze/math/functors/Conj.h>
#include <blaze/math/functors/Cos.h>
#include <blaze/math/functors/Cosh.h>
#include <blaze/math/functors/CTrans.h>
#include <blaze/math/functors/DeclDiag.h>
#include <blaze/math/functors/DeclHerm.h>
#include <blaze/math/functors/DeclId.h>
#include <blaze/math/functors/DeclLow.h>
#include <blaze/math/functors/DeclStrLow.h>
#include <blaze/math/functors/DeclStrUpp.h>
#include <blaze/math/functors/DeclSym.h>
#include <blaze/math/functors/DeclUniLow.h>
#include <blaze/math/functors/DeclUniUpp.h>
#include <blaze/math/functors/DeclUpp.h>
#include <blaze/math/functors/DeclZero.h>
#include <blaze/math/functors/Div.h>
#include <blaze/math/functors/DivAssign.h>
#include <blaze/math/functors/Erf.h>
#include <blaze/math/functors/Erfc.h>
#include <blaze/math/functors/Eval.h>
#include <blaze/math/functors/Exp.h>
#include <blaze/math/functors/Exp2.h>
#include <blaze/math/functors/Exp10.h>
#include <blaze/math/functors/Floor.h>
#include <blaze/math/functors/Greater.h>
#include <blaze/math/functors/Hypot.h>
#include <blaze/math/functors/Imag.h>
#include <blaze/math/functors/Inv.h>
#include <blaze/math/functors/InvAdd.h>
#include <blaze/math/functors/InvCbrt.h>
#include <blaze/math/functors/InvSqrt.h>
#include <blaze/math/functors/Join.h>
#include <blaze/math/functors/Kron.h>
#include <blaze/math/functors/L1Norm.h>
#include <blaze/math/functors/L2Norm.h>
#include <blaze/math/functors/L3Norm.h>
#include <blaze/math/functors/L4Norm.h>
#include <blaze/math/functors/LeftShiftAssign.h>
#include <blaze/math/functors/Less.h>
#include <blaze/math/functors/LGamma.h>
#include <blaze/math/functors/Log.h>
#include <blaze/math/functors/Log1p.h>
#include <blaze/math/functors/Log2.h>
#include <blaze/math/functors/Log10.h>
#include <blaze/math/functors/LpNorm.h>
#include <blaze/math/functors/MAC.h>
#include <blaze/math/functors/MakePair.h>
#include <blaze/math/functors/Max.h>
#include <blaze/math/functors/Min.h>
#include <blaze/math/functors/Minmax.h>
#include <blaze/math/functors/ModuloAssign.h>
#include <blaze/math/functors/Mult.h>
#include <blaze/math/functors/MultAssign.h>
#include <blaze/math/functors/NoAlias.h>
#include <blaze/math/functors/Noop.h>
#include <blaze/math/functors/NoSIMD.h>
#include <blaze/math/functors/Not.h>
#include <blaze/math/functors/Or.h>
#include <blaze/math/functors/OrAssign.h>
#include <blaze/math/functors/Pow.h>
#include <blaze/math/functors/Pow2.h>
#include <blaze/math/functors/Pow3.h>
#include <blaze/math/functors/Pow4.h>
#include <blaze/math/functors/Qdrt.h>
#include <blaze/math/functors/Real.h>
#include <blaze/math/functors/Reset.h>
#include <blaze/math/functors/RightShiftAssign.h>
#include <blaze/math/functors/Round.h>
#include <blaze/math/functors/Schur.h>
#include <blaze/math/functors/Serial.h>
#include <blaze/math/functors/ShiftLI.h>
#include <blaze/math/functors/ShiftLV.h>
#include <blaze/math/functors/ShiftRI.h>
#include <blaze/math/functors/ShiftRV.h>
#include <blaze/math/functors/Sign.h>
#include <blaze/math/functors/Sin.h>
#include <blaze/math/functors/Sinh.h>
#include <blaze/math/functors/SqrAbs.h>
#include <blaze/math/functors/Sqrt.h>
#include <blaze/math/functors/Sub.h>
#include <blaze/math/functors/SubAssign.h>
#include <blaze/math/functors/Tan.h>
#include <blaze/math/functors/Tanh.h>
#include <blaze/math/functors/Trans.h>
#include <blaze/math/functors/Trunc.h>
#include <blaze/math/functors/XorAssign.h>

#endif
