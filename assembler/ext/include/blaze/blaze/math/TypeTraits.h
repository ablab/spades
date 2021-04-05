//=================================================================================================
/*!
//  \file blaze/math/TypeTraits.h
//  \brief Header file for all type traits
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

#ifndef _BLAZE_MATH_TYPETRAITS_H_
#define _BLAZE_MATH_TYPETRAITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/CustomOppositeType.h>
#include <blaze/math/typetraits/CustomTransposeType.h>
#include <blaze/math/typetraits/DynamicAllocator.h>
#include <blaze/math/typetraits/GetAllocator.h>
#include <blaze/math/typetraits/HasAbs.h>
#include <blaze/math/typetraits/HasAcos.h>
#include <blaze/math/typetraits/HasAcosh.h>
#include <blaze/math/typetraits/HasAdd.h>
#include <blaze/math/typetraits/HasAsin.h>
#include <blaze/math/typetraits/HasAsinh.h>
#include <blaze/math/typetraits/HasAtan.h>
#include <blaze/math/typetraits/HasAtan2.h>
#include <blaze/math/typetraits/HasAtanh.h>
#include <blaze/math/typetraits/HasCbrt.h>
#include <blaze/math/typetraits/HasCeil.h>
#include <blaze/math/typetraits/HasClamp.h>
#include <blaze/math/typetraits/HasCompositeType.h>
#include <blaze/math/typetraits/HasConj.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasCos.h>
#include <blaze/math/typetraits/HasCosh.h>
#include <blaze/math/typetraits/HasDiv.h>
#include <blaze/math/typetraits/HasErf.h>
#include <blaze/math/typetraits/HasErfc.h>
#include <blaze/math/typetraits/HasExp.h>
#include <blaze/math/typetraits/HasExp2.h>
#include <blaze/math/typetraits/HasExp10.h>
#include <blaze/math/typetraits/HasFloor.h>
#include <blaze/math/typetraits/HasHypot.h>
#include <blaze/math/typetraits/HasImag.h>
#include <blaze/math/typetraits/HasInvCbrt.h>
#include <blaze/math/typetraits/HasInvSqrt.h>
#include <blaze/math/typetraits/HasLGamma.h>
#include <blaze/math/typetraits/HasLoad.h>
#include <blaze/math/typetraits/HasLog.h>
#include <blaze/math/typetraits/HasLog1p.h>
#include <blaze/math/typetraits/HasLog2.h>
#include <blaze/math/typetraits/HasLog10.h>
#include <blaze/math/typetraits/HasMax.h>
#include <blaze/math/typetraits/HasMin.h>
#include <blaze/math/typetraits/HasMult.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasPow.h>
#include <blaze/math/typetraits/HasReal.h>
#include <blaze/math/typetraits/HasResultType.h>
#include <blaze/math/typetraits/HasRound.h>
#include <blaze/math/typetraits/HasSign.h>
#include <blaze/math/typetraits/HasSIMDAbs.h>
#include <blaze/math/typetraits/HasSIMDAcos.h>
#include <blaze/math/typetraits/HasSIMDAcosh.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDAsin.h>
#include <blaze/math/typetraits/HasSIMDAsinh.h>
#include <blaze/math/typetraits/HasSIMDAtan.h>
#include <blaze/math/typetraits/HasSIMDAtan2.h>
#include <blaze/math/typetraits/HasSIMDAtanh.h>
#include <blaze/math/typetraits/HasSIMDBitand.h>
#include <blaze/math/typetraits/HasSIMDBitor.h>
#include <blaze/math/typetraits/HasSIMDBitxor.h>
#include <blaze/math/typetraits/HasSIMDCbrt.h>
#include <blaze/math/typetraits/HasSIMDCeil.h>
#include <blaze/math/typetraits/HasSIMDConj.h>
#include <blaze/math/typetraits/HasSIMDCos.h>
#include <blaze/math/typetraits/HasSIMDCosh.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDEqual.h>
#include <blaze/math/typetraits/HasSIMDErf.h>
#include <blaze/math/typetraits/HasSIMDErfc.h>
#include <blaze/math/typetraits/HasSIMDExp.h>
#include <blaze/math/typetraits/HasSIMDExp2.h>
#include <blaze/math/typetraits/HasSIMDExp10.h>
#include <blaze/math/typetraits/HasSIMDFloor.h>
#include <blaze/math/typetraits/HasSIMDHypot.h>
#include <blaze/math/typetraits/HasSIMDInvCbrt.h>
#include <blaze/math/typetraits/HasSIMDInvSqrt.h>
#include <blaze/math/typetraits/HasSIMDLGamma.h>
#include <blaze/math/typetraits/HasSIMDLog.h>
#include <blaze/math/typetraits/HasSIMDLog1p.h>
#include <blaze/math/typetraits/HasSIMDLog2.h>
#include <blaze/math/typetraits/HasSIMDLog10.h>
#include <blaze/math/typetraits/HasSIMDMax.h>
#include <blaze/math/typetraits/HasSIMDMin.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDPow.h>
#include <blaze/math/typetraits/HasSIMDRound.h>
#include <blaze/math/typetraits/HasSIMDShiftLI.h>
#include <blaze/math/typetraits/HasSIMDShiftLV.h>
#include <blaze/math/typetraits/HasSIMDShiftRI.h>
#include <blaze/math/typetraits/HasSIMDShiftRV.h>
#include <blaze/math/typetraits/HasSIMDSign.h>
#include <blaze/math/typetraits/HasSIMDSin.h>
#include <blaze/math/typetraits/HasSIMDSinh.h>
#include <blaze/math/typetraits/HasSIMDSqrt.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/HasSIMDTan.h>
#include <blaze/math/typetraits/HasSIMDTanh.h>
#include <blaze/math/typetraits/HasSIMDTrunc.h>
#include <blaze/math/typetraits/HasSin.h>
#include <blaze/math/typetraits/HasSinh.h>
#include <blaze/math/typetraits/HasSqrt.h>
#include <blaze/math/typetraits/HasSub.h>
#include <blaze/math/typetraits/HasTan.h>
#include <blaze/math/typetraits/HasTanh.h>
#include <blaze/math/typetraits/HasTrunc.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsAddExpr.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsBand.h>
#include <blaze/math/typetraits/IsBinaryMapExpr.h>
#include <blaze/math/typetraits/IsBLASCompatible.h>
#include <blaze/math/typetraits/IsColumn.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsColumns.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsCommutative.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsCrossExpr.h>
#include <blaze/math/typetraits/IsCUDAAssignable.h>
#include <blaze/math/typetraits/IsCustom.h>
#include <blaze/math/typetraits/IsDeclaration.h>
#include <blaze/math/typetraits/IsDeclDiagExpr.h>
#include <blaze/math/typetraits/IsDeclExpr.h>
#include <blaze/math/typetraits/IsDeclHermExpr.h>
#include <blaze/math/typetraits/IsDeclLowExpr.h>
#include <blaze/math/typetraits/IsDeclStrLowExpr.h>
#include <blaze/math/typetraits/IsDeclStrUppExpr.h>
#include <blaze/math/typetraits/IsDeclSymExpr.h>
#include <blaze/math/typetraits/IsDeclUniLowExpr.h>
#include <blaze/math/typetraits/IsDeclUniUppExpr.h>
#include <blaze/math/typetraits/IsDeclUppExpr.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsDivExpr.h>
#include <blaze/math/typetraits/IsEigenExpr.h>
#include <blaze/math/typetraits/IsElements.h>
#include <blaze/math/typetraits/IsEvalExpr.h>
#include <blaze/math/typetraits/IsExpandExpr.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsGeneral.h>
#include <blaze/math/typetraits/IsGenExpr.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsInitializer.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsKronExpr.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatEvalExpr.h>
#include <blaze/math/typetraits/IsMatExpExpr.h>
#include <blaze/math/typetraits/IsMatGenExpr.h>
#include <blaze/math/typetraits/IsMatInvExpr.h>
#include <blaze/math/typetraits/IsMatMapExpr.h>
#include <blaze/math/typetraits/IsMatMatAddExpr.h>
#include <blaze/math/typetraits/IsMatMatKronExpr.h>
#include <blaze/math/typetraits/IsMatMatMapExpr.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsMatMatSolveExpr.h>
#include <blaze/math/typetraits/IsMatMatSubExpr.h>
#include <blaze/math/typetraits/IsMatNoAliasExpr.h>
#include <blaze/math/typetraits/IsMatNoSIMDExpr.h>
#include <blaze/math/typetraits/IsMatReduceExpr.h>
#include <blaze/math/typetraits/IsMatRepeatExpr.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsMatScalarDivExpr.h>
#include <blaze/math/typetraits/IsMatScalarMultExpr.h>
#include <blaze/math/typetraits/IsMatSerialExpr.h>
#include <blaze/math/typetraits/IsMatTransExpr.h>
#include <blaze/math/typetraits/IsMatVecMultExpr.h>
#include <blaze/math/typetraits/IsMatVecSolveExpr.h>
#include <blaze/math/typetraits/IsModification.h>
#include <blaze/math/typetraits/IsMultExpr.h>
#include <blaze/math/typetraits/IsNoAliasExpr.h>
#include <blaze/math/typetraits/IsNoSIMDExpr.h>
#include <blaze/math/typetraits/IsOperation.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsPaddingEnabled.h>
#include <blaze/math/typetraits/IsProxy.h>
#include <blaze/math/typetraits/IsReduceExpr.h>
#include <blaze/math/typetraits/IsRepeatExpr.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRow.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRows.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSchurExpr.h>
#include <blaze/math/typetraits/IsSerialExpr.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/IsSIMDPack.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSolveExpr.h>
#include <blaze/math/typetraits/IsSparseElement.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStatic.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyTriangular.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSubExpr.h>
#include <blaze/math/typetraits/IsSubmatrix.h>
#include <blaze/math/typetraits/IsSubvector.h>
#include <blaze/math/typetraits/IsSVDExpr.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsTransformation.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsTVecMatMultExpr.h>
#include <blaze/math/typetraits/IsUnaryMapExpr.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsVecEvalExpr.h>
#include <blaze/math/typetraits/IsVecExpandExpr.h>
#include <blaze/math/typetraits/IsVecGenExpr.h>
#include <blaze/math/typetraits/IsVecMapExpr.h>
#include <blaze/math/typetraits/IsVecNoAliasExpr.h>
#include <blaze/math/typetraits/IsVecNoSIMDExpr.h>
#include <blaze/math/typetraits/IsVecRepeatExpr.h>
#include <blaze/math/typetraits/IsVecScalarDivExpr.h>
#include <blaze/math/typetraits/IsVecScalarMultExpr.h>
#include <blaze/math/typetraits/IsVecSerialExpr.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/IsVecTransExpr.h>
#include <blaze/math/typetraits/IsVecTVecMapExpr.h>
#include <blaze/math/typetraits/IsVecTVecMultExpr.h>
#include <blaze/math/typetraits/IsVecVecAddExpr.h>
#include <blaze/math/typetraits/IsVecVecDivExpr.h>
#include <blaze/math/typetraits/IsVecVecKronExpr.h>
#include <blaze/math/typetraits/IsVecVecMapExpr.h>
#include <blaze/math/typetraits/IsVecVecMultExpr.h>
#include <blaze/math/typetraits/IsVecVecSubExpr.h>
#include <blaze/math/typetraits/IsView.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/MakeComplex.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/typetraits/YieldsDiagonal.h>
#include <blaze/math/typetraits/YieldsHermitian.h>
#include <blaze/math/typetraits/YieldsIdentity.h>
#include <blaze/math/typetraits/YieldsLower.h>
#include <blaze/math/typetraits/YieldsStrictlyLower.h>
#include <blaze/math/typetraits/YieldsStrictlyTriangular.h>
#include <blaze/math/typetraits/YieldsStrictlyUpper.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsTriangular.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsUniLower.h>
#include <blaze/math/typetraits/YieldsUniTriangular.h>
#include <blaze/math/typetraits/YieldsUniUpper.h>
#include <blaze/math/typetraits/YieldsUpper.h>
#include <blaze/math/typetraits/YieldsZero.h>

#endif
