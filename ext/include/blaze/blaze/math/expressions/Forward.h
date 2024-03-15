//=================================================================================================
/*!
//  \file blaze/math/expressions/Forward.h
//  \brief Header file for all forward declarations for expression class templates
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

#ifndef _BLAZE_MATH_EXPRESSIONS_FORWARD_H_
#define _BLAZE_MATH_EXPRESSIONS_FORWARD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/ReductionFlag.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename > struct AddExpr;
template< typename > struct BinaryMapExpr;
template< typename > struct CrossExpr;
template< typename > struct DeclDiagExpr;
template< typename > struct DeclExpr;
template< typename > struct DeclHermExpr;
template< typename > struct DeclLowExpr;
template< typename > struct DeclStrLowExpr;
template< typename > struct DeclStrUppExpr;
template< typename > struct DeclSymExpr;
template< typename > struct DeclUniLowExpr;
template< typename > struct DeclUniUppExpr;
template< typename > struct DeclUppExpr;
template< typename > struct DivExpr;
template< typename, bool > class DenseMatrix;
template< typename, bool > class DenseVector;
template< typename, bool > class DMatDeclDiagExpr;
template< typename, bool > class DMatDeclHermExpr;
template< typename, bool > class DMatDeclLowExpr;
template< typename, bool > class DMatDeclStrLowExpr;
template< typename, bool > class DMatDeclStrUppExpr;
template< typename, bool > class DMatDeclSymExpr;
template< typename, bool > class DMatDeclUniLowExpr;
template< typename, bool > class DMatDeclUniUppExpr;
template< typename, bool > class DMatDeclUppExpr;
template< typename, typename, bool > class DMatDMatAddExpr;
template< typename, typename, bool > class DMatDMatKronExpr;
template< typename, typename, typename, bool > class DMatDMatMapExpr;
template< typename, typename, bool, bool, bool, bool > class DMatDMatMultExpr;
template< typename, typename, bool > class DMatDMatSchurExpr;
template< typename, typename, bool > class DMatDMatSolveExpr;
template< typename, typename, bool > class DMatDMatSubExpr;
template< typename, typename > class DMatDVecMultExpr;
template< typename, typename, bool > class DMatDVecSolveExpr;
template< typename, bool > class DMatEigenExpr;
template< typename, bool > class DMatEvalExpr;
template< typename, bool > class DMatExpExpr;
template< typename, typename, bool > class DMatGenExpr;
template< typename, bool > class DMatInvExpr;
template< typename, typename, bool > class DMatMapExpr;
template< typename, bool > class DMatNoAliasExpr;
template< typename, bool > class DMatNoSIMDExpr;
template< typename, typename, ReductionFlag > class DMatReduceExpr;
template< typename, bool, size_t... > class DMatRepeatExpr;
template< typename, typename, bool > class DMatScalarDivExpr;
template< typename, typename, bool > class DMatScalarMultExpr;
template< typename, bool > class DMatSerialExpr;
template< typename, typename, bool > class DMatSMatAddExpr;
template< typename, typename, bool > class DMatSMatKronExpr;
template< typename, typename, bool, bool, bool, bool > class DMatSMatMultExpr;
template< typename, typename > class DMatSMatSchurExpr;
template< typename, typename, bool > class DMatSMatSubExpr;
template< typename, bool > class DMatSVDExpr;
template< typename, typename > class DMatSVecMultExpr;
template< typename, typename > class DMatTDMatAddExpr;
template< typename, typename, typename > class DMatTDMatMapExpr;
template< typename, typename, bool, bool, bool, bool > class DMatTDMatMultExpr;
template< typename, typename > class DMatTDMatSchurExpr;
template< typename, typename > class DMatTDMatSubExpr;
template< typename, bool > class DMatTransExpr;
template< typename, bool > class DMatTransposer;
template< typename, typename > class DMatTSMatAddExpr;
template< typename, typename, bool, bool, bool, bool > class DMatTSMatMultExpr;
template< typename, typename > class DMatTSMatSchurExpr;
template< typename, typename > class DMatTSMatSubExpr;
template< typename, typename, bool > class DVecDVecAddExpr;
template< typename, typename, bool > class DVecDVecCrossExpr;
template< typename, typename, bool > class DVecDVecDivExpr;
template< typename, typename, bool > class DVecDVecKronExpr;
template< typename, typename, typename, bool > class DVecDVecMapExpr;
template< typename, typename, bool > class DVecDVecMultExpr;
template< typename, typename,typename > class DVecDVecOuterExpr;
template< typename, typename, bool > class DVecDVecSubExpr;
template< typename, bool > class DVecEvalExpr;
template< typename, bool, size_t... > class DVecExpandExpr;
template< typename, typename, bool > class DVecGenExpr;
template< typename, typename, bool > class DVecMapExpr;
template< typename, bool > class DVecNoAliasExpr;
template< typename, bool > class DVecNoSIMDExpr;
template< typename, bool, size_t... > class DVecRepeatExpr;
template< typename, typename, bool > class DVecScalarDivExpr;
template< typename, typename, bool > class DVecScalarMultExpr;
template< typename, bool > class DVecSerialExpr;
template< typename, typename, bool > class DVecSVecAddExpr;
template< typename, typename, bool > class DVecSVecCrossExpr;
template< typename, typename, bool > class DVecSVecKronExpr;
template< typename, typename, bool > class DVecSVecMultExpr;
template< typename, typename > class DVecSVecOuterExpr;
template< typename, typename, bool > class DVecSVecSubExpr;
template< typename, bool > class DVecTransExpr;
template< typename, bool > class DVecTransposer;
template< typename > struct EigenExpr;
template< typename > struct EvalExpr;
template< typename > struct ExpandExpr;
template< typename > struct Expression;
template< typename > struct GenExpr;
template< typename > struct KronExpr;
template< typename > struct MatEvalExpr;
template< typename > struct MatExpExpr;
template< typename > struct MatGenExpr;
template< typename > struct MatInvExpr;
template< typename > struct MatMapExpr;
template< typename > struct MatMatAddExpr;
template< typename > struct MatMatKronExpr;
template< typename > struct MatMatMapExpr;
template< typename > struct MatMatMultExpr;
template< typename > struct MatMatSolveExpr;
template< typename > struct MatMatSubExpr;
template< typename > struct MatNoAliasExpr;
template< typename > struct MatNoSIMDExpr;
template< typename, ReductionFlag > struct MatReduceExpr;
template< typename, size_t... > struct MatRepeatExpr;
template< typename, bool > class Matrix;
template< typename > struct MatScalarDivExpr;
template< typename > struct MatScalarMultExpr;
template< typename > struct MatSerialExpr;
template< typename > struct MatTransExpr;
template< typename > struct MatVecMultExpr;
template< typename > struct MatVecSolveExpr;
template< typename > struct MultExpr;
template< typename > struct NoAliasExpr;
template< typename > struct NoSIMDExpr;
template< typename > struct ReduceExpr;
template< typename > struct RepeatExpr;
template< typename > struct SchurExpr;
template< typename > struct SerialExpr;
template< typename, bool > class SMatDeclDiagExpr;
template< typename, bool > class SMatDeclHermExpr;
template< typename, bool > class SMatDeclLowExpr;
template< typename, bool > class SMatDeclStrLowExpr;
template< typename, bool > class SMatDeclStrUppExpr;
template< typename, bool > class SMatDeclSymExpr;
template< typename, bool > class SMatDeclUniLowExpr;
template< typename, bool > class SMatDeclUniUppExpr;
template< typename, bool > class SMatDeclUppExpr;
template< typename, typename, bool > class SMatDMatKronExpr;
template< typename, typename, bool, bool, bool, bool > class SMatDMatMultExpr;
template< typename, typename > class SMatDMatSchurExpr;
template< typename, typename, bool > class SMatDMatSubExpr;
template< typename, typename > class SMatDVecMultExpr;
template< typename, bool > class SMatEvalExpr;
template< typename, typename, bool > class SMatMapExpr;
template< typename, bool > class SMatNoAliasExpr;
template< typename, typename, ReductionFlag > class SMatReduceExpr;
template< typename, bool, size_t... > class SMatRepeatExpr;
template< typename, typename, bool > class SMatScalarDivExpr;
template< typename, typename, bool > class SMatScalarMultExpr;
template< typename, bool > class SMatSerialExpr;
template< typename, typename > class SMatSMatAddExpr;
template< typename, typename > class SMatSMatKronExpr;
template< typename, typename > class SMatSMatMultExpr;
template< typename, typename > class SMatSMatSchurExpr;
template< typename, typename > class SMatSMatSubExpr;
template< typename, typename > class SMatSVecMultExpr;
template< typename, typename, bool, bool, bool, bool > class SMatTDMatMultExpr;
template< typename, typename > class SMatTDMatSubExpr;
template< typename, bool > class SMatTransExpr;
template< typename, bool > class SMatTransposer;
template< typename, typename > class SMatTSMatAddExpr;
template< typename, typename > class SMatTSMatKronExpr;
template< typename, typename > class SMatTSMatMultExpr;
template< typename, typename > class SMatTSMatSchurExpr;
template< typename, typename > class SMatTSMatSubExpr;
template< typename, ReductionFlag > class SMatVarExpr;
template< typename > struct SolveExpr;
template< typename, bool > class SparseMatrix;
template< typename, bool > class SparseVector;
template< typename > struct SubExpr;
template< typename > struct SVDExpr;
template< typename, typename, bool > class SVecDVecCrossExpr;
template< typename, typename, bool > class SVecDVecDivExpr;
template< typename, typename, bool > class SVecDVecKronExpr;
template< typename, typename, bool > class SVecDVecMultExpr;
template< typename, typename > class SVecDVecOuterExpr;
template< typename, typename, bool > class SVecDVecSubExpr;
template< typename, bool > class SVecEvalExpr;
template< typename, bool, size_t... > class SVecExpandExpr;
template< typename, typename, bool > class SVecMapExpr;
template< typename, bool > class SVecNoAliasExpr;
template< typename, bool, size_t... > class SVecRepeatExpr;
template< typename, typename, bool > class SVecScalarDivExpr;
template< typename, typename, bool > class SVecScalarMultExpr;
template< typename, bool > class SVecSerialExpr;
template< typename, typename, bool > class SVecSVecAddExpr;
template< typename, typename, bool > class SVecSVecCrossExpr;
template< typename, typename, bool > class SVecSVecKronExpr;
template< typename, typename, bool > class SVecSVecMultExpr;
template< typename, typename > class SVecSVecOuterExpr;
template< typename, typename, bool > class SVecSVecSubExpr;
template< typename, bool > class SVecTransExpr;
template< typename, bool > class SVecTransposer;
template< typename, typename, bool, bool, bool, bool > class TDMatDMatMultExpr;
template< typename, typename > class TDMatDVecMultExpr;
template< typename, typename > class TDMatSMatAddExpr;
template< typename, typename, bool, bool, bool, bool > class TDMatSMatMultExpr;
template< typename, typename > class TDMatSMatSubExpr;
template< typename, typename > class TDMatSVecMultExpr;
template< typename, typename, bool, bool, bool, bool > class TDMatTDMatMultExpr;
template< typename, typename, bool, bool, bool, bool > class TDMatTSMatMultExpr;
template< typename, typename > class TDVecDMatMultExpr;
template< typename, typename > class TDVecSMatMultExpr;
template< typename, typename > class TDVecTDMatMultExpr;
template< typename, typename > class TDVecTSMatMultExpr;
template< typename > struct TransExpr;
template< typename, typename, bool, bool, bool, bool > class TSMatDMatMultExpr;
template< typename, typename > class TSMatDMatSchurExpr;
template< typename, typename > class TSMatDMatSubExpr;
template< typename, typename > class TSMatDVecMultExpr;
template< typename, typename > class TSMatSMatKronExpr;
template< typename, typename > class TSMatSMatMultExpr;
template< typename, typename > class TSMatSMatSchurExpr;
template< typename, typename > class TSMatSMatSubExpr;
template< typename, typename > class TSMatSVecMultExpr;
template< typename, typename, bool, bool, bool, bool > class TSMatTDMatMultExpr;
template< typename, typename > class TSMatTSMatAddExpr;
template< typename, typename > class TSMatTSMatKronExpr;
template< typename, typename > class TSMatTSMatMultExpr;
template< typename, typename > class TSMatTSMatSchurExpr;
template< typename, typename > class TSMatTSMatSubExpr;
template< typename, typename > class TSVecDMatMultExpr;
template< typename, typename > class TSVecSMatMultExpr;
template< typename, typename > class TSVecTDMatMultExpr;
template< typename, typename > class TSVecTSMatMultExpr;
template< typename > struct TVecMatMultExpr;
template< typename > struct UnaryMapExpr;
template< typename > struct VecEvalExpr;
template< typename, size_t... > struct VecExpandExpr;
template< typename > struct VecGenExpr;
template< typename > struct VecMapExpr;
template< typename > struct VecNoAliasExpr;
template< typename > struct VecNoSIMDExpr;
template< typename, size_t... > struct VecRepeatExpr;
template< typename > struct VecScalarDivExpr;
template< typename > struct VecScalarMultExpr;
template< typename > struct VecSerialExpr;
template< typename, bool > class Vector;
template< typename > struct VecTransExpr;
template< typename > struct VecTVecMapExpr;
template< typename > struct VecTVecMultExpr;
template< typename > struct VecVecAddExpr;
template< typename > struct VecVecDivExpr;
template< typename > struct VecVecKronExpr;
template< typename > struct VecVecMapExpr;
template< typename > struct VecVecMultExpr;
template< typename > struct VecVecSubExpr;
template< typename > struct View;


template< typename VT1, typename VT2, bool TF >
decltype(auto) operator+( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator+( const DenseVector<VT1,false>&, const DenseVector<VT2,true>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator+( const DenseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator+( const SparseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator+( const SparseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename MT1, typename MT2, bool SO >
decltype(auto) operator+( const DenseMatrix<MT1,SO>&, const DenseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const DenseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const DenseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2, bool SO >
decltype(auto) operator+( const DenseMatrix<MT1,SO>&, const SparseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const DenseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const DenseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2, bool SO >
decltype(auto) operator+( const SparseMatrix<MT1,SO>&, const DenseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const SparseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const SparseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator+( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );


template< typename VT1, typename VT2, bool TF >
decltype(auto) operator-( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator-( const DenseVector<VT1,false>&, const DenseVector<VT2,true>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator-( const DenseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator-( const SparseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator-( const SparseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename MT1, typename MT2, bool SO >
decltype(auto) operator-( const DenseMatrix<MT1,SO>&, const DenseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const DenseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const DenseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2, bool SO >
decltype(auto) operator-( const DenseMatrix<MT1,SO>&, const SparseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const DenseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const DenseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2, bool SO >
decltype(auto) operator-( const SparseMatrix<MT1,SO>&, const DenseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const SparseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const SparseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator-( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );


template< typename VT1, typename VT2, bool TF >
decltype(auto) operator*( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const DenseVector<VT1,true>&, const DenseVector<VT2,false>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const DenseVector<VT1,false>&, const DenseVector<VT2,true>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator*( const DenseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const DenseVector<VT1,true>&, const SparseVector<VT2,false>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const DenseVector<VT1,false>&, const SparseVector<VT2,true>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator*( const SparseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const SparseVector<VT1,true>&, const DenseVector<VT2,false>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const SparseVector<VT1,false>&, const DenseVector<VT2,true>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator*( const SparseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const SparseVector<VT1,true>&, const SparseVector<VT2,false>& );

template< typename VT1, typename VT2 >
decltype(auto) operator*( const SparseVector<VT1,false>&, const SparseVector<VT2,true>& );

template< typename MT, typename VT >
decltype(auto) operator*( const DenseMatrix<MT,false>&, const DenseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const DenseMatrix<MT,true>&, const DenseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const DenseMatrix<MT,false>&, const SparseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const DenseMatrix<MT,true>&, const SparseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const SparseMatrix<MT,false>&, const DenseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const SparseMatrix<MT,true>&, const DenseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const SparseMatrix<MT,false>&, const SparseVector<VT,false>& );

template< typename MT, typename VT >
decltype(auto) operator*( const SparseMatrix<MT,true>&, const SparseVector<VT,false>& );

template< typename VT, typename MT >
decltype(auto) operator*( const DenseVector<VT,true>&, const DenseMatrix<MT,false>& );

template< typename VT, typename MT >
decltype(auto) operator*( const DenseVector<VT,true>&, const DenseMatrix<MT,true>& );

template< typename VT, typename MT >
decltype(auto) operator*( const DenseVector<VT,true>&, const SparseMatrix<MT,false>& );

template< typename VT, typename MT >
decltype(auto) operator*( const DenseVector<VT,true>&, const SparseMatrix<MT,true>& );

template< typename VT, typename MT >
decltype(auto) operator*( const SparseVector<VT,true>&, const DenseMatrix<MT,false>& );

template< typename VT, typename MT >
decltype(auto) operator*( const SparseVector<VT,true>&, const DenseMatrix<MT,true>& );

template< typename VT, typename MT >
decltype(auto) operator*( const SparseVector<VT,true>&, const SparseMatrix<MT,false>& );

template< typename VT, typename MT >
decltype(auto) operator*( const SparseVector<VT,true>&, const SparseMatrix<MT,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,false>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,true>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const DenseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,false>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,true>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator*( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );


template< typename VT1, typename VT2, bool TF >
decltype(auto) operator/( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2 >
decltype(auto) operator/( const DenseVector<VT1,false>&, const DenseVector<VT2,true>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator/( const SparseVector<VT1,TF>&, const DenseVector<VT2,TF>& );


template< typename VT1, typename VT2, bool TF >
decltype(auto) operator%( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator%( const DenseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator%( const SparseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) operator%( const SparseVector<VT1,TF>&, const SparseVector<VT2,TF>& );


template< typename MT1, typename MT2, bool SO >
decltype(auto) operator%( const DenseMatrix<MT1,SO>&, const DenseMatrix<MT2,SO>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const DenseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const DenseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const DenseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const DenseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const DenseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const DenseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,false>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,false>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,true>&, const DenseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,true>&, const DenseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) operator%( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );


template< typename VT1, typename VT2, bool TF >
decltype(auto) kron( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) kron( const DenseVector<VT1,TF>&, const SparseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) kron( const SparseVector<VT1,TF>&, const DenseVector<VT2,TF>& );

template< typename VT1, typename VT2, bool TF >
decltype(auto) kron( const SparseVector<VT1,TF>&, const SparseVector<VT2,TF>& );


template< typename MT1, bool SO1, typename MT2, bool SO2 >
decltype(auto) kron( const DenseMatrix<MT1,SO1>&, const DenseMatrix<MT2,SO2>& );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
decltype(auto) kron( const DenseMatrix<MT1,SO1>&, const SparseMatrix<MT2,SO2>& );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
decltype(auto) kron( const SparseMatrix<MT1,SO1>&, const DenseMatrix<MT2,SO2>& );

template< typename MT1, typename MT2 >
decltype(auto) kron( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) kron( const SparseMatrix<MT1,false>&, const SparseMatrix<MT2,true>& );

template< typename MT1, typename MT2 >
decltype(auto) kron( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,false>& );

template< typename MT1, typename MT2 >
decltype(auto) kron( const SparseMatrix<MT1,true>&, const SparseMatrix<MT2,true>& );


template< typename VT, bool TF >
decltype(auto) trans( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) trans( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) trans( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) trans( const SparseMatrix<MT,SO>& );


template< bool TTF, typename VT, bool TF >
decltype(auto) transTo( const DenseVector<VT,TF>& );

template< bool TTF, typename VT, bool TF >
decltype(auto) transTo( const SparseVector<VT,TF>& );


template< bool B, typename MT, bool SO >
decltype(auto) transIf( const DenseMatrix<MT,SO>& );

template< bool B, typename MT, bool SO >
decltype(auto) transIf( const SparseMatrix<MT,SO>& );


template< typename VT, bool TF >
decltype(auto) eval( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) eval( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) eval( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) eval( const SparseMatrix<MT,SO>& );


template< typename VT, bool TF >
decltype(auto) serial( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) serial( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) serial( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) serial( const SparseMatrix<MT,SO>& );


template< typename VT, bool TF >
decltype(auto) noalias( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) noalias( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) noalias( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) noalias( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) inv( const DenseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) matexp( const DenseMatrix<MT,SO>& );


template< typename MT, bool SO, typename VT, bool TF >
decltype(auto) solve( const DenseMatrix<MT,SO>&, const DenseVector<VT,TF>& );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
decltype(auto) solve( const DenseMatrix<MT1,SO1>&, const DenseMatrix<MT2,SO2>& );


template< typename MT, bool SO >
decltype(auto) eigen( const DenseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) svd( const DenseMatrix<MT,SO>& );


template< typename VT, bool TF, typename OP >
decltype(auto) map( const DenseVector<VT,TF>&, OP );

template< typename VT, bool TF, typename OP >
decltype(auto) map( const SparseVector<VT,TF>&, OP );

template< typename MT, bool SO, typename OP >
decltype(auto) map( const DenseMatrix<MT,SO>&, OP );

template< typename MT, bool SO, typename OP >
decltype(auto) map( const SparseMatrix<MT,SO>&, OP );

template< typename VT1, typename VT2, bool TF, typename OP >
decltype(auto) map( const DenseVector<VT1,TF>&, const DenseVector<VT2,TF>&, OP );

template< typename VT1, typename VT2, bool TF, typename OP >
decltype(auto) map( const DenseVector<VT1,false>&, const DenseVector<VT2,true>&, OP );

template< typename MT1, typename MT2, bool SO, typename OP >
decltype(auto) map( const DenseMatrix<MT1,SO>&, const DenseMatrix<MT2,SO>&, OP );

template< typename MT1, typename MT2, typename OP >
decltype(auto) map( const DenseMatrix<MT1,false>&, const DenseMatrix<MT2,true>&, OP );

template< typename MT1, typename MT2, typename OP >
decltype(auto) map( const DenseMatrix<MT1,true>&, const DenseMatrix<MT2,false>&, OP );


template< typename VT, bool TF, typename OP >
decltype(auto) reduce( const DenseVector<VT,TF>&, OP );

template< typename VT, bool TF, typename OP >
decltype(auto) reduce( const SparseVector<VT,TF>&, OP );

template< typename MT, bool SO, typename OP >
decltype(auto) reduce( const DenseMatrix<MT,SO>&, OP );

template< ReductionFlag RF, typename MT, bool SO, typename OP >
decltype(auto) reduce( const DenseMatrix<MT,SO>&, OP );

template< typename MT, bool SO, typename OP >
decltype(auto) reduce( const SparseMatrix<MT,SO>&, OP );

template< ReductionFlag, typename MT, bool SO, typename OP >
decltype(auto) reduce( const SparseMatrix<MT,SO>&, OP );


template< typename VT, bool TF >
decltype(auto) expand( const DenseVector<VT,TF>&, size_t );

template< size_t E, typename VT, bool TF >
decltype(auto) expand( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) expand( const SparseVector<VT,TF>&, size_t );

template< size_t E, typename VT, bool TF >
decltype(auto) expand( const SparseVector<VT,TF>& );


template< typename VT, bool TF >
decltype(auto) repeat( const DenseVector<VT,TF>&, size_t );

template< size_t R0, typename VT, bool TF >
decltype(auto) repeat( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) repeat( const SparseVector<VT,TF>&, size_t );

template< size_t R0, typename VT, bool TF >
decltype(auto) repeat( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) repeat( const DenseMatrix<MT,SO>&, size_t, size_t );

template< size_t R0, size_t R1, typename MT, bool SO >
decltype(auto) repeat( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) repeat( const SparseMatrix<MT,SO>&, size_t, size_t );

template< size_t R0, size_t R1, typename MT, bool SO >
decltype(auto) repeat( const SparseMatrix<MT,SO>& );


template< typename VT, bool TF >
decltype(auto) mean( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) mean( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) mean( const DenseMatrix<MT,SO>& );

template< ReductionFlag, typename MT, bool SO >
decltype(auto) mean( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) mean( const SparseMatrix<MT,SO>& );

template< ReductionFlag, typename MT, bool SO >
decltype(auto) mean( const SparseMatrix<MT,SO>& );


template< typename VT, bool TF >
decltype(auto) var( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) var( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) var( const DenseMatrix<MT,SO>& );

template< ReductionFlag, typename MT, bool SO >
decltype(auto) var( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) var( const SparseMatrix<MT,SO>& );

template< ReductionFlag, typename MT, bool SO >
decltype(auto) var( const SparseMatrix<MT,SO>& );


template< typename VT, bool TF >
decltype(auto) stddev( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
decltype(auto) stddev( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
decltype(auto) stddev( const DenseMatrix<MT,SO>& );

template< ReductionFlag, typename MT, bool SO >
decltype(auto) stddev( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) stddev( const SparseMatrix<MT,SO>& );

template< ReductionFlag, typename MT, bool SO >
decltype(auto) stddev( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) declsym( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) declsym( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) declherm( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) declherm( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) decllow( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) decllow( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) declunilow( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) declunilow( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) declstrlow( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) declstrlow( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) declupp( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) declupp( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) decluniupp( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) decluniupp( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) declstrupp( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) declstrupp( const SparseMatrix<MT,SO>& );


template< typename MT, bool SO >
decltype(auto) decldiag( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) decldiag( const SparseMatrix<MT,SO>& );

} // namespace blaze

#endif
