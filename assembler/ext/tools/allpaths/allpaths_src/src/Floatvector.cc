///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Floatvector.cc
 * \author tsharpe
 * \date Sep 9, 2009
 *
 * \brief Feudal vectors of floats.
 */
#include "Floatvector.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec< float, MempoolAllocator<float> >;
template class OuterVec<FloatVec>;

template class SmallVec< double, MempoolAllocator<double> >;
template class OuterVec<DoubleVec>;
