///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Floatvector.h
 * \author tsharpe
 * \date Sep 9, 2009
 *
 * \brief Feudal vectors of floats.
 */
#ifndef FLOATVECTOR_H_
#define FLOATVECTOR_H_

#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

typedef SerfVec<float> FloatVec;
typedef MasterVec< FloatVec > VecFloatVec;

typedef SerfVec<double> DoubleVec;
typedef MasterVec< FloatVec > VecDoubleVec;

#endif /* FLOATVECTOR_H_ */
