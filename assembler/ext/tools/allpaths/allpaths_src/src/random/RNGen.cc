///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file RNGen.cc
 * \author tsharpe
 * \date Nov 17, 2011
 *
 * \brief
 */
#include "random/RNGen.h"

RNGen RNGen::gDflt(1u);
RNGen RNGen::gSystem;
SpinLockedData RNGen::gLock;
