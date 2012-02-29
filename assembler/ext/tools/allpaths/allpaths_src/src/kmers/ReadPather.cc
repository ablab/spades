///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadPather.cc
 * \author tsharpe
 * \date Dec 13, 2011
 *
 * \brief
 */
#include "kmers/ReadPather.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

CanonicalForm EdgeDesc::gCompStatus[3] = { REV, FWD, PALINDROME };
unsigned const EdgeID::NULLVAL;
unsigned const EdgeOffset::NULLVAL;
template class SmallVec< UnipathEvidence, MempoolAllocator<UnipathEvidence> >;
template class OuterVec<UnipathEvidenceVec>;
