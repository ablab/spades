///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadHBVs.h
 * \author tsharpe
 * \date Sep 28, 2011
 *
 * \brief
 */
#ifndef PATHS_READHBVS_H_
#define PATHS_READHBVS_H_

#include "String.h"
#include "Vec.h"
#include "paths/HyperBasevector.h"

/// Reads the localized.hbvs file in the specified directory, and loads
/// the selected HBVs into pHBVs. Warning messages can be turned off
/// setting verbose to false.
void readHBVs( String const& sub_dir,
	       vec<bool> const& selected,
	       vec<HyperBasevector>* pHBVs,
	       bool verbose = true );

#endif /* PATHS_READHBVS_H_ */
