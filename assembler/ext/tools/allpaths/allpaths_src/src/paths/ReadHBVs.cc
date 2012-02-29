///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadHBVs.cc
 * \author tsharpe
 * \date Sep 28, 2011
 *
 * \brief
 */
#include "paths/ReadHBVs.h"
#include "feudal/BinaryStream.h"
#include <cstddef>

void readHBVs( String const& sub_dir,
	       vec<bool> const& selected,
	       vec<HyperBasevector>* pHBVs,
	       bool verbose )
{
    String hbvs_file = sub_dir + "/localized.hbvs";
    BinaryReader rdr(hbvs_file.c_str());

    size_t nSeeds;
    rdr.read(&nSeeds);
    ForceAssertEq(nSeeds,selected.size());

    pHBVs->clear();
    HyperBasevector exemplar(0);
    pHBVs->resize(nSeeds,exemplar);

    size_t seedId;
    while ( rdr.read(&seedId) != ~0UL )
    {
        ForceAssertLt(seedId,nSeeds);
        if ( selected[seedId] )
            rdr.read(&pHBVs[0][seedId]);
        else
            rdr.read(&exemplar); // move the file pointer, but discard the HBV
    }
    ForceAssert(rdr.atEOF());

    if ( ! verbose ) return;

    for ( size_t idx = 0; idx < nSeeds; ++idx )
        if ( selected[idx] && !pHBVs[0][idx].K() )
            std::cout << "WARNING: Could not find HyperBasevector for seed #"
                      << idx << std::endl;
}
