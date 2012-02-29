///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FakeQualb.cc
 * \author tsharpe
 * \date Feb 7, 2012
 *
 * \brief Generate qualb file of the same shape as a fastb.
 * The scores are set to a constant value.
 */

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include <algorithm>
#include <cstddef>

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String(FASTB);
    CommandArgument_Int_OrDefault(SCORE,0);
    EndCommandArguments;

    File fastbFile(FASTB);
    VirtualMasterVec<bvec> fastb(fastbFile);
    IncrementalWriter<qvec> qualb(fastbFile.changeExtension(".qualb"));

    qvec qv;
    typedef VirtualMasterVec<bvec>::const_iterator Itr;
    for ( Itr itr(fastb.begin()), end(fastb.end()); itr != end; ++itr )
    {
        qv.resize(itr->size(),SCORE);
        qualb.add(qv);
    }
    qualb.close();
}
