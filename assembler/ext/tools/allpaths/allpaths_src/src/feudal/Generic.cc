///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Generic.cc
 * \author tsharpe
 * \date Mar 10, 2010
 *
 * \brief
 */
#include "feudal/Generic.h"
#include "String.h"
#include <cstdlib>
#include <unistd.h>

UtilizationReporter UtilizationReporter::gInstance;

UtilizationReporter::UtilizationReporter() : mpFetcher(0)
{
    if ( getenv("FEUDAL_USAGE") )
    {
        String fileName = "feudal_utilization_" + ToString(getpid()) + ".txt";
        mLog.open(fileName.c_str());
        mpFetcher = new ThreadsafeOStreamFetcher(mLog);
    }
}

UtilizationReporter::~UtilizationReporter()
{
    delete mpFetcher;
    if ( mLog.is_open() )
        mLog.close();
}

void UtilizationReporter::writeLine( void const* addr, size_t bits, size_t siz,
                                        size_t cap, char const* type )
{
    mpFetcher->get() << '(' << addr << ")\t[" << type << "] [Bits: " << bits
            << "] [Sz: " << siz << "] [Cap: " << cap << ']' << std::endl;
}
