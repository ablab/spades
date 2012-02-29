///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

///\file Assert.cc
/// Assert.cc provides some of the guts of the Assert macros
#include "system/Assert.h"
#include <iostream>

void Assert::reportVals( char const* loc, char const* func, char const* vals )
{
    std::cout << loc << " failed in function\n" << func << "\n";

    if ( vals )
        std::cout << "with values " << vals;

    std::cout << std::endl;
}

void Assert::reportValsAndDie( char const* loc, char const* func, char const* vals )
{
    reportVals(loc,func,vals);
    CRD::exit(1);
}
