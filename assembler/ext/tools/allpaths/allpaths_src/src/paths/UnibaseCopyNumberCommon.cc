///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file UnibaseCopyNumberCommon.cc
 * \author tsharpe
 * \date Dec 16, 2011
 *
 * \brief
 */
#include "paths/UnibaseCopyNumberCommon.h"

void GetMinLength( const vecbvec& unibases, int& min_length, int& nlongest )
{
    vec<unsigned> ulens(unibases.size(), 0);
    if ( unibases.size() == 0 ) // hope this can't ever happen
    {
        min_length = 0;
        nlongest = 0;
        return;
    }

    for ( size_t i = 0; i < unibases.size(); i++ )
        ulens[i] = unibases[i].size();

    ReverseSort(ulens);

    // compute kmer starts at top 100 longest unipaths, but don't go any lower
    // than 10% of the longest unipath
    int top = 100;
    double min_frac = 0.1;
    nlongest = Min(top, ulens.isize());
    while ( nlongest >= 2 && ulens[nlongest - 1] < min_frac * double(ulens[0]) )
    {
        nlongest--;
    }
    min_length = ulens[nlongest - 1];
}

void ComputeBias(
     // INPUTS:
     const int K,
     double occCnPloidy,
     const vec<longlong>& biasOccs,
     const vec<longlong>& biasInst,
     // OUTPUT:
     vec<double>& biasCurveLoc,
     // LOGGING:
     const Bool VERBOSE,
     ostream& logout )
{
    logout << Date() << ": Performing bias computation..." << endl;
    if ( VERBOSE )
    {
        PRINT_TO( logout, occCnPloidy );
        biasOccs.Print(logout);
        logout << endl;
        biasInst.Print(logout);
        logout << endl;
    }
    ForceAssert( occCnPloidy > 0.0 );

    int CTHRESHOLD = 5; // at least 5 instances to perform bias computation
    for ( int i = 0; i <= K; i++ )
        if ( biasInst[i] >= CTHRESHOLD )
            biasCurveLoc[i] = (double) biasOccs[i] / (double) biasInst[i];

    if ( VERBOSE )
    {
        biasCurveLoc.Print(logout);
        logout << endl;
    }

    for ( int i = 0; i <= K; i++ )
        if ( biasInst[i] >= CTHRESHOLD )
            biasCurveLoc[i] = biasCurveLoc[i] / occCnPloidy;

    double sumOld = 0, sumNew = 0;
    for ( int i = 0; i <= K; i++ )
    {
        sumOld += biasOccs[i];
        sumNew += biasOccs[i] / biasCurveLoc[i];
    }

    for ( int i = 0; i <= K; i++ )
        if ( biasInst[i] >= CTHRESHOLD )
            biasCurveLoc[i] *= sumNew / sumOld;

    {
        double sumOld = 0, sumNew = 0;
        for ( int i = 0; i <= K; i++ )
        {
            sumOld += biasOccs[i];
            sumNew += biasOccs[i] / biasCurveLoc[i];
        }
        if ( VERBOSE )
        {
            PRINT2_TO( logout, sumOld, sumNew );
            biasCurveLoc.Print(logout);
        }
    }
}
