///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PowerOf2.cc
 * \author tsharpe
 * \date Feb 4, 2011
 *
 * \brief
 */
#include "math/PowerOf2.h"

unsigned long PowerOf2::getNearbyPrime( unsigned powerOf2 )
{
    // table is from http://primes.utm.edu/lists/2small/
    static size_t gTab[] = { 1, 2, 3, 7, 13, 31, 61, 127, 251,
                             (1ul<<9)-3ul, (1ul<<10)-3ul, (1ul<<11)-9ul,
                             (1ul<<12)-3ul, (1ul<<13)-1ul, (1ul<<14)-3ul,
                             (1ul<<15)-19ul, (1ul<<16)-15ul, (1ul<<17)-1ul,
                             (1ul<<18)-5ul, (1ul<<19)-1ul, (1ul<<20)-3ul,
                             (1ul<<21)-9ul, (1ul<<22)-3ul, (1ul<<23)-15ul,
                             (1ul<<24)-3ul, (1ul<<25)-39ul, (1ul<<26)-5ul,
                             (1ul<<27)-39ul, (1ul<<28)-57ul, (1ul<<29)-3ul,
                             (1ul<<30)-35ul, (1ul<<31)-1ul, (1ul<<32)-5ul,
                             (1ul<<33)-9ul, (1ul<<34)-41ul, (1ul<<35)-31ul,
                             (1ul<<36)-5ul, (1ul<<37)-25ul, (1ul<<38)-45ul,
                             (1ul<<39)-7ul, (1ul<<40)-87ul, (1ul<<41)-21ul,
                             (1ul<<42)-11ul, (1ul<<43)-57ul, (1ul<<44)-17ul,
                             (1ul<<45)-55ul, (1ul<<46)-21ul, (1ul<<47)-115ul,
                             (1ul<<48)-59ul, (1ul<<49)-81ul, (1ul<<50)-27ul,
                             (1ul<<51)-129ul, (1ul<<52)-47ul, (1ul<<53)-111ul,
                             (1ul<<54)-33ul, (1ul<<55)-55ul, (1ul<<56)-5ul,
                             (1ul<<57)-13ul, (1ul<<58)-27ul, (1ul<<59)-55ul,
                             (1ul<<60)-93ul, (1ul<<61)-1ul, (1ul<<62)-57ul,
                             (1ul<<63)-25ul, 0ul-59ul };
  AssertLe(powerOf2,64u);
  return gTab[powerOf2];
}
