///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalString.h
 * \author ghall
 * \author tsharpe
 * \date Oct 5, 2009
 *
 * \brief Exactly like std::string, but only 2^32 elements, uses SmallVec for containment
 */

#ifndef FEUDAL_TOOLS_H
#define FEUDAL_TOOLS_H

#include "String.h"
#include "Vec.h"
#include <cstddef>
#include <vector>

void MergeMastervecs(const String& dstFile, std::vector<String> const& inputs);
void MergeMastervecs(const String& file1, const String& file2, const String& dest_file);

size_t MastervecFileObjectCount(const String& filename);
size_t MastervecFileRawCount( const String& filename, size_t dataSize=0 );
bool IsGoodFeudalFile( const String& filename, bool verbose=false );
void RemoveMastervecFiles(const String& filename);

template<class SwappableVec>
void PermuteSwappableVec(SwappableVec & v,
const vec<int> & permutation)
{
    AssertEq((longlong) v.size(), (longlong) permutation.size());
    vec<int> o = permutation;
    for (int i = 0; i != (longlong) v.size(); ++i)
    {
        while (o[i] != i && o[i] != -1)
        {
            v.SwapElements(i, o[i]);
            std::swap(o[i], o[o[i]]);
        }
    }
}

template<class Mastervec1, class Mastervec2>
bool SameSizes(const Mastervec1 & v1, const Mastervec2 & v2)
{
    bool result = v1.size() == v2.size();
    if ( result )
    {
        typedef typename Mastervec1::const_iterator Itr1;
        typedef typename Mastervec2::const_iterator Itr2;
        Itr1 end(v1.end());
        Itr2 itr2(v2.begin());
        for ( Itr1 itr(v1.begin()); result && itr != end; ++itr, ++itr2 )
            result = itr->size() == itr2->size();
    }
    return result;
}


#endif // FEUDAL_TOOLS_H
