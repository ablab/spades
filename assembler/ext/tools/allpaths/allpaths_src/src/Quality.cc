// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "system/Assert.h"
#include "math/Functions.h"
#include "Quality.h"
#include "Qualvector.h"
#include "String.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

void ReadQualityScores( String score_file, const vec<int>& lengths,
     vecqualvector& Q, vecqualvector& Qrc, Bool noQrc, Bool append )
{
    ForceAssert( IsRegularFile(score_file) ^ IsRegularFile( score_file + ".gz" ) );
    if ( IsRegularFile(score_file + ".gz") )
        System("gzip -d " + score_file + ".gz");

    if ( !append )
        Q.clear();
    Q.reserve(Q.size() + lengths.size());

    vecqvec::size_type orig_size = Q.size();

    Q.ReadRange(score_file, 0, lengths.size());

    if ( !noQrc )
    {
        if ( !append )
            Qrc.clear();
        Qrc.reserve(Qrc.size() + lengths.size());

        for ( vecqvec::size_type i = orig_size; i < Q.size(); i++ )
        {
            Qrc.push_back(Q[i]);
            Qrc[i].ReverseMe();
        }
    }
}

void ReadSubsetOfQualityScores( String score_file, const vec<int>& lengths,
     vecqualvector& Q, vec<int> use, int extra_space, Bool append )
{    ForceAssert( IsRegularFile(score_file) );
     if ( !append ) Q.clear( );
     Q.reserve( use.size() + Q.size() );
     Q.SparseRead( score_file, use, extra_space, True );    }
