// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#ifndef QUALITY
#define QUALITY

#include "Qualvector.h"
#include "String.h"
#include "system/Types.h"
#include "Vec.h"

//
// Writes the quality scores in Q, into file score_file, and gzips if gzip==True
//

inline void WriteQualityScores( String score_file, const vecqualvector& Q,
                                Bool ignored )
{ Q.WriteAll(score_file); }

//
// Reads the quality scores Q from score_file, where read lengths are given in
// 'lengths' vec.  Qrc is also initialized (reverse complement qualities) if
// noQrc==False.
//

void ReadQualityScores( String score_file, const vec<int>& lengths,
     vecqualvector& Q, vecqualvector& Qrc,
     Bool noQrc = False, Bool append = False );

inline void ReadQualityScores( String score_file, const vec<int>& lengths,
     vecqualvector& Q, Bool append = False )
{    vecqualvector Qrc;
     ReadQualityScores( score_file, lengths, Q, Qrc, True, append );   }

void ReadSubsetOfQualityScores( String score_file, const vec<int>& lengths,
     vecqualvector& Q, vec<int> use, int extra_space = 0, Bool append = False );

#endif
