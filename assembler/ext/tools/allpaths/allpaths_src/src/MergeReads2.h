// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef MERGEREADS2
#define MERGEREADS2

#include "Basevector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "system/Types.h"
#include "Vec.h"

const unsigned char GapBase = 4;
const unsigned char NoBase = 5;

void PickWinner

     (  /* inputs */  const vec<unsigned char>& base, 
                      const qualvector& score, const vec<int>& lids,
                      const vec<Bool>& rc, Bool forward, int contig_count,
        /* outputs */ vec<unsigned char>& consensus, qualvector& quality,
                      vec<int>& start_on_contig, vec<int>& stop_on_contig,
                      vec<Bool>& RC, vec<int>& contig_no,
                      int scoring );

Bool CenterMobileGaps( align& a, const basevector& rd1, const basevector& rd2,
     bool verbose, ostream& log );

#endif
