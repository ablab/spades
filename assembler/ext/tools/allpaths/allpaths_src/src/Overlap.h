// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef OVERLAP
#define OVERLAP

#include "Alignment.h"
#include "Basevector.h"
#include "PackAlign.h"

int EstimatedOverlap( const alignment& a,
     const basevector& rd1, const basevector& rd2 );

int EstimatedOverlap( const align& a,
     const basevector& rd1, const basevector& rd2 );

#endif
