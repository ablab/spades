/***************************************************************************
 * Title:          IntegralTupleStatic.h
 * Author:         Glenn Tesler
 * Created:        2010
 * Last modified:  02/25/2010
 *
 * Copyright (c) 2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

// For using in libraries, or when included from other .h files:
//     #include "IntegralTuple.h"
// For using in executables:
//     #include "IntegralTupleStatic.h"

#ifndef INTEGRAL_TUPLE_STATIC_H_
#define INTEGRAL_TUPLE_STATIC_H_

#include "IntegralTuple.h"

// Allocate all the static variables of the IntegralTuple class

//LongTupleWord IntegralTuple::MASK_OFF_FIRST = 0;
LongTupleWord IntegralTuple::VERTEX_MASK    = 0;
int IntegralTuple::tupleSize = 0;
int IntegralTuple::numWords = 0;
int IntegralTuple::numWords1 = 0;

#endif
