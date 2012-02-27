/***************************************************************************
 * Title:          modify_tuple.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _TUPLE_H_
#define _TUPLE_H_

#include "TupleDef.h"

Tuple snp(Tuple index, ssize_t n, ssize_t v, ssize_t k);
Tuple insert(Tuple index, ssize_t n, ssize_t v, ssize_t k);
ssize_t getNuc(Tuple index, ssize_t n, ssize_t k);
Tuple subset(Tuple index, ssize_t lb, ssize_t ub, ssize_t k);
Tuple del(Tuple orig, ssize_t pos);
Tuple ReverseComplement(Tuple orig, ssize_t tuple_len);
#endif
