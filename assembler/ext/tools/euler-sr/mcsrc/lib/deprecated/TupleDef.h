/***************************************************************************
 * Title:          TupleDef.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _TUPLE_DEF_H
#define _TUPLE_DEF_H

#ifdef _INTEL_ 
typedef long long ssize_t Tuple;
#else
#ifdef _3_3_
typedef long long ssize_t Tuple;
#else
#ifdef _powerpc_
typedef long long ssize_t Tuple;

#else
typedef long  Tuple;
#endif

#endif
 
#endif

#endif
