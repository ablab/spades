/***************************************************************************
 * Title:          word.h 
 * Author:         Mark Chaisson, Glenn Tesler
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/*
** Copyright 2004 Ross Lippert
** 
** This file is part of bbbwt.
** 
** bbbwt is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** bbbwt is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with bbbwt; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
*/
#ifndef WORD_H
#define WORD_H

// 64 bit uint stuff goes here
// This stuff is somewhat architecture dependent

// Glenn Tesler-- changing all macros to assume C99 compliance
// g++ needs the following flags to bring in certain definitions

// needed for UINT64_C(x)
#ifndef __STDC_CONSTANT_MACROS
#  define __STDC_CONSTANT_MACROS
#endif

// needed for PRIu64
#ifndef __STDC_FORMAT_MACROS
#  define __STDC_FORMAT_MACROS
#endif

#include <inttypes.h>

typedef uint64_t word_t;
#define word_LIT(x) UINT64_C(x)
#define word_FMT PRIu64 

// BEGIN DISABLING ALL ORIGINAL MACROS
#if 0

// One day, when everything is C99 compliant, this won't be necessary

#if(__i386__)
#define _WORD_IS_LONGLONG
#endif
#if(__ppc__)
#define _WORD_IS_LONGLONG
#endif
#if(__ia64__)
#define _WORD_IS_LONG
#endif

// 64 bit words are either long int or long long int

#ifdef _WORD_IS_LONGLONG

typedef long long unsigned word_t;
#define word_LIT(x) (x ## ULL)
#define word_FMT "%llu"

#endif

#ifdef _WORD_IS_LONG

typedef long unsigned word_t;
#define word_LIT(x) (x ## UL)
#define word_FMT "%lu"

#endif

#endif
// END DISABLING ALL ORIGINAL MACROS


#endif
