/*===-- include/Support/DataTypes.h - Define fixed size types -----*- C -*-===*\
|*                                                                            *|
|*                     The LLVM Compiler Infrastructure                       *|
|*                                                                            *|
|* This file is distributed under the University of Illinois Open Source      *|
|* License. See LICENSE.TXT for details.                                      *|
|*                                                                            *|
|*===----------------------------------------------------------------------===*|
|*                                                                            *|
|* This file contains definitions to figure out the size of _HOST_ data types.*|
|* This file is important because different host OS's define different macros,*|
|* which makes portability tough.  This file exports the following            *|
|* definitions:                                                               *|
|*                                                                            *|
|*   [u]int(32|64)_t : typedefs for signed and unsigned 32/64 bit system types*|
|*   [U]INT(8|16|32|64)_(MIN|MAX) : Constants for the min and max values.     *|
|*                                                                            *|
|* No library is required when using these functions.                         *|
|*                                                                            *|
|*===----------------------------------------------------------------------===*/

/* Please leave this file C-compatible. */

/* Please keep this file in sync with DataTypes.h.in */

#ifndef SUPPORT_DATATYPES_H
#define SUPPORT_DATATYPES_H

#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UINT64_T 1
#define HAVE_U_INT64_T 1

#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif

#include <inttypes.h>

/* Note that <inttypes.h> includes <stdint.h>, if this is a C99 system. */
#include <sys/types.h>

/* Set defaults for constants which we cannot find. */
#if !defined(INT64_MAX)
# define INT64_MAX 9223372036854775807LL
#endif
#if !defined(INT64_MIN)
# define INT64_MIN ((-INT64_MAX)-1)
#endif
#if !defined(UINT64_MAX)
# define UINT64_MAX 0xffffffffffffffffULL
#endif

#endif  /* SUPPORT_DATATYPES_H */
