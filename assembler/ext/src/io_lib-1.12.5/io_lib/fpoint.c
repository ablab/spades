/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <math.h>
#include "io_lib/fpoint.h"
/*
extern double log ( double x ) ;
extern double exp ( double x ) ;
*/
#define IEEE

float int_to_float(int in)
/*
** interpret the integer in as a
** floating point number in IEEE format
*/
{
   /*
  Assume `in' is stored as a float according to the 
  ANSI IEEE 754-1985 standard. See the tables below:

  s = sign ( 1 bit)
  e = biased exponent (8 bits)
  f = fraction (23 bits)

  floating point number =  (-1)^s 2^(e-127) 1.f

     Bits  Name      Content
      31   Sign      1 iff number is negative
    23-30  Exponent  Eight-Bit exponent, biased by 127
     0-22  Fraction  23-bit fraction component of normalised significant.
		     The "one" bit is "hidden"

  If IEEE floating point format is supported on your machine...
  ensure there is a #define IEEE somewhere. 
  */

#ifdef IEEE
  union {
    int i;
    float f;
  } cvt;
  cvt.i = in;
  return cvt.f;
#else
  int fraction;
  int exponent;
  int sign;

  fraction = in & ( (1<<23)-1 );
  exponent = (in >> 23) & ( (1<<8)-1 );
  sign = (in >> 31);

  return
    (float) (
      (sign?-1.0:1.0) *
      exp ( log ( (double) 2.0) * (double) (exponent - 127 - 23) ) *
      (double) ((1<<23)+fraction)) ;
#endif
}
