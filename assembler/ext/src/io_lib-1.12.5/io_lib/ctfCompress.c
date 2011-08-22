/*  File: ctftrace.c
    %W% %G%
*/

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "io_lib/seqIOCTF.h"
#include "io_lib/stdio_hack.h"

#define ACEDB4

#define A_ 1
#define T_ 2
#define G_ 4
#define C_ 8

/*#include "regular.h"
#include "keyset.h"
#include "dna.h"
*/
#include "io_lib/array.h"
#include "io_lib/Read.h"

/*
 * #defines to modify this file to compile easily as part of the Staden
 * Package io_lib.
 */
#define freeSeq(_s) (read_deallocate(_s), (_s) = 0)
#define seqMax(_seq)  ((_seq)->NPoints )
#define seqMaxBase(_seq)  ((_seq)->NBases)

#define xmalloc malloc
#define xcalloc calloc
#define xrealloc realloc
#define xfree free

#define arrayCreate(s,t) ArrayCreate(sizeof(t),(s))
#define array(a,n,t) ARR(t,a,n)
#define arrayMax ArrayMax
#define arrayDestroy ArrayDestroy
#undef arrp
#define arrp(a,n,t) \
    &((t*)((a)->base))[n]
#define arrayReCreate(a,n,t) arrayCreate(n,t)
#define arrayExists(a) ((a)->base != NULL)

#define BOOL int
#define mysize_t size_t
#define FALSE 0
#define TRUE 1


#define MAGIC 523747
#define PREDICTIONMODE 3  /* predictor degree */
#define COMPRESSIONMODE 3 /* compressor version */

static char *ctfType = 0 ;

/**********************************************************/
/**********************************************************/

static int ctfTracePeakValue ( Read *seq)
{
  int i, ii, max = 0 ;
  TRACE *bp[4], *u ;
  
  bp[0] = seq->traceA ;
  bp[1] = seq->traceC ;
  bp[2] = seq->traceG ;
  bp[3] = seq->traceT ;
  
  ii = 4 ; 
  while (ii--)
    {
      u = bp[ii] - 1 ;
      i = seqMax(seq) ;
      while (u++, i--)
	if (*u > max) max = *u ;
    }
  
  seq->maxTraceVal = max ;
  return max ;
}

/**********************************************************/
/* get/store in sex independant way */
static int ctfGetInt (unsigned char *cp)
{
  int n = 0, i = 4 ;
  
  while (i--) 
    { n <<= 8 ; n |= *cp++ ; }
  return n ;
}

/**********/

static void ctfStoreInt (unsigned char *cp, int n)
{
  int i = 4 ;
  
  cp += 4 ;
  while (cp--, i--) 
    { *cp = n & 0xff ; n >>= 8 ; }
}

/**********************************************************/
/**********************************************************/
/* Method zero
   store the shorts as a pair of char
*/
static Array ctfCompress0 (Array a)
{ 
  int i = arrayMax (a) ;
  Array b = arrayCreate (2 * i, unsigned char) ; /* always true */
  unsigned char *cp, *cp0 ;
  short *sp, z ;

  array (b, 2*i - 1 , unsigned char) = 0 ;  /* make room */
  cp = cp0 = arrp (b , 0, unsigned char) ;
  sp = arrp (a, 0, short) ;

  while (i--)
    {
      z = *sp++ ;
      *cp++ = (z >> 8) & 0xff ;
      *cp++ = z &  0xff ;
    }
  arrayMax(b) = cp - cp0 ;
  return b ;
}

/**********************************************************/

static Array ctfDecompress0 (int dataMax, int shMax,
			     unsigned char *cp)
{ 
  int i = shMax ;
  Array b = 0 ;
  short *sp ;

  if (dataMax != 2 * shMax)
    return 0 ;

  b = arrayCreate (shMax, short) ;
  array (b, shMax - 1, short) = 0 ;  /* make room */
  sp = arrp (b, 0, short) ;

  while (i--)
    {
      *sp++ = ((*cp)<< 8) | (*(cp + 1)) ;
      cp += 2 ;
    }
  return b ;
}

/**********************************************************/
/**********************************************************/
/* Method one
   converts short s to unsigned char cc:
   if s in range -126, +126, cc = s + 128
   else transmit 0xFF then value on 2 bytes 
*/
static Array ctfCompress1 (Array a)
{
  int i = arrayMax (a), j ;
  Array b = arrayCreate (3 *i, char) ; /* unreliable size, use arrayp */
  short *sp, z ;

  array (b, 3*i - 1 , unsigned char) = 0 ;  /* make room */
  j = 0 ;
  sp = arrp (a, 0, short) ;

  while (i--)
    {
      z = *sp++ + 128 ;
      while (z >= 254) { z -= 252 ; array (b , j++, unsigned char) = 254 ;}
      while (z <= 1) { z += 252 ; array (b , j++, unsigned char) = 1 ; }
      array (b , j++, unsigned char) = z ;
    }
  arrayMax (b) = j ;
  return b ;
}

/**********************************************************/

static Array ctfDecompress1 (int dataMax, int shMax,
			    unsigned char *cp)
{
  int i = dataMax, dz = 0 ;
  short *sp, *spMax ;
  Array b = arrayCreate (shMax, short) ;

  array (b, shMax  - 1, short) = 0 ;  /* make room */
  sp = arrp (b, 0, short) ;
  spMax = sp + shMax ;

  cp-- ;
  while (cp++, i-- && sp < spMax)
    switch (*cp)
      {
      case 1: dz -= 252 ; break ;
      case 254: dz += 252 ; break ;
      default: *sp++ = dz + *cp - 128 ; dz = 0 ; break ;
      }
  if (i != -1 || sp != spMax)
    arrayDestroy (b) ;
  return b ;
}

/**********************************************************/
/**********************************************************/
/* Method two
   convert strings of zeroes as single chars
   small values as 7 bit codes
   rest as shorts in next 2 bytes
*/

static Array ctfCompress2 (Array a)
{
  int n0, n1, n2, n3, n4 ;
  int i = arrayMax (a), j = 0 ;
  Array b = arrayCreate (3 *i, char) ; /* worst case */
  unsigned char *cp, *cp0 ;
  short *sp, z ;

  array (b, 3*i - 1 , unsigned char) = 0 ;  /* make room */
  cp = cp0 = arrp (b , 0, unsigned char) ;
  sp = arrp (a, 0, short) ;

  n0 = n1 = n2 = n3 = n4 = 0 ;
  while (i--)
    {
      z = *sp++ ;
      if (!z)  /* string of zeroes */
	{
	  j =1 ;
	  while (i > 0 && j < 126 && !(*sp)) { j++ ; sp++; i-- ; } ;
	  *cp++ =  (j & 0x7f) ; /* bit 1 = 0 */
	  n0 += j ; n1++ ;
	}
      else if ( z < 63 && z > - 63)
	{
	  j = z + 63 ; /* range 1 ... 125 */
	  *cp++ = 0x80 | (j & 0x7f) ;
	  n2++ ;
	}
      else if ( z < 128 && z > -129)
	{
	  j = z + 128 ; /* range 0 ... 255 */
	  *cp++ = 0x80 | 126 ;
	  *cp++ = j & 0xff ;
	  n3++ ;
	}
      else
	{
	  j = z ; 
	  *cp++ = 0x80 | 127 ;
	  *cp++ = (j >> 8) & 0xff ;
	  *cp++ = j & 0xff ;	
	  n4++ ;
	}
    }
  arrayMax(b) = cp - cp0 ;
  printf ("compress2 : %d zeros in %d bytes, %d < 7 , %d byte, %d short total %d char for %d shrt\n",
	  n0, n1, n2, n3, n4, arrayMax(b), arrayMax(a)) ;
  return b ;
}
/*
compress2 : 
  12338 zeros in 3328 bytes, 14616 < 7 , 23 byte, 3 short 
          total 17999 char for 26980 shrt
*/

/**********************************************************/

static Array ctfDecompress2 (int dataMax, int shMax,
			    unsigned char *cp)
{
  int i = dataMax, mode, arg ;
  unsigned char cc, cc1, cc2 ;
  short *sp, *spMax ;
  Array b = arrayCreate (shMax, short) ;

  array (b, shMax  - 1, short) = 0 ;  /* make room */
  sp = arrp (b, 0, short) ;
  spMax = sp + shMax ;

  while (i-- && sp < spMax)
    {
      cc = *cp++ ;
      mode = cc & 0x80 ; arg = cc & 0x7f ;
      switch (mode)
	{
	case 0: /* initial zero = string of zero */
	  while (arg-- && sp < spMax) *sp++ = 0 ; 
	  break ;
	case 0x80:
	  switch (arg)
	    {
	    case 127:   /* next 2 bytes is a short */
	      i -= 2 ;
	      cc1 = *cp++ ; cc2 = *cp++ ;
	      *sp++ = (cc1 << 8) | cc2 ;
	      break ;
	    case 126:   /* next byte is a byte */
	      i-- ;
	      cc1 = *cp++ ;
	      *sp++ = cc1 - 128 ;
	      break ;
	    default:   /* 7 bytes is sufficient */
	      *sp++ = arg - 63 ;
	      break ;
	    }
	}
    }
  if (i != -1 || sp != spMax)
    arrayDestroy (b) ;
  return b ;
}

/**********************************************************/
/**********************************************************/
/* Method two
   convert strings of zeroes as single chars
   small values as 7 bit codes
   rest as shorts in next 2 bytes
*/

/**********************************************************/
/* create a code for the 125 most frequent words */
static void ctfCompress3Init (Array *aap, int **lp, int **mp, int *maxCodep)
{
  short *sp ; int  i, j, k ;
  static int lng[128], mark[128], maxCode = 0 ;
  static Array aa = 0 ; 


  *aap = aa ; *lp = lng ; *mp = mark ; *maxCodep = maxCode ;
  if (aa) return ;
  *aap = aa = arrayCreate (512, short) ;
  array (aa, 511, short) = 0 ; /* make room */
  sp = arrp (aa, 0, short) ;
  j = 0 ;

  i = 0 ;  /* empty word */
  mark[i] = j ; j += lng[i] ; 
  
  /* single values up to +- 8 */
  for (k = 1 ; i < 126 && k < 12 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ;
      lng [i] = 1 ; j++ ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ;
      lng [i] = 1 ; j++ ;
    }
  /* double values up to 50 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ; *sp++ = 0 ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ; *sp++ = 0 ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to 51 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ; *sp++ = 1 ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ; *sp++ = 1 ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to 5-1 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ; *sp++ = -1 ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ; *sp++ = -1 ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to 15 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = 1 ; *sp++ = k ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = 1 ; *sp++ = -k ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to -15 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = -1 ; *sp++ = k ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -1 ; *sp++ = -k ;
      lng [i] = 2 ; j += 2 ;
    }
  /* triple values up to 111 */

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;

  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;

  /* quadruple values up to 1111 */

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++ =- 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = 1 ;
  lng [i] = 5 ; j += 5 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 5 ; j += 5 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = 1 ;
  lng [i] = 5 ; j += 5 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 5 ; j += 5 ;

  if (i >= 126) 
    { fprintf (stderr, "FATAL ERROR, ctfCompress3Init i=%d > 126", i) ;
      exit (1) ;
    }
  if (j > 511)  
    { 
      fprintf (stderr, "FATAL ERROR, ctfCompress3Init j=%d > 511", j) ; 
      exit (1) ;
    }
  *maxCodep = maxCode = i ;
}

/**********************************************************/

static Array ctfCompress3 (Array a)
{
  int n0, n1, n10, n2, n3, n4 ;
  int i = arrayMax (a), j = 0, n ;
  Array b = arrayCreate (3 *i, char) ; /* worst case */
  unsigned char *cp, *cp0 ;
  short *sp, *sp1, *wp, z ;
  int w, *lng, *mark, maxCode ;
  Array aa = 0 ; 
  BOOL debug = FALSE ;

  ctfCompress3Init (&aa, &lng, &mark, &maxCode) ;

  array (b, 3*i - 1 , unsigned char) = 0 ;  /* make room */
  cp = cp0 = arrp (b , 0, unsigned char) ;
  sp = arrp (a, 0, short) ;

  n0 = n1 = n10 = n2 = n3 = n4 = 0 ;
  while (i--)
    {
      z = *sp++ ;
      if (!z)  /* string of zeroes */
	{
	  j =1 ;
	  while (i > 0 && j < 126 && !(*sp)) { j++ ; sp++; i-- ; } ;
	  *cp++ =  j & 0x7f ;  /* bit 1 = 0 */
	  n0 += j ; n1++ ;
	  continue ;
	}
       /* search for code word */
      for (w = maxCode ; w > 1 ; w--)  /* w-- favors long code words */
	{
	  n = lng[w] ; wp = arrp (aa, mark[w], short) ; sp1 = sp - 1 ;
	  while (n-- && *wp++ == *sp1++) ;
	  if (n == -1) break ;
	}
      if (w > 1 && lng[w] < i) /* code word found */
	{
	  n2++ ; n10 += lng[w] ;  *cp++ = 0x80 | (w & 0x7f) ;
	  sp += lng[w] - 1 ; i -= lng[w] - 1 ;
	  if (lng[w] <= 0)
	    {
	      fprintf (stderr, "FATAL ERROR in ctfCompress3 bad coding lng[w]") ;
	      exit (1) ;
	    }
	}
      else if ( z < 128 && z > -129) /* transmit byte */
	{
	  j = z + 128 ; /* range 0 ... 255 */
	  *cp++ = 0x80 | 126 ;
	  *cp++ = j & 0xff ;
	  n3++ ;
	}
      else     /* transmit short */
	{
	  j = z ; 
	  *cp++ = 0x80 | 127 ;
	  *cp++ = (j >> 8) & 0xff ;
	  *cp++ = j & 0xff ;	
	  n4++ ;
	}
    }
  arrayMax(b) = cp - cp0 ;
  if (debug) 
    printf (" // compress3:\n//  %d zeros in %d bytes, %d values coded in %d byte, %d bytes, %d short. \n// Total %d char for %d shrt\n",
	  n0, n1, n10, n2, n3, n4, arrayMax(b), arrayMax(a)) ;
  if (arrayMax(a) != n0 + n10 + n3 + n4)
    { 
      fprintf (stderr, "FATAL ERROR in ctfCompress3, codind error in compress 3") ;
      exit (1) ;
    }
  return b ;
}
/*
compress3 : 
//found   10829 zeros in 1865 bytes, 16524 values coded in 9114 byte, 183 bytes, 0 short. 
// Total 11162 char for 27536 shrt
*/

/**********************************************************/

static Array ctfDecompress3 (int dataMax, int shMax,
			    unsigned char *cp)
{
  int i = dataMax, mode, arg, n ;
  unsigned char cc, cc1, cc2 ;
  short *sp, *spMax, *wp ;
  Array b = arrayCreate (shMax, short) ;
  int *lng, *mark, maxCode ;
  Array aa = 0 ; 
  int n0, n1, n10, n2, n3, n4 ;
  BOOL debug = FALSE ;

  ctfCompress3Init (&aa, &lng, &mark, &maxCode) ;

  array (b, shMax - 1, short) = 0 ;  /* make room */
  sp = arrp (b, 0, short) ;
  spMax = sp + shMax ;
  n0 = n1 = n10 = n2 = n3 = n4 = 0 ;

  while (i-- && sp < spMax)
    {
      cc = *cp++ ;
      mode = cc & 0x80 ; arg = cc & 0x7f ;
      switch (mode)
	{
	case 0: /* initial zero = string of zero */
	  if (arg <= 0) /* should not happen */
	    { 
	      fprintf (stderr,"bad decompress3") ; 
	      goto abort ;
	    }  
	  n1++ ; n0 += arg ;
	  while (arg-- && sp < spMax) *sp++ = 0 ; 
	  break ;
	case 0x80:
	  switch (arg)
	    {
	    case 127:   /* next 2 bytes is a short */
	      i -= 2 ;  /* I need 3 bytes to code a short */
	      cc1 = *cp++ ; cc2 = *cp++ ;
	      *sp++ = (cc1 << 8) | cc2 ;
	      n4++ ;
	      break ;
	    case 126:   /* next byte is a byte */
	      i-- ;     /* I need 2 bytes to code a char */
	      cc1 = *cp++ ;
	      *sp++ = cc1 - 128 ;
	      n3++ ;
	      break ;
	    default:   /* 7 bytes is a code */
	      n = lng[arg] ; 
	      n2++ ; n10 += n ;
	      wp = arrp (aa, mark[arg], short) ;
	      while (n-- && sp < spMax) *sp++ = *wp++ ;
	      break ;
	    }
	}
    }
  if (debug)
    printf (" // compress3:\n//found   %d zeros in %d bytes, %d values coded in %d byte, %d bytes, %d short. \n// Total %d char for %d shrt\n",
	    n0, n1, n10, n2, n3, n4, n1 + n2 + n3 + n4, n0 + n10 + n3 + n4) ;
  
  if (i != -1 || sp != spMax)
    goto abort ;
  return b ;

abort:
  arrayDestroy (b) ;
  return 0 ;
}

/**********************************************************/
/**********************************************************/

static Array ctfDecompress (int compressionMode, 
			    int dataMax, int traceMax, 
			    unsigned char **cpp)
{ 
  Array a = 0 ;

  switch (compressionMode)
    {
    case 0:
      a = ctfDecompress0 (dataMax, 4*traceMax, *cpp) ;
      break ;
    case 1:
      a = ctfDecompress1 (dataMax, 4*traceMax, *cpp) ;
      break ;
    case 2:
      a = ctfDecompress2 (dataMax, 4*traceMax, *cpp) ;
      break ;
    case 3:
      a = ctfDecompress3 (dataMax, 4*traceMax, *cpp) ;
      break ;
    default:  /* unknown compression mode */
      break ;
    }

  *cpp += dataMax ;
  return a ;
}

/**********************************************************/

static Array ctfCompress (int compressionMode, Array a)
{
  switch (compressionMode)
    {
    case 0:
      return ctfCompress0 (a) ;
    case 1:
      return ctfCompress1 (a) ;
    case 2:
      return ctfCompress2 (a) ;
    case 3:
      return ctfCompress3 (a) ;
    default:
      fprintf (stderr,"FATAL ERROR in ctfCompress, Non existing compression mode") ;
      exit (1) ;
      return 0 ; /* for compiler happiness */
    }
}

/**********************************************************/
/**********************************************************/
/* called by saucisse fill, a system to test the efficiency of the system */

static Array ctfDecorrelate (Read *read, int predictionMode)
{
  int j, j1, u1, u2, u3, u4 ;
  short *zp, z = 0 ;
  TRACE *tt[4], *sp ;
  int traceMax = read->NPoints ;
  Array a = arrayCreate (4 * traceMax, short) ;

  if (predictionMode == -1)
    predictionMode = PREDICTIONMODE ;

  tt[0] = read->traceA ;
  tt[1] = read->traceG ;
  tt[2] = read->traceC ;
  tt[3] = read->traceT ;

  array (a, 4 * traceMax - 1 , short) = 0 ;  /* make room */
  zp = arrp (a, 0, short) ;
  for (j1 = 0 ; j1 < 4 ; j1++)
    { 
      sp = tt[j1] ;

      u1 = u2 = u3 = u4 = 0 ;
      for (j=0 ; j < traceMax ; zp++, sp++, j++)
	{ 
	  switch (predictionMode)
	    {
	    case 1: z = u1 ; break ; /* predict flat, transmit derivative */
	    case 2: z = 2*u1 - u2 ; break ; /* predict line trans dd2 */
	    case 3: z = 3*u1 - 3*u2 + u3 ; break ; /* predict parabole */
	    case 4: z = 4*u1 - 6*u2 + 4*u3 - u4; break ; /* overpredict ! */
	    case 0: 
	    default: z = 0 ; break ; /* predict zero, transmit value */
	    }
	  u4 = u3  ; u3 = u2 ; u2 = u1 ; u1 = *sp ;
	  *zp = u1 - z ;
	}
    }
  return a ;
}

/**********************************************************/

static BOOL ctfRecorrelate (Read *read, int predictionMode, Array a)
{ 
  int j, j1, u1, u2, u3, u4, z = 0 ;
  short *zp ;
  TRACE *sp, *tt[4] ;
  int traceMax = read->NPoints ;

  if (!a || arrayMax(a) != 4 * traceMax)
    return FALSE ;

  for (j1 = 0 ; j1 < 4 ; j1++) tt[j1] = 0 ; /* to allow harmless abort */
  zp = arrp (a, 0, short) ;
  for (j1 = 0 ; j1 < 4 ; j1++)
    { 
      /* staden's allocation system */
      sp = tt[j1] = (TRACE *) malloc(traceMax * sizeof (TRACE)) ;
      memset (sp, 0, traceMax * sizeof (TRACE)) ;
      u1 = u2 = u3 = u4 = 0 ;
      for (j=0 ; j < traceMax ; zp++, sp++, j++)
	{ 
	  switch (predictionMode)
	    {
	    case 1: z = u1 ; break ; /* predict flat, transmit derivative */
	    case 2: z = 2*u1 - u2 ; break ; /* predict line trans dd2 */
	    case 3: z = 3*u1 - 3*u2 + u3 ; break ; /* predict parabole */
	    case 4: z = 4*u1 - 6*u2 + 4*u3 - u4; break ; /* overpredict ! */
	    case 0: 
	    default: z = 0 ; break ; /* predict zero, transmit value */
	    }
	  u4 = u3  ; u3 = u2 ; u2 = u1 ; u1 = *sp = z + *zp ;
	}
    }
  read->traceA  = tt[0];
  read->traceG  = tt[1];
  read->traceC  = tt[2];
  read->traceT  = tt[3];

  return TRUE ;
}

/**********************************************************/
/**********************************************************/
/* returns 0: no probability available
           1: single base proba or equal proba on all bases,
              in this case, fill mixProba
           2: independant proba for the various bases
	   */

static int ctfProbInfoLevel  (Read *read, unsigned char *mixProb)
{
  int i, a, t, g, c, n, probInfoLevel = 0 ;

  n = read->NBases ;
  probInfoLevel = 0 ; i = 0 ;	  
  if (read->prob_A && read->prob_C && read->prob_G && read->prob_T)
    while (probInfoLevel < 4 && n--)
      {
	i = 0 ; a = t = g = c = 0 ;
	if ((a = read->prob_A [n])) { mixProb[n] = read->prob_A [n] ; i++ ; }
	if ((c = read->prob_C [n])) { mixProb[n] = read->prob_C [n] ; i++ ; }
	if ((g = read->prob_G [n])) { mixProb[n] = read->prob_G [n] ; i++ ; }
	if ((t = read->prob_T [n])) { mixProb[n] = read->prob_T [n] ; i++ ; }
	
	switch (i)
	  {
	  case 0: break ;
	  case 1: probInfoLevel = 1 ; break ;
	  case 4: if (a == c && a == g && a == t) { probInfoLevel = 1 ; break ; }
	              /* else fall to default */
	  default: probInfoLevel = 4 ; break ;
	  } 
	/*
	  if (p++ < 12 ) fprintf (stderr, "probInfoLevel %d %d %d %d %d \n", 
		 probInfoLevel, a, c, g, t) ;
		 */
      }
  /*  fprintf (stderr, "probInfoLevel %d\n", probInfoLevel) ; */
  return probInfoLevel ;  
}

/**********************************************************/

static Array  ctfPackTraces (Read *read)
{ 
  signed int x, dx;
  int section, sectionLength ;
  int i, n, dataMax, probInfoLevel, 
    traceMax = read->NPoints, baseMax = read->NBases,
    safe0 = 0, safe = 50 ;
  Array a = 0 , a1 = 0, a2 = 0 ;
  unsigned char *cp, *b ;
  unsigned short *bp ;
  unsigned char * mixProb = 0 ; 
  unsigned char *cq ;
  Array bb = 0 ;
  TRACE *ap, *bbp ;
  
  ctfTracePeakValue (read) ; /* sets read-> maxTraceVal */
  mixProb = (unsigned char *) malloc (read->NBases) ;

  probInfoLevel = ctfProbInfoLevel (read, mixProb) ;

lao:  /* the idea is that i will never have to loop */
  safe = 12 * traceMax + 2 * baseMax + probInfoLevel * baseMax + safe0 * traceMax + 64 ;
  if (a) arrayDestroy (a) ;
  a = arrayCreate (safe, unsigned char) ;
  array (a, safe + 150, unsigned char) = 0 ; /* make room */
  cp = arrp (a, 0, unsigned char) ;
  

  ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;

  /* SECTION 11 bases  */
  if (!ctfType || strchr(ctfType, '1'))
    {
      section = 11 ;
      sectionLength = baseMax + 8 ;
      n = baseMax ; 
      ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
      
      ctfStoreInt (cp, baseMax) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
      
      n = baseMax ; safe -= n ;
      b = (unsigned char *)read->base ; 
      memcpy (cp, b, n) ; cp += n ;
      
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
    }

  /* SECTION 12 base positions, upgrades old section 1 */
  if (!ctfType || strchr(ctfType, '1'))
    {
      section = 12 ;
      sectionLength = baseMax + 8 ;
      n = baseMax ; x = 0 ;
      bp = read->basePos ;
      while (n--)
	{
	  dx = *bp++ - x + 32 ;
	  if (dx < 0 || dx > 255)
	    { sectionLength += 2 ; safe -= 2 ; }
	}
      ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
      
      ctfStoreInt (cp, baseMax) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
      
      n = baseMax ; safe -= n ;
      x = 0 ;
      bp = read->basePos ;
      
      while (n--)
	{ 
	  dx = *bp++ - x + 32 ;
	  /* these will smoothe away after a few steps */
	  if (dx < 0)
	    {
	      *cp++ = 254 ; dx = -dx ;
	      *cp++ = (dx >> 8) & 0xff;
	      *cp++ = (dx >> 0) & 0xff;
	    }
	  else if (dx < 254)
	    *cp++ = dx;
	  else /* dx >= 254 */
	    {
	      *cp++ = 255 ;
	      *cp++ = (dx >> 8) & 0xff;
	      *cp++ = (dx >> 0) & 0xff;
	    }
	  x += dx - 32; 
	}
      
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
    }

  /* SECTION 2 the traces */
  if (!ctfType || strchr(ctfType, '2'))
    {
      a1 = ctfDecorrelate (read, PREDICTIONMODE) ;
      a2 = ctfCompress (COMPRESSIONMODE, a1) ;
      dataMax = arrayMax (a2) ;
      
      if (1)  /* debugging */
	{
	  cq = arrp (a2, 0, unsigned char) ;
	  bb = ctfDecompress (COMPRESSIONMODE, dataMax, traceMax, &cq) ;
	  i = 4*traceMax ;
	  
	  ap = arrp (a1, 0, TRACE) ;
	  bbp = arrp (bb, 0, TRACE) ;
	  while (i--)
	    if (*ap++ != *bbp++)
	      {
		fprintf (stderr, 
			 "FATAL ERROR bad compress decompress at i = %d\n", i) ;
		exit (1) ;
	      }
	  arrayDestroy (bb) ;
	}
      
      section = 2 ;
      sectionLength = 16 + dataMax ;
      ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
      
      ctfStoreInt (cp, PREDICTIONMODE) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, COMPRESSIONMODE) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, traceMax) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, dataMax) ; cp += 4 ; safe -= 4 ;
      
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
      
      if (safe > dataMax)
	memcpy (cp, arrp(a2, 0, unsigned char), (mysize_t) dataMax) ;
      cp += dataMax ; safe -= dataMax ;
      arrayDestroy (a1) ;
      arrayDestroy (a2) ;
      if (safe < 0)
	{ safe0++ ; goto lao ; }
      
      /* end section */
      if (safe < 12)
	{ safe0++ ; goto lao ; }
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
    }

  
  /* SECTION 3 miscelaneous info */
  if (!ctfType || strchr(ctfType, '3'))
    {
      section = 3 ;
      
      if (read->info)
	n = strlen(read->info) ;
      else
	n = 0 ;
      if (!read->rightCutoff)
	read->rightCutoff = read->NBases + 1 ;
      
      sectionLength = 20 + n ;
      ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
      
      ctfStoreInt (cp, read->leftCutoff) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, read->rightCutoff) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, read->orig_trace_format) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, read->maxTraceVal) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, n) ; cp += 4 ; safe -= 4 ;
      if (n > 0 && safe > n) 
	{ strncpy ((char *)cp, read->info, n) ; cp += n ; } ; 
      safe -= n ;
      /* end section */
      if (safe < 12)
	{ safe0++ ; goto lao ; }
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
    }

  /* SECTION 4 Probability information */  
  if (!ctfType || strchr(ctfType, '4'))
    {
      switch (probInfoLevel)
	{
	case 4:
	  section = 4 ;
	  sectionLength = 4 * baseMax + 4 ; n = baseMax ;
	  ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
	  ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
	  
	  ctfStoreInt (cp, baseMax) ; cp += 4 ; safe -= 4 ;
	  if (safe > n) { memcpy (cp, read->prob_A, n) ; cp += n ;} 
	  safe -= n ; 
	  if (safe > n) { memcpy (cp, read->prob_C, n) ; cp += n ; }
	  safe -= n ; 
	  if (safe > n) { memcpy (cp, read->prob_G, n) ; cp += n ; }
	  safe -= n ; 
	  if (safe > n) { memcpy (cp, read->prob_T, n) ; cp += n ; }
	  safe -= n ; 
	  if (safe < 12)
	    { safe0++ ; goto lao ; }

	  ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
	  break ;
	case 1:
	  section = 5 ;
	  sectionLength = baseMax + 4 ; n = baseMax ;
	  ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
	  ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
	  
	  ctfStoreInt (cp, baseMax) ; cp += 4 ; safe -= 4 ;
	  if (safe > n) { memcpy (cp, mixProb, n) ; cp += n ; safe -= n ; }
	  if (safe < 12)
	    { safe0++ ; goto lao ; }
	  
	  ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
	  break ;
	}
    }

  /* SECTION 999 the end */
  if (TRUE)
    {
      section = 999 ;
      sectionLength = 0 ;
      ctfStoreInt (cp, section) ; cp += 4 ; safe -= 4 ;
      ctfStoreInt (cp, sectionLength) ; cp += 4 ; safe -= 4 ;
      
      ctfStoreInt (cp, MAGIC) ; cp += 4 ; safe -= 4 ;
      
      if (safe < 0)
	{ safe0++ ; goto lao ; }
    }

  arrayMax (a) = cp - arrp (a, 0, unsigned char) ;

  free (mixProb) ;
  return a ;  
}

/**********************************************************/

static void ctfUnmixProb (Read *read, int n, unsigned char *mixProb) 
{
  unsigned char *cp, *cq, *ca, *cg, *ct, *cc ;

  ca = (unsigned char *)(read->prob_A = (char *) malloc (n)) ;
  cc = (unsigned char *)(read->prob_C = (char *) malloc (n)) ;
  cg = (unsigned char *)(read->prob_G = (char *) malloc (n)) ;
  ct = (unsigned char *)(read->prob_T = (char *) malloc (n)) ;
  cp = (unsigned char *)read->base ;
  cq = mixProb ;
  while (n--)
    switch (*cp++)
      {
      case 'A': case 'a':
	*ca++ = *cq++ ; *cc++ = *cg++ = *ct++ = 0 ;
	break ;
      case 'C': case 'c':
	*cc++ = *cq++ ; *cg++ = *ct++ = *ca++ = 0 ;
	break ;
      case 'G': case 'g':	
	*cg++ = *cq++ ; *ct++ = *ca++ = *cc++ = 0 ;
	break ;
      case 'T': case 't':	
	*ct++ = *cq++ ; *ca++ = *cc++ = *cg++ = 0 ;
	break ;
      default:
	*ca++ = *cc++ = *cg++ = *ct++ = *cq++ ;
	break ;
      }
}

/**********************************************************/

#define CHECKMAGIC   \
magic = ctfGetInt (cp) ; \
cp += 4 ; nn -= 4 ; \
if (magic != MAGIC) \
{ \
  fprintf (stderr, "Error reading compressed trace file, sorry\n") ;  \
  goto abort ; \
}


static BOOL ctfUnPackTraces (Read *read, Array a)
{ 
  int compressionMode, predictionMode, magic ;
  int section, sectionLength ;
  int n, nn, traceMax, baseMax, dataMax, nMixProb = 0 ;
  signed int x, dx;
  Array  decompressedData = 0 ;
  unsigned char *cp, *b ;
  unsigned short *bp ;
  unsigned char *ucp, *mixProb = 0 ;

  if (!arrayExists (a))
    return FALSE ;

  cp = arrp (a, 0, unsigned char) ; nn = arrayMax (a) ;

  while (TRUE)
    {
      CHECKMAGIC ; 
      section = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
      sectionLength = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
      switch (section)
	{
	case 999:  /* end of record */
	  goto done ;
	case 1:  /* read the bases  and positions, problem if negative dx */
	  baseMax = read->NBases = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  CHECKMAGIC ;

	  n = baseMax ; nn -= 2 * n ;

	  read->base = (char *)malloc (baseMax) ; /* staden's alloc */
	  b = (unsigned char *)read->base;
	  bp = read->basePos = (unsigned short *) malloc (baseMax * sizeof (unsigned short)) ;
	  
	  x = 0 ;
	  
	  while (n--)
	    { x += *cp++ ;
	    *bp++ = x ;
	    switch ((*cp++) & 0xf)
	     {
	     case A_: *b++ = 'A' ; break ;
	     case T_: *b++ = 'T' ; break ;
	     case G_: *b++ = 'G' ; break ;
	     case C_: *b++ = 'C' ; break ;
	     default: *b++ = 'N' ; break ;
	     }
	    }
	  break ;
	  
	case 11:  /* read the basecall, fullproof method */
	  baseMax = read->NBases = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  CHECKMAGIC ;

	  n = baseMax ; nn -= n ;
	  b = (unsigned char *)(read->base = (char *)malloc (baseMax)) ; 
	  memcpy (b, cp, n) ;
	  cp += n ;
	  break ;
	  
	case 12:  /* read the baspos, fullproof method */
	  baseMax = read->NBases = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  CHECKMAGIC ;

	  n = baseMax ; nn -= n ;

	  bp = read->basePos = (unsigned short *) malloc (baseMax * sizeof (unsigned short)) ;
	  
	  x = 0 ;
	  while (n--)
	    { 
	      dx = (unsigned char) *cp++;
	      if (dx == 254) 
		{
		  dx = (cp[0] << 8) | cp[1];
		  cp += 2; nn -= 2 ;
		  dx = -dx ;
		}
	      else if (dx == 255) 
		{
		  dx = (cp[0] << 8) | cp[1];
		  cp += 2; nn -= 2 ;
		}
	      dx -= 32 ;
	      x += dx;
	      *bp++ = x ;
	    }
	  break ;
	  
	case 2:  /* read the traces */
	  predictionMode = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  compressionMode = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  traceMax = read->NPoints = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  dataMax = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  CHECKMAGIC ;
	  decompressedData = ctfDecompress (compressionMode, dataMax, traceMax, &cp) ;
	  if (!ctfRecorrelate (read, predictionMode, decompressedData))
	    goto abort ;	  
	  arrayDestroy (decompressedData) ;
	  break ;
	  
	case 3:  /* read miscelaneous info */
	  read->leftCutoff = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  read->rightCutoff = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  read->orig_trace_format = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ; 
	  read->maxTraceVal = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  n = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  if (n > 0)
	    {
	      read->info = (char *) malloc (n+1) ; 
	      strncpy (read->info, (char *)cp, n) ; cp += n ; nn -= n ;
	      read->info[n] = 0 ; /* zero terminate the string */
	    }
	  break ;

	case 4:  /*  Probability information */  
	  n = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  read->prob_A = (char *) malloc (n) ;
	  read->prob_C = (char *) malloc (n) ;
	  read->prob_G = (char *) malloc (n) ;
	  read->prob_T = (char *) malloc (n) ;
	  memcpy (read->prob_A, cp, n) ; cp += n ; nn -= n ;
	  memcpy (read->prob_C, cp, n) ; cp += n ; nn -= n ;
	  memcpy (read->prob_G, cp, n) ; cp += n ; nn -= n ;
	  memcpy (read->prob_T, cp, n) ; cp += n ; nn -= n ; 
	  break ;

	case 5:  /*  Mixed Probability information */  
	  n = ctfGetInt (cp) ; cp += 4 ; nn -= 4 ;
	  mixProb = (unsigned char*) malloc (n) ; nMixProb = n ;
	  ucp = cp ; /* to please memcpy */
	  memcpy (mixProb, ucp, n) ; cp += n ; nn -= n ; 
	  break ;

	case 6:
	  /* The original input format data, used in staden's
	     int orig_trace_format;
	     void *orig_trace;
	     
	     i do not yet support that because i am not sure it is used
	     so i
	     */

	  /*   fallthru  to default */

	default: /* not yet defined stuff */
	  cp += sectionLength ; 
	  break ;
	}
    }

 done:
  if (mixProb)
    {
      if (read->NBases == nMixProb)
	ctfUnmixProb (read, nMixProb, mixProb) ; 
      else
	fprintf (stderr, "mixProb problem, read->NBases = %d nMixProb = %d",
		   read->NBases, nMixProb ) ;
    }
  if (!read->rightCutoff)
    read->rightCutoff = read->NBases + 1 ;
  

  if (mixProb) free (mixProb) ;
  CHECKMAGIC ;   /* terminal CHECK */
  return TRUE ;

 abort:
  arrayDestroy (decompressedData) ;

  freeSeq (read) ;
  return FALSE ;
}

/**********************************************************/
/**********************************************************/
/**********************************************************/

static BOOL ctfWriteTrace (FILE *ff, Array a)
{
  int n ; char *cp ;

  n = arrayMax (a) ;  cp = arrp (a, 0, char) ;
  fwrite (cp, n, 1, ff) ;
  return TRUE ;
} 

/**********************************************************/

static Array ctfReadTrace (FILE *fil)
{
  unsigned int i = 0, nr, nb = 100000, size = 1 ; 
  unsigned char *cp ;
  Array a = arrayCreate (nb, unsigned char) ;

  do
    {
      array(a,(++i)*nb,unsigned char) = 0 ; /* to create space */
      cp = arrp(a,nb*(i-1),unsigned char) ; /* possible relocation */
    }
  while ((nr = fread (cp, size, nb, fil)) == nb) ;
  
  arrayMax(a) -= nb - nr; /* artificial space removed */
 
  if(!arrayMax(a))
    arrayDestroy(a) ;

  return a ;
} 

/**********************************************************/
/**********************************************************/
/*****   interaction with staden 's makeSCF ***************/
/**********************************************************/

int ctfFWrite (FILE *ff, Read *read) 
{
  int result = -1 ; /* assume error */

  if (read &&  read->NBases && read->NPoints && ff)
    {
      Array a = ctfPackTraces (read);
      
      ctfWriteTrace (ff, a) ;
      arrayDestroy (a) ;
  
      result = 0 ; /* success */
    }

  return result ;
} /* ctfFWrite */

/**********************************************************/

Read *ctfFRead (FILE *ff)
{
  Read * read = 0 ;
  int NBases = 0 ;
  Array a = 0 ;

  if ((a = ctfReadTrace (ff)) &&
      (read =  (Read *) malloc (sizeof(Read))))
    {
      memset (read, 0, sizeof(Read)) ;
      if (ctfUnPackTraces (read, a))
	{
	  read->format = TT_CTF ;
	  NBases = read->NBases ;
	  
	  if (!read->prob_A) 
	    {
	      read->prob_A = (char *) malloc (NBases) ;
	      if (!read->prob_C) read->prob_C = (char *) malloc (NBases) ;
	      if (!read->prob_G) read->prob_G = (char *) malloc (NBases) ;
	      if (!read->prob_T) read->prob_T = (char *) malloc (NBases) ;
	      memset (read->prob_A, 0, NBases) ;
	      memset (read->prob_C, 0, NBases) ;
	      memset (read->prob_G, 0, NBases) ;
	      memset (read->prob_T, 0, NBases) ;
	    }
	  
	  read->orig_trace = 0x0 ;
	}
      else
	read = 0 ;
    }

 arrayDestroy (a) ;
 return read ;
}  /* ctfFRead */

/* Examples
   run -s -any /users/mieg/CTFtest/tt/a.scf -ctfout a.scf.ctf
   run -s -any /users/mieg/CTFtest/tt/a.scf.ctf -output a.scf.ctf.scf
   */

/**********************************************************/
/**********************************************************/

