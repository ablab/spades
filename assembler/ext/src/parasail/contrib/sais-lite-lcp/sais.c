/*
 * sais.c for sais-lite
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sais.h"

#ifndef UCHAR_SIZE
# define UCHAR_SIZE 256
#endif
#ifndef MINBUCKETSIZE
# define MINBUCKETSIZE 256
#endif

#define sais_index_type int
#define sais_bool_type  int
#define SAIS_LMSSORT2_LIMIT 0x3fffffff

#define SAIS_MYMALLOC(_num, _type) ((_type *)malloc((_num) * sizeof(_type)))
#define SAIS_MYFREE(_ptr, _num, _type) free((_ptr))
#define chr(_a) (cs == sizeof(sais_index_type) ? ((sais_index_type *)T)[(_a)] : ((unsigned char *)T)[(_a)])

/* qsort int comparison function */ 
int int_cmp(const void *a, const void *b) 
{ 
  const int *ia = (const int *)a; /* casting pointer types  */
  const int *ib = (const int *)b;
  return *ia  - *ib; 
} 

/* find the start or end of each bucket */
static
void
getCounts(const void *T, sais_index_type *C, sais_index_type n, sais_index_type k, int cs) {
  sais_index_type i;
  for(i = 0; i < k; ++i) { C[i] = 0; }
  for(i = 0; i < n; ++i) { ++C[chr(i)]; }
}
static
void
getBuckets(const sais_index_type *C, sais_index_type *B, sais_index_type k, sais_bool_type end) {
  sais_index_type i, sum = 0;
  if(end) { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum; } }
  else { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum - C[i]; } }
}

/* sort all type LMS suffixes */
static
void
LMSsort1(const void *T, sais_index_type *SA,
         sais_index_type *C, sais_index_type *B,
         sais_index_type n, sais_index_type k, int cs) {
  sais_index_type bb, i, j;
  sais_index_type c0, c1;

  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  bb = B[c1 = chr(j)];
  --j;
  SA[bb++] = (chr(j) < c1) ? ~j : j;
  for(i = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      assert(chr(j) >= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = bb; bb = B[c1 = c0]; }
      assert(i < bb);
      --j;
      SA[bb] = (chr(j) < c1) ? ~j : j;
      ++bb;
      SA[i] = 0;
    } else if(j < 0) {
      SA[i] = ~j;
    }
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, bb = B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      assert(chr(j) <= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = bb; bb = B[c1 = c0]; }
      assert((bb) <= i);
      --j;
      SA[--bb] = (chr(j) > c1) ? ~(j + 1) : j;
      SA[i] = 0;
    }
  }
}
static
sais_index_type
LMSpostproc1(const void *T, sais_index_type *SA,
             sais_index_type n, sais_index_type m, int cs) {
  sais_index_type i, j, p, q, plen, qlen, name;
  sais_index_type c0, c1;
  sais_bool_type diff;

  /* compact all the sorted substrings into the first m items of SA
      2*m must be not larger than n (proveable) */
  assert(0 < n);
  for(i = 0; (p = SA[i]) < 0; ++i) { SA[i] = ~p; assert((i + 1) < n); }
  if(i < m) {
    for(j = i, ++i;; ++i) {
      assert(i < n);
      if((p = SA[i]) < 0) {
        SA[j++] = ~p; SA[i] = 0;
        if(j == m) { break; }
      }
    }
  }

  /* store the length of all substrings */
  i = n - 1; j = n - 1; c0 = chr(n - 1);
  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
  for(; 0 <= i;) {
    do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
    if(0 <= i) {
      SA[m + ((i + 1) >> 1)] = j - i; j = i + 1;
      do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
    }
  }

  /* find the lexicographic names of all substrings */
  for(i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
    p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
    if((plen == qlen) && ((q + plen) < n)) {
      for(j = 0; (j < plen) && (chr(p + j) == chr(q + j)); ++j) { }
      if(j == plen) { diff = 0; }
    }
    if(diff != 0) { ++name, q = p, qlen = plen; }
    SA[m + (p >> 1)] = name;
  }

  return name;
}
static
void
LMSsort2(const void *T, sais_index_type *SA,
         sais_index_type *C, sais_index_type *B, sais_index_type *D,
         sais_index_type n, sais_index_type k, int cs) {
  sais_index_type *b, i, j, t, d;
  sais_index_type c0, c1;
  assert(C != B);

  /* compute SAl */
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  --j;
  t = (chr(j) < c1);
  j += n;
  *b++ = (t & 1) ? ~j : j;
  for(i = 0, d = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      if(n <= j) { d += 1; j -= n; }
      assert(chr(j) >= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      --j;
      t = c0; t = (t << 1) | (chr(j) < c1);
      if(D[t] != d) { j += n; D[t] = d; }
      *b++ = (t & 1) ? ~j : j;
      SA[i] = 0;
    } else if(j < 0) {
      SA[i] = ~j;
    }
  }
  for(i = n - 1; 0 <= i; --i) {
    if(0 < SA[i]) {
      if(SA[i] < n) {
        SA[i] += n;
        for(j = i - 1; SA[j] < n; --j) { }
        SA[j] -= n;
        i = j;
      }
    }
  }

  /* compute SAs */
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, d += 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      if(n <= j) { d += 1; j -= n; }
      assert(chr(j) <= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      --j;
      t = c0; t = (t << 1) | (chr(j) > c1);
      if(D[t] != d) { j += n; D[t] = d; }
      *--b = (t & 1) ? ~(j + 1) : j;
      SA[i] = 0;
    }
  }
}
static
sais_index_type
LMSpostproc2(sais_index_type *SA, sais_index_type n, sais_index_type m) {
  sais_index_type i, j, d, name;

  /* compact all the sorted LMS substrings into the first m items of SA */
  assert(0 < n);
  for(i = 0, name = 0; (j = SA[i]) < 0; ++i) {
    j = ~j;
    if(n <= j) { name += 1; }
    SA[i] = j;
    assert((i + 1) < n);
  }
  if(i < m) {
    for(d = i, ++i;; ++i) {
      assert(i < n);
      if((j = SA[i]) < 0) {
        j = ~j;
        if(n <= j) { name += 1; }
        SA[d++] = j; SA[i] = 0;
        if(d == m) { break; }
      }
    }
  }
  if(name < m) {
    /* store the lexicographic names */
    for(i = m - 1, d = name + 1; 0 <= i; --i) {
      if(n <= (j = SA[i])) { j -= n; --d; }
      SA[m + (j >> 1)] = d;
    }
  } else {
    /* unset flags */
    for(i = 0; i < m; ++i) {
      if(n <= (j = SA[i])) { j -= n; SA[i] = j; }
    }
  }

  return name;
}

/* compute SA and BWT */
static void induceSA(const void *T, sais_index_type *SA,
		     sais_index_type *C, sais_index_type *B,
		     sais_index_type n, sais_index_type k, int cs) {
  sais_index_type i, j;
  sais_index_type bb;
  sais_index_type c0, c1;
  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  bb = B[c1 = chr(j)];
  SA[bb++] = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
    j = SA[i], SA[i] = ~j;
    if(0 < j) {
      --j;
      assert(chr(j) >= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = bb; bb = B[c1 = c0]; }
      assert(i < bb);
      SA[bb] = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
      ++bb;
    }
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, bb = B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      --j;
      assert(chr(j) <= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = bb; bb = B[c1 = c0]; }
      assert(bb <= i);
      SA[--bb] = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
    } else {
      SA[i] = ~j;
    }
  }
}

static void induceSAandLCP(const void *T, sais_index_type *SA,
			   sais_index_type *LCP,
			   sais_index_type *C, sais_index_type *B,
			   sais_index_type n, sais_index_type k, int cs) {
  /*
    When entering this procedure, we are in the following situation:
    all S*-suffixes have been sorted and put at the end of their
    corresponding buckets in SA. Further, all their LCP-values have been
    computed (in LCP). A value of -1 in LCP means that "this is the first
    S*-suffix in its bucket." As in in the inducing step an L-suffix
    can be placed before the leftmost S*-suffix, this means that the actual
    LCP-value at this "L/S-seam" has to be recomputed. All other uncomputed
    LCP-entries are marked -2.
  */

  sais_index_type i, bb;         /* indices in SA/LCP (origin/target) */
  sais_index_type j;             /* position in text */
  sais_index_type c0, c1;        /* characters (new/last) */
  
  sais_index_type lcp;           /* LCP-value */
  sais_index_type l;             /* for finding LCP at L/S-seam */
  sais_index_type start, end, stack_end; /* for inducing the LCP-values */

  sais_index_type *D; /* store beginnings of buckets (not CURRENT beginnings!) */
  sais_index_type *LastW; /* store last written L or S-suffix for every bucket */
  sais_index_type sigma;       /* (true) alphabet size */
  sais_index_type *TranslateSigma; /* general to effective alphabet ([0..k-1] |--> [0..sigma-1]) */
  sais_index_type *LastOcc; /* store last occurrences of characters */
  const sais_index_type stack_extra_space = 1024; /* this is the same Min-Stack as in Gog's sdsl */
  sais_index_type stack_size;
  sais_index_type *MinStack; /* Min-Stack */

  if ((D = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { exit(-1); } /* TODO: check if D is necessary!!! (first write to bucket=>0) */
  if ((LastW = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { exit(-1); }
  for (i = 0; i < k; ++i) LastW[i] = n-1; /* point to $ */
  /* todo: move memory management to sais_main */

  /* compute SAl and LCPl*/
  if(C == B) { getCounts(T, C, n, k, cs); }  /* re-calculate character counts */
  getBuckets(C, B, k, 0);                    /* find starts of buckets */
  memcpy(D,B,k*sizeof(sais_index_type));     /* store starts of buckets */
  j = n - 1;                                 /* go to last character $ */
  bb = B[c1 = chr(j)];        /* bb = position in induced bucket */
  LCP[bb]  = 0;               /* set LCP-value of $ (first value in bucket => 0) */
  SA[bb++] = (chr(j - 1) < c1) ? ~j : j; /* put last character $ into its bucket */
                        /* negative values mean "don't induce from here anymore" */

  /* in case LCP[0], which is always 0, hadn't been set yet */
  LCP[0] = 0;

  /* Variant 3: stack */
  sigma = 0;       /* (true) alphabet size */
  if ((TranslateSigma = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { exit(-1); }
  for (i = 0; i < k; ++i)  {   /* calculate effective alphabet size */
    TranslateSigma[i] = sigma; /* (also stores values for unused characters) */
    if (C[i] > 0) ++sigma;     /* count characters */
  }
  if ((LastOcc = SAIS_MYMALLOC(sigma, sais_index_type)) == NULL) { exit(-1); }
  for (i = 0; i < sigma; ++i) LastOcc[i] = -1; /* init with impossible values */

  stack_size = 2 * (stack_extra_space + sigma + 4);
  if ((MinStack = SAIS_MYMALLOC(stack_size+4, sais_index_type)) == NULL) { exit(-1); }
  MinStack[0] = -1; MinStack[1] = -1; /* (pos, LCP-value) */
  stack_end = 1;

  for(i = 0; i < n; ++i) {
    j = SA[i], SA[i] = ~j;
    if(0 < j) { /* otherwise don't induce anymore from j */
      lcp = LCP[i];
      if (lcp == -1) {
	/* here we are at the seam between L and S in the same bucket */
	c0 = chr(j);  /* i's bucket */
	lcp = 0;
	while (chr(j+lcp) == chr(LastW[c0]+lcp)) lcp++; /* naive LCP-computation (overall linear!) */
	/* no need to store LCP[i]=lcp (will be re-calculated in right-to-left scan!) */
      }
      --j;      /* move to suffix T[SA[i]-1] */
      assert(chr(j) >= chr(j + 1));   /* induced suffix must be type L */
      if((c0 = chr(j)) != c1) {       /* induced SA-value in new bucket c0 */
	B[c1] = bb;                   /* store current end in old bucket */
	bb = B[c1 = c0];              /* go to position in new bucket */
      }
      assert(i < bb);                 /* can only induce to the right */
      LastW[c0] = j;                  /* store last written L-suffix for every bucket */
      SA[bb] = ((0 < j) && (chr(j - 1) < c0)) ? ~j : j;

      /* Variant 3: use stack: */
      assert(lcp >= 0);                                 /* lcp already computed */
      while (lcp <= MinStack[stack_end]) stack_end -= 2; /* pop from stack */
      MinStack[++stack_end] = i;   /* push position on stack */
      MinStack[++stack_end] = lcp; /* push lcp-value */

      start = LastOcc[TranslateSigma[c0]] + 1; /* start of query */
      assert(stack_end-3 >= 0);                /* stopper (-1) and last (i) are on stack */
      end = stack_end - 3;
      while (start <= MinStack[end]) end -= 2; /* search until smaller element found */
      if (bb == D[c0]) LCP[bb] = 0;            /* store 0 at bucket beginnings */
      else LCP[bb] = MinStack[end+3] + 1;      /* induce LCP-value! */

      LastOcc[TranslateSigma[c0]] = i; /* store origin of last occurrence of c0 */
      ++bb;                            /* advance in bucket */
    }
    else { /* don't induce, but update stack with LCP[i] */
      lcp = LCP[i];      /* get current LCP-value */
      assert(lcp != -1); /* -1 only for S*, but we induce from S* */
      if (lcp >= 0) {    /* check if already computed */
	while (lcp <= MinStack[stack_end]) stack_end -= 2; /* pop from stack */
	MinStack[++stack_end] = i;   /* push position on stack */
	MinStack[++stack_end] = lcp; /* push lcp-value */
      }
    }
    if (stack_end > stack_size) { /* re-adjust stack: */
      sais_index_type *LastOccCopy; /* Copy of LastOcc */
      if ((LastOccCopy = SAIS_MYMALLOC(sigma, sais_index_type)) == NULL) { exit(-1); }
      memcpy(LastOccCopy, LastOcc, sigma*sizeof(sais_index_type));
      qsort(LastOccCopy, sigma, sizeof(sais_index_type), int_cmp);

      end = 1;
      for (j = 0, l=2; j < sigma; ++j) {
	start = LastOccCopy[j] + 1; /* start of next largest query */
	if (start > MinStack[end-1]) {         /* otherwise correct element already taken */
	  while (l < stack_end && start > MinStack[l]) l += 2;
	  if (l > stack_end) break;
	  assert(l < stack_end);
	  MinStack[++end] = MinStack[l];       /* take first element >= start */
	  MinStack[++end] = MinStack[l+1];
	}
      }
      stack_end = end;
      free(LastOccCopy);
    }
  }

  /* compute SAs and LCPl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for (i = 0; i < sigma; ++i) LastOcc[i] = n-1; /* init with impossible values */
  MinStack[0] = n; MinStack[1] = -1; /* (pos, LCP-value) */
  stack_end = 1;

  for(i = n - 1, bb = B[c1 = 0]; 0 <= i; --i) {
    lcp = LCP[i];
    if (0 < i && lcp < 0 && LCP[i-1] >= 0) { /* calculate LCP at L/S-seam */
      j = SA[i];         /* go to suffix */
      if (j < 0) j = ~j; /* entry in SA could be negative => adjust */
      l = SA[i-1];       /* go to lex. previous suffix */
      if (l < 0) l = ~l; /* entry in SA could be negative => adjust */
      lcp = 0;
      while (chr(j+lcp) == chr(l+lcp)) lcp++; /* naive LCP-computation (overall linear!) */
      LCP[i] = lcp;                           /* this time set LCP-value at seam */
    }
    if(0 < (j = SA[i])) { /* induce SA and LCP */
      --j;                           /* go to suffix T[SA[i]-1] (to be induced) */
      assert(chr(j) <= chr(j + 1));  /* must be type S */
      if((c0 = chr(j)) != c1) {
	B[c1] = bb; bb = B[c1 = c0]; /* switch bucket */
      }
      assert(bb <= i); /* induced suffix must be placed to the left of i */
      SA[--bb] = ((j == 0) || (chr(j - 1) > c0)) ? ~j : j; /* continue if type L */
      assert(c0+1<k);        /* we cannot induce into the last bucket */

      /* search MinStack: */
      start = LastOcc[TranslateSigma[c0]];     /* end of query */
      assert(stack_end-1 >= 0);                /* stopper (-1) is on stack */
      end = stack_end - 1;
      while (start >= MinStack[end]) end -= 2; /* search until smaller element found */
      if (bb+1 == D[c0+1]) LCP[bb+1] = 0;      /* store 0 at bucket beginnings */
      else LCP[bb+1] = MinStack[end+3] + 1;    /* induce LCP-value! */
      if (bb+1 == i) lcp = LCP[i];             /* update if inducing changed my LCP-value */
      LastOcc[TranslateSigma[c0]] = i;         /* store origin of last occurrence of c0 */
    } else { /* don't induce */
      SA[i] = ~j;
    }
    /* update MinStack: */
    assert(lcp >= 0); /* LCP must already have been computed */
    while ((lcp <= MinStack[stack_end]) && (stack_end >= 0))  stack_end -= 2; /* pop from stack */
    MinStack[++stack_end] = i;   /* push position on stack */
    MinStack[++stack_end] = lcp; /* push lcp-value */

    if (stack_end > stack_size) { /* re-adjust stack: */
      sais_index_type *LastOccCopy; /* Copy of LastOcc */
      if ((LastOccCopy = SAIS_MYMALLOC(sigma, sais_index_type)) == NULL) { exit(-1); }
      memcpy(LastOccCopy, LastOcc, sigma*sizeof(sais_index_type));
      qsort(LastOccCopy, sigma, sizeof(sais_index_type), int_cmp);

      end = 1;
      for (j = sigma-1, l=2; j >= 0; --j) {
	start = LastOccCopy[j]; /* start of next largest query */
	if (start < MinStack[end-1]) {         /* otherwise correct element already taken */
	  while (l < stack_end && start < MinStack[l]) l += 2;
	  if (l > stack_end) break;
	  assert(l < stack_end);
	  MinStack[++end] = MinStack[l];       /* take first element >= start */
	  MinStack[++end] = MinStack[l+1];
	}
      }
      stack_end = end;
      free(LastOccCopy);
    }
  }

  free(D);
  free(LastOcc);
  free(MinStack);
  free(TranslateSigma);
  free(LastW);
}

static
sais_index_type
computeBWT(const void *T, sais_index_type *SA,
           sais_index_type *C, sais_index_type *B,
           sais_index_type n, sais_index_type k, int cs) {
  sais_index_type *b, i, j, pidx = -1;
  sais_index_type c0, c1;
  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      --j;
      assert(chr(j) >= chr(j + 1));
      SA[i] = ~((sais_index_type)(c0 = chr(j)));
      if(c0 != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
    } else if(j != 0) {
      SA[i] = ~j;
    }
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      --j;
      assert(chr(j) <= chr(j + 1));
      SA[i] = (c0 = chr(j));
      if(c0 != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      *--b = ((0 < j) && (chr(j - 1) > c1)) ? ~((sais_index_type)chr(j - 1)) : j;
    } else if(j != 0) {
      SA[i] = ~j;
    } else {
      pidx = i;
    }
  }
  return pidx;
}

/* find the suffix array SA of T[0..n-1] in {0..255}^n */
static sais_index_type sais_main(const void *T, sais_index_type *SA,
				 sais_index_type *LCP,
				 sais_index_type fs, sais_index_type n, sais_index_type k, int cs,
				 sais_bool_type isbwt,
				 sais_bool_type level0) { /* level0 = 1 iff recursion depth is 0 */
  sais_index_type *C, *B, *D, *RA, *PLCP, *PHI, *DELTA, *b;
  sais_index_type i, j, m, /* m: number of S*-suffixes */
    p, q, t, name, pidx = 0, newfs;
  sais_index_type c0, c1;
  unsigned int flags;

  assert((T != NULL) && (SA != NULL));
  assert((0 <= fs) && (0 < n) && (1 <= k));

  if(k <= MINBUCKETSIZE) {
    if((C = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { return -2; }
    if(k <= fs) {
      B = SA + (n + fs - k);
      flags = 1;
    } else {
      if((B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { SAIS_MYFREE(C, k, sais_index_type); return -2; }
      flags = 3;
    }
  } else if(k <= fs) {
    C = SA + (n + fs - k);
    if(k <= (fs - k)) {
      B = C - k;
      flags = 0;
    } else if(k <= (MINBUCKETSIZE * 4)) {
      if((B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { return -2; }
      flags = 2;
    } else {
      B = C;
      flags = 8;
    }
  } else {
    if((C = B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { return -2; }
    flags = 4 | 8;
  }
  if((n <= SAIS_LMSSORT2_LIMIT) && (2 <= (n / k))) {
    if(flags & 1) { flags |= ((k * 2) <= (fs - k)) ? 32 : 16; }
    else if((flags == 0) && ((k * 2) <= (fs - k * 2))) { flags |= 32; }
  }

  /* stage 1: reduce the problem by at least 1/2
     sort all the LMS-substrings */
  getCounts(T, C, n, k, cs); getBuckets(C, B, k, 1); /* find ends of buckets */

  for(i = 0; i < n; ++i) { SA[i] = 0; }
  b = &t; i = n - 1; j = n; m = 0; c0 = chr(n - 1);
  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
  for(; 0 <= i;) {
    do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
    if(0 <= i) {
      *b = j;
      b = SA + --B[c1]; j = i; ++m;
      do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
    }
  }

  if(1 < m) {
    if(flags & (16 | 32)) {
      if(flags & 16) {
        if((D = SAIS_MYMALLOC(k * 2, sais_index_type)) == NULL) {
          if(flags & (1 | 4)) { SAIS_MYFREE(C, k, sais_index_type); }
          if(flags & 2) { SAIS_MYFREE(B, k, sais_index_type); }
          return -2;
        }
      } else {
        D = B - k * 2;
      }
      assert((j + 1) < n);
      ++B[chr(j + 1)];
      for(i = 0, j = 0; i < k; ++i) {
        j += C[i];
        if(B[i] != j) { assert(SA[B[i]] != 0); SA[B[i]] += n; }
        D[i] = D[i + k] = 0;
      }
      LMSsort2(T, SA, C, B, D, n, k, cs);
      name = LMSpostproc2(SA, n, m);
      if(flags & 16) { SAIS_MYFREE(D, k * 2, sais_index_type); }
    } else {
      LMSsort1(T, SA, C, B, n, k, cs);
      name = LMSpostproc1(T, SA, n, m, cs);
    }
  } else if (m == 1) { /* only one S*-suffix => set immediately */
    *b = j + 1;        /* set entry in SA */
    if (level0) { LCP[b-SA] = -1; } /* mark first (=only) S*-suffix in bucket */
    name = 1;
  } else {
    name = 0;
  }

  /* stage 2: solve the reduced problem
     recurse if names are not yet unique */
  if(name < m) {
    if(flags & 4) { SAIS_MYFREE(C, k, sais_index_type); }
    if(flags & 2) { SAIS_MYFREE(B, k, sais_index_type); }
    newfs = (n + fs) - (m * 2);
    if((flags & (1 | 4 | 8)) == 0) {
      if((k + name) <= newfs) { newfs -= k; }
      else { flags |= 8; }
    }
    assert((n >> 1) <= (newfs + m));
    RA = SA + m + newfs;
    for(i = m + (n >> 1) - 1, j = m - 1; m <= i; --i) {
      if(SA[i] != 0) {
        RA[j--] = SA[i] - 1;
      }
    }
    
    if(sais_main(RA, SA, NULL, newfs, m, name, sizeof(sais_index_type), 0, 0) != 0) {
      if(flags & 1) { SAIS_MYFREE(C, k, sais_index_type); }
      return -2;
    }

    /* (re)compute starting positions of S*-suffixes (stored in RA): */
    i = n - 1; j = m - 1; c0 = chr(n - 1);
    do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
    for(; 0 <= i;) {
      do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
      if(0 <= i) {
        RA[j--] = i + 1;
        do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
      }
    }

    /* construct LCP for S*-suffixes: */
    /* PHI: "to whom I want to be compared" (pos. in T) */
    /* DELTA: "distance (in T) to next S*" (in PHI-order) */
    if (level0) {
      if (m < n/3) { /* hence we can store PHI and DELTA interleaved */
	PHI = LCP+m;       /* use space in LCP-array for PHI and DELTA */
	RA[m] = n;         /* stopper */
	j = SA[0];         /* j stores SA[i-1] in the following loop */
	PHI[j<<1] = n-1;   /* set PHI[SA[0]] to $ (causes mismatch in char. comp.) */
	PHI[(j<<1)+1] = 0; /* set DELTA */
	for (i = 1; i < m; ++i) {
	  q = SA[i];                /* text position */
	  p = q<<1;                 /* for interleaving */
	  PHI[p]=RA[j];             /* set PHI-array */
	  PHI[p+1]=RA[j+1]-RA[j];   /* set DELTA */
	  j = q;                    /* store for next loop iteration */
	}

	PLCP = PHI; /* overwrite DELTA in following loop */
	p = 0; /* guaranteed LCP-value */
	j = 0; /* position in PLCP and RA */
	for (i = 0; i < n; ++i) {
	  if (i == RA[j]) {
	    sais_index_type twoj = j << 1;
	    if (p < 0) p = 0;
	    while (chr(i+p) == chr(PHI[twoj]+p)) ++p;
	    t = PHI[twoj+1];      /* accesses DELTA-value */
	    q = RA[j+1]-RA[j];    /* length difference */
	    PLCP[twoj] = p;       /* overwrite PHI with PLCP */
	    ++j;
	    p -= (t > q) ? t : q; /* decrease p by larger of t and q */
	  }
	}

	/* translate PLCP-values to SA-order: */
	for (j = 0; j < m; ++j) LCP[j] = PLCP[SA[j]<<1];
      }
      else { /* non-interleaved */
	PHI = LCP;     /* use space in LCP-array for PHI */
	DELTA = LCP+m; /* because we compute only m < n/2 values, this is valid */
	RA[m] = n;     /* stopper */
	j = SA[0];     /* j stores SA[i-1] in the following loop */
	PHI[j] = n-1;  /* set PHI[SA[0]] to $ (causes mismatch in char. comp.) */
	DELTA[j] = 0;
	for (i = 1; i < m; ++i) {
	  q = SA[i];              /* text position */
	  PHI[q]=RA[j];           /* set PHI-array */
	  DELTA[q]=RA[j+1]-RA[j]; /* set DELTA */
	  j = q;                  /* store for next loop iteration */
	}

	PLCP = DELTA; /* overwrite DELTA in following loop */
	p = 0; /* guaranteed LCP-value */
	j = 0; /* position in PLCP and RA */
	for (i = 0; i < n; ++i) {
	  if (i == RA[j]) {
	    if (p < 0) p = 0;
	    while (chr(i+p) == chr(PHI[j]+p)) ++p;
	    t = PLCP[j];          /* accesses DELTA-value */
	    q = RA[j+1]-RA[j];    /* length difference */
	    PLCP[j++] = p;
	    p -= (t > q) ? t : q; /* decrease p by larger of t and q */
	  }
	}

	/* translate PLCP-values to SA-order: */
	for (j = 0; j < m; ++j) LCP[j] = PLCP[SA[j]];
      }
    }

    /* translate indices in RA to indices in T: */
    for(i = 0; i < m; ++i) SA[i] = RA[SA[i]];

    if(flags & 4) {
      if((C = B = SAIS_MYMALLOC(k, int)) == NULL) { return -2; }
    }
    if(flags & 2) {
      if((B = SAIS_MYMALLOC(k, int)) == NULL) {
        if(flags & 1) { SAIS_MYFREE(C, k, sais_index_type); }
        return -2;
      }
    }
  } /* endif (name < m) */
  else if (level0) { /* this should only occur for small or pathetic inputs */
    /* all names unique => computing LCP for *S naively takes linear time */
    /*printf("*** computing LCP naively..."); */
    j = SA[0]; /* j = SA[i-1] in the following loop */
    for (i = 1; i < m; ++i) {
      p = 0;
      while (chr(SA[i]+p) == chr(j+p)) p++;
      LCP[i] = p;
      j = SA[i];
    }
    /*printf("done.\n"); */
  }

  /* stage 3: induce the result for the original problem */
  if(flags & 8) { getCounts(T, C, n, k, cs); }
  /* put all S*-suffixes (and their LCP-values) into their buckets */
  if(1 < m) { /* otherwise SA (and LCP) is already correct */
    getBuckets(C, B, k, 1); /* find ends of buckets */
    i = m - 1, j = n, p = SA[m - 1], c1 = chr(p);
    if (level0) {
      newfs = LCP[m-1]; /* newfs stores LCP[i] in the following loop */
      do {
	q = B[c0 = c1];
	while(q < j) {
	  SA[--j] = 0; LCP[j] = -2; /* set remaining entries in old bucket to 0/-2 */
	}
	
	do { /* step through bucket c0 and write S*-suffixes to SA: */
	  SA[--j] = p; LCP[j] = newfs;
	  if(--i < 0) break;
	  newfs = LCP[i]; p = SA[i];
	} while((c1 = chr(p)) == c0);
	/*assert(LCP[j]==0); *//* first S*-suffix in bucket must have LCP-value 0 */
	LCP[j] = -1;       /* mark first S*-suffix in every bucket */
      } while(0 <= i);
      while(0 < j) {
	SA[--j] = 0; LCP[j] = -2; /* set remaining entries in smallest buckets to 0/-2 */
      }
    }
    else {
      do {
	q = B[c0 = c1];
	while(q < j) SA[--j] = 0; /* set remaining entries in old bucket to 0 */
	do { /* step through bucket c0 */
	  SA[--j] = p;
	  if(--i < 0) break;
	  p = SA[i];
	} while((c1 = chr(p)) == c0);
      } while(0 <= i);
      while(0 < j) SA[--j] = 0; /* set remaining entries in 1st bucket to 0 */
    }
  }

  if(isbwt == 0) {
    if (level0) induceSAandLCP(T, SA, LCP, C, B, n, k, cs);
    else induceSA(T, SA, C, B, n, k, cs);
  }
  else { pidx = computeBWT(T, SA, C, B, n, k, cs); }
  if(flags & (1 | 4)) { SAIS_MYFREE(C, k, sais_index_type); }
  if(flags & 2) { SAIS_MYFREE(B, k, sais_index_type); }

  return pidx;
}

/*---------------------------------------------------------------------------*/

int
sais(const unsigned char *T, int *SA, int* LCP, int n) {
  if((T == NULL) || (SA == NULL) || (LCP == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; LCP[0] = 0; } return 0; }
  return sais_main(T, SA, LCP, 0, n, UCHAR_SIZE, sizeof(unsigned char), 0,1);
}

int
sais_int(const int *T, int *SA, int n, int k) {
  if((T == NULL) || (SA == NULL) || (n < 0) || (k <= 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; } return 0; }
  return sais_main(T, SA, NULL, 0, n, k, sizeof(int), 0, 1);
}

int
sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n) {
  int i, pidx;
  if((T == NULL) || (U == NULL) || (A == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }
  pidx = sais_main(T, A, NULL, 0, n, UCHAR_SIZE, sizeof(unsigned char), 1,1);
  if(pidx < 0) { return pidx; }
  U[0] = T[n - 1];
  for(i = 0; i < pidx; ++i) { U[i + 1] = (unsigned char)A[i]; }
  for(i += 1; i < n; ++i) { U[i] = (unsigned char)A[i]; }
  pidx += 1;
  return pidx;
}

int
sais_int_bwt(const int *T, int *U, int *A, int n, int k) {
  int i, pidx;
  if((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (k <= 0)) { return -1; }
  if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }
  pidx = sais_main(T, A, NULL, 0, n, k, sizeof(int), 1,1);
  if(pidx < 0) { return pidx; }
  U[0] = T[n - 1];
  for(i = 0; i < pidx; ++i) { U[i + 1] = A[i]; }
  for(i += 1; i < n; ++i) { U[i] = A[i]; }
  pidx += 1;
  return pidx;
}
