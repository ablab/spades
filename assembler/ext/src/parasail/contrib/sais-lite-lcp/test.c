/*
 * test.c for sais-lite
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sais.h"


static
int
cmp_suf(const unsigned char *T, int n, int p1, int p2) {
  int r, s = (p1 < p2) ? 1 : ((p1 > p2) ? -1 : 0);
  for(r = 0; (p1 < n) && (p2 < n) && ((r = T[p1] - T[p2]) == 0); ++p1, ++p2) { }
  return (r != 0) ? r : s;
}

int
main(int argc, const char *argv[]) {
  unsigned char *T1;
  int *T3;
  unsigned char *T1BWT;
  int *T3BWT;
  int *SA1;
  int *LCP1;
  int *SA3;
  int *A;
  int i, j, n, p1, p3;
  unsigned int bits;

  fprintf(stdout, "start test...\n");
  for(n = 1; n <= 24; ++n) {
    T1 = malloc(n * sizeof(unsigned char));
    T1BWT = malloc(n * sizeof(unsigned char));
    T3 = malloc(n * sizeof(int));
    T3BWT = malloc(n * sizeof(int));
    SA1 = malloc((n+1) * sizeof(int));
    LCP1 = malloc(n * sizeof(int));
    SA3 = malloc(n * sizeof(int));
    A = malloc(n * sizeof(int));
    for(bits = 0; bits < (1U << n); ++bits) {
      if((bits & 4095) == 0) {
        fprintf(stderr, "  n=%2d : %3d%%\r", n, (int)((double)bits / (double)((1U << n) - 1) * 100.0));
      }
      for(i = 0; i < n; ++i) {
        T1[i] = (bits >> i) & 1;
        T3[i] = T1[i] * 511;
      }

      /* construct sa and bwt */
      if(sais(T1, SA1, LCP1, n) != 0) {
        fprintf(stderr, "  n=%2d, bits=%u : failure - sais\n", n, bits);
        exit(EXIT_FAILURE);
      }
      if((p1 = sais_bwt(T1, T1BWT, A, n)) < 0) {
        fprintf(stderr, "  n=%2d, bits=%u : failure - sais_bwt\n", n, bits);
        exit(EXIT_FAILURE);
      }
      if(sais_int(T3, SA3, n, 512) != 0) {
        fprintf(stderr, "  n=%2d, bits=%u : failure - sais_int\n", n, bits);
        exit(EXIT_FAILURE);
      }
      if((p3 = sais_int_bwt(T3, T3BWT, A, n, 512)) < 0) {
        fprintf(stderr, "  n=%2d, bits=%u : failure - sais_int_bwt\n", n, bits);
        exit(EXIT_FAILURE);
      }

      /* check SA1 */
      for(i = 1; i < n; ++i) {
        if(0 <= cmp_suf(T1, n, SA1[i - 1], SA1[i])) {
          fprintf(stderr, "  n=%2d, bits=%u : failure - SA1\n", n, bits);
          for(i = 0; i < n; ++i) {
            fprintf(stderr, "    SA[%d]=%d: ", i, SA1[i]);
            for(j = SA1[i]; j < n; ++j) { fprintf(stderr, "%d", T1[j]); }
            fprintf(stderr, "\n");
          }
          exit(EXIT_FAILURE);
        }
      }

      /* check SA3 */
      for(i = 0; i < n; ++i) {
        if(SA1[i] != SA3[i]) {
          fprintf(stderr, "  n=%2d, bits=%u : failure - SA3\n", n, bits);
          for(i = 0; i < n; ++i) {
            fprintf(stderr, "    SA1[%d]=%d, SA3[%d]=%d: ", i, SA1[i], i, SA3[i]);
            for(j = SA3[i]; j < n; ++j) { fprintf(stderr, "%d", T3[j] / 511); }
            fprintf(stderr, "\n");
          }
          exit(EXIT_FAILURE);
        }
      }

      /* check T1BWT */
      for(i = 0, j = 0; i <= n; ++i) {
        if(i != 0) {
          if(SA1[i - 1] == 0) { if(p1 != i) { break; } }
          else if(n <= j) { break; }
          else { if(T1BWT[j++] != T1[SA1[i - 1] - 1]) { break; } }
        } else {
          if(T1BWT[j++] != T1[n - 1]) { break; }
        }
      }
      if((i != (n + 1)) || (j != n)) {
        fprintf(stderr, "  n=%2d, bits=%u : failure - T1BWT\n", n, bits);
        fprintf(stderr, "    T1BWT=");
        for(i = 0; i < n; ++i) { fprintf(stderr, "%d", T1BWT[i]); }
        fprintf(stderr, ", p1=%d\n", p1);
        for(i = 0; i < n; ++i) {
          fprintf(stderr, "    SA[%d]=%d: ", i, SA1[i]);
          for(j = SA1[i]; j < n; ++j) { fprintf(stderr, "%d", T1[j]); }
          fprintf(stderr, " ");
          for(j = 0; j < SA1[i]; ++j) { fprintf(stderr, "%d", T1[j]); }
          fprintf(stderr, "\n");
        }
        exit(EXIT_FAILURE);
      }

      /* check T3BWT */
      for(i = 0; i < n; ++i) { if(T1BWT[i] != T3BWT[i] / 511) { break; } }
      if((i != n) || (p1 != p3)) {
        fprintf(stderr, "  n=%2d, bits=%u : failure - T3BWT\n", n, bits);
        fprintf(stderr, "    T3BWT=");
        for(i = 0; i < n; ++i) { fprintf(stderr, "%d", T3BWT[i]); }
        fprintf(stderr, ", p3=%d\n", p3);
        for(i = 0; i < n; ++i) {
          fprintf(stderr, "    SA[%d]=%d: ", i, SA3[i]);
          for(j = SA3[i]; j < n; ++j) { fprintf(stderr, "%d", T3[j] / 511); }
          fprintf(stderr, " ");
          for(j = 0; j < SA3[i]; ++j) { fprintf(stderr, "%d", T3[j] / 511); }
          fprintf(stderr, "\n");
        }
        exit(EXIT_FAILURE);
      }

    }
    fprintf(stderr, "  n=%2d : success\n", n);
    free(T1);
    free(T1BWT);
    free(T3);
    free(T3BWT);
    free(SA1);
    free(SA3);
    free(A);
  }
  fprintf(stderr, "finish test\n");

  return 0;
}
