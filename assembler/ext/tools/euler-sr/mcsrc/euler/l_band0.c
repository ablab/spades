
/* A PACKAGE FOR LOCALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call LOCAL_ALIGN(A,B,M,N,L,U,W,G,H,S,&SI,&SJ,&EI,&EJ,MW).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	*SI : starting position of sequence A in the optimal local alignment
	*SJ : starting position of sequence B in the optimal local alignment
	*EI : ending position of sequence A in the optimal local alignment
	*EJ : ending position of sequence B in the optimal local alignment
	MW  : maximum window size
*/

#include <stdinc.h>
#include <extfunc.h>

#define MININT -999999999

int LOCAL_ALIGN0(char *A, char *B, int M, int N, int low, int up, int W[][15],
		int G, int H, int *psi, int *psj, int *pei, int *pej, int MW);

int LOCAL_ALIGN0(char *A, char *B, int M, int N, int low, int up, int W[][15],
		int G, int H, int *psi, int *psj, int *pei, int *pej, int MW)
{ 
  int band;
  int i, j, si, ei;
  int c, d, e, t, m;
  int leftd, rightd;
  int best_score, starti, startj, endi, endj;
  int *wa, curd;
  int ib;
  char flag;
  int *CC, *DD;
  
  m = G+H;
  low = max(-M, low);
  up = min(N, up);
  
  if (N <= 0) { 
    *psi = *psj = *pei = *pej;
    return 0;
  }
  if (M <= 0) {
    *psi = *psj = *pei = *pej;
    return 0;
  }
  band = up-low+1;
  if (band < 1) {
    printf("low %d up %d\n", low, up);
    printf("low > up is unacceptable!\n");
    exit(1);
  }
  j = (MW+1+2) * sizeof(int);
  CC = (int *) ckalloc(j);
  DD = (int *) ckalloc(j);
  
  if (low > 0) leftd = 1;
  else if (up < 0) leftd = band;
  else leftd = 1-low;
  rightd = band;
  si = max(0,-up);
  ei = min(M,N-low);
  CC[leftd] = 0;
  for (j = leftd+1; j <= rightd; j++) {
    CC[j] = 0;
    DD[j] = -G;
  }
  CC[rightd+1] = MININT;
  DD[rightd+1] = MININT;
/*
  best_score = 0;
*/
/*	By Haixu, best_score are set up according to at least k matches,
	where k = length of the k-mer	*/
  best_score = 0;

  endi = si;
  endj = si+low;
  CC[leftd-1] = MININT;
  DD[leftd] = -G;
  for (i = si+1; i <= ei; i++) {
    if (i > N-up) rightd--;
    if (leftd > 1) leftd--;
    wa = W[A[i]];
    if ((c = CC[leftd+1]-m) > (d = DD[leftd+1]-H)) d = c;
    if ((ib = leftd+low-1+i ) > 0) c = CC[leftd]+wa[B[ib]];
/*
    if (ib > N) fprintf(stderr,"B[%d] out of range %d\n",ib,N);
*/
    if (d > c) c = d;
    if (c < 0) c = 0;
    e = c-G;
    DD[leftd] = d;
    CC[leftd] = c;
    if (c > best_score) {
      best_score = c;
      endi = i;
      endj = ib;
    }
    for (curd=leftd+1; curd <= rightd; curd++) {
      if ((c = c-m) > (e = e-H)) e = c;
      if ((c = CC[curd+1]-m) > (d = DD[curd+1]-H)) d = c;
/*
      if ((ib=curd+low-1+i) <= 0 || ib > N)
	fprintf(stderr,"B[%d]:%d\n",ib,B[ib]);
*/
      c = CC[curd] + wa[B[curd+low-1+i]];
      if (e > c) c = e;
      if (d > c) c = d;
      if (c < 0) c = 0;
      CC[curd] = c;
      DD[curd] = d;
      if (c > best_score) {
	best_score = c;
	endi = i;
	endj = curd+low-1+i;
      }
    }
  }
  
  leftd = max(1,-endi-low+1);
  rightd = band-(up-(endj-endi));
  CC[rightd] = 0;
  t = -G;
  for (j = rightd-1; j >= leftd; j--) {
    CC[j] = t = t-H;
    DD[j] = t-G;
  }
  for (j = rightd+1; j <= band; ++j) CC[j] = MININT;
  CC[leftd-1] = DD[leftd-1] = MININT;
  DD[rightd] = -G;
  flag = 0;
  for (i = endi; i >= 1; i--) {
    if (i+low <= 0) leftd++;
    if (rightd < band) rightd++;
    wa = W[A[i]];
    if ((c = CC[rightd-1]-m) > (d = DD[rightd-1]-H)) d = c;
    if ((ib = rightd+low-1+i) <= N) c = CC[rightd]+wa[B[ib]];

/*
    if (ib <= 0) fprintf(stderr,"rB[%d] <1\n",ib);
*/
    if (d > c) c = d;
    e = c-G;
    DD[rightd] = d;
    CC[rightd] = c;
    if (c == best_score) {
      starti = i;
      startj = ib;
      flag = 1;
      break;
    }
    for (curd=rightd-1; curd >= leftd; curd--) {
      if ((c = c-m) > (e = e-H)) e = c;
      if ((c = CC[curd-1]-m) > (d = DD[curd-1]-H)) d = c;

/*
      if ((ib=curd+low-1+i) <= 0 || ib > N)
	fprintf(stderr,"i: %d, B[%d]:%d\n",i,ib,B[ib]);
*/
      c = CC[curd] + wa[B[curd+low-1+i]];
      if (e > c) c = e;
      if (d > c) c = d;
      CC[curd] = c;
      DD[curd] = d;
      if (c == best_score) {
	starti = i;
	startj = curd+low-1+i;
	flag = 1;
	break;
      }
    }
    if (flag == 1) break;
  }
  
  if (starti < 0 || starti > M || startj < 0 || startj > N) {
    printf("starti=%d, startj=%d\n",starti,startj);
    *psi = *psj = *pei = *pej;
    exit(1);
  }
  *psi = starti;
  *psj = startj;
  *pei = endi;
  *pej = endj;

  free((void *) CC);
  free((void *) DD);

  return best_score;
}
