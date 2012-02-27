/***************************************************************************
 * Title:          extvab.h
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
extern int total_nuc;
extern char na_name[17];
extern int LOW_COV, LOW_COV_PATH, MID_COV, SHORT_D;
extern int MIN_END, END_LEG, INT_LEG;
extern int MIN_INT;
extern int idum;
extern int max_seq;
extern int MAX_TMP_LEG;
extern int MAX_BRA;
extern int g, h;
extern int match_score, mov_qual;
extern int W[15][15];
extern double MIN_PERC, MIN_IDENTITY;
extern int MIN_OVERLAP;
extern int start_ent, end_ent, band;
extern int SHORTEST_OVERLAP;
extern long n_ban;
extern int overlaplen, min_numcov;
extern int word_len;
extern int gap_k;
extern char END_MERGE;
extern int MAX_NODES, MAX_EDGE;
extern int MAX_DIF;
extern int SMALL_CYCLE;
extern int MIS_SCORE, LPAT;
extern int MIN_END2, MIN_OVERLAP2;
extern char librule[500][100];
extern char namerule[2][10];
extern char pairrule[500][2];
extern int ntyperule[2], ntypepair;
extern int pairrange[500][2];
extern int platerule[500][2];
extern int matetype;
extern int VERTEX_SIZE;
extern int EndLength;
extern int LINK_COV, LINK_MIN_LEN;
extern int BulgeLength, BulgeCoverage, WhirlLength,
  ChimericTerm, ChimericCoverage, ErosionLength,
  SecondChimericCoverage;
extern int SMALL_EDGE;
extern int FILTER_THRESH_ABSOLUTE;
extern double FILTER_THRESH_PERC;
void RemovePathFromVertex(NODES *vertex, int path, int pos);
