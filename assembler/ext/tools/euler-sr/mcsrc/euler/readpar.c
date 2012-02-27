/***************************************************************************
 * Title:          readpar.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>
void set_score();

void readpar()
{
	int	i;

	LOW_COV = 15;
	LOW_COV = 10;
	LOW_COV_PATH = 2;
	MID_COV = 2;
	MIN_END = 100;
	MIN_OVERLAP = 20;
	END_LEG = 100;
	INT_LEG = 500;
	MIN_INT = 20;
	SHORT_D = 20;
	idum = 10;
	start_ent = 100;
	end_ent = 450;
	band = 15;
	max_seq = 75000;
	MIN_IDENTITY = 0.90;
	MIN_PERC = 0.6;
	MIS_SCORE = 4;
	MAX_BRA = 200;
	MAX_NODES = 40000;
	MAX_EDGE = 60000;
	MAX_DIF = 50;
	MAX_TMP_LEG = 5000000;
	MIN_OVERLAP2 = 50;
	MIN_END2 = 500;
	SHORTEST_OVERLAP = 10;
	SMALL_EDGE = 500;
	/*SMALL_EDGE = 20;*/
	LPAT = 1;
	g = 1;
	h = 4;
	overlaplen = 25;
	EndLength = 100;
	BulgeLength = 400;
	BulgeCoverage = 30;
	WhirlLength = 100;
	ChimericTerm = 200;
	ChimericCoverage = 1;
	SecondChimericCoverage = 2;
	LINK_COV = 1;
	LINK_MIN_LEN = 1000;
	VERTEX_SIZE = 20;
	ErosionLength = 60;

	FILTER_THRESH_PERC = 0.6;
	FILTER_THRESH_ABSOLUTE = 200;

	min_numcov = 10;
	word_len = 10;
	match_score = 6;
	n_ban = 1;
	for(i = 0; i < word_len; i ++)	n_ban *= 4;

	random1(&idum);
}

void set_score()
{
	int	i, j;

	for(i = 0; i < 15; i ++)	{
		for(j = 0; j < 15; j ++)	{
			if(W[i][j] < 0)		W[i][j] = -match_score;
		}
	}
}
