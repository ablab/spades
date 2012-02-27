/***************************************************************************
 * Title:          ChopReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("usage: chopReads infile readlen outfile\n");
		exit(0);
	}
	FILE* in, *out;
	in = fopen(argv[1], "r");
	ssize_t newLen = atoi(argv[2]);
	out = fopen(argv[3], "w");
	//UNUSED+// char read2[200],title2[200];
	char title[200], read1[200]  ;
	//UNUSED// char c;
	//UNUSED// ssize_t r2len;
	while(in) {
		if (fgets(title, 199, in) == NULL) break;
		if (fgets(read1, 199, in) == NULL) break;
		/*		if (fgets(title2, 199, in) == NULL) break;
			if (fgets(read2, 199, in) == NULL) break;*/
		
		fprintf(out, "%s", title);
		read1[newLen] = '\0';
		fprintf(out, "%s\n", read1);
		/*		printf("%s", title2);
		r2len = strlen(read2);
		read2[r2len-1] = '\0';
		r2len--;
		printf("%s\n", read2+r2len-newLen);
		*/
	}
}
