/***************************************************************************
 * Title:          html.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extfunc.h>

extern char htmlout, ***content, caption[2000];

void print_text_line(FILE *fp, int length);
void print_emptyline(FILE *fp);
void print_chimtable(FILE *fp, READTABLE *RT);
void print_hl(FILE *fp);
void print_header(FILE *fp, char *ptitle);
void print_line(FILE *fp, char *line);
void print_tailor(FILE *fp);
void print_table(FILE *fp, int nrow, int ncol, char ***content, char *caption);
char ***allocate_content(int nrow, int ncol, int len);
char ***free_content(char ***content, int nrow, int ncol);

void print_emptyline(FILE *fp)
{
	fprintf(fp, "<BR>\n");
}


void print_hl(FILE *fp)
{
	fprintf(fp, "<hr>\n");
}

void print_header(FILE *fp, char *ptitle)
{
	fprintf(fp, "<html>\n");
	fprintf(fp, "<head>\n");
	fprintf(fp, "<title>%s</title>\n", ptitle);
	fprintf(fp,"<link REL=\"STYLESHEET\" HREF=\"euler.css\" TYPE=\"text/css\">\n");
	fprintf(fp, "<meta http-equiv=\"Content-Language\" content=\"en-us\">\n");
	fprintf(fp, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=windows-1252\">\n");
	fprintf(fp, "</head>\n");
	fprintf(fp, "<body>\n");
	fprintf(fp, "<H4>%s</H4>", ptitle);
}

void print_line(FILE *fp, char *line)
{
	fprintf(fp, "%s<BR>\n", line);
}

/* print the section header */
void print_section_head(FILE *fp, char *line)
{
	fprintf(fp, "<H5>%s</H5>\n", line);
}


void print_tailor(FILE *fp)
{
	fprintf(fp, "</body>\n");
	fprintf(fp, "</html>\n");
}

/* The first row is the header */

void print_table(FILE *fp, int nrow, int ncol, char ***content, char *caption)
{
	int	i, j, k;

	print_emptyline(fp);
	fprintf(fp, "<TABLE cellspacing = 1 cellpadding = 5 border = 1 bordercolor = \"#AAAAAA\">\n");
	fprintf(fp, "<CAPTION align = \"left\">%s</CAPTION>\n", caption);
	fprintf(fp, "<TR>\n");
	for(j = 0; j < ncol; j ++)	{
	  //fprintf(fp, "<TH>%s</TH>\n", content[0][j]);
		fprintf(fp, "<TD class = \"header\"><DIV class = \"tbl_hd_top\">%s</DIV></TD>\n", content[0][j]);
	}
	fprintf(fp, "</TR>\n");
	for(i = 1; i < nrow; i ++)	{
		fprintf(fp, "<TR>\n");
		for(j = 0; j < ncol; j ++)	{
		  //fprintf(fp, "<TD>%s</TD>\n", content[i][j]);
			fprintf(fp, "<TD class = \"data\"><DIV class = \"tbl_data\">%s</DIV></TD>\n", content[i][j]);
		}
		fprintf(fp, "</TR>\n");
	}
	fprintf(fp, "</TABLE>\n");
	print_emptyline(fp);
}

char ***allocate_content(int nrow, int ncol, int len)
{
	int	i, j, k, l;
	char	***content;

	content = (char ***) ckalloc(nrow * sizeof(char **));
	for(i = 0; i < nrow; i ++)	{
		content[i] = (char **) ckalloc(ncol * sizeof(char *));
		for(j = 0; j < ncol; j ++)	{
			content[i][j] = (char *) ckalloc(len * sizeof(char));
		}
	}
	return(content);
}

char ***free_content(char ***content, int nrow, int ncol)
{
	int	i, j, k, l;

	for(i = 0; i < nrow; i ++)	{
		for(j = 0; j < ncol; j ++)	{
			free((void **) content[i][j]);
		}
		free((void **) content[i]);
	}
	free((void ***) content);
	
	free((char ***) NULL);
}

void print_chimtable(FILE *fp, READTABLE *RT)
{
	int	i, j, k, n;
	char	temp[1000];

	if(RT -> num_chim == 0)	{
		return;
	}

/*	Count # of chimeric reads	*/
	n = 0;
	for(i = 0; i < RT -> num_chim; i ++)	{
		if(RT -> chim[i] < RT -> num_seq)	{
			n ++;
		}
	}
	k = 0;
	if(!htmlout)	{
		fprintf(fp, "# chimeric reads found: %d.\n", n);
		print_text_line(fp, LINE_LENGTH);
		fprintf(fp, "Reads     Index     Name\n");
		print_text_line(fp, LINE_LENGTH);
		for(i = 0; i < RT -> num_chim; i ++)	{
			if(RT -> chim[i] < RT -> num_seq)	{
				k ++;
				fprintf(fp, "%-10d%-10d%-20s\n", k, RT -> chim[i], RT -> src_name[RT -> chim[i]]);
			}
		}
	} else	{
		sprintf(temp, "# chimeric reads found: %d.", n);
		print_line(fp, temp);
		sprintf(caption, "List of chimeric reads");
		content = allocate_content(RT -> num_chim + 1, 3, 50);
		strcpy(content[0][0], "Reads");
		strcpy(content[0][1], "Index");
		strcpy(content[0][2], "Name");
		for(i = 0; i < RT -> num_chim; i ++)	{
			if(RT -> chim[i] < RT -> num_seq)	{
				k ++;
				sprintf(content[k][0], "%d", k);
				sprintf(content[k][1], "%d", RT -> chim[i]);
				sprintf(content[k][2], "%s", RT -> src_name[RT -> chim[i]]);
			}
		}
		print_table(fp, k + 1, 3, content, caption);
		content = free_content(content, RT -> num_chim + 1, 3);
	}
	fflush(fp);
}

void print_gaps(FILE *fp, int *dist)
{
	int	i, j, k;

	if(!htmlout)	{
		printf("------------------------------------------------------------------\n");
		printf("Distribution of # gaps in overlap.\n");
		printf("------------------------------------------------------------------\n");
		for(i = 2; i < 6; i ++)	{
			printf("%7d ", i);
		}
		printf("     >5\n");
		printf("------------------------------------------------------------------\n");
		for(i = 2; i < 7; i ++)	{
			printf("%7d ", dist[i]);
		}
		printf("\n");
		printf("------------------------------------------------------------------\n");
	} else	{
		sprintf(caption, "Distribution of # gaps in overlap");
		content = allocate_content(2, 6, 50);
		strcpy(content[0][0], "# gaps in alignment");
		strcpy(content[1][0], "# alignments");
		for(i = 2; i < 6; i ++)	{
			sprintf(content[0][i - 1], "%7d ", i);
			sprintf(content[1][i - 1], "%7d ", dist[i]);
		}
		sprintf(content[0][5], ">5");
		sprintf(content[1][5], "%d", dist[6]);
		print_table(fp, 2, 6, content, caption);
		content = free_content(content, 2, 6);
	}
}

void print_text_line(FILE *fp, int length)
{
	int	i;
	for(i = 0; i < length; i ++)	{
		fprintf(fp, "-");
	}
	fprintf(fp, "\n");
}

void print_text(FILE *fp, char *text)
{
	fprintf(fp, text);
}
