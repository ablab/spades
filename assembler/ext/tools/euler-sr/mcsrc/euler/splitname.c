/***************************************************************************
 * Title:          splitname.c
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

char splitname(char *src_name, int *plate, char *library);
char extractname(char *src_name);
int loc_pair(char *name, char ext, char **src_name, char *ext_name, int num_seq);

char splitname(char *src_name, int *plate, char *library)
{
	int	i, j, k, l;
	char	c;

	l = strlen(src_name);
	for(i = 0; i < l; i ++)	{
		if(src_name[i] == '_' || src_name[i] == '-')	{
			strncpy(library, src_name, i);
			library[i + 1] = '\0';
			sscanf(&src_name[i + 1], "%d", plate);
			break;
		}
	}
	for(i = 0; i < l; i ++)	{
		if(src_name[i] == '.')	{
			c = src_name[i + 1];
			src_name[i] = '\0';
			return(c);
		}
	}

/*
	printf("Strange read name %s\n", src_name);
*/
	return('s');
}

char extractname(char *src_name)
{
	int	i, j, k, l;
	char	c;

	l = strlen(src_name);
	for(i = l - 2; i >= 0; i --)	{
		if(src_name[i] == '.')	{
			c = src_name[i + 1];
			return(c);
		}
	}

/*
	printf("Strange read name %s\n", src_name);
*/
	return('s');
}

int loc_pair(char *name, char ext, char **src_name, char *ext_name, int num_seq)
{
	int	i, j, k, l;

	for(i = 0; i < num_seq; i ++)	{
		if(ext_name[i] == ext && !strcmp(name, src_name[i]))	{
			return(i);
		}
	}

	return(-1);
}
