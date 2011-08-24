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

#include "io_lib/misc.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>

int fstrlen(char *f, int max_f)
{
    for (; max_f > 0 && (isspace(f[max_f-1]) || f[max_f-1]=='\0'); max_f--);
    return max_f;
}





void f2cstr(char *f, int max_f, char *c, int max_c)
{
    int i;

    i = min(fstrlen(f,max_f),max_c);
    strncpy(c,f,i);
    c[i]='\0';
}


void c2fstr(char *c, int max_c, char *f, int max_f)
{
    int i;
    i = min((int)strlen(c),max_f);
    strncpy(f,c,i);
    for( ; i<max_f; i++) f[i]=' ';

}




char *mystrtok(char *s, char *ct)
/*
** When strtok isn't good enough
*/
{
    char *this;
    static char *look;
    static int last;

    if (s == NULL) {
	if (last) return NULL;
    } else {
	look = s;
	last = 0;
    }
    this = look;

    for ( ; *look && strchr(ct,*look)==NULL; look++ ) ;
    last = (! *look);
    *look++ = '\0';
    
    return this;
}


void str_tolower (char *s)
/*
** Convert string to lower case
*/
{
    if (!s) return;
    for ( ; *s ; s++ )
	if (isupper(*s))
	    *s = tolower(*s);
}

void str_toupper (char *s)
/*
** Convert string to upper case
*/
{
    if (!s) return;
    for ( ; *s ; s++ )
	if (islower(*s))
	    *s = toupper(*s);
}

#ifdef NOSTRSTR
/*
** My routines for nice sun ones.
*/
char *strstr(char *cs, char *ct)
/*
** ANSI C has the function strstr().
**
**     strstr() returns a pointer to the first  occurrence  of  the
**     pattern  string  s2  in  s1.   For example, if s1 is "string
**     thing" and s2 is "ing", strstr() returns "ing thing".  If s2
**     does not occur in s1, strstr() returns NULL.
**
** It's not always implemented. Here's my cludge:
*/
{
    int i;
    int len_ct;
    int end;
    len_ct = strlen(ct);
    end = strlen(cs) - len_ct;
    for (i=0;i<=end;i++)
      if (strncmp(&cs[i],ct,len_ct)==0)
	return &cs[i];

    return NULL;
}
#endif

#ifdef NOSTRDUP
char *strdup(const char *str)
/*
** SunOS has a nice strdup() function.
**
**     strdup() returns a pointer to a new string which is a dupli-
**     cate  of the string pointed to by s1.  The space for the new
**     string is obtained using malloc(3V).  If the new string can-
**     not be created, a NULL pointer is returned.
**
** Other ANSI C libraries don't have this. Here is my kludge:
*/
{
    char *newstr;
    int i = strlen(str);

    if ((newstr = (char *)malloc((unsigned int)(i+1))) == NULL)
        return NULL;

    for (; i>=0; i--)
        newstr[i] = str[i];

    return newstr;
}
#endif

#ifdef NOSTRCASECMP
int strcasecmp(const char *s1, const char *s2) {
    while (tolower(*s1) == tolower(*s2)) {
        /* If at the end of the string, then they're equal */
        if (0 == *s1)
	    return 0;
        s1++;
	s2++;
    }
  
    /* One ended before the other, so return 1 or -1 */
    return (*(unsigned char *)s1) < (*(unsigned char *)s2) ? -1 : 1;
}
#endif
