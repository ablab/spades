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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 19/3/99 johnt - added Corba support */
#ifdef USE_CORBA
#include "stcorba.h"
#endif

#ifdef USE_BIOLIMS
#include "spBiolims.h"
#endif

char *myfind(char *file, char* searchpath, int (*found) (char *) )
{
    static char wholePath[1024];
    char *path;
    char *f;

    f = NULL;
    if (found(file)) {
	strcpy(wholePath,file);
	f = wholePath;
    } else if (searchpath != NULL) {
	char *paths;
	char *next;

	paths = (char *) malloc(strlen(searchpath)+1);
	strcpy(paths,searchpath);

	path = paths;
	next = strchr(path,':');
        while( next && (*(next+1) == ':' )){
    	   /* 26/03/99 johnt - allow : to be entered into path by using :: */
 	   memmove(next,next+1,strlen(next+1)+1); /* shuffle up data [including \0]*/
	   next = strchr(next+1,':');
	}
        if(next)
	    *next = '\0';

	while (path!= NULL) {

#ifdef USE_CORBA	
	  /* 19/03/99 johnt - if it is a corba path - look there */
          if( !strncmp( CORBATAG,path,strlen(CORBATAG))){
	    if(corba_found(wholePath,path+strlen(CORBATAG),file)){
	      f = wholePath;
	      break;
	    }
	  }
	  else
#endif
#ifdef USE_BIOLIMS
	  if( !strncmp( BIOLIMS_TAG,path,strlen(BIOLIMS_TAG))){
	    if(biolims_found(wholePath,path+strlen(BIOLIMS_TAG),file)){
	      f = wholePath;
	      break;
	    }
	  }
	  else
#endif    
	  {
	    (void) strcpy(wholePath,path);
	    (void) strcat(wholePath,"/");
	    (void) strcat(wholePath,file);
	    if (found(wholePath)) {
	      f = wholePath;
	      break;
	    }
	  }
	  path = next;
	  if( path ){
	      path++;
	      next = strchr(path,':');
	      while( next && (*(next+1) == ':' )){
    		 /* 26/03/99 johnt - allow : to be entered into path by using :: */
 		 memmove(next,next+1,strlen(next+1)+1); /* shuffle up data */
		 next = strchr(next+1,':');
	      }
	      if(next)
		*next='\0';
	  }
	}
	free(paths);
    }

    return f;
}
