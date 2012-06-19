/*
Copyright 2009 John Marshall (jm18@sanger.ac.uk) 

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>

#include "globals.h"
#include "utility.h"

static void allocExitError(const char *function, unsigned long long count,
                          unsigned long long size, const char *name)
{
       if (size == 1)
               exitErrorf(EXIT_FAILURE, true,
                          "Can't %s %llu %ss",
                          function, count, name);
       else
               exitErrorf(EXIT_FAILURE, true,
                          "Can't %s %llu %ss totalling %llu bytes",
                          function, count, name, count * size);
}

void *mallocOrExit3(size_t count, size_t size, const char *name)
{
       void *p = malloc(count * size);
       if (p == NULL && count > 0)
               allocExitError("malloc", count, size, name);

       return p;
}

void *callocOrExit3(size_t count, size_t size, const char *name)
{
       void *p = calloc(count, size);
       if (p == NULL && count > 0)
               allocExitError("calloc", count, size, name);

       return p;
}

void *reallocOrExit4(void *ptr, size_t count, size_t size, const char *name)
{
       void *p = realloc(ptr, count * size);
       if (p == NULL && count > 0)
               allocExitError("realloc", count, size, name);

       return p;
}


static const char *programName = NULL;

void setProgramName(const char *name)
{
       programName = name;
}

void exitErrorf(int exitStatus, boolean showErrno, const char *format, ...)
{
       int savedErrno = errno;
       va_list args;
       va_start(args, format);
       if (programName)
               fprintf(stderr, "%s: ", programName);
       vfprintf(stderr, format, args);
       if (showErrno)
               fprintf(stderr, ": %s", strerror(savedErrno));
       fprintf(stderr, "\n");
       va_end(args);

       exit(exitStatus);
}
