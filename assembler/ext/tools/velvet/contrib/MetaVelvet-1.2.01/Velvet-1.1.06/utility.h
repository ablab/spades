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
#ifndef UTILITY_H_
#define UTILITY_H_

#ifdef __GNUC__
#define ATTRIBUTE(list)  __attribute__ (list)
#else
#define ATTRIBUTE(list)
#endif

#include <stdio.h>

// Wrappers for malloc(), calloc(), and realloc() that always succeed.
// These functions print an error message and exit on failure, rather than
// requiring the calling function to check for NULL.  The arguments contain
// the type itself -- mallocOrExit(n, Foo) rather than malloc(n * sizeof Foo)
// -- to enable type checking and so that it can be shown in error messages.
#define mallocOrExit(count, type) \
               ((type *) mallocOrExit3((count), sizeof(type), #type))
#define callocOrExit(count, type) \
               ((type *) callocOrExit3((count), sizeof(type), #type))
#define reallocOrExit(ptr, count, type) \
               ((type *) reallocOrExit4((ptr), (count), sizeof(type), #type))

// However there are types for which just appending a '*' produces the
// wrong type or a syntax error, rather than a pointer-to-<type>.  These
// less type-safe wrappers are provided for use in these unusual cases.
#define mallocOrExitWithoutCast(count, type) \
               (mallocOrExit3((count), sizeof(type), #type))
#define callocOrExitWithoutCast(count, type) \
               (callocOrExit3((count), sizeof(type), #type))
#define reallocOrExitWithoutCast(ptr, count, type) \
               (reallocOrExit4((ptr), (count), sizeof(type), #type))

// (Implementation functions -- use the macro wrappers above.)
void *mallocOrExit3(size_t count, size_t size, const char *name);
void *callocOrExit3(size_t count, size_t size, const char *name);
void *reallocOrExit4(void *ptr, size_t count, size_t size, const char *name);

// Sets the program name to be prepended to error messages.
void setProgramName(const char *name);

// Prints an error message to standard error (with printf-style formatting
// and optionally appending a perror-style description of errno), and calls
// exit() with the specified exit status.
void exitErrorf(int exitStatus, boolean showErrno, const char *format, ...)
       ATTRIBUTE((format(printf, 3, 4), noreturn));

// Velvet-specific logging utility
void velvetLog(const char *format, ...)
	ATTRIBUTE((format(printf, 1, 2)));

// fprintf wrapper which exits upon failure (e.g. disk full)
void velvetFprintf(FILE * file, const char * format, ...)
	ATTRIBUTE((format(printf, 2, 3)));

// String Buffer

typedef struct
{
	char *str;
	size_t length;
	size_t allocated;
}
StringBuffer;

StringBuffer *newStringBuffer(size_t size);
void destroyStringBuffer(StringBuffer *buffer, boolean freeString);
void appendStringBuffer(StringBuffer *buffer, char *str);
void resetStringBuffer(StringBuffer *buffer);

// Solaris 10 and earlier versions lack titmersub()
#if defined (__SVR4) && defined (__sun) && !defined(timersub)
#define timersub(a, b, res) \
       *res.tv_sec = *a.tv_sec - *b.tv_sec; \
       *res.tv_usec = 0;
#endif

#endif
