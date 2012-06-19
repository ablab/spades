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

// Wrappers for malloc(), calloc(), and realloc() that always succeed.
// These functions print an error message and exit on failure, rather than
// requiring the calling function to check for NULL.  The arguments contain
// the type itself -- mallocOrExit(n, Foo) rather than malloc(n * sizeof Foo)
// -- mainly so that it can be shown in error messages.
#define mallocOrExit(count, type)  (mallocOrExit3((count), sizeof(type), #type))
#define callocOrExit(count, type)  (callocOrExit3((count), sizeof(type), #type))
#define reallocOrExit(ptr, count, type) \
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

#endif
