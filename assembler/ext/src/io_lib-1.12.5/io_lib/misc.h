#ifndef _misc_h
#define _misc_h

#include "io_lib/os.h"

#include <stdio.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */
#include <sys/types.h>

#include "io_lib/xalloc.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int is_directory(char * fn);
extern int is_file(char * fn);
extern int file_exists(char * fn);
extern int compressed_file_exists(char *fname);
extern int file_size(char * fn);
extern FILE *open_fofn(char *files);
extern char *read_fofn(FILE *fp);
extern void close_fofn(FILE *fp);
extern int fstrlen(char *f, int max_f);
extern void f2cstr(char *f, int max_f, char *c, int max_c);
extern void c2fstr(char *c, int max_c, char *f, int max_f);
extern char *mystrtok(char *s, char *ct);
extern char *myfind(char *file, char* searchpath, int (*found) (char *) );
extern void crash (char* format,...);
extern void str_tolower (char *s);
extern void str_toupper (char *s);
extern char *fn_tail (char *s);
extern void fn_tolower (char *s);
extern void fn_toupper (char *s);
extern void shell_call(char *command, char *output, int len);
extern char *date_str(void);
#ifdef NOSTRDUP
extern char *strdup(const char *s);
#endif
#ifdef NOSTRSTR
extern char *strstr(char *cs, char *ct);
#endif

#ifdef NOMEMMOVE
#define memmove(d,s,l) bcopy(s,d,l)
#endif
extern int myusleep(unsigned int useconds);

extern void errout(char *fmt, ...);
extern void messout(char *fmt, ...);

/*
 * Useful macros
 */
#define findfile(F,S) myfind((F),(S),file_exists)
/*is_file fails for symbolic links*/
/*#define findfile(F,S) myfind((F),(S),is_file)*/

#if defined(min)
#undef min
#undef max
#endif

#define min(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define max(A,B) ( ( (A) > (B) ) ? (A) : (B) )
#define sgn(A) ( (A) ? ( ( (A) < 0 ) ? -1 : 1 ) : 0 )

#ifdef MIN
#undef MIN
#endif
#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#ifdef MAX
#undef MAX
#endif
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )
#define SGN(A) ( (A) ? ( ( (A) < 0 ) ? -1 : 1 ) : 0 )
#define ABS(A) ( (A) < 0 ? -(A) : (A) )

/* Number of elements in array */
#define Number(A) ( sizeof(A) / sizeof((A)[0]) )

/*
 * Things taken from the new gap text_output.h. They'll be used globally
 * across all the programs in the end.
 */

/*
 * Usage: verror(priority, format, args...);
 * NB: don't pass more than 8K per call
 */
#define ERR_WARN 0
#define ERR_FATAL 1
void verror(int priority, char *name, char *fmt, ...);

/*
 * Usage: vmessage(format, args...);
 * NB: don't pass more than 8K per call
 */
void vmessage(char *fmt, ...);

/*
 * Adds a new header to the text output window.
 */
void vfuncheader(char *fmt, ...);

/*
 * As vfuncheader, but only outputting when necessary.
 */
void vfuncgroup(int group, char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /*_misc_h*/
