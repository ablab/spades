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

#include <sys/types.h>
#include <sys/stat.h>
/* Alliant's Concentrix <sys/stat.h> is hugely deficient */
/* Define things we require in this program              */
/* Methinks S_IFMT and S_IFDIR aren't defined in POSIX   */
#ifndef S_ISDIR
#define S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#endif /*!S_ISDIR*/
#ifndef S_ISREG
#define S_ISREG(m)      (((m)&S_IFMT) == S_IFREG)
#endif /*!S_ISREG*/

int is_directory(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISDIR(buf.st_mode);
}

int is_file(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISREG(buf.st_mode);
}

int file_exists(char * fn)
{
    struct stat buf;
    return ( stat(fn,&buf) == 0);
}

int compressed_file_exists(char *fname)
{
    struct stat buf;
    char fn[2048];

    if (stat(fname, &buf) == 0) return 1;

    sprintf(fn, "%s.gz", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.bz", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.bz2", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.Z", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.z", fname);
    if (stat(fn, &buf) == 0) return 1;

    return 0;
}

int file_size(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) != 0) return 0;
    return buf.st_size;
}

/*
 * ---------------------------------------------------------------------------
 * File of filename management
 */

FILE *open_fofn(char *files) {
  return fopen(files, "r");
}

char *read_fofn(FILE *fp) {
  char line[256];
  static char name[256];
  
  while (fgets(line, 254, fp)) {
    if (1 == sscanf(line, "%s", name))
      return name;
  }

  return NULL;
}

void close_fofn(FILE *fp) {
  fclose(fp);
}


/*
 * ---------------------------------------------------------------------------
 * Temporary file handling.
 */

#ifdef _WIN32
/*
 * On UNIX systems we use tmpfile().
 *
 * On windows this is broken because it always attempts to create files in
 * the root directory of the current drive, and fails if the user does not
 * have write permission.
 *
 * We can't wrap up mkstemp() either as that doesn't exist under windows.
 * Instead we roll our own tmpfile() using the native windows API.
 */
#include <windows.h>
#include <winbase.h>
#include <stdio.h>
#include <io.h>
#include <fcntl.h>
#define _POSIX_ /* needed to get PATH_MAX */
#include <limits.h>

static void display_win_error(char *msg) {
    LPVOID lpMsgBuf;
    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |
    	          FORMAT_MESSAGE_FROM_SYSTEM |
    	          FORMAT_MESSAGE_IGNORE_INSERTS,
    	          NULL,
    	          GetLastError(),
    	          0,                 /* Default language */
    	          (LPTSTR)&lpMsgBuf, /* Got to love void* to str casts! */
    	          0, NULL);
    fprintf(stderr, "%s: error #%d, %s", msg, GetLastError(), lpMsgBuf);
    LocalFree(lpMsgBuf);
}

/*
 * Creates a temporary file and returns a FILE pointer to it.
 * The file will be automatically deleted when it is closed or the
 * applicaton exits.
 *
 * Returns NULL on failure.
 */
FILE *tmpfile(void) {
    DWORD ret;
    char tmp_path[PATH_MAX], shrt_path[PATH_MAX];
    int fd;
    FILE *fp;

    /* The Windows Way: get the temp directory and a file within it */
    ret = GetTempPath(PATH_MAX, tmp_path);
    if (ret == 0 || ret > PATH_MAX) {
	display_win_error("GetTempPath()");
        return NULL;
    }

    if (0 == GetTempFileName(tmp_path, "fubar", 0, shrt_path)) {
	display_win_error("GetTempFileName()");
	return NULL;
    }

    /*
     * O_TRUNC incase anyone has managed to inject data in the newly created
     * file already via race-conditions.
     *
     * O_EXCL to (in theory) stop anyone else opening it and to die if someone
     * beat us to it - although this appears to not actually work on Windows.
     *
     * O_TEMPORARY so that the file is removed on close.
     */
    if (-1 == (fd = _open(shrt_path, O_RDWR | O_TRUNC | O_EXCL |
			  O_BINARY | O_TEMPORARY, 0600))) {
	perror(shrt_path);
    }

    /* Replace fd with FILE*. No need to close fd */
    if (NULL == (fp = _fdopen(fd, "r+b"))) {
	perror(shrt_path);
    }

    return fp;
}

#endif /* _WIN32 */
