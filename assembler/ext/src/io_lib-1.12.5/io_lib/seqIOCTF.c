/* 
 * Title:       seqIOCTF
 * 
 * File: 	seqIOCTF.c
 * Purpose:	Reading/writing of CTF sequences
 * Last update: March 2000
 *
 * Change log:
 * Created mieg, march 2000, importing code from wabi/ctftrace.c
 */


/* ---- Imports ---- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "io_lib/Read.h"
#include "io_lib/seqIOCTF.h"
#include "io_lib/abi.h"
#include "io_lib/mach-io.h"
#include "io_lib/xalloc.h"
#include "io_lib/stdio_hack.h"

/* ---- Constants ---- */

/*
 * Read the CTF format sequence from FILE *fp into a Read structure.

 * Returns:
 *   Read *	- Success, the Read structure read.
 *   NULLRead	- Failure.
 */
Read *fread_ctf (FILE *fp) {
    Read *read = ctfFRead (fp) ;

    return read ;
}

/*
 * Read the CTF format sequence from file 'fn' into a Read structure.
 */

Read *read_ctf (char *fn) {
    Read *read;
    FILE *fp;

    /* Open file */
    if ((fp = fopen(fn, "rb")) == NULL)
	return NULLRead;

    read = fread_ctf(fp);
    fclose(fp);

    if (read && (read->trace_name = (char *)xmalloc(strlen(fn)+1)))
	strcpy(read->trace_name, fn);

    return read;
}
    

/*
 * Write to an CTF file - unsupported.
 */
/* ARGSUSED */
int fwrite_ctf (FILE *fp, Read *read) {
  return ctfFWrite (fp, read) ;
}

/*
 * Write to an CTF file 
 */
/* ARGSUSED */
int write_ctf(char *fn, Read *read) {
  FILE *fp;
  
  /* Open file */
  if ((fp = fopen(fn, "wb")) == NULL)
    return -1 ;
  
  fwrite_ctf (fp, read) ;
  fclose(fp);

  return 0 ;
}

