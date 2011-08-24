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

/* 
 * File: 	Read.c
 * Purpose:	Performs read/write IO on the Read data stucture.
 * Last update: 01/09/94
 */


/*
    The Read data type is designed so that it can hold a varying degree
    of information about sequences, yet have a single set of calls
    to access the data.

    There are plenty of assumptions around that both the number of
    bases and the number of points will fit into an int_2, a short.

*/

/* ---- Includes ---- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h> /* Only need on windows for _O_BINARY */
#include <unistd.h>

#include "io_lib/Read.h"
#include "io_lib/mFILE.h"

#ifdef IOLIB_ABI
# include "io_lib/abi.h"
#endif
#ifdef IOLIB_SCF
# include "io_lib/scf.h"
#endif
#ifdef IOLIB_ALF
# include "io_lib/alf.h"
#endif
#ifdef IOLIB_PLN
# include "io_lib/plain.h"
#endif
#ifdef IOLIB_ZTR
# include "io_lib/ztr.h"
#endif
#ifdef IOLIB_CTF
# include "io_lib/seqIOCTF.h"
#endif
#ifdef IOLIB_SFF
# include "io_lib/sff.h"
#endif
#ifdef IOLIB_EXP
# include "io_lib/expFileIO.h"
#endif
#ifdef USE_BIOLIMS
# include "spBiolims.h"
#endif

#include "io_lib/xalloc.h"
#include "io_lib/translate.h"
#include "io_lib/traceType.h"
#include "io_lib/misc.h"
#include "io_lib/open_trace_file.h"

/*
 * Read a sequence from a file "fnin" of format "format". If "format" is 0
 * (ANY_FORMAT), we automatically determine the correct format.
 * Returns:
 *   Read *   for success
 *   NULLRead for failure
 */
Read *read_reading(char *fn, int format) {
    Read *read;
    mFILE *fp;

#ifdef USE_BIOLIMS
    if( !strncmp(fn,BIOLIMS_TAG,strlen(BIOLIMS_TAG))){
	return spReadBiolimsReading(fn);
   }
#endif

    /*
     * If we're asking for an Experiment file, read it.
     * If the format is ANY then attempt EXP first following by trace.
     * Otherwise use the trace search mechanism.
     *
     * Note this is purely for locating files and not for forcing the file
     * format. It's here so that experiment files and trace files may be
     * given identical names but accessed through different search paths
     * (as is the case with the trace server).
     */
    if (format == TT_EXP) {
	if (NULL == (fp = open_exp_mfile(fn, NULL))) {
	    errout("'%s': couldn't open\n", fn);
	    return NULL;
	}
    } else {
	fp = NULL;
	if (format == TT_ANY)
	    fp = open_exp_mfile(fn, NULL);

	if (!fp && NULL == (fp = open_trace_mfile(fn, NULL))) {
	    errout("'%s': couldn't open\n", fn);
	    return NULL;
	}
    }

    read = mfread_reading(fp, fn, format);
    mfclose(fp);

    return read;
}

/*
 * Read a sequence from a FILE *fp of format "format". If "format" is 0
 * (ANY_FORMAT), we automatically determine the correct format.
 * We still pass a filename 'fn' although this isn't used other than for
 * filling in the read->trace_name field.
 *
 * NB this function should NOT be used when Biolims support is required
 * (as biolims readings are not stored in a file)
 *
 * Returns:
 *   Read *   for success
 *   NULLRead for failure
 */
Read *mfread_reading(mFILE *fp, char *fn, int format) {
    Read *read;
    mFILE *newfp;

    if (!fn)
	fn = "(unknown)";

    newfp = freopen_compressed(fp, NULL);
    if (newfp != fp) {
	fp = newfp;
    } else {
	newfp = NULL;
    }

#ifdef _WIN32
    /*
     * jkb 16/05/00 comment below
     *
     * On windows "prog < file.abi" will work wrongly (compared to
     * "prog file.abi") because windows is rather stupid. It treats ascii
     * and binary streams differently, it considers stdin to be ascii unless
     * told otherwise, and it can only be told otherwise by using non-ansi
     * windows-specific function calls.
     */
    if (format != TT_EXP && format != TT_PLN && fp->fp)
	_setmode(_fileno(fp->fp), _O_BINARY);
#endif

    if (format == TT_ANY || format == TT_ANYTR) {
	format = fdetermine_trace_type(fp);
	mrewind(fp);
    }

    switch (format) {
    case TT_UNK:
    case TT_ERR:
	errout("File '%s' has unknown trace type\n", fn);
	read = NULLRead;
	break;

#ifdef IOLIB_SCF
    case TT_SCF: {
        Scf *scf;
	scf = mfread_scf(fp);

	if (scf) {
	    read = scf2read(scf);
	    scf_deallocate(scf);
	} else
	    read = NULLRead;

	break;
    }
#endif

#ifdef IOLIB_CTF
    case TT_CTF:
	read = mfread_ctf(fp);
	break;
#endif

#ifdef IOLIB_SFF
    case TT_SFF:
	read = mfread_sff(fp);
	break;
#endif

#ifdef IOLIB_ZTR
    case TT_ZTR:
    case TT_ZTR1:
    case TT_ZTR2:
    case TT_ZTR3: {
        ztr_t *ztr;

	if ((ztr = mfread_ztr(fp))) {
	    uncompress_ztr(ztr);
	    read = ztr2read(ztr);
	    delete_ztr(ztr);
	} else {
	    read = NULLRead;
	}
	break;
    }
#endif

#ifdef IOLIB_ABI
    case TT_ABI:
	read = mfread_abi(fp);
	break;
#endif

#ifdef IOLIB_ALF
    case TT_ALF:
	read = mfread_alf(fp);
	break;
#endif

#ifdef IOLIB_EXP
    case TT_EXP: {
	/* FIXME: we shouldn't redirect like this */
	Exp_info *e = exp_mfread_info(fp);
	
	read = e ? exp2read(e,fn) : NULLRead;
	break;
    }
#endif

#ifdef IOLIB_PLN
    case TT_PLN:
	read = mfread_pln(fp);
	break;
#endif

    default:
	errout("Unknown format %d specified to read_reading()\n", format);
	read = NULLRead;
    }

    if (read != NULLRead && (read->trace_name = (char *)xmalloc(strlen(fn)+1)))
	strcpy(read->trace_name, fn);

    if (newfp) mfclose(newfp);

    return read;
}

Read *fread_reading(FILE *fp, char *fn, int format) {
    return mfread_reading(mfreopen(fn, "rb", fp), fn, format);
}

/*
 * Write a sequence to a FILE *fp of format "format". If "format" is 0,
 * we choose our favourite - SCF.
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int mfwrite_reading(mFILE *fp, Read *read, int format) {
    int r = -1;
    int no_compress = 0;

#ifdef _WIN32
    /*
     * jkb 09/06/00 comment below
     *
     * On windows "prog > file.scf" will work wrongly (compared to
     * "prog file.scf") because windows is rather stupid. It treats ascii
     * and binary streams differently, it considers stdout to be ascii unless
     * told otherwise, and it can only be told otherwise by using non-ansi
     * windows-specific function calls.
     */
    if (format != TT_EXP && format != TT_PLN && fp->fp)
	_setmode(_fileno(fp->fp), _O_BINARY);
#endif

    switch (format) {
    default:
	/* Defaults to ZTR type */

#ifdef IOLIB_ZTR
    case TT_ZTR:
    case TT_ZTR2: {
        ztr_t *ztr;
	ztr = read2ztr(read);
	compress_ztr(ztr, 2);
	r = mfwrite_ztr(fp, ztr); 
	delete_ztr(ztr);
	no_compress = 1;
	break;
    }
    case TT_ZTR1: {
        ztr_t *ztr;
	ztr = read2ztr(read);
	compress_ztr(ztr, 1);
	r = mfwrite_ztr(fp, ztr); 
	delete_ztr(ztr);
	break;
    }
    case TT_ZTR3: {
        ztr_t *ztr;
	ztr = read2ztr(read);
	compress_ztr(ztr, 3);
	r = mfwrite_ztr(fp, ztr); 
	delete_ztr(ztr);
	no_compress = 1;
	break;
    }
#endif

#ifdef IOLIB_SCF
    case TT_SCF: {
        Scf *scf;
	scf = read2scf(read);
	r = mfwrite_scf(scf, fp);
	scf_deallocate(scf);
	break;
    }
#endif

#ifdef IOLIB_CTF
    case TT_CTF:
	r = mfwrite_ctf(fp, read); 
	break;
#endif

#ifdef IOLIB_ABI
    case TT_ABI:
	/*return mfwrite_abi(fp, read); */
	break;
#endif

#ifdef IOLIB_SFF
    case TT_SFF:
	/*return mfwrite_sff(fp, read); */
	break;
#endif

#ifdef IOLIB_ALF
    case TT_ALF:
	/* return mfwrite_alf(fp, read); */
	break;
#endif

#ifdef IOLIB_EXP
    case TT_EXP: {
	Exp_info *e = read2exp(read, read->ident ? read->ident : "unknown");
	
	if (NULL == e) {
	    fprintf(stderr, "Failed to create experiment file.\n");
	    r = -1;
	} else {
	    exp_print_mfile(fp, e);
	    exp_destroy_info(e);
	    r = 0;
	}
	break;
    }
#endif

#ifdef IOLIB_PLN
    case TT_PLN:
	r = mfwrite_pln(fp, read);
	break;
#endif
    }

    mftruncate(fp, -1);
    if (r == 0 && !no_compress) {
	fcompress_file(fp);
    }
    mfflush(fp);

    return r;
}

int fwrite_reading(FILE *fp, Read *read, int format) {
    int ret;
    mFILE *mf = mfreopen(NULL, "wbx", fp);
    if (mf) {
	ret = mfwrite_reading(mf, read, format);
	mfflush(mf);
	mf->fp = NULL; /* Don't want this closed here */
	mfclose(mf);
    } else {
	return -1;
    }

    return ret;
}

/*
 * Write a sequence to a file "fn" of format "format". If "format" is 0,
 * we choose our favourite - SCF.
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int write_reading(char *fn, Read *read, int format) {
    int ret;
    mFILE *fp = mfopen(fn, "wb");
    if (!fp)
	return -1;
    
    ret = mfwrite_reading(fp, read, format);
    mfclose(fp);
    return ret;
}

/*
 * Old style stub interfaces implemented simply as redirection through
 * fread_reading and frwrite_reading.
 */
#ifdef IOLIB_ABI
Read *fread_abi(FILE *fp) {
    return fread_reading(fp, NULL, TT_ABI);
}

int fwrite_abi(FILE *fp, Read *read) {
    return fwrite_reading(fp, read, TT_ABI);
}
#endif

#ifdef IOLIB_ALF
Read *fread_alf(FILE *fp) {
    return fread_reading(fp, NULL, TT_ALF);
}

int fwrite_alf(FILE *fp, Read *read) {
    return fwrite_reading(fp, read, TT_ALF);
}
#endif

#ifdef IOLIB_CTF
Read *fread_ctf(FILE *fp) {
    return fread_reading(fp, NULL, TT_CTF);
}

int fwrite_ctf(FILE *fp, Read *read) {
    return fwrite_reading(fp, read, TT_CTF);
}
#endif

#ifdef IOLIB_PLN
Read *fread_pln(FILE *fp) {
    return fread_reading(fp, NULL, TT_PLN);
}

int fwrite_pln(FILE *fp, Read *read) {
    return fwrite_reading(fp, read, TT_PLN);
}
#endif

#ifdef IOLIB_ZTR
ztr_t *fread_ztr(FILE *fp) {
    ztr_t *z;
    mFILE *mf;

    if (NULL == (mf = mfreopen(NULL, "rb", fp)))
	return NULL;

    z = mfread_ztr(mf);
    mfclose(mf);
    return z;
}

int fwrite_ztr(FILE *fp, ztr_t *z) {
    mFILE *mf;
    int r;

    if (NULL == (mf = mfreopen(NULL, "wbx", fp)))
	return -1;

    r = mfwrite_ztr(mf, z);
    mfflush(mf);
    mf->fp = NULL; /* Don't want this closed here */
    mfclose(mf);
    return r;
}
#endif

#ifdef IOLIB_SCF
Scf *fread_scf(FILE *fp) {
    Scf *s;
    mFILE *mf;

    if (NULL == (mf = mfreopen(NULL, "rb", fp)))
	return NULL;

    s = mfread_scf(mf);
    mf->fp = NULL; /* Don't want this closed here */
    mfclose(mf);
    return s;
}

int fwrite_scf(Scf *s, FILE *fp) {
    mFILE *mf;
    int r;

    if (NULL == (mf = mfreopen(NULL, "wbx", fp)))
	return -1;

    r = mfwrite_scf(s, mf);
    mfflush(mf);
    mf->fp = NULL; /* Don't want this closed here */
    mfclose(mf);
    return r;
}
#endif

#ifdef IOLIB_EXP
Exp_info *exp_fread_info(FILE *fp) {
    Exp_info *e;
    mFILE *mf;

    if (NULL == (mf = mfreopen(NULL, "rb", fp)))
	return NULL;

    e = exp_mfread_info(mf);
    mf->fp = NULL; /* Don't want this closed here */
    mfclose(mf);
    return e;
}

void exp_print_file(FILE *fp, Exp_info *e) {
    mFILE *mf;

    if (NULL == (mf = mfreopen(NULL, "wbx", fp)))
	return;

    exp_print_mfile(mf, e);
    mfflush(mf);
    mf->fp = NULL; /* Don't want this closed here */
    mfclose(mf);
}
#endif
