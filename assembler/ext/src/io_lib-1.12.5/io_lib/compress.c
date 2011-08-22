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
 * Handles compression and decompression.
 * Two functions are available. One compresses files, and the other opens
 * (read only) a compressed file and returns a FILE pointer.
 * Neither of these two are likely to work under Windows or MacOS.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_SYS_WAIT_H
#    include <sys/wait.h>
#    define DO_PIPE2
#endif
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#include "io_lib/os.h" /* for ftruncate() under WINNT */
#include "io_lib/compress.h"

#ifdef HAVE_ZLIB
#include <zlib.h>

/* ------------------------------------------------------------------------- */
/* GZIP reading and writing code via ZLIB. */

#define BS 8192

#define MAX_WBITS 15

#define FTEXT    (1<<0)
#define FHCRC    (1<<1)
#define FEXTRA   (1<<2)
#define FNAME    (1<<3)
#define FCOMMENT (1<<4)
/* Given a gzip file this returns the size of the gzip header */
static int gzheadersize(unsigned char *data) {
    int offset = 10;
    int flags = data[3];

    if (flags & FEXTRA)
	offset += 2 + data[offset] + data[offset+1]*256;

    if (flags & FNAME)
	while (data[offset++]);

    if (flags & FCOMMENT)
	while (data[offset++]);

    if (flags & FHCRC)
	offset += 2;

    return offset;
}

char *memgunzip(char *data, size_t size, size_t *udata_size) {
    int gzheader;
    z_stream s;
    char *udata = NULL;
    int udata_alloc = 0;
    int udata_pos = 0;

    /* Compute gzip header size */
    gzheader = gzheadersize((unsigned char *)data);

    /* Initialise zlib stream starting after the header */
    s.zalloc = (alloc_func)0;
    s.zfree = (free_func)0;
    s.opaque = (voidpf)0;
    s.next_in = (unsigned char *)data + gzheader;
    s.avail_in = size - gzheader;
    inflateInit2(&s, -MAX_WBITS);

    /* Decode to 'udata' array */
    for (;;) {
	int err;
	if (udata_alloc - udata_pos < 1) {
	    udata_alloc = udata_alloc ? udata_alloc * 2 : 256;
	    udata = realloc(udata, udata_alloc);
	}
	s.next_out = (unsigned char *)&udata[udata_pos];
	s.avail_out = udata_alloc - udata_pos;
	err = inflate(&s, Z_NO_FLUSH);
	udata_pos = udata_alloc - s.avail_out;
	if (err) {
	    if (err == Z_STREAM_END) {
		break;
	    } else {
		inflateEnd(&s);
		return NULL;
	    }
	}
    }

    inflateEnd(&s);
    *udata_size = udata_pos;

    return udata;
}

char *memgzip(char *data, size_t size, size_t *cdata_size) {
    z_stream s;
    char *cdata = NULL;
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;
    uint32_t i32;

    /* Create a minimal gzip header */
    cdata = malloc(cdata_alloc = size*1.02+10+8);
    memcpy(cdata, "\037\213\010\000\000\000\000\000\000\377", 10);
    cdata_pos = 10;

    /* Initialise zlib stream starting after the header */
    s.zalloc = (alloc_func)0;
    s.zfree = (free_func)0;
    s.opaque = (voidpf)0;
    s.next_in = (unsigned char *)data;
    s.avail_in = size;

    err = deflateInit2(&s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -MAX_WBITS,
		       9 /* DEF_MEM_LEVEL */, Z_DEFAULT_STRATEGY);

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = (unsigned char *)&cdata[cdata_pos];
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Gzip produced larger output than expected. Abort\n"); 
	    return NULL;
	}
	err = deflate(&s, Z_NO_FLUSH);
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK)
	    break;
    }
    deflate(&s, Z_FINISH);
    cdata_pos = 10 + s.total_out;

    i32 = crc32(0L, (unsigned char *)data, size);
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;

    i32 = size;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;
    cdata[cdata_pos++] = (int)(i32 & 0xff); i32 >>= 8;

    deflateEnd(&s);
    *cdata_size = cdata_pos;

    return cdata;
}
#endif

#ifdef DO_PIPE2
/* ------------------------------------------------------------------------- */
/* pipe_into - for piping via compression and decompression tools */

/*
 * This pipes 'input' data of length 'size' into a unix 'command'.
 * The output is then returned as an allocated block of memory. It is the
 * caller's responsibility to free this data.
 *
 * Returns malloc()ed data on success
 *         NULL on failure
 */
#define PIPEBS 8192
static char *pipe_into(const char *command, char *input, size_t insize,
		       size_t *outsize) {
    char *output = NULL;
    int output_alloc = 0;
    int output_used = 0;
    int fdp[2][2];
    fd_set rdfds, wrfds;
    int n = 0;
    pid_t pid;
    char buf[PIPEBS];
    int len, status;
    int eof_rd = 0, eof_wr = 0;

    /*
     * Make the connections:
     *
     * fdp[0] is stdin for the child
     * fdp[1] is stdout for the child
     * fdp[x][0] is the read end, and fdp[x][1] is the write end.
     * Hence:
     *     fdp[0][1] (parent's output) ->  (child's stdin) fdp[0][0]
     *     fdp[1][1] (child's stdout)  -> (parent's input) fdp[1][0]
     */

    if (-1 == pipe(fdp[0]))
        return NULL;

    if (-1 == pipe(fdp[1])) {
        close(fdp[0][0]);
        close(fdp[0][1]);
        return NULL;
    }

    if (n < fdp[1][0] + 1)
	n = fdp[1][0] + 1;
    if (n < fdp[0][1] + 1)
	n = fdp[0][1] + 1;

    switch(pid = fork()) {
    case 0: /* child */
        dup2(fdp[0][0], 0);
        dup2(fdp[1][1], 1);
        close(fdp[0][1]);
        close(fdp[1][0]);

        execlp("sh", "sh", "-c", command, NULL);
        exit(1);

    default: /* parent */
        close(fdp[0][0]);
        close(fdp[1][1]);
        break;

    case -1: /* error */
	return NULL;
    }

    /*
     * Set both parent ends to be non blocking. Deadlock can not be
     * completely avoided in a double pipe, so if something's going to
     * break we want to make sure it'll be the child, not the parent.
     */
    (void)fcntl(fdp[0][1], F_SETFL, O_NONBLOCK);
    (void)fcntl(fdp[1][0], F_SETFL, O_NONBLOCK);

    do {
	struct timeval tv;
	int r;

	FD_ZERO(&rdfds);
	FD_ZERO(&wrfds);
	if (!eof_wr)
	    FD_SET(fdp[0][1], &wrfds);
	if (!eof_rd)
	    FD_SET(fdp[1][0], &rdfds);

	tv.tv_sec = 1;
	tv.tv_usec = 0;
	if (-1 == (r = select(n, &rdfds, &wrfds, NULL, &tv)))
	    /* Handle EINTR etc... */
	    break;

	if (r) {
	    if (FD_ISSET(fdp[1][0], &rdfds)) {
		len = read(fdp[1][0], buf, PIPEBS);
		if (len > 0) {
		    while (output_used + len > output_alloc) {
			output_alloc = output_alloc ? output_alloc*2 : PIPEBS;
			output = realloc(output, output_alloc);
		    }
		    memcpy(&output[output_used], buf, len);
		    output_used += len;
		} else {
		    close(fdp[1][0]);
		    eof_rd = 1;
		}
	    }
	    if (FD_ISSET(fdp[0][1], &wrfds)) {
		len = write(fdp[0][1], input, insize>PIPEBS ? PIPEBS : insize);
		if (len > 0) {
		    input += len;
		    insize -= len;

		    if (insize == 0) {
			close(fdp[0][1]);
			eof_wr = 1;
		    }
		}
	    }
	}
	
    } while(!eof_rd || !eof_wr);

    close(fdp[0][1]);   /* should be closed already, but being doubly- */
    close(fdp[1][0]);   /* sure in case of error */
    waitpid(pid, &status, 0);

    *outsize = output_used;
    return output;
}
#endif /* DO_PIPE2 */

/* ------------------------------------------------------------------------- */
/* The main external routines for io_lib */

/*
 * This contains the last used compression method.
 */
static int compression_used = 0;

typedef struct {
    unsigned char magic[3];
    int magicl;
    char *compress;
    char *uncompress;
    char *extension;
} Magics;

/*
 * The list of magic numbers. The attempted order for compression is the
 * order of entries in this file.
 *
 * NB: bzip gives very good (better than gzip) results, is sometimes faster for
 * compression, but unfortunately much slower (4x?) for decompression. Most
 * people won't have it anyway.
 *
 * szip is definitely the best in compression ratios, and is faster than bzip.
 * However it's still slower than gzip. For comparable ratios, but much faster,
 * see the ztr format.
 */
static Magics magics[] = {
    {{'B',   'Z',    '0'},	3,	"bzip",		"bzip -d",   ".bz"},
    {{'\037', 0213, '\0'},	2,	"gzip",		"gzip -d",   ".gz"},
    {{'\037', 0235, '\0'},	2,	"compress",	"uncompress",".Z"},
    {{'B',   'Z',    'h'},	3,	"bzip2",	"bzip2 -d",  ".bz2"},
    {{'S',   'Z',   '\n'},	3,	"szip",	        "szip -d",   ".sz"},
};

void set_compression_method(int method) {
    compression_used = method;
}

int get_compression_method(void) {
    return compression_used;
}

/*
 * Converts compress mode strings (eg "gzip") to numbers.
 */
int compress_str2int(char *mode) {
    if (strcmp(mode, "bzip") == 0)
	return COMP_METHOD_BZIP;
    else if (strcmp(mode, "bzip2") == 0)
	return COMP_METHOD_BZIP2;
    else if (strcmp(mode, "gzip") == 0)
	return COMP_METHOD_GZIP;
    else if (strcmp(mode, "compress") == 0)
	return COMP_METHOD_COMPRESS;
    else if (strcmp(mode, "szip") == 0)
	return COMP_METHOD_SZIP;
    else return 0;

}
/*
 * Converts compress mode numbers to strings (eg "gzip").
 */
char *compress_int2str(int mode) {
    switch (mode) {
    case COMP_METHOD_BZIP:      return "bzip";
    case COMP_METHOD_GZIP:      return "gzip";
    case COMP_METHOD_BZIP2:     return "bzip2";
    case COMP_METHOD_COMPRESS:  return "compress";
    case COMP_METHOD_SZIP:      return "szip";
    }
    return "none";
}

/*
 * Compress a file using the method set in the compression_used value
 * (set by set_compression_method and fopen_compressed).
 *
 * If compression succeeds, we rename the file back its original name.
 *
 * When compression_used is 0 no compression is done.
 */
int compress_file(char *file) {
    char fname[2048];
    mFILE *mf;
    FILE *fp;

    /* Do nothing unless requested */
    if (compression_used == 0)
	return 0;

    mf = mfopen(file, "r");
    fcompress_file(mf);

    sprintf(fname, "%s%s", file, magics[compression_used-1].extension);
    if (NULL == (fp = fopen(fname, "wb")))
	return -1;

    fwrite(mf->data, 1, mf->size, fp);
    fclose(fp);
    mfclose(mf);


    return 0;
}

/*
 * Compress an mFILE using the method set in the compression_used value
 * (set by set_compression_method and fopen_compressed). This is done
 * in-memory by using a pipe to and from the compression program, or zlib
 * if we want to use gzip.
 *
 * When compression_used is 0 no compression is done.
 */
int fcompress_file(mFILE *fp) {
    size_t size;
    char *data;

    /* Do nothing unless requested */
    if (compression_used == 0)
	return 0;

#ifdef HAVE_ZLIB
    /*
     * If zlib is used then we use it to implement gzip internally, thus
     * saving starting up a separate process. This is substantially faster.
     */
    if (compression_used == 2) {
	data = memgzip(fp->data, fp->size, &size);
    } else
#endif
    {
#ifdef DO_PIPE2
	/*
	 * We have to pipe the data via an external tool, avoiding temporary
	 * files for speed.
	 */
	data = pipe_into(magics[compression_used-1].compress,
			 fp->data, fp->size, &size);
#else
	return -1;
#endif
    }

    mfrecreate(fp, data, size);
    mfseek(fp, size, SEEK_SET);

    return 0;
}


/*
 * Returns a file pointer of an uncompressed copy of 'file'.
 * 'file' need not exist if 'file'.ext (eg file.gz)
 * exists and can be uncompressed.
 *
 * NO LONGER SUPPORTED:-
 * If ofp is non NULL then the original file pointer will also be returned
 * (opened for update) to allow writing back to the original file. In cases
 * of uncompressed data this is the same as the returned file pointer.
 */
mFILE *fopen_compressed(char *file, mFILE **ofp) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    int i;
    char fext[1024];

    if (ofp) {
	fprintf(stderr, "ofp not supported in fopen_compressed() yet\n");
	*ofp = NULL;
    }

    /*
     * Try opening the file and reading the magic number.
     * If this doesn't work, then don't worry - the filename may be
     * the original name which has been renamed due to compression.
     * (eg file.gz).
     */
    for (i = -1; i < num_magics; i++) {
	mFILE *fp, *newfp;

	if (i == -1) {
	    if (NULL == (fp = mfopen(file, "rb")))
		continue;
	} else {
	    sprintf(fext, "%s%s", file, magics[i].extension);
	    if (NULL == (fp = mfopen(fext, "rb")))
		continue;
	}

	newfp = freopen_compressed(fp, NULL);
	if (fp != newfp)
	    /* Was compressed, so free compressed copy & return uncompressed */
	    mfclose(fp);

	if (newfp) {
	    return newfp;
	}
    }

    return NULL;
}

/*
 * Returns a file pointer of an uncompressed copy of 'fp'.
 * This may be the input fp or it may be a new fp.
 * The input fp is not modified and is left open. Therefore it is left up
 * to the caller to close the input fp and to check whether the returned fp
 * differs, and if so to close that too.
 */
mFILE *freopen_compressed(mFILE *fp, mFILE **ofp) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    unsigned char mg[3];
    int i;
    char *udata;
    size_t usize;

    if (ofp) {
	fprintf(stderr, "ofp not supported in fopen_compressed() yet\n");
	*ofp = NULL;
    }

    /* Test that it's compressed with full magic number */
    mfread(mg, 1, 3, fp);
    mrewind(fp);
    for (i = 0; i < num_magics; i++) {
	if (0 == memcmp(mg, magics[i].magic, magics[i].magicl))
	    break;
    }
    if (i == num_magics) {
	compression_used = 0;
	return fp;
    }

#ifdef HAVE_ZLIB
    if (i == 1) {
	udata = memgunzip(fp->data, fp->size, &usize);
    } else
#endif
    {
#ifdef DO_PIPE2
	udata = pipe_into(magics[i].uncompress, fp->data, fp->size, &usize);
#else
	return NULL;
#endif
    }

    compression_used = i+1;

    return mfcreate(udata, usize);
}

/*
 * Given a filename remove a known compression extension
 *
 * Returns: None
 */
void remove_extension(char *file) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    int i;
    for (i=0;i<num_magics;i++) {
      char *cp = file+strlen(file)-strlen(magics[i].extension);
      if (strcmp(cp, magics[i].extension) == 0) {
            *cp = '\0';
            return;
        }
    }
    return;
}
