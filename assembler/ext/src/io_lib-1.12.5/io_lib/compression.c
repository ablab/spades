#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#include "io_lib/ztr.h"
#include "io_lib/os.h"
#include "io_lib/compression.h"
#include "io_lib/xalloc.h"

#ifndef NDEBUG 
#  define NDEBUG 
#endif 

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_ZLIB
 * ---------------------------------------------------------------------------
 */

/*
 * Some comments on zlib usage.
 *
 * - Ideally for trace data, after decorrelation, we should use Z_FILTERED.
 *   Empirical studies show that this gives the best compression ratio, but
 *   it is slow (about the same speed as normal gzip). MUCH faster is huffman
 *   only, and it doesn't give radically different compression ratios.
 *
 * - When compressing using Z_HUFFMAN_ONLY we used compression level
 *   '1' as this invokes the deflate_fast() algorithm. It makes no
 *   difference to the compression level, but it seems to be quicker still.
 *
 */

/*
 * zlib_huff()
 *
 * Compresses data using huffman encoding, as implemented by zlib.
 *
 * Arguments:
 *	uncomp		Uncompressed input data
 *	uncomp_len	Length of uncomp data
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *zlib_huff(char *uncomp, int uncomp_len, int strategy, int *comp_len) {
    z_stream zstr;
    int err;
    int comp_len_tmp = (int)(uncomp_len * 1.001 + 12); /* Maximum expansion size */
    char *comp = (char *)xmalloc(comp_len_tmp+5);
    int c_len;

    /* Initialise zlib */
    zstr.zalloc = (alloc_func)0;
    zstr.zfree = (free_func)0;
    zstr.opaque = (voidpf)0;

    if ((err = deflateInit2(&zstr, 1, Z_DEFLATED, 15, 8, strategy)) != Z_OK) {
	fprintf(stderr, "zlib errror in deflateInit2(): %d\n", err);
	return NULL;
    }

    /* Set up input and output buffers */
    zstr.next_in = (unsigned char *)uncomp;
    zstr.avail_in = uncomp_len;
    zstr.next_out = (unsigned char *)comp+5;
    zstr.avail_out = comp_len_tmp;
    
    /* Do the compression */
    if ((err = deflate(&zstr, Z_FINISH)) != Z_STREAM_END) {
	fprintf(stderr, "zlib errror in deflate(): %d\n", err);
	return NULL;
    }

    /* Tidy up */
    deflateEnd(&zstr);

    c_len = zstr.total_out;

    /* Return */
    comp[0] = ZTR_FORM_ZLIB;
    comp[1] = (uncomp_len >>  0) & 0xff;
    comp[2] = (uncomp_len >>  8) & 0xff;
    comp[3] = (uncomp_len >> 16) & 0xff;
    comp[4] = (uncomp_len >> 24) & 0xff;

    if (comp_len)
	*comp_len = c_len+5;

    return comp;
}

/*
 * zlib_dehuff()
 *
 * Uncompresses data using huffman encoding, as implemented by zlib.
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *zlib_dehuff(char *comp, int comp_len, int *uncomp_len) {
    z_stream zstr;
    int err;
    char *uncomp;
    int ulen;

    /* Allocate */
    ulen =
	((unsigned char)comp[1] <<  0) +
	((unsigned char)comp[2] <<  8) +
	((unsigned char)comp[3] << 16) +
	((unsigned char)comp[4] << 24);
    uncomp = (char *)xmalloc(ulen);

    /* Initialise zlib */
    zstr.zalloc = (alloc_func)0;
    zstr.zfree = (free_func)0;
    zstr.opaque = (voidpf)0;

    if ((err = inflateInit(&zstr)) != Z_OK) {
	fprintf(stderr, "zlib errror in inflateInit(): %d\n", err);
	return NULL;
    }

    /* Set up input and output buffers */
    zstr.next_in = (unsigned char *)comp+5;
    zstr.avail_in = comp_len-5;
    zstr.next_out = (unsigned char *)uncomp;
    zstr.avail_out = ulen;
    
    /* Do the decompression */
    if ((err = inflate(&zstr, Z_FINISH)) != Z_STREAM_END) {
	fprintf(stderr, "zlib errror in deflate(): %d\n", err);
	return NULL;
    }

    /* Tidy up */
    inflateEnd(&zstr);

    if (uncomp_len)
	*uncomp_len = ulen;

    return uncomp;
}


/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_RLE
 * ---------------------------------------------------------------------------
 */

/*
 * Run length encoding.
 *
 * Any run of 3 or more identical characters (up to 255 in a row) are replaced
 * by a 'guard' byte followed by the number of characters followed by
 * the character value itself.
 * Any single guard value in the input is escaped using 'guard 0'.
 *
 * Specifying guard as -1 will automatically pick one of the least used
 * characters in the input as the guard.
 *
 * Arguments:
 *	uncomp		Input data
 *	uncomp_len	Length of input data 'uncomp'
 *	guard		Guard byte - used to encode "N" copies of data
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *rle(char *uncomp, int uncomp_len, int guard, int *comp_len) {
    int i, k, c_len = 0;
    char *comp = xmalloc(2 * uncomp_len + 6);
    char *out = comp + 6;

    /* A guard of -1 implies to search for the least frequent symbol */
    if (guard == -1) {
	int cnt[256];
	int bestcnt = uncomp_len + 1;

	for (i = 0; i < 256; i++)
	    cnt[i] = 0;

	for (i = 0; i < uncomp_len; i++) {
	    cnt[(unsigned char)uncomp[i]]++;
	}

	for (i = 0; i < 256; i++) {
	    if (cnt[i] < bestcnt) {
		bestcnt = cnt[i];
		guard = i;
	    }
	}
    }

    for (i = 0; i < uncomp_len; i=k) {
	/*
	 * Detect blocks of up identical bytes up to 255 bytes long.
	 */
	for (k = i; k < uncomp_len && uncomp[i] == uncomp[k]; k++)
	    if (k-i == 255)
		break;

	/* 1, 2 or 3 bytes are best stored "as is" */
	if (k-i < 4) {
	    do {
		/*
		 * If we find 'guard' in our sequence, escape it by
		 * outputting 'guard' <zero>. (We know that we'll never
		 * write out zero copies of a token in our rle compression
		 * algorithm.)
		 */
		if ((unsigned char)(uncomp[i]) == guard) {
		    out[c_len++] = guard;
		    out[c_len++] = 0;
		} else {
		    out[c_len++] = uncomp[i];
		}
		i++;
	    } while (k >= i+1);
	} else {
	    /* More than 3 bytes: store as ('guard', length, byte value) */
	    out[c_len++] = guard;
	    out[c_len++] = k-i;
	    out[c_len++] = uncomp[i];
	}
    }

    /* Return */
    comp[0] = ZTR_FORM_RLE;
    comp[1] = (uncomp_len >>  0) & 0xff;
    comp[2] = (uncomp_len >>  8) & 0xff;
    comp[3] = (uncomp_len >> 16) & 0xff;
    comp[4] = (uncomp_len >> 24) & 0xff;
    comp[5] = guard;

    if (comp_len)
	*comp_len = c_len+6;

    return comp;
}

/*
 * Reverses run length encoding.
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
/* ARGSUSED */
char *unrle(char *comp, int comp_len, int *uncomp_len) {
    int in_i, out_i, i, val, count, out_len;
    char *uncomp;
    unsigned char *in = (unsigned char *)comp+6;
    char *out;
    int guard = (unsigned char)comp[5];

    /* Allocate */
    out_len =
	((unsigned char)comp[1] <<  0) +
	((unsigned char)comp[2] <<  8) +
	((unsigned char)comp[3] << 16) +
	((unsigned char)comp[4] << 24);
    out = uncomp = (char *)xmalloc(out_len);
    
    for (in_i = out_i = 0; out_i < out_len; in_i++) {
	if (in[in_i] != guard) {
	    /* When not 'guard' it's easy - just output this token */
	    assert(out_i >= 0 && out_i < out_len);
	    out[out_i++] = in[in_i];

	} else {
	    /*
	     * Found an 'guard' token. If next token is zero, then
	     * we were simply escaping a real 'guard' token in the input
	     * data, otherwise output a string of bytes.
	     */
	    count = in[++in_i];
	    if (count != 0) {
		val = in[++in_i];
		for (i = 0; i < count; i++) {
		    assert(out_i >= 0 && out_i < out_len);
		    out[out_i++] = val;
		}
	    } else {
		assert(out_i >= 0 && out_i < out_len);
		out[out_i++] = guard;
	    }
	}
    }

    if (uncomp_len)
	*uncomp_len = out_len;

    return uncomp;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_XRLE
 * ---------------------------------------------------------------------------
 */

/*
 * Mutli-byte run length encoding.
 *
 * Any run of 3 or more identical characters (up to 255 in a row) are replaced
 * by a 'guard' byte followed by the number of characters followed by
 * the character value itself.
 * Any single guard value in the input is escaped using 'guard 0'.
 *
 * Specifying guard as -1 will automatically pick one of the least used
 * characters in the input as the guard.
 *
 * Arguments:
 *	uncomp		Input data
 *	uncomp_len	Length of input data 'uncomp'
 *	guard		Guard byte - used to encode "N" copies of data
 *      rsz             Size of blocks to compare for run checking.
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *xrle(char *uncomp, int uncomp_len, int guard, int rsz, int *comp_len) {
    char *comp = (char *)malloc(2 * uncomp_len + 3);
    char *out = comp;
    int i, j, k;

    /* A guard of -1 implies to search for the least frequent symbol */
    if (guard == -1) {
	int cnt[256];
	int bestcnt = uncomp_len + 1;

	for (i = 0; i < 256; i++)
	    cnt[i] = 0;

	for (i = 0; i < uncomp_len; i++) {
	    cnt[(unsigned char)uncomp[i]]++;
	}

	for (i = 0; i < 256; i++) {
	    if (cnt[i] < bestcnt) {
		bestcnt = cnt[i];
		guard = i;
	    }
	}
    }

    *out++ = ZTR_FORM_XRLE;
    *out++ = rsz;
    *out++ = guard;

    for (i = 0; i < uncomp_len; i = k) {
	/* Count repeats */
	k = i + rsz;
	while (k <= uncomp_len - rsz && !memcmp(&uncomp[i], &uncomp[k], rsz)) {
	       k += rsz;
	       if ((k-i)/rsz == 255)
		   break;
	}

	if (k-i > rsz) {
	    /* Duplicates, so RLE */
	    *out++ = guard;
	    *out++ = (k-i)/rsz;
	    for (j = 0; j < rsz; j++)
		*out++ = uncomp[i+j];
	} else {
	    /* No dups, store as is escaping guarding as appropriate */
	    if ((unsigned char)(uncomp[i]) == guard) {
		*out++ = guard;
		*out++ = 0;
	    } else {
		*out++ = uncomp[i];
	    }
	    k = i+1;
	}
    }

    *comp_len = out-comp;

    return comp;
}


/*
 * Reverses multi-byte run length encoding.
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *unxrle(char *comp, int comp_len, int *uncomp_len) {
    char *uncomp;
    char *out;
    int rsz = (unsigned char)comp[1];
    int guard = (unsigned char)comp[2];
    unsigned char *in;
    int unclen = 0, cpos, len;
    
    /* Calculate uncompressed length */
    for (in = (unsigned char *)comp, cpos = 3; cpos < comp_len; unclen++) {
	if (in[cpos++] == guard) {
	    if ((len = in[cpos++])) {
		cpos += rsz;
		unclen += len*rsz -1;
	    }
	}
    }
    *uncomp_len = unclen;

    /* Expand */
    uncomp = out = (char *)malloc(unclen+1);
    for (in = (unsigned char *)comp, cpos = 3; cpos < comp_len;) {
	char c;
	if ((c = in[cpos++]) != guard) {
	    *out++ = c;
	} else {
	    int len = in[cpos++];
	    if (len) {
		int i, j;
		for (i = 0; i < len; i++) {
		    for (j = 0; j < rsz; j++) {
			*out++ = in[cpos+j];
		    }
		}
		cpos += rsz;
	    } else {
		*out++ = guard;
	    }
	}
    }
    *out++ = 0;

    return uncomp;
}

/*
 * Mutli-byte run length encoding.
 *
 * Steps along in words of size 'rsz'. Unlike XRLE above this does run-length
 * encoding by writing out an additional "length" word every time 2 or more
 * words in a row are spotted. This removes the need for a guard byte.
 *
 * Additionally this method ensures that both input and output formats remain
 * aligned on words of size 'rsz'.
 *
 * Arguments:
 *	uncomp		Input data
 *	uncomp_len	Length of input data 'uncomp'
 *      rsz             Size of blocks to compare for run checking.
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *xrle2(char *uncomp, int uncomp_len, int rsz, int *comp_len) {
    char *comp = (char *)malloc(1.4*uncomp_len + rsz);
    char *last = uncomp;
    char *out = comp;
    int i, j, k, run_len = 0;

    *out++ = ZTR_FORM_XRLE2;
    *out++ = rsz;
    for (i = 2; i < rsz; i++)
	*out++ = -40;

    /* FIXME: how to deal with uncomp_len not being a multiple of rsz */
    for (i = 0; i < uncomp_len; i += rsz) {
	/* FIXME: use inline #def versions of memcmp/memcpy for speed? */
	memcpy(out,  &uncomp[i], rsz);
	out += rsz;

	if (memcmp(last, &uncomp[i], rsz) == 0) {
	    run_len++;
	} else {
	    run_len = 1;
	    last = &uncomp[i];
	}

	/* NB: >= 3 is more optimal in many cases */
	if (run_len >= 2) {
	    /* Count remaining copies */
	    for (k = i+rsz; k < uncomp_len && run_len < 257; k += rsz) {
		if (memcmp(last, &uncomp[k], rsz) != 0)
		    break;
		run_len++;
	    }
	    run_len -= 2;

	    *out++ = run_len;
	    for (j = 1; j < rsz; j++) {
		/* Padding with last reduces entropy compared to padding
		 * with copies of run_len
		 */
		*out++ = last[j];
	    }
	    i = k-rsz;

	    last = out-rsz;
	    run_len = 0;
	}
    }

    *comp_len = out-comp;

    return comp;
}

/*
 * Reverses multi-byte run length encoding (xrle_new).
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *unxrle2(char *comp, int comp_len, int *uncomp_len) {
    char *out, *last;
    int out_len, out_alloc, rsz, i, j, run_len;

    out_alloc = comp_len*2; /* just an estimate */
    out_len = 0;
    if (NULL == (out = (char *)malloc(out_alloc)))
	return NULL;

    if (*comp++ != ZTR_FORM_XRLE2)
	return NULL;

    /* Read rsz and swallow padding */
    rsz = *comp++;
    comp_len -= 2;
    for (i = 2; i < rsz; i++) {
	comp++;
	comp_len--;
    }
    
    /* Uncompress */
    run_len = 0;
    last = comp;
    for (i = 0; i < comp_len;) {
	while (out_len + rsz > out_alloc) {
	    out_alloc *= 2;
	    if (NULL == (out = (char *)realloc(out, out_alloc)))
		return NULL;
	}
	memcpy(&out[out_len], &comp[i], rsz);

	if (memcmp(&out[out_len], last, rsz) == 0) {
	    run_len++;
	} else {
	    run_len = 1;
	}

	i += rsz;
	out_len += rsz;

	/* NB: >= 3 is more optimal in many cases */
	if (run_len >= 2) {
	    /* Count remaining copies */
	    run_len = (unsigned char)comp[i];
	    i += rsz;

	    while (out_len + run_len * rsz > out_alloc) {
		out_alloc *= 2;
		if (NULL == (out = (char *)realloc(out, out_alloc)))
		    return NULL;
	    }

	    for (j = 0; j < run_len; j++) {
		memcpy(&out[out_len], last, rsz);
		out_len += rsz;
	    }
	    run_len = 0;
	}

	last = &comp[i-rsz];
    }

    /* Shrink back down to avoid excessive memory usage */
    out = realloc(out, out_len);
    *uncomp_len = out_len;

    return out;
}


/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_DELTA1
 * ---------------------------------------------------------------------------
 */

/*
 * This replaces 'samples' with successive differences between samples.
 * These implementations support 'level's of 1, 2 and 3.
 *
 * NB: This is analogous to our SCF delta_samples1 (etc) function, except that
 * this function about 40% faster.
 *
 * Implementation details taken from Jean Thierry-Mieg's CTF code.
 */

/*
 * decorrelate1()
 *
 * Produce successive deltas from a 1-byte array.
 *
 * Arguments:
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data
 *	level		Differencing level (must be 1, 2 or 3)
 *	comp_len	Return: where to store new compressed length
 *
 * Returns:
 *	Success: A decorrelated buffer (malloced)
 *	Failure: NULL
 */
char *decorrelate1(char *x_uncomp,
		   int uncomp_len,
		   int level,
		   int *comp_len) {
    int i, z;
    int u1 = 0, u2 = 0, u3 = 0;
    char *comp = (char *)xmalloc(uncomp_len + 2);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;

    if (!comp)
	return NULL;

    comp+=2;
    switch (level) {
    case 1:
	for (i = 0; i < uncomp_len; i++) {
	    z = u1;
	    u1 = u_uncomp[i];
	    comp[i] = u_uncomp[i] - z;
	}
	break;
	
    case 2:
	for (i = 0; i < uncomp_len; i++) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 = u_uncomp[i];
	    comp[i] = u_uncomp[i] - z;
	}
	break;

    case 3:
	for (i = 0; i < uncomp_len; i++) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 = u_uncomp[i];
	    comp[i] = u_uncomp[i] - z;
	}
	break;

    default:
	return NULL;
    }
    comp-=2;
    comp[0] = ZTR_FORM_DELTA1;
    comp[1] = level;

    *comp_len = uncomp_len+2;

    return comp;
}

#define ABS(a) ((a) > 0 ? (a) : -(a))

/* ZTR_FORM_DDELTA1 - experimental */
char *decorrelate1dyn(char *x_uncomp,
		      int uncomp_len,
		      int *comp_len) {
    int i, j, z[4];
    int u1 = 0, u2 = 0, u3 = 0;
    char *comp = (char *)xmalloc(uncomp_len + 2);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;
    int level = 3; /* default level */
    int last_level = level;
    int best;

    if (!comp)
	return NULL;

    comp+=2;
    for (i = 0; i < uncomp_len; i++) {
	z[1] = u1;
	z[2] = 2*u1 - u2;
	z[3] = 3*u1 - 3*u2 + u3;
	comp[i] = u_uncomp[i] - z[last_level];
	best = 10000;
	for (j = 1; j < 3; j++) {
	    if (ABS(u_uncomp[i] - z[j]) < best) {
		best = ABS(u_uncomp[i] - z[j]);
		last_level = j;
	    }
	}
	u3 = u2;
	u2 = u1;
	u1 = u_uncomp[i];
    }
    comp-=2;
    comp[0] = ZTR_FORM_DDELTA1;
    comp[1] = level;

    *comp_len = uncomp_len+2;

    return comp;
}

/*
 * recorrelate1()
 *
 * The reverse of decorrelate1()
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Success: uncompressed data
 *	Failure: NULL
 */
char *recorrelate1(char *x_comp,
		   int comp_len,
		   int *uncomp_len) {
    int i, z;
    int u1 = 0, u2 = 0, u3 = 0;
    int level = x_comp[1];
    char *uncomp;

    uncomp = (char *)xmalloc(comp_len-2);
    if (!uncomp)
	return NULL;

    x_comp+=2;
    comp_len-=2;
    *uncomp_len = comp_len;

    switch (level) {
    case 1:
	for (i = 0; i < comp_len; i++) {
	    z = u1;
	    u1 = uncomp[i] = x_comp[i] + z;
	}
	break;

    case 2:
	for (i = 0; i < comp_len; i++) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 = uncomp[i] = x_comp[i] + z;
	}
	break;
	
    case 3:
	for (i = 0; i < comp_len; i++) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 = uncomp[i] = x_comp[i] + z;
	}
	break;
    }

    return uncomp;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_DELTA2
 * ---------------------------------------------------------------------------
 */

/*
 * decorrelate2()
 *
 * Produce successive deltas from a 2-byte array (big endian)
 *
 * Arguments:
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data
 *	level		Differencing level (must be 1, 2 or 3)
 *	comp_len	Return: where to store new compressed length
 *
 * Returns:
 *	Success: A decorrelated buffer (malloced)
 *	Failure: NULL
 */
char *decorrelate2(char *x_uncomp,
		   int uncomp_len,
		   int level,
		   int *comp_len) {
    int i, z, delta;
    int u1 = 0, u2 = 0, u3 = 0;
    char *comp = (char *)xmalloc(uncomp_len + 2);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;

    if (!comp)
	return NULL;

    comp+=2;
    switch (level) {
    case 1:
	for (i = 0; i < uncomp_len; i+=2) {
	    z = u1;
	    u1 = (u_uncomp[i] << 8) + u_uncomp[i+1];
	    delta = u1 - z;
	    comp[i  ] = (delta >> 8) & 0xff;
	    comp[i+1] = (delta >> 0) & 0xff;
	}
	break;
	
    case 2:
	for (i = 0; i < uncomp_len; i+=2) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 = (u_uncomp[i] << 8) + u_uncomp[i+1];
	    delta = u1 - z;
	    comp[i  ] = (delta >> 8) & 0xff;
	    comp[i+1] = (delta >> 0) & 0xff;
	}
	break;

    case 3:
	for (i = 0; i < uncomp_len; i+=2) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 = (u_uncomp[i] << 8) + u_uncomp[i+1];
	    delta = u1 - z;
	    comp[i  ] = (delta >> 8) & 0xff;
	    comp[i+1] = (delta >> 0) & 0xff;
	}
	break;

    default:
	return NULL;
    }
    comp-=2;
    comp[0] = ZTR_FORM_DELTA2;
    comp[1] = level;

    *comp_len = uncomp_len+2;

    return comp;
}

char *decorrelate2dyn(char *x_uncomp,
		      int uncomp_len,
		      int *comp_len) {
    int i, j, z[4];
    int u1 = 0, u2 = 0, u3 = 0;
    char *comp = (char *)xmalloc(uncomp_len + 2);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;
    int level = 2; /* minimum level */
    int last_level = level;
    int best;

    if (!comp)
	return NULL;

    comp+=2;
    for (i = 0; i < uncomp_len; i+=2) {
	z[0] = 0;
	z[1] = u1;
	z[2] = 2*u1 - u2;
	z[3] = 3*u1 - 3*u2 + u3;
	u3 = u2;
	u2 = u1;
	u1 = (u_uncomp[i]<<8) + u_uncomp[i+1];
	comp[i  ] = ((u1 - z[last_level]) >> 8) & 0xff;
	comp[i+1] = ((u1 - z[last_level]) >> 0) & 0xff;
	best = 10000;
	for (j = level; j < 4; j++) {
	    if (ABS(u1 - z[j]) < best) {
		best = ABS(u1 - z[j]);
		last_level = j;
	    }
	}
    }
    comp-=2;
    comp[0] = ZTR_FORM_DDELTA2;
    comp[1] = level;

    *comp_len = uncomp_len+2;

    return comp;
}

/*
 * recorrelate2()
 *
 * The reverse of decorrelate2()
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Success: uncompressed data
 *	Failure: NULL
 */
char *recorrelate2(char *x_comp,
		   int comp_len,
		   int *uncomp_len) {
    int i, z;
    int u1 = 0, u2 = 0, u3 = 0;
    int level = x_comp[1];
    char *uncomp;
    unsigned char *u_comp = (unsigned char *)x_comp;

    uncomp = (char *)xmalloc(comp_len-2);
    if (!uncomp)
	return NULL;

    u_comp+=2;
    comp_len-=2;
    *uncomp_len = comp_len;

    switch (level) {
    case 1:
	for (i = 0; i < comp_len; i+=2) {
	    z = u1;
	    u1 = ((u_comp[i] << 8) | u_comp[i+1]) + z;
	    uncomp[i  ] = (u1 >> 8) & 0xff;
	    uncomp[i+1] = (u1 >> 0) & 0xff;
	}
	break;

    case 2:
	for (i = 0; i < comp_len; i+=2) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 = ((u_comp[i] << 8) | u_comp[i+1]) + z;
	    uncomp[i  ] = (u1 >> 8) & 0xff;
	    uncomp[i+1] = (u1 >> 0) & 0xff;
	}
	break;
	
    case 3:
	for (i = 0; i < comp_len; i+=2) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 = ((u_comp[i] << 8) | u_comp[i+1]) + z;
	    uncomp[i  ] = (u1 >> 8) & 0xff;
	    uncomp[i+1] = (u1 >> 0) & 0xff;
	}
	break;
    }

    return uncomp;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_DELTA4
 * ---------------------------------------------------------------------------
 */

/*
 * decorrelate4()
 *
 * Produce successive deltas from a 4-byte array (big endian)
 *
 * Arguments:
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data
 *	level		Differencing level (must be 1, 2 or 3)
 *	comp_len	Return: where to store new compressed length
 *
 * Returns:
 *	Success: A decorrelated buffer (malloced)
 *	Failure: NULL
 */
char *decorrelate4(char *x_uncomp,
		   int uncomp_len,
		   int level,
		   int *comp_len) {
    int i, z, delta;
    int u1 = 0, u2 = 0, u3 = 0;
    char *comp = (char *)xmalloc(uncomp_len + 4);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;

    if (!comp)
	return NULL;

    comp+=4;
    switch (level) {
    case 1:
	for (i = 0; i < uncomp_len; i+=4) {
	    z = u1;
	    u1 =(u_uncomp[i  ] << 24) +
		(u_uncomp[i+1] << 16) +
		(u_uncomp[i+2] <<  8) +
		(u_uncomp[i+3] <<  0);
	    delta = u1 - z;
	    comp[i  ] = (delta >> 24) & 0xff;
	    comp[i+1] = (delta >> 16) & 0xff;
	    comp[i+2] = (delta >>  8) & 0xff;
	    comp[i+3] = (delta >>  0) & 0xff;
	}
	break;
	
    case 2:
	for (i = 0; i < uncomp_len; i+=4) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 =(u_uncomp[i  ] << 24) +
		(u_uncomp[i+1] << 16) +
		(u_uncomp[i+2] <<  8) +
		(u_uncomp[i+3] <<  0);
	    delta = u1 - z;
	    comp[i  ] = (delta >> 24) & 0xff;
	    comp[i+1] = (delta >> 16) & 0xff;
	    comp[i+2] = (delta >>  8) & 0xff;
	    comp[i+3] = (delta >>  0) & 0xff;
	}
	break;

    case 3:
	for (i = 0; i < uncomp_len; i+=4) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 =(u_uncomp[i  ] << 24) +
		(u_uncomp[i+1] << 16) +
		(u_uncomp[i+2] <<  8) +
		(u_uncomp[i+3] <<  0);
	    delta = u1 - z;
	    comp[i  ] = (delta >> 24) & 0xff;
	    comp[i+1] = (delta >> 16) & 0xff;
	    comp[i+2] = (delta >>  8) & 0xff;
	    comp[i+3] = (delta >>  0) & 0xff;
	}
	break;

    default:
	return NULL;
    }
    comp-=4;
    comp[0] = ZTR_FORM_DELTA4;
    comp[1] = level;
    comp[2] = 0; /* dummy - to align on 4-byte boundary */
    comp[3] = 0; /* dummy - to align on 4-byte boundary */

    *comp_len = uncomp_len+4;

    return comp;
}

/*
 * recorrelate4()
 *
 * The reverse of decorrelate4()
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Success: uncompressed data
 *	Failure: NULL
 */
char *recorrelate4(char *x_comp,
		   int comp_len,
		   int *uncomp_len) {
    int i, z;
    int u1 = 0, u2 = 0, u3 = 0;
    int level = x_comp[1];
    char *uncomp;
    unsigned char *u_comp = (unsigned char *)x_comp;

    uncomp = (char *)xmalloc(comp_len-4);
    if (!uncomp)
	return NULL;

    u_comp+=4;
    comp_len-=4;
    *uncomp_len = comp_len;

    switch (level) {
    case 1:
	for (i = 0; i < comp_len; i+=4) {
	    z = u1;
	    u1 = z +
		((u_comp[i  ] << 24) |
		 (u_comp[i+1] << 16) |
		 (u_comp[i+2] <<  8) |
		  u_comp[i+3]);
	    uncomp[i  ] = (u1 >> 24) & 0xff;
	    uncomp[i+1] = (u1 >> 16) & 0xff;
	    uncomp[i+2] = (u1 >>  8) & 0xff;
	    uncomp[i+3] = (u1 >>  0) & 0xff;
	}
	break;

    case 2:
	for (i = 0; i < comp_len; i+=4) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 = z +
		((u_comp[i  ] << 24) |
		 (u_comp[i+1] << 16) |
		 (u_comp[i+2] <<  8) |
		  u_comp[i+3]);
	    uncomp[i  ] = (u1 >> 24) & 0xff;
	    uncomp[i+1] = (u1 >> 16) & 0xff;
	    uncomp[i+2] = (u1 >>  8) & 0xff;
	    uncomp[i+3] = (u1 >>  0) & 0xff;
	}
	break;
	
    case 3:
	for (i = 0; i < comp_len; i+=4) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 = z +
		((u_comp[i  ] << 24) |
		 (u_comp[i+1] << 16) |
		 (u_comp[i+2] <<  8) |
		  u_comp[i+3]);
	    uncomp[i  ] = (u1 >> 24) & 0xff;
	    uncomp[i+1] = (u1 >> 16) & 0xff;
	    uncomp[i+2] = (u1 >>  8) & 0xff;
	    uncomp[i+3] = (u1 >>  0) & 0xff;
	}
	break;
    }

    return uncomp;
}


/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_16TO8
 * ---------------------------------------------------------------------------
 */

/*
 * shrink_16to8()
 *
 * Stores an array of 16-bit (big endian) array elements in an 8-bit array.
 * We assume that most 16-bit elements encode numbers that fit in an 8-bit
 * value. When not possible, we store a marker followed by the 16-bit value
 * stored as multiple 8-bit values.
 *
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data (in bytes)
 *	comp_len	Return: where to store new compressed length
 *	
 * Returns:
 *	Success: An 8-bit array (malloced)
 *	Failure: NULL
 */
char *shrink_16to8(char *x_uncomp, int uncomp_len, int *comp_len) {
    char *comp;
    int i, j, i16;
    signed char *s_uncomp = (signed char *)x_uncomp;

    /* Allocation - worst case is 3 * (uncomp_len/2) + 1 */
    if (NULL == (comp = (char *)xmalloc(3 * (uncomp_len/2) + 1)))
	return NULL;

    comp[0] = ZTR_FORM_16TO8;
    for (i = 0, j = 1; i < uncomp_len; i+=2) {
	i16 = (s_uncomp[i] << 8) | (unsigned char)s_uncomp[i+1];
	if (i16 >= -127 && i16 <= 127) {
	    comp[j++] = i16;
	} else {
	    comp[j++] = -128;
	    comp[j++] = s_uncomp[i];
	    comp[j++] = s_uncomp[i+1];
	}
    }

    /* Reclaim unneeded memory */
    comp = xrealloc(comp, j);
    
    *comp_len = j;
    return comp;
}


/*
 * expand_8to16()
 *
 * The opposite of the shrink_16to8() function.
 *
 *	comp		Compressed input data
 *	comp_len	Length of comp data (in bytes)
 *	uncomp_len	Output: length of uncompressed data (in bytes)
 *	
 * Returns:
 *	Success: Uncompressed data (char *)
 *	Failure: NULL
 */
char *expand_8to16(char *x_comp, int comp_len, int *uncomp_len) {
    int i, j;
    char *uncomp;
    signed char *s_comp = (signed char *)x_comp;

    /* Allocation - worst case is twice comp_len */
    if (NULL == (uncomp = (char *)xmalloc(comp_len*2)))
	return NULL;

#if 0
    for (i = 0, j = 1; j < comp_len; i+=2) {
	if (s_comp[j] != -128) {
	    uncomp[i  ] = s_comp[j] < 0 ? -1 : 0;
	    uncomp[i+1] = s_comp[j++];
	} else {
	    j++;
	    uncomp[i  ] = s_comp[j++];
	    uncomp[i+1] = s_comp[j++];
	}
    }
#endif

    for (i = 0, j = 1; j < comp_len; i+=2) {
	if (s_comp[j] >= 0) {
	    uncomp[i  ] = 0;
	    uncomp[i+1] = s_comp[j++];
	} else {
	    if (s_comp[j] != -128) {
		uncomp[i+1] = s_comp[j++];
		uncomp[i  ] = -1;
	    } else {
		j++;
		uncomp[i  ] = s_comp[j++];
		uncomp[i+1] = s_comp[j++];
	    }
	}
    }

    /* Reclaim unneeded memory */
    uncomp = xrealloc(uncomp, i);
    
    *uncomp_len = i;
    return uncomp;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_32TO8
 * ---------------------------------------------------------------------------
 */

/*
 * shrink_32to8()
 *
 * Stores an array of 32-bit (big endian) array elements in an 8-bit array.
 * We assume that most 32-bit elements encode numbers that fit in an 8-bit
 * value. When not possible, we store a marker followed by the 32-bit value
 * stored as multiple 8-bit values.
 *
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data (in bytes)
 *	comp_len	Return: where to store new compressed length
 *	
 * Returns:
 *	Success: An 8-bit array (malloced)
 *	Failure: NULL
 */
char *shrink_32to8(char *x_uncomp, int uncomp_len, int *comp_len) {
    char *comp;
    int i, j, i32;
    signed char *s_uncomp = (signed char *)x_uncomp;

    /* Allocation - worst case is 5 * (uncomp_len/4) + 1 */
    if (NULL == (comp = (char *)xmalloc(5 * (uncomp_len/4) + 1)))
	return NULL;

    comp[0] = ZTR_FORM_32TO8;
    for (i = 0, j = 1; i < uncomp_len; i+=4) {
	i32 = (s_uncomp[i] << 24) |
	    (s_uncomp[i+1] << 16) |
	    (s_uncomp[i+2] <<  8) |
	    (unsigned char)s_uncomp[i+3];
	if (i32 >= -127 && i32 <= 127) {
	    comp[j++] = i32;
	} else {
	    comp[j++] = -128;
	    comp[j++] = s_uncomp[i];
	    comp[j++] = s_uncomp[i+1];
	    comp[j++] = s_uncomp[i+2];
	    comp[j++] = s_uncomp[i+3];
	}
    }

    /* Reclaim unneeded memory */
    comp = xrealloc(comp, j);
    
    *comp_len = j;
    return comp;
}

/*
 * expand_8to32()
 *
 * The opposite of the shrink_32to8() function.
 *
 *	comp		Compressed input data
 *	comp_len	Length of comp data (in bytes)
 *	uncomp_len	Output: length of uncompressed data (in bytes)
 *	
 * Returns:
 *	Success: Uncompressed data (char *)
 *	Failure: NULL
 */
char *expand_8to32(char *comp, int comp_len, int *uncomp_len) {
    int i, j;
    char *uncomp;
    signed char *s_comp = (signed char *)comp; 

    /* Allocation - worst case is four times comp_len */
    if (NULL == (uncomp = (char *)xmalloc(comp_len*4)))
	return NULL;

    for (i = 0, j = 1; j < comp_len; i+=4) {
	if (s_comp[j] != -128) {
	    uncomp[i  ] = s_comp[j] < 0 ? -1 : 0;
	    uncomp[i+1] = s_comp[j] < 0 ? -1 : 0;
	    uncomp[i+2] = s_comp[j] < 0 ? -1 : 0;
	    uncomp[i+3] = s_comp[j++];
	} else {
	    j++;
	    uncomp[i  ] = s_comp[j++];
	    uncomp[i+1] = s_comp[j++];
	    uncomp[i+2] = s_comp[j++];
	    uncomp[i+3] = s_comp[j++];
	}
    }

    /* Reclaim unneeded memory */
    uncomp = xrealloc(uncomp, i);
    
    *uncomp_len = i;
    return uncomp;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_FOLLOW1
 * ---------------------------------------------------------------------------
 */
static int follow_tab[256][256];
char *follow1(char *x_uncomp,
	      int uncomp_len,
	      int *comp_len) {
    char *comp = (char *)xmalloc(uncomp_len + 256 + 1);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;
    signed char *s_uncomp = ((signed char *)x_uncomp);
    int i, j;
    char next[256];
    int count[256];

    if (!comp)
	return NULL;

    /* Count di-freqs */
    memset(follow_tab, 0, 256*256*sizeof(int));
#if 0
    for (i = 0; i < uncomp_len-1; i++)
	follow_tab[u_uncomp[i]][u_uncomp[i+1]]++;

    /* Pick the most frequent next byte from the preceeding byte */
    for (i = 0; i < 256; i++) {
	int bestval, bestind;

	bestval = bestind = 0;
	for (j = 0; j < 256; j++) {
	    if (follow_tab[i][j] > bestval) {
		bestval = follow_tab[i][j];
		bestind = j;
	    }
	}
	next[i] = bestind;
    }
#endif
    memset(next, 0, 256*sizeof(*next));
    memset(count, 0, 256*sizeof(*count));
    /* Pick the most frequent next byte from the preceeding byte */
    for (i = 0; i < uncomp_len-1; ) {
	int cur = u_uncomp[i];
	int nxt = u_uncomp[++i];
	int folcnt = ++follow_tab[cur][nxt];
	if (folcnt > count[cur]) {
	    count[cur] = folcnt;
	    next[cur] = nxt;
	}
    }

    j = 0;
    comp[j++] = ZTR_FORM_FOLLOW1;

    /* Output 'next' array */
    for (i = 0; i < 256; i++, j++)
	comp[j] = next[i];

    /* Output new 'uncomp' as next['uncomp'] */
    comp[j++] = u_uncomp[0];
    for (i = 1; i < uncomp_len; i++, j++) {
	comp[j] = next[u_uncomp[i-1]] - s_uncomp[i];
    }
    *comp_len = j;

    return comp;
}

char *unfollow1(char *x_comp,
		int comp_len,
		int *uncomp_len) {
    unsigned char *u_uncomp;
    int i, j;
    char next[256];
    unsigned char *u_comp = (unsigned char *)x_comp;
    signed   char *s_comp = (signed   char *)x_comp;

    u_uncomp = (unsigned char *)xmalloc(comp_len-256-1);
    if (!u_uncomp)
	return NULL;

    /* Load next[] array */
    j = 1;
    for (i = 0; i < 256; i++, j++)
	next[i] = u_comp[j];

    /* Replace comp[x] with next[comp[x-1]] - comp[x]*/
    u_uncomp[0] = u_comp[j++];

    comp_len -= 257;
    s_comp   += 257;
    for (i = 1; i < comp_len; i++) {
	u_uncomp[i] = next[u_uncomp[i-1]] - s_comp[i];
    }

    *uncomp_len = i;

    return (char *)u_uncomp;
}


/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_CHEB445
 * ---------------------------------------------------------------------------
 */

#if 0
/*
 * Compresses using chebyshev polynomials to predict the next peak.
 * Based on around 96 modern ABI files it compresses by around 5% (varied
 * between 3.9 and 6.6).
 */

/*
 * For now this has been disabled in favour of the integer version below
 * as we cannot guarantee all systems to have the same floating point
 * roundings, especially with the final conversion to integer.
 * (Also, for some unknown reason, the integer version produces smaller
 * files.)
 */

char *cheb445comp(char *uncomp, int uncomp_len, int *data_len)
{
    int i, k, l, z;
    int datap;
    float frac[5];
    float fz[25];
    signed short *d16 = (signed short *)uncomp;
    int nwords = uncomp_len / 2;
    short *dptr = d16;
    signed short *data = (signed short *)malloc((nwords+1)*sizeof(short));

    data[0] = le_int2(ZTR_FORM_CHEB445);
    /* Check for boundary cases */
    if (nwords <= 4) {
	switch (nwords) {
	case 4:
	    data[4] = be_int2(be_int2(d16[3])-be_int2(d16[2]));
	case 3:
	    data[3] = be_int2(be_int2(d16[2])-be_int2(d16[1]));
	case 2:
	    data[2] = be_int2(be_int2(d16[1])-be_int2(d16[0]));
	case 1:
	    data[1] = be_int2(d16[0]);
	}
	*data_len = nwords*2;
	return (char *)data;
    }

    /* First 4 values are just direct deltas */
    data[1] = be_int2(d16[0]);
    data[2] = be_int2(be_int2(d16[1])-be_int2(d16[0]));
    data[3] = be_int2(be_int2(d16[2])-be_int2(d16[1]));
    data[4] = be_int2(be_int2(d16[3])-be_int2(d16[2]));
    datap = 5;

    /* Initialise - speeds up loop */
    for (k = 0; k < 5; k++) {
	float kx = cos(M_PI*(k+0.5)/5)*1.5;
	frac[k] = (kx + 1.5) - (int)(kx + 1.5);
    }
    for (z = l = 0; l < 5; l++) {
	for (k = 0; k < 5; k++, z++) {
	    fz[z] = 0.4 * cos(l * M_PI*(k+0.5)/5);
	}
    }

    /* Loop */
    for (i = 0; i < nwords-4; i++) {
	float dd, y = 10/3.0;
	float f[5], coef[5];
	signed short diff;
	int p;
	    
	f[0] = be_int2(dptr[2])*frac[4] + be_int2(dptr[3])*frac[0];
	f[1] = be_int2(dptr[2])*frac[3] + be_int2(dptr[3])*frac[1];
	f[2] = be_int2(dptr[1])*frac[2] + be_int2(dptr[2])*frac[2];
	f[3] = be_int2(dptr[0])*frac[1] + be_int2(dptr[1])*frac[3];
	f[4] = be_int2(dptr[0])*frac[0] + be_int2(dptr[1])*frac[4];
    
	for (z = l = 0; l < 5; l++, z+=5)
	    coef[l] = f[0] * fz[z+0] +
		      f[1] * fz[z+1] +
		      f[2] * fz[z+2] +
		      f[3] * fz[z+3] +
		      f[4] * fz[z+4];

	dd = y*coef[3]+coef[2];
	p = 0.5 + 5/3.0*(y*dd-coef[3]+coef[1])-dd+coef[0]/2.0;

	if (p < 0) p = 0;
	diff = be_int2(dptr[4]) - p;
	data[datap++] = be_int2(diff);

	dptr++;
    }

    *data_len = datap*2;

    return (char *)data;
}

char *cheb445uncomp(char *comp, int comp_len, int *uncomp_len)
{
    int i, k, l, z;
    float frac[5];
    float fz[25];
    signed short *d16 = (signed short *)comp;
    int nwords = comp_len / 2;
    signed short *data = (signed short *)xmalloc(comp_len);
    short *dptr = data, *dptr2 = d16;

    /* Check for boundary cases */
    if (nwords <= 3) {
	switch (nwords) {
	case 3:
	    data[0] = be_int2(d16[1]);
	    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
	    data[2] = be_int2(be_int2(d16[3])+be_int2(data[1]));
	    break;
	case 2:
	    data[0] = be_int2(d16[1]);
	    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
	    break;
	case 1:
	    data[0] = be_int2(d16[1]);
	    break;
	}
	*uncomp_len = (nwords-1)*2;
	return (char *)data;
    }

    /* First 3 values are just direct deltas */
    data[0] = be_int2(d16[1]);
    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
    data[2] = be_int2(be_int2(d16[3])+be_int2(data[1]));
    data[3] = be_int2(be_int2(d16[4])+be_int2(data[2]));
    dptr2 += 5;

    /* Initialise - speeds up loop */
    for (k = 0; k < 5; k++) {
	float kx = cos(M_PI*(k+0.5)/5)*1.5;
	frac[k] = (kx + 1.5) - (int)(kx + 1.5);
    }
    for (z = l = 0; l < 5; l++) {
	for (k = 0; k < 5; k++, z++) {
	    fz[z] = 0.4 * cos(l * M_PI*(k+0.5)/5);
	}
    }

    /* Loop */
    for (i = 0; i < nwords-3; i++) {
	float dd, y = 10/3.0;
	float f[5], coef[5];
	signed short diff;
	int p;
	    
	f[0] = be_int2(dptr[2])*frac[4] + be_int2(dptr[3])*frac[0];
	f[1] = be_int2(dptr[2])*frac[3] + be_int2(dptr[3])*frac[1];
	f[2] = be_int2(dptr[1])*frac[2] + be_int2(dptr[2])*frac[2];
	f[3] = be_int2(dptr[0])*frac[1] + be_int2(dptr[1])*frac[3];
	f[4] = be_int2(dptr[0])*frac[0] + be_int2(dptr[1])*frac[4];
    
	for (z = l = 0; l < 5; l++, z+=5)
	    coef[l] = f[0] * fz[z+0] +
		      f[1] * fz[z+1] +
		      f[2] * fz[z+2] +
		      f[3] * fz[z+3] +
		      f[4] * fz[z+4];

	dd = y*coef[3]+coef[2];
	p = 0.5 + 5/3.0*(y*dd-coef[3]+coef[1])-dd+coef[0]/2.0;

	if (p < 0) p = 0;
	diff = be_int2(*dptr2) + p;
	dptr[4] = be_int2(diff);

	dptr++;
	dptr2++;
    }

    *uncomp_len = (nwords-1)*2;
    return (char *)data;
}
#endif

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_ICHEB
 * ---------------------------------------------------------------------------
 */

/*
 * Integer versions of the chebyshev polynomial compressor. This uses
 * the polynomials to predict the next peak from the preceeding 3.
 * Tested on 100 ABI-3700, Megabace and Licor files it compressed by
 * around 7-8%. (Oddly this is slightly more than the floating point
 * version.)
 *
 * These require 32-bit integers and have code to make sure that arithmetic
 * does not overflow this.
 */
#define CH1 150
#define CH2 105
char *ichebcomp(char *uncomp, int uncomp_len, int *data_len)
{
    int i, l, z;
    int datap;
    int frac[5] = {139,57,75,93,11};
    int fz[20] = {42, 42, 42, 42, 42,
		  39, 24,  0,-24,-39,
		  33,-12,-42,-12, 33,
		  24,-39,  0, 39,-24};
    int dfac;
    signed short *d16 = (signed short *)uncomp;
    int nwords = uncomp_len / 2;
    short *dptr = d16;
    signed short *data = (signed short *)malloc((nwords+1)*sizeof(short));

    data[0] = le_int2(ZTR_FORM_ICHEB);
    /* Check for boundary cases */
    if (nwords <= 4) {
	switch (nwords) {
	case 4:
	    data[4] = be_int2(be_int2(d16[3])-be_int2(d16[2]));
	case 3:
	    data[3] = be_int2(be_int2(d16[2])-be_int2(d16[1]));
	case 2:
	    data[2] = be_int2(be_int2(d16[1])-be_int2(d16[0]));
	case 1:
	    data[1] = be_int2(d16[0]);
	}
	*data_len = nwords*2;
	return (char *)data;
    }

    /* First 4 values are just direct deltas */
    data[1] = be_int2(d16[0]);
    data[2] = be_int2(be_int2(d16[1])-be_int2(d16[0]));
    data[3] = be_int2(be_int2(d16[2])-be_int2(d16[1]));
    data[4] = be_int2(be_int2(d16[3])-be_int2(d16[2]));
    datap = 5;

    /* Loop */
    for (i = 4; i < nwords; i++) {
	int dd, f[5];
	signed int coef[4];
	signed short diff;
	int p;

	/*
	 * FIXME: As an alternative to the range checking below, if we
	 * scale dptr[X] such that it's never higher than 2800 then
	 * the 32-bit arithmetic will never overflow. Practically speaking,
	 * all observed ABI and Megabace files have vales up to 1600 only.
	 */

	/*
	 * frac[N] is always paired with frac[4-N], summing to 1.0 - or
	 * 150 when scaled.
	 * Hence f[0] has range 0 to 65536*150.
	 */
	f[0] = ((unsigned short)be_int2(dptr[2]))*frac[4] +
	       ((unsigned short)be_int2(dptr[3]))*frac[0];
	f[1] = ((unsigned short)be_int2(dptr[2]))*frac[3] +
	       ((unsigned short)be_int2(dptr[3]))*frac[1];
	f[2] = ((unsigned short)be_int2(dptr[1]))*frac[2] +
	       ((unsigned short)be_int2(dptr[2]))*frac[2];
	f[3] = ((unsigned short)be_int2(dptr[0]))*frac[1] +
	       ((unsigned short)be_int2(dptr[1]))*frac[3];
	f[4] = ((unsigned short)be_int2(dptr[0]))*frac[0] +
	       ((unsigned short)be_int2(dptr[1]))*frac[4];

	/*
	 * fz[z+0..5] sums to no more than 210 (5*42) and no less than 0.
	 * Therefore coef[l] has range 0 to 65536*150*210, which (just)
	 * fits in 31-bits, plus 1 for the sign.
	 */
	for (z = l = 0; l < 4; l++, z+=5)
	    coef[l] = (f[0] * fz[z+0] +
		       f[1] * fz[z+1] +
		       f[2] * fz[z+2] +
		       f[3] * fz[z+3] +
		       f[4] * fz[z+4]);

	/*
	 * computing p requires at most a temporary variable of 
	 * 24.1 * coef, but coef may be a full 32-bit integer.
	 * If coef is sufficiently close to cause an integer overflow then
	 * we scale it down.
	 */
	{
	    int max = 0;

	    for (l = 0; l < 4; l++) {
		if (max < ABS(coef[l]))
		    max = ABS(coef[l]);
	    }

    
	    if (max > 1<<26) {
		dfac = max / (1<<26) + 1;
		for (l = 0; l < 4; l++)
		    coef[l] /= dfac;
	    } else {
		dfac = 1;
	    }
	}

	dd = (coef[3]/3)*10+coef[2];
	p = ((((dd/3)*10-coef[3]+coef[1])/3)*5-dd+coef[0]/2)/(CH1*CH2);
	p *= dfac;

	if (p < 0) p = 0;
	diff = be_int2(dptr[4]) - p;
	data[datap++] = be_int2(diff);

	dptr++;
    }

    *data_len = datap*2;

    return (char *)data;
}

char *ichebuncomp(char *comp, int comp_len, int *uncomp_len)
{
    int i, l, z;
    int frac[5] = {139,57,75,93,11};
    int fz[20] = {42, 42, 42, 42, 42,
		  39, 24,  0,-24,-39,
		  33,-12,-42,-12, 33,
		  24,-39,  0, 39,-24};
    signed short *d16 = (signed short *)comp;
    int nwords = comp_len / 2 - 1;
    signed short *data = (signed short *)xmalloc(comp_len);
    short *dptr = data, *dptr2 = d16;
    int dfac;

    /* Check for boundary cases */
    if (nwords <= 4) {
	switch (nwords) {
	case 4:
	    data[0] = be_int2(d16[1]);
	    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
	    data[2] = be_int2(be_int2(d16[3])+be_int2(data[1]));
	    data[3] = be_int2(be_int2(d16[4])+be_int2(data[2]));
	    break;
	case 3:
	    data[0] = be_int2(d16[1]);
	    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
	    data[2] = be_int2(be_int2(d16[3])+be_int2(data[1]));
	    break;
	case 2:
	    data[0] = be_int2(d16[1]);
	    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
	    break;
	case 1:
	    data[0] = be_int2(d16[1]);
	    break;
	}
	*uncomp_len = nwords*2;
	return (char *)data;
    }

    /* First 4 values are just direct deltas */
    data[0] = be_int2(d16[1]);
    data[1] = be_int2(be_int2(d16[2])+be_int2(data[0]));
    data[2] = be_int2(be_int2(d16[3])+be_int2(data[1]));
    data[3] = be_int2(be_int2(d16[4])+be_int2(data[2]));
    dptr2 += 5;

    /* Loop */
    for (i = 4; i < nwords; i++) {
	int dd, coef[5], f[5];
	signed short diff;
	int p;
	    
	f[0] = ((unsigned short)be_int2(dptr[2]))*frac[4] +
	       ((unsigned short)be_int2(dptr[3]))*frac[0];
	f[1] = ((unsigned short)be_int2(dptr[2]))*frac[3] +
	       ((unsigned short)be_int2(dptr[3]))*frac[1];
	f[2] = ((unsigned short)be_int2(dptr[1]))*frac[2] +
	       ((unsigned short)be_int2(dptr[2]))*frac[2];
	f[3] = ((unsigned short)be_int2(dptr[0]))*frac[1] +
	       ((unsigned short)be_int2(dptr[1]))*frac[3];
	f[4] = ((unsigned short)be_int2(dptr[0]))*frac[0] +
	       ((unsigned short)be_int2(dptr[1]))*frac[4];
    
	for (z = l = 0; l < 4; l++, z+=5)
	    coef[l] = f[0] * fz[z+0] +
		      f[1] * fz[z+1] +
		      f[2] * fz[z+2] +
		      f[3] * fz[z+3] +
		      f[4] * fz[z+4];

	/*
	 * computing p requires at most a temporary variable of 
	 * 24.1 * coef, but coef may be a full 32-bit integer.
	 * If coef is sufficiently close to cause an integer overflow then
	 * we scale it down.
	 */
	{
	    int max = 0;

	    for (l = 0; l < 4; l++) {
		if (max < ABS(coef[l]))
		    max = ABS(coef[l]);
	    }

    
	    if (max > 1<<26) {
		dfac = max / (1<<26) + 1;
		for (l = 0; l < 4; l++)
		    coef[l] /= dfac;
	    } else {
		dfac = 1;
	    }
	}

	dd = (coef[3]/3)*10+coef[2];
	p = ((((dd/3)*10-coef[3]+coef[1])/3)*5-dd+coef[0]/2)/(CH1*CH2);
	p *= dfac;

	if (p < 0) p = 0;
	diff = be_int2(*dptr2) + p;
	dptr[4] = be_int2(diff);

	dptr++;
	dptr2++;
    }

    *uncomp_len = nwords*2;
    return (char *)data;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_LOG
 * ---------------------------------------------------------------------------
 */

/*
 * This is a LOSSY compression. It replaces N with 10 * log2(N).
 * (Just an idea, and not a great one it seems.)
 */
char *log2_data(char *x_uncomp,
		int uncomp_len,
		int *comp_len) {
    int i, u1, l1;
    char *comp = (char *)xmalloc(uncomp_len + 2);
    unsigned char *u_uncomp = (unsigned char *)x_uncomp;
    
    if (!comp)
	return NULL;

    comp+=2;
    for (i = 0; i < uncomp_len; i+=2) {
	u1 =(u_uncomp[i  ] <<  8) +
	    (u_uncomp[i+1] <<  0);
	l1 = (int)(10 * log(u1+1) / log(2));
	comp[i  ] = (l1 >> 8) & 0xff;
	comp[i+1] = (l1 >> 0) & 0xff;
    }

    comp-=2;
    comp[0] = ZTR_FORM_LOG2;
    comp[1] = 0; /* dummy - to align on 2-byte boundary */

    *comp_len = uncomp_len+2;

    return comp;
}

char *unlog2_data(char *x_comp,
		  int comp_len,
		  int *uncomp_len) {
    int i, u1, l1;
    char *uncomp;
    unsigned char *u_comp = (unsigned char *)x_comp;

    uncomp = (char *)xmalloc(comp_len-2);
    if (!uncomp)
	return NULL;

    u_comp+=2;
    comp_len-=2;
    *uncomp_len = comp_len;

    for (i = 0; i < comp_len; i+=2) {
	l1 = ((u_comp[i  ] << 8) |
	      (u_comp[i+1] << 0));
	u1 = (int)pow(2.0, l1/10.0)-1;
	uncomp[i  ] = (u1 >> 8) & 0xff;
	uncomp[i+1] = (u1 >> 0) & 0xff;
    }

    return uncomp;
}

/*
 * ---------------------------------------------------------------------------
 * ZTR_FORM_STHUFF
 * ---------------------------------------------------------------------------
 */

/*
 * Implements compression using a set of static huffman codes stored using
 * the Deflate algorithm (and so in this respect it's similar to zlib).
 *
 * The huffman codes though can be previously stored in the ztr object
 * using ztr_add_hcode(). "cset" indicates which numbered stored huffman
 * code set is to be used, or passing zero will use inline codes (ie they
 * are stored in the data stream itself, just as in standard deflate).
 *
 * Arguments:
 *	ztr		ztr_t pointer; used to find stored code-sets
 *	uncomp		The uncompressed input data
 *	uncomp_len	Length of uncomp
 *	cset		Stored code-set number, zero for inline
 *	recsz		Record size - only used when cset == 0.
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data stream if successful + comp_len
 *      NULL on failure
 */
char *sthuff(ztr_t *ztr, char *uncomp, int uncomp_len, 
	     int cset, int recsz, int *comp_len) {
    block_t *blk = block_create(NULL, 2);
    unsigned char bytes[2];
    huffman_codeset_t *c = NULL;
    unsigned char *comp = NULL;
    ztr_hcode_t *hc = NULL;

    if (cset >= CODE_USER) {
	if (NULL == (hc = ztr_find_hcode(ztr, cset)))
	    return NULL;
	c = hc->codes;
    } else if (cset != CODE_INLINE) {
	c = generate_code_set(cset, 1, NULL, 0, 1, MAX_CODE_LEN, 0);
    }

    if (!c) {
	/* No cached ones found, so inline some instead */
	cset = 0;
	c = generate_code_set(0, recsz, (unsigned char *)uncomp, uncomp_len, 1,
			      MAX_CODE_LEN, 0);
    }

    bytes[0] = ZTR_FORM_STHUFF;
    bytes[1] = cset;
    store_bytes(blk, bytes, 2);

    if (hc) {
	if (!c->blk) {
	    c->blk = block_create(NULL, 2);
	    store_codes(c->blk, c, 1);
	}
	blk->bit = c->blk->bit;
    } else {
	store_codes(blk, c, 1);
    }

    /*
    {int i;
    for (i = 0; i < c->ncodes; i++) {
	output_code_set(stderr, c->codes[i]);
    }}
    */


    /*
     * Unless CODE_INLINE, all we wanted to know is what bit number
     * to start on. The above is therefore somewhat inefficient.
     */
    if (cset != 0) {
	blk->byte = 2;
	memset(&blk->data[2], 0, blk->alloc - 2);
    }

    if (0 == huffman_multi_encode(blk, c, cset,
			    (unsigned char *)uncomp, uncomp_len)) {
	comp = blk->data;
	*comp_len = blk->byte + (blk->bit != 0);
	block_destroy(blk, 1);
    }

    if (cset == 0)
	huffman_codeset_destroy(c);

    return (char *)comp;
}

char *unsthuff(ztr_t *ztr, char *comp, int comp_len, int *uncomp_len) {
    int cset = (unsigned char)(comp[1]);
    huffman_codeset_t *cs = NULL, *cs_free = NULL;
    block_t *blk_in = block_create(NULL, comp_len),
	*blk_out = block_create(NULL, 1000);
    int bfinal = 1, bit_num = 0;
    char *uncomp;

    if (cset >= CODE_USER) {
	/* Scans through HUFF chunks */
	ztr_hcode_t *hc = ztr_find_hcode(ztr, cset);
	if (!hc)
	    return NULL;

	cs = hc->codes;
	bit_num = cs->bit_num;
	blk_in->bit = 0;
    } else if (cset > 0) {
	/* Create some temporary huffman_codes to stringify */
	cs_free = cs = generate_code_set(cset, 1, NULL, 0, 1, MAX_CODE_LEN, 0);
	if (!cs)
	    return NULL;

	bit_num = cs->bit_num;
	blk_in->bit = 0;
    } /* else inline codes */

    /*
     * We need to know at what bit the huffman codes would have ended on
     * so we can store our huffman encoded symbols immediately following it.
     * For speed though this bit-number is cached.
     */
    blk_in->data[blk_in->byte++] |= *(comp+2);
    store_bytes(blk_in, (unsigned char *)comp+3, comp_len-3);


    /* Rewind */
    blk_in->byte = 0;
    blk_in->bit = bit_num;

    do {
	block_t *out;

	/*
	 * We're either at the start of a block with codes to restore
	 * (cset == INLINE or the 2nd onwards block) or we've already
	 * got some codes in cs and we're at the position where huffman
	 * encoded symbols are stored. 
	 */
	if (!cs)
	    if (NULL == (cs = cs_free = restore_codes(blk_in, &bfinal)))
		return NULL;

	/*
	{int i;
	for (i = 0; i < cs->ncodes; i++) {
	    output_code_set(stderr, cs->codes[i]);
	}}
	*/

	if (NULL == (out = huffman_multi_decode(blk_in, cs))) {
	    huffman_codeset_destroy(cs);
	    return NULL;
	}

	/* Could optimise this for the common case of only 1 block */
	store_bytes(blk_out, out->data, out->byte);
	block_destroy(out, 0);
	if (cs_free)
	    huffman_codeset_destroy(cs_free);
	cs = cs_free = NULL;
    } while (!bfinal);

    *uncomp_len = blk_out->byte;
    uncomp = (char *)blk_out->data;

    block_destroy(blk_in, 0);
    block_destroy(blk_out, 1);

    return uncomp;
}


#ifndef NDEBUG
#define SYM_EOF 256
static void output_code_set(FILE *fp, huffman_codes_t *cds) {
    int i, j;
    int nbits_in = 0, nbits_out = 0;
    huffman_code_t *codes = cds->codes;
    int ncodes = cds->ncodes;

    fprintf(fp, "static huffman_code_t codes_FIXME[] = {\n");
    for (i = j = 0; i < ncodes; i++) {
	nbits_out += codes[i].nbits * codes[i].freq;
	nbits_in  += 8*codes[i].freq;
	if (j == 0)
	    fprintf(fp, "    ");
	if (codes[i].symbol == SYM_EOF) {
	    fprintf(fp, "{SYM_EOF,%3d}, ", codes[i].nbits);
	    j = 10;
	} else {
	    if (isalnum(codes[i].symbol)) {
		fprintf(fp, "{'%c',%3d}, ", codes[i].symbol, codes[i].nbits);
	    } else {
		fprintf(fp, "{%3d,%3d}, ", codes[i].symbol, codes[i].nbits);
	    }
	}
	j++;
	
	if (j >= 6) {
	    fputc('\n', fp);
	    j = 0;
	}
    }
    if (j)
	fputc('\n', fp);
    fprintf(fp, "};\n");
    fprintf(fp, "/* Expected compression to %f of input */\n",
	    (double)nbits_out/nbits_in);
}
#endif

/*
 * Reorders quality data from its RAW format to an interleaved 4-byte
 * aligned format.
 *
 * Starting with sequence A1 C2 G3 the raw format is quality of called
 * bases followed by quality of remaining bases in triplets per base:
 * 0 (RAW format)
 * Q(A1) Q(C2) Q(G3)
 * Q(C1) Q(G1) Q(T1) 
 * Q(A2) Q(G2) Q(T2) 
 * Q(A3) Q(C3) Q(T3) 
 *
 * We reorder it to:
 * ZTR_FORM_QSHIFT <any> <any> 0(raw)
 * Q(A1) Q(C1) Q(G1) Q(T1)
 * Q(C2) Q(A2) Q(G2) Q(T2)
 * Q(G3) Q(A3) Q(C3) Q(T3)
 * 
 * Returns shifted data on success
 *         NULL on failure
 */
char *qshift(char *qold, int qlen, int *new_len) {
    int i, j, k;
    char *qnew;
    int nbases;

    /*
     * Correct input is raw encoding + 4x nbases bytes
     */
    if ((qlen-1)%4 != 0 || *qold != 0)
	return NULL;

    nbases = (qlen-1)/4;
    qnew = (char *)malloc((nbases+1)*4);
    qnew[0] = ZTR_FORM_QSHIFT; /* reorder code */
    qnew[1] = -40;  /* pad */
    qnew[2] = -40;  /* pad */
    qnew[3] = qold[0];

    for (i = 0, j = 4, k = nbases; i < nbases; i++, j+=4, k+=3) {
	qnew[j  ] = qold[1+i  ];
	qnew[j+1] = qold[1+k  ];
	qnew[j+2] = qold[1+k+1];
	qnew[j+3] = qold[1+k+2];
    }

    *new_len = (nbases+1)*4;
    return qnew;
}

/*
 * The opposite transform from qshift() above.
 *
 * Returns unshifted data on success
 *         NULL on failure.
 */
char *unqshift(char *qold, int qlen, int *new_len) {
    int i, j, k;
    char *qnew;
    int nbases;

    /*
     * Correct input is 4x (nbases+1) bytes
     */
    if (qlen%4 != 0 || *qold != ZTR_FORM_QSHIFT)
	return NULL;

    nbases = qlen/4-1;
    qnew = (char *)malloc(nbases*4+1);
    qnew[0] = 0; /* raw byte */

    for (i = 0, j = 4, k = nbases; i < nbases; i++, j+=4, k+=3) {
	qnew[1+i  ] = qold[j  ];
	qnew[1+k  ] = qold[j+1];
	qnew[1+k+1] = qold[j+2];
	qnew[1+k+2] = qold[j+3];
    }

    *new_len = nbases*4+1;
    return qnew;
}

/*
 * Given a sequence ACTG this shifts trace data from the order:
 *
 *     A1A2A3A4 C1C2C3C4 G1G2G3G4 T1T2T3T4
 *
 * to
 *
 *     A1C1G1T1 C2A2G2T2 T3A3C3G3 G4C4C4T4
 *
 * Ie for each base it ouputs the signal for the called base first
 * followed by the remaining 3 signals in A,C,G,T order (minus the
 * called signal already output).
 */

/*
 * NCBI uses a rotation mechanism. Thus instead of converting acGt (G
 * called) to Gact they produce Gtac.
 *
 * Given 4 columns of data the two schemes produce value orderings as:
 *
 * Rotate: Acgt       Shift: Acgt
 *         Cgta              Cagt
 *         Gtac              Gact
 *         Tacg              Tacg
 *
 * As a consequence any channel bias for a/c/g/t is better preserved
 * by shifting than it is by rotating as each column (except #1) are
 * only populated by two base types and 2 of those columns are on a
 * 3:1 ratio.
 */
/* #define ROTATE_INSTEAD */
char *tshift(ztr_t *ztr, char *told_c, int tlen, int *new_len) {
    int nc, i;
    ztr_chunk_t **base = ztr_find_chunks(ztr, ZTR_TYPE_BASE, &nc);
    char *bases;
    int nbases;
    unsigned short *tnew, *told;

    if (nc == 0)
	return NULL;

    if (*told_c != 0)
	return NULL; /* assume RAW format trace input */

    /* Use last BASE chunk if multiple are present */
    /* FIXME: ensure uncompressed first */
    bases  = base[nc-1]->data+1;
    nbases = base[nc-1]->dlength-1;

    if (nbases != (tlen-2)/8) {
	fprintf(stderr, "Mismatch in number of base calls to samples\n");
	return NULL;
    }

    /* Allocate and initialise header */
    told = ((unsigned short *)told_c) + 1;
    *new_len = (nbases*4+4) * sizeof(*tnew);
    tnew = (unsigned short *)malloc(*new_len);
    for (i = 0; i < 4; i++) {
	tnew[i] = 0;
    }
    ((char *)tnew)[0] = ZTR_FORM_TSHIFT;

#ifdef ROTATE_INSTEAD
    /* Reorder */
    for (i = 0; i < nbases; i++) {
	switch(bases[i]) {
	case 'T':
	    tnew[4+i*4+0] = told[3*nbases+i]; /* TACG */
	    tnew[4+i*4+1] = told[0*nbases+i];
	    tnew[4+i*4+2] = told[1*nbases+i];
	    tnew[4+i*4+3] = told[2*nbases+i];
	    break;
	case 'G':
	    tnew[4+i*4+0] = told[2*nbases+i]; /* GTAC */
	    tnew[4+i*4+1] = told[3*nbases+i];
	    tnew[4+i*4+2] = told[0*nbases+i];
	    tnew[4+i*4+3] = told[1*nbases+i];
	    break;
	case 'C':
	    tnew[4+i*4+0] = told[1*nbases+i]; /* CGTA */
	    tnew[4+i*4+1] = told[2*nbases+i];
	    tnew[4+i*4+2] = told[3*nbases+i];
	    tnew[4+i*4+3] = told[0*nbases+i];
	    break;
	default:
	    tnew[4+i*4+0] = told[0*nbases+i]; /* ACGT */
	    tnew[4+i*4+1] = told[1*nbases+i];
	    tnew[4+i*4+2] = told[2*nbases+i];
	    tnew[4+i*4+3] = told[3*nbases+i];
	    break;
	}
    }
#else
    /* Reorder */
    for (i = 0; i < nbases; i++) {
	switch(bases[i]) {
	case 'T':
	    tnew[4+i*4+0] = told[3*nbases+i]; /* TACG */
	    tnew[4+i*4+1] = told[0*nbases+i];
	    tnew[4+i*4+2] = told[1*nbases+i];
	    tnew[4+i*4+3] = told[2*nbases+i];
	    break;
	case 'G':
	    tnew[4+i*4+0] = told[2*nbases+i]; /* GACT */
	    tnew[4+i*4+1] = told[0*nbases+i];
	    tnew[4+i*4+2] = told[1*nbases+i];
	    tnew[4+i*4+3] = told[3*nbases+i];
	    break;
	case 'C':
	    tnew[4+i*4+0] = told[1*nbases+i]; /* CAGT */
	    tnew[4+i*4+1] = told[0*nbases+i];
	    tnew[4+i*4+2] = told[2*nbases+i];
	    tnew[4+i*4+3] = told[3*nbases+i];
	    break;
	default:
	    tnew[4+i*4+0] = told[0*nbases+i]; /* ACGT */
	    tnew[4+i*4+1] = told[1*nbases+i];
	    tnew[4+i*4+2] = told[2*nbases+i];
	    tnew[4+i*4+3] = told[3*nbases+i];
	    break;
	}
    }
#endif

    xfree(base);

    return (char *)tnew;
}

char *untshift(ztr_t *ztr, char *told_c, int tlen, int *new_len) {
    unsigned short *tnew, *told = (unsigned short *)told_c;
    int nc, nbases, i;
    char *bases;
    ztr_chunk_t **base = ztr_find_chunks(ztr, ZTR_TYPE_BASE, &nc);

    if (nc == 0)
	return NULL;

    /* Use last BASE chunk if multiple are present */
    uncompress_chunk(ztr, base[nc-1]);
    bases  = base[nc-1]->data+1;
    nbases = base[nc-1]->dlength-1;

    *new_len = 2 + nbases*4 * sizeof(*tnew);
    tnew = (unsigned short *)malloc(*new_len);
    tnew[0] = 0;
    told += 4;

#ifdef ROTATE_INSTEAD
    /* Reorder */
    for (i = 0; i < nbases; i++) {
	switch(bases[i]) {
	case 'T':
	    tnew[1+3*nbases+i] = *told++; /* TACG */
	    tnew[1+0*nbases+i] = *told++;
	    tnew[1+1*nbases+i] = *told++;
	    tnew[1+2*nbases+i] = *told++;
	    break;
	case 'G':
	    tnew[1+2*nbases+i] = *told++; /* GTAC */
	    tnew[1+3*nbases+i] = *told++;
	    tnew[1+0*nbases+i] = *told++;
	    tnew[1+1*nbases+i] = *told++;
	    break;
	case 'C':
	    tnew[1+1*nbases+i] = *told++; /* CGTA */
	    tnew[1+2*nbases+i] = *told++;
	    tnew[1+3*nbases+i] = *told++;
	    tnew[1+0*nbases+i] = *told++;
	    break;
	default:
	    tnew[1+0*nbases+i] = *told++; /* ACGT */
	    tnew[1+1*nbases+i] = *told++;
	    tnew[1+2*nbases+i] = *told++;
	    tnew[1+3*nbases+i] = *told++;
	    break;
	}
    }
#else
    /* Reorder */
    for (i = 0; i < nbases; i++) {
	switch(bases[i]) {
	case 'T':
	    tnew[1+3*nbases+i] = *told++; /* TACG */
	    tnew[1+0*nbases+i] = *told++;
	    tnew[1+1*nbases+i] = *told++;
	    tnew[1+2*nbases+i] = *told++;
	    break;
	case 'G':
	    tnew[1+2*nbases+i] = *told++; /* GACT */
	    tnew[1+0*nbases+i] = *told++;
	    tnew[1+1*nbases+i] = *told++;
	    tnew[1+3*nbases+i] = *told++;
	    break;
	case 'C':
	    tnew[1+1*nbases+i] = *told++; /* CAGT */
	    tnew[1+0*nbases+i] = *told++;
	    tnew[1+2*nbases+i] = *told++;
	    tnew[1+3*nbases+i] = *told++;
	    break;
	default:
	    tnew[1+0*nbases+i] = *told++; /* ACGT */
	    tnew[1+1*nbases+i] = *told++;
	    tnew[1+2*nbases+i] = *told++;
	    tnew[1+3*nbases+i] = *told++;
	    break;
	}
    }
#endif

    xfree(base);

    return (char *)tnew;
}
