#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <io_lib/ztr.h>
#include <io_lib/compression.h>
#include <io_lib/xalloc.h>

static char *format2str(int format) {
    static char unk[100];

    switch (format) {
    case ZTR_FORM_RAW:     return "raw";
    case ZTR_FORM_RLE:     return "rle";
    case ZTR_FORM_XRLE:    return "xrle";
    case ZTR_FORM_XRLE2:   return "xrle2";
    case ZTR_FORM_ZLIB:    return "zlib";
    case ZTR_FORM_DELTA1:  return "delta1";
    case ZTR_FORM_DELTA2:  return "delta2";
    case ZTR_FORM_DELTA4:  return "delta4";
    case ZTR_FORM_DDELTA1: return "ddelta1";
    case ZTR_FORM_DDELTA2: return "ddelta2";
    case ZTR_FORM_DDELTA4: return "ddelta4";
    case ZTR_FORM_16TO8:   return "16to8";
    case ZTR_FORM_32TO8:   return "32to8";
    case ZTR_FORM_FOLLOW1: return "follow1";
    case ZTR_FORM_CHEB445: return "cheb445";
    case ZTR_FORM_ICHEB:   return "icheb";
    case ZTR_FORM_LOG2:    return "log2";
    case ZTR_FORM_STHUFF:  return "sthuff";
    case ZTR_FORM_QSHIFT:  return "qshift";
    case ZTR_FORM_TSHIFT:  return "tshift";
    }

    sprintf(unk, "?%d?\n", format);
    return unk;
}

/*
 * Shannon showed that for storage in base 'b' with alphabet symbols 'a' having
 * a probability of ocurring in any context of 'Pa' we should encode
 * symbol 'a' to have a storage width of -logb(Pa).
 *
 * Eg. b = 26, P(e) = .22. => width .4647277.
 *
 * We use this to calculate the entropy of a signal by summing over all letters
 * in the signal. In this case, our storage has base 256.
 */
#define EBASE 256
double entropy(unsigned char *data, int len) {
    double E[EBASE];
    double P[EBASE];
    double e;
    int i;
    
    for (i = 0; i < EBASE; i++)
        P[i] = 0;

    for (i = 0; i < len; i++)
        P[data[i]]++;

    for (i = 0; i < EBASE; i++) {
        if (P[i]) {
            P[i] /= len;
            E[i] = -(log(P[i])/log(EBASE));
        } else {
            E[i] = 0;
        }
    }

    for (e = i = 0; i < len; i++)
        e += E[data[i]];

    return e;
}

/* Debug version of the ztr.c uncompress_chunk function. */
static int explode_chunk(ztr_t *ztr, ztr_chunk_t *chunk) {
    char *new_data = NULL;
    int new_len;

    while (chunk->dlength > 0 && chunk->data[0] != ZTR_FORM_RAW) {
	double ent = entropy((unsigned char *)chunk->data, chunk->dlength);

	switch (chunk->data[0]) {
	case ZTR_FORM_RLE:
	    new_data = unrle(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_XRLE:
	    new_data = unxrle(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_XRLE2:
	    new_data = unxrle2(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_ZLIB:
	    new_data = zlib_dehuff(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_DELTA1:
	    new_data = recorrelate1(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_DELTA2:
	    new_data = recorrelate2(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_DELTA4:
	    new_data = recorrelate4(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_16TO8:
	    new_data = expand_8to16(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_32TO8:
	    new_data = expand_8to32(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_FOLLOW1:
	    new_data = unfollow1(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_ICHEB:
	    new_data = ichebuncomp(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_LOG2:
	    new_data = unlog2_data(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_STHUFF:
	    new_data = unsthuff(ztr, chunk->data, chunk->dlength, &new_len);
	    break;
	    
	case ZTR_FORM_QSHIFT:
	    new_data = unqshift(chunk->data, chunk->dlength, &new_len);
	    break;

	case ZTR_FORM_TSHIFT:
	    new_data = untshift(ztr, chunk->data, chunk->dlength, &new_len);
	    break;

	default:
	    fprintf(stderr, "Unknown encoding format %d\n", chunk->data[0]);
	    return -1;
	}
	    
	if (!new_data) {
	    fprintf(stderr, "Failed to decode chunk with format %s\n",
		    format2str(chunk->data[0]));
	    return -1;
	}

	printf("    format %8s => %6d to %6d, entropy %8.1f to %8.1f\n",
	       format2str(chunk->data[0]), chunk->dlength, new_len,
	       ent, entropy((unsigned char *)new_data, new_len));

	chunk->dlength = new_len;
	xfree(chunk->data);
	chunk->data = new_data;
    }

    return 0;
}

int main(int argc, char **argv) {
    ztr_t *ztr;
    mFILE *fp;
    int i;

    if (argc >= 2) {
	if (NULL == (fp = mfopen(argv[1], "rb"))) {
	    perror(argv[1]);
	    return 1;
	}
    } else {
	fp = mstdin();
    }

    if (NULL == (ztr = mfread_ztr(fp))) {
	perror("fread_ztr");
	return 1;
    }

    printf("Nchunks = %d\n", ztr->nchunks);
    for (i = 0; i < ztr->nchunks; i++) {
	char str[5];
	int complen;
	int rawlen;
	char *val;

	(void)ZTR_BE2STR(ztr->chunk[i].type, str);
	complen = ztr->chunk[i].dlength;
	val = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "TYPE");
	if (val)
	    printf("-- %s (%s) --\n", str, val);
	else
	    printf("-- %s --\n", str);
	explode_chunk(ztr, &ztr->chunk[i]);
	rawlen = ztr->chunk[i].dlength;
	printf("SUMMARY %s  mlen %3d, dlen %6d, rawlen %6d, ratio %f\n",
	       str, ztr->chunk[i].mdlength,
	       complen, rawlen, (double)complen/rawlen);
#if 0
	fflush(stdout);
	puts("\n========================================");
	write(1, ztr->chunk[i].data, ztr->chunk[i].dlength);
	puts("\n========================================");
#endif
    }

    delete_ztr(ztr);
    
    return 0;
}
