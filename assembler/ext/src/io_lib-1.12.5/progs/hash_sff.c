/*
 * This adds a hash table index (".hsh" v1.00 format) to an SFF archive.
 * It does this either inline on the file itself (provided it doesn't already
 * have an index) or by producing a new indexed SFF archive.
 *
 * It has been coded to require only the memory needed to store the index
 * and so does quite a lot of I/O but with minimised memory. For a 460,000
 * SFF archive it took about 22 seconds real time on a 1.7GHz P4 when copying
 * or 10 seconds when updating inline.
 */

/* ---------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <io_lib/hash_table.h>
#include <io_lib/sff.h>
#include <io_lib/os.h>
#include <io_lib/mFILE.h>

/*
 * Override the sff.c functions to use FILE pointers instead. This means
 * we don't have to load the entire archive into memory, which is optimal when
 * dealing with a single file (ie in sff/sff.c), but not when indexing it.
 *
 * Done with minimal error checking I'll admit...
 */
static sff_read_header *fread_sff_read_header(FILE *fp) {
    sff_read_header *h;
    unsigned char rhdr[16];

    if (16 != fread(rhdr, 1, 16, fp))
	return NULL;
    h = decode_sff_read_header(rhdr);

    if (h->name_len != fread(h->name, 1, h->name_len, fp))
	return free_sff_read_header(h), NULL;
    
    /* Pad to 8 chars */
    fseek(fp, (ftell(fp) + 7)& ~7, SEEK_SET);

    return h;
}

static sff_common_header *fread_sff_common_header(FILE *fp) {
    sff_common_header *h;
    unsigned char chdr[31];

    if (31 != fread(chdr, 1, 31, fp))
	return NULL;
    h = decode_sff_common_header(chdr);
    if (h->flow_len != fread(h->flow, 1, h->flow_len, fp))
	return free_sff_common_header(h), NULL;
    if (h->key_len != fread(h->key , 1, h->key_len,  fp))
	return free_sff_common_header(h), NULL;

    /* Pad to 8 chars */
    fseek(fp, (ftell(fp) + 7)& ~7, SEEK_SET);

    return h;
}

void usage(void) {
    fprintf(stderr, "Usage: hash_sff [-o outfile] [-t] sff_file ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    HashFile *hf;
    sff_common_header *ch;
    sff_read_header *rh;
    int i, dot, arg;
    char *sff;
    char hdr[31];
    uint64_t index_offset = 0;
    uint32_t index_size, index_skipped;
    FILE *fp, *fpout = NULL;
    int copy_archive = 1;
    

    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-o") == 0 && argc > 1) {
	    if (NULL == (fpout = fopen(argv[1], "wb+"))) {
		perror(argv[1]);
		return 1;
	    }
	    argv++;
	    argc--;

	} else if (strcmp(*argv, "-t") == 0) {
	    copy_archive = 0;

	} else if (**argv == '-') {
	    usage();
	}

    }

    if (argc < 1)
	usage();

    if (copy_archive == 0 && argc != 1) {
	fprintf(stderr, "-t option only supported with a single sff argument\n");
	return 1;
    }

    /* Create the hash table */
    hf = HashFileCreate(0, HASH_DYNAMIC_SIZE);
    hf->nheaders = 0;
    hf->headers = NULL;

    for (arg = 0; arg < argc; arg++) {
	/* open (and read) the entire sff file */
	sff = argv[arg];

	printf("Indexing %s:\n", sff);
	if (fpout) {
	    if (NULL == (fp = fopen(sff, "rb"))) {
		perror(sff);
		return 1;
	    }
	} else { 
	    if (NULL == (fp = fopen(sff, "rb+"))) {
		perror(sff);
		return 1;
	    }
	}

	/* Read the common header */
	ch = fread_sff_common_header(fp);

	if (ch->index_len && !fpout) {
	    fprintf(stderr, "Archive already contains index.\nReplacing the"
		    " index requires the \"-o outfile\" option.\n");
	    return 1;
	}

	/* Add the SFF common header as a hash file-header */
	hf->nheaders++;
	hf->headers = (HashFileSection *)realloc(hf->headers, hf->nheaders *
						 sizeof(*hf->headers));
	hf->headers[hf->nheaders-1].pos = 0;
	hf->headers[hf->nheaders-1].size = ch->header_len;
	hf->headers[hf->nheaders-1].cached_data = NULL;

	/* Read the index items, adding to the hash */
	index_skipped = 0;
	dot = 0;
	printf("                                                                       |\r|");
	for (i = 0; i < ch->nreads; i++) {
	    int dlen;
	    uint32_t offset;
	    HashData hd;
	    HashFileItem *hfi;
	    
	    if (i >= dot * (ch->nreads/69)) {
		putchar('.');
		fflush(stdout);
		dot++;
	    }

	    /* Skip old index if present */
	    offset = ftell(fp);
	    if (offset == ch->index_offset) {
		fseek(fp, ch->index_len, SEEK_CUR);
		index_skipped = ch->index_len;
		continue;
	    }

	    hfi = (HashFileItem *)calloc(1, sizeof(*hfi));
	    rh = fread_sff_read_header(fp);
	    dlen = (2*ch->flow_len + 3*rh->nbases + 7) & ~7;
	    fseek(fp, dlen, SEEK_CUR);
	
	    hfi->header = hf->nheaders;
	    hfi->footer = 0;
	    hfi->pos = offset - index_skipped;
	    hfi->size = (ftell(fp) - index_skipped) - hfi->pos;
	    hd.p = hfi;

	    HashTableAdd(hf->h, rh->name, rh->name_len, hd, NULL);
	}
	printf("\n");
	HashTableStats(hf->h, stdout);

	index_offset = ftell(fp) - index_skipped;

	/* Copy the archive if needed, minus the old index */
	if (fpout && copy_archive) {
	    char block[8192];
	    size_t len;
	    uint64_t pos = 0;

	    printf("\nCopying archive\n");

	    fseek(fp, 0, SEEK_SET);
	    while (len = fread(block, 1, 8192, fp)) {
		/* Skip previous index */
		if (pos < ch->index_offset && pos+len > ch->index_offset) {
		    len = ch->index_offset - pos;
		    fseek(fp, ch->index_offset + ch->index_len, SEEK_SET);
		}
		if (len && len != fwrite(block, 1, len, fpout)) {
		    fprintf(stderr, "Failed to output new archive\n");
		    return 1;
		}
		pos += len;
	    }
	}
	
	if (!fpout) {
	    /* Save the hash */
	    printf("Saving index\n");
	    fseek(fp, 0, SEEK_END);
	    index_size = HashFileSave(hf, fp, 0);
	    HashFileDestroy(hf);

	    /* Update the common header */
	    fseek(fp, 0, SEEK_SET);
	    fread(hdr, 1, 31, fp);
	    *(uint64_t *)(hdr+8)  = be_int8(index_offset);
	    *(uint32_t *)(hdr+16) = be_int4(index_size);
	    fseek(fp, 0, SEEK_SET);
	    fwrite(hdr, 1, 31, fp);
	}

	fclose(fp);
    }

    if (fpout) {
	/* Save the hash */
	printf("Saving index\n");

	if (!copy_archive) {
	    hf->narchives = 1;
	    hf->archives = malloc(sizeof(char **));
	    hf->archives[0] = strdup(argv[0]);
	    index_offset = 0;
	}

	fseek(fpout, 0, SEEK_END);
	index_size = HashFileSave(hf, fpout, 0);
	HashFileDestroy(hf);

	/* Update the common header to indicate index location */
	if (copy_archive) {
	    fseek(fpout, 0, SEEK_SET);
	    fread(hdr, 1, 31, fpout);
	    *(uint64_t *)(hdr+8)  = be_int8(index_offset);
	    *(uint32_t *)(hdr+16) = be_int4(index_size);
	    fseek(fpout, 0, SEEK_SET);
	    fwrite(hdr, 1, 31, fpout);
	}
	fclose(fpout);
    }
    
    return 0;
}
