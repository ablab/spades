#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <inttypes.h>

#include <io_lib/os.h>
#include <io_lib/hash_table.h>
#include <io_lib/srf.h>

/* Command line options */
typedef struct {
    int long_format;
    int count_only;
    int verbose;
} opts;

/*
 * Lists the contents of an SRF file.
 *
 * Returns num_reads for success
 *        -1 for failure
 */
int64_t list_file(char *fname, opts *opts) {
    srf_t *srf;
    char name[512];
    int64_t count = 0;
    int type;
    uint64_t pos;

    if (NULL == (srf = srf_open(fname, "r"))) {
	perror(fname);
	return -1;
    }

    /* Scan through file gathering the details to index in memory */
    while ((type = srf_next_block_details(srf, &pos, name)) >= 0) {
	if (type == SRFB_TRACE_BODY) {
	    count++;
	    if (!opts->count_only) {
		if (opts->long_format)
		    printf("%-30s %10"PRId64" + %4d + %5d\n",
			   name, pos,
			   srf->tb.trace_size,
			   srf->th.trace_hdr_size);
		else
		    puts(name);
	    }
	}
    }

    srf_destroy(srf, 1);

    return count;
}

/*
 * Counts the contents of an SRF file.
 * If the hash index exists it uses this instead.
 *
 * Returns num_reads for success
 *        -1 for failure
 */
int64_t count_file(char *fname, opts *opts) {
    srf_t *srf;
    srf_index_hdr_t hdr;
    int64_t count = 0;
    off_t ipos, skip;
    int item_sz = 9, i;
    unsigned char data[16];

    if (NULL == (srf = srf_open(fname, "r"))) {
	perror(fname);
	return -1;
    }

    /* Read the index header */
    if (0 != srf_read_index_hdr(srf, &hdr, 0)) {
	srf_destroy(srf, 1);
	return list_file(fname, opts);
    }
    ipos = ftello(srf->fp);

    /* Compute the remaining size of the index and divide by item_sz */
    if (hdr.dbh_pos_stored_sep)
	item_sz += 4;
    skip = hdr.index_hdr_sz 
	 + hdr.n_container * 8
	 + hdr.n_data_block_hdr * 8
	 + hdr.n_buckets * 8;

    srf_destroy(srf, 1);

    return (hdr.size - skip - 16/* footer*/) / item_sz;
}

void usage(int error) {
    printf("Usage: srf_list [options] srf_file ...\n");
    printf("Options:  -c\tCount only - do not list filenames\n");
    printf("          -v\tVerbose - gives summary data per file too\n");
    printf("          -l\tList in long format. Lines contain:\n");
    printf("            \t    name position body-size header-size\n");

    exit(error);
}

/*
 * Lists the contents of a .hash file
 */
int main(int argc, char **argv) {
    opts opts;
    int i, c;
    int64_t count = 0;

    opts.long_format = 0;
    opts.count_only = 0;
    opts.verbose = 0;

    while ((c = getopt(argc, argv, "lcvh")) != -1) {
	switch (c) {
	case 'l':
	    opts.long_format = 1;
	    break;

	case 'c':
	    opts.count_only = 1;
	    break;

	case 'v':
	    opts.verbose = 1;
	    break;

	case 'h':
	    usage(0);

	default:
	    usage(1);
	}
    }

    for (i = optind; i < argc; i++) {
	int64_t c;

	if (opts.count_only)
	    c = count_file(argv[i], &opts);
	else
	    c = list_file(argv[i], &opts);

	if (c < 0)
	    return 1;

	if (opts.verbose)
	    printf("%s: %d sequences\n", argv[i], c);
	count += c;
    }
    
    if (opts.count_only)
	printf("%"PRId64"\n", count);

    return 0;
}
