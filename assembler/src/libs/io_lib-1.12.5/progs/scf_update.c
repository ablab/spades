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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <io_lib/scf.h>

/*
 * Update a v2.x SCF files to a v3.x SCF files.
 */

int main(int argc, char **argv) {
    Scf *scf;
    int version = 3;

    if (argc > 2 && strcmp(argv[1], "-v") == 0) {
	version = atoi(argv[2]);
	if (version != 2 && version != 3) {
	    fprintf(stderr, "Please specify version 2 or version 3\n");
	    return 1;
	}

	argc-=2;
	argv+=2;
    }

    if (argc != 3) {
	fprintf(stderr, "Usage: scf_update [-v version] source destination\n");
	return 1;
    }

    if (NULL == (scf = read_scf(argv[1]))) {
	perror(argv[1]);
	return 1;
    }

    set_scf_version(version);

    if (-1 == write_scf(scf, argv[2])) {
	perror(argv[2]);
	return 1;
    }

    scf_deallocate(scf);
    return 0;
}
