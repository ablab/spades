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
#include <io_lib/scf.h>

int main(int argc, char **argv) {
    Scf *scf;

    if (argc != 2) {
	fprintf(stderr, "Usage: scf_info scf_filename\n");
	return 1;
    }

    scf = read_scf(argv[1]);

    if (scf == NULL) {
	fprintf(stderr, "read_scf failed\n");
	return 1;
    }

    printf("Version_number       %.4s\n", scf->header.version);
    printf("Number_of_samples    %d\n", scf->header.samples);
    printf("Samples_offset       %d\n", scf->header.samples_offset);
    printf("Samples_size         %d\n", scf->header.sample_size);
    printf("Number_of_bases      %d\n", scf->header.bases);
    printf("Bases_offset         %d\n", scf->header.bases_offset);
    printf("Comments_size        %d\n", scf->header.comments_size);
    printf("Comments_offset      %d\n", scf->header.comments_offset);
    printf("Left_clip            %d\n", scf->header.bases_left_clip);
    printf("Right_clip           %d\n", scf->header.bases_right_clip);
    printf("Code set             %d\n", scf->header.code_set);

    scf_deallocate(scf);

    return 0;
}
