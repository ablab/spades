#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <io_lib/hash_table.h>

/*
 * Copies a single named file to stdout.
 * Returns 0 on success
 *         1 on failure
 */
int extract(HashFile *hf, char *file) {
    size_t len;
    char *data;

    if (data = HashFileExtract(hf, file, &len)) {
	fwrite(data, len, 1, stdout);
	free(data);
	return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    char *fofn = NULL;
    char *hash;
    HashFile *hf;
    int ret = 0;

    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-I") == 0) {
	    argv++;
	    fofn = *argv;
	    argc--;
	}
    }

    if (argc < 2 && !fofn) {
	fprintf(stderr, "Usage: hash_extract [-I fofn] hashfile [name ...]\n");
	return 1;
    }
    hash = argv[0];
    argc--;
    argv++;

    if (NULL == (hf = HashFileOpen(hash))) {
	perror(hash);
	return 1;
    }

    if (fofn) {
	FILE *fofnfp;
	char file[256];

	if (strcmp(fofn, "-") == 0) {
	    fofnfp = stdin;
	} else {
	    if (NULL == (fofnfp = fopen(fofn, "r"))) {
		perror(fofn);
		return 1;
	    }
	}

	while (fgets(file, 255, fofnfp)) {
	    char *c;
	    if (c = strchr(file, '\n'))
		*c = 0;

	    ret |= extract(hf, file);
	}

	fclose(fofnfp);
    }

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    for (; argc; argc--, argv++) {
	ret |= extract(hf, *argv);
    }

    HashFileDestroy(hf);

    return ret;
}
