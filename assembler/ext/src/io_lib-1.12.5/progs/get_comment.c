/*
 * Fetches all or several specific comment(s) from a trace file TEXT section.
 *
 * Usage:
 *	get_comment [options] [field ...] < infile
 *
 * Options:
 *	-c	Suppresses display of field-ID
 *	-h	Help
 *
 * Return codes:
 *	0	Success
 *	1	At least one field was not found
 *	2	Failed to read file, or usage message displayed
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/Read.h>

/* Nasty Microsoft bits */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif

/*
 * From unistd.h
extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;
*/

void usage(void) {
    puts("Usage:");
    puts("    get_comment [options] [field ...] < infile");
    puts("\nOptions:");
    puts("    -c	Suppresses display of field-ID");
    puts("    -h	Help");
    puts("\nReturn codes:");
    puts("    0	Success");
    puts("    1	At least one field was not found");
    puts("    2	Failed to read file, or usage message displayed");
    exit(2);
}

int main(int argc, char **argv) {
    Read *r;
    char *ident, *value;
    int ident_len, value_len;
    enum state_t {
	NAME, EQUALS, VALUE, NL
    } state;
    size_t len;
    int i, j, found;
    int suppress = 0;
    int *found_args = NULL;

    /* Parse arguments */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-')
	    break;

	if (strcmp(*argv, "-c") == 0) {
	    suppress = 1;

	} else {
	    usage();
	}
    }

    /* Read the file */
    read_sections(READ_COMMENTS);
    if (NULL == (r = fread_reading(stdin, "(stdin)", TT_ANY))) {
	fprintf(stderr, "failed to read trace from stdin\n");
	return 2;
    }

    if (!r->info)
	return 1;

    if (argc == 0) {
	/* Display all of them */
	puts(r->info);

    } else {
	/* Display only the ones listed on the command line */

	found_args = (int *)calloc(argc, sizeof(int));

	len = strlen(r->info);
	state = NAME;
	ident = r->info;
	found = 0;
	/* Not needed, but avoids "might be used uninitialized" message */
	value = NULL;
	ident_len = value_len = 0;
	for (i = 0; i <= len; i++) {
	    switch (state) {
	    case NAME:
		if (r->info[i] == '=') {
		    state = EQUALS;
		    ident_len = &r->info[i] - ident;
		    value_len = 0;
		}
		break;

	    case EQUALS:
		for (j = 0; j < argc; j++) {
		    if (strncmp(ident, argv[j], ident_len) == 0) {
			found = 1;
			found_args[j] = 1;
		    }
		}
		state = VALUE;
		value = &r->info[i];

		/* DELIBERATE FLOW THROUGH */

	    case VALUE:
		if (r->info[i] == '\n' || r->info[i] == 0) {
		    value_len = &r->info[i] - value;
		    state = NL;
		}
		break;

	    case NL:
		if (found) {
		    if (suppress) {
			printf("%.*s\n",
			       value_len, value);
		    } else {
			printf("%.*s=%.*s\n",
			       ident_len, ident, value_len, value);
		    }
		}
		state = NAME;
		ident = &r->info[i];
		found = 0;
		break;
	    }
	}
    }
    
    read_deallocate(r);

    if (found_args) {
	for (j = 0; j < argc; j++) {
	    if (found_args[j] == 0)
		return 1;
	}
    }

    return 0;
}
