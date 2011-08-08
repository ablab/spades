#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <fcntl.h>
#include <ctype.h>
#include <io_lib/tar_format.h>
#include <io_lib/hash_table.h>

typedef struct {
    int   directories;
    int   verbose;
    int   append_mode;
    int   prepend_mode;
    int   basename;
    char *header;
    char *footer;
    char *archive; /* when reading from stdin */
    HashTable *map;
} options_t;

typedef struct {
    char member[256];
    unsigned char archive;
    uint64_t pos;
    uint32_t size;
} tar_file;

static tar_file *files = NULL;
static int files_alloc = 0;
static int nfiles = 0;

void seek_forward(FILE *fp, int size) {
    if (fp != stdin) {
	fseek(fp, size, SEEK_CUR);
    } else {
	/* Seeking on a pipe isn't supported, even for fwd seeks */
	char buf[8192];
	while (size) {
	    size -= fread(buf, 1, size > 8192 ? 8192 : size, fp);
	}
    }
}

HashTable *load_map(char *fn) {
    HashTable *h = HashTableCreate(65536, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    FILE *fp;
    char line[8192];

    if (NULL == (fp = fopen(fn, "r"))) {
	perror(fn);
	return NULL;
    }

    while(fgets(line, 8192, fp)) {
	char *cp, *from, *to;
	HashData hd;

	for (from = cp = line; *cp && !isspace(*cp); cp++);
	if (!*cp) {
	    fprintf(stderr, "Malformed line '%s'\n", line);
	    return NULL;
	}
	*cp++ = 0;
	
	for (to = cp; isprint(*cp); cp++);
	*cp++ = 0;

	hd.p = strdup(to);
	if (!HashTableAdd(h, from, strlen(from), hd, NULL))
	    return NULL;
    }

    close(fp);

    return h;
}

/*
 * Adds tar index data to the global files[] array.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int accumulate(HashFile *hf, FILE *fp, char *archive, options_t *opt) {
    tar_block blk;
    char member[256];
    int LongLink = 0;
    size_t size, extra;
    size_t offset = 0;

    /* Add to HashFile archives list */
    if (archive) {
	hf->narchives++;
	hf->archives = realloc(hf->archives, hf->narchives * sizeof(char *));
	hf->archives[hf->narchives-1] = strdup(archive);	
    }

    /* Fill out the files[] array with the offsets, size and names */
    while(fread(&blk, sizeof(blk), 1, fp) == 1) {
	/*
	 * If a directory is too large to fit in the name (>100) but short
	 * enough to fit in the prefix the name field will be empty, this is
	 * not the cas for ordinary files where the name field is always
	 * non-empty
	 */
	if (!blk.header.name[0] && !blk.header.prefix[0])
	    break;

        /* get size of member, rounded to a multiple of TBLOCK */
	size = strtoul(blk.header.size, NULL, 8);
        extra = TBLOCK*((size+TBLOCK-1)/TBLOCK) - size;

        /* skip directories unless requested */
        if (opt->directories || blk.header.typeflag != DIRTYPE) {

            /*
	     * extract member name (prefix + name), unless last member
	     * was ././@LongLink
	     */
            if (LongLink == 0) {
                (void) strncpy(member, blk.header.prefix, 155);
	        if (strlen(blk.header.prefix) > 0 && blk.header.name[0])
		    (void) strcat(member, "/");
    	        (void) strncat(member, blk.header.name, 100);
            }
            
            /* account for gtar ././@LongLink */
            if (strcmp(member, "././@LongLink") == 0) {
                /* still expect filenames to fit into 256 bytes */
                if (size > 256) {
                    fread(member, 1, size > 256 ? 256 : size, fp);
                    fprintf(stderr,"././@LongLink too long size=%ld\n",
			    (long)size);
                    fprintf(stderr,"%s...\n", member);
                    exit(1);
                }
                /*
		 * extract full name of next member then rewind to start
		 * of header
		 */
                fread(member, 1, size > 256 ? 256 : size, fp);
                fseek(fp, -size, SEEK_CUR);
                LongLink = 1;
            } else {
                /* output offset, member name */
                /* printf("%lu %.256s\n", (long)offset, member); */
                LongLink = 0;

		if (nfiles >= files_alloc) {
		    if (files_alloc)
			files_alloc *= 2;
		    else
			files_alloc = 1024;
		    files = (tar_file *)realloc(files,
						files_alloc*sizeof(tar_file));
		}
		if (opt->basename) {
		    char *cp = strrchr(member, '/');
		    if (cp)
			memmove(member, cp+1, strlen(cp+1)+1);
		}

		if (opt->map) {
		    HashItem *hi = HashTableSearch(opt->map,
						   member,
						   strlen(member));
		    if (hi) {
			//fprintf(stderr, "Mapped %s to %s\n",
			//	member, hi->data.p);
			strcpy(files[nfiles].member, hi->data.p);
		    } else {
			//fprintf(stderr, "No map for %s\n",
			//	member);
			strcpy(files[nfiles].member, member);
		    }
		} else {
		    strcpy(files[nfiles].member, member);
		}

		files[nfiles].archive = hf->narchives-1;
		files[nfiles].pos = offset+sizeof(blk);
		files[nfiles].size = size;
		if (opt->verbose)
		    fprintf(stderr, "File %d: pos %010ld+%06d: %s\n",
			    nfiles,
			    (long)files[nfiles].pos,
			    files[nfiles].size,
			    files[nfiles].member);

		nfiles++;
            }
        }

        /* increment offset */
        size += extra;
	seek_forward(fp, size);
        offset += sizeof(blk) + size;
    }
}

void link_footers(HashFile *hf, options_t *opt) {
    int found_header = 0, found_footer = 0;
    int i;

    for (i = 0; i < nfiles; i++) {
	if (opt->header && strncmp(opt->header, files[i].member, 256) == 0) {
	    hf->headers[0].archive_no = 0; /* hard-coded, sorry */
	    hf->headers[0].pos  = files[i].pos;
	    hf->headers[0].size = files[i].size;
	    hf->headers[0].cached_data = NULL;
	    found_header++;
	}
	if (opt->footer && strncmp(opt->footer, files[i].member, 256) == 0) {
	    hf->footers[0].archive_no = 0; /* hard-coded, sorry */
	    hf->footers[0].pos  = files[i].pos;
	    hf->footers[0].size = files[i].size;
	    hf->footers[0].cached_data = NULL;
	    found_footer++;
	}
    }

    if (opt->header && !found_header) {
	fprintf(stderr, "Warning: could not find header '%s' in file\n",
		opt->header);
	hf->nheaders = 0;
    }

    if (opt->footer && !found_footer) {
	fprintf(stderr, "Warning: could not find footer '%s' in file\n",
		opt->footer);
	hf->nfooters = 0;
    }
}

void construct_hash(HashFile *hf) {
    int i;

    for (i = 0; i < nfiles; i++) {
	HashData hd;
	HashFileItem *hfi = (HashFileItem *)calloc(1, sizeof(*hfi));

	/* Just use the last head/foot defined as we only allow 1 at the mo. */
	hfi->header  = hf->nheaders;
	hfi->footer  = hf->nfooters;
	hfi->pos     = files[i].pos;
	hfi->size    = files[i].size;
	hfi->archive = files[i].archive;
	hd.p = hfi;
	HashTableAdd(hf->h, files[i].member, strlen(files[i].member),
		     hd, NULL);
    }
}


void save_hash(HashFile *hf, options_t *opt) {
    HashTableStats(hf->h, stderr);
	
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    HashFileSave(hf, stdout, opt->prepend_mode ? HASHFILE_PREPEND : 0);
    HashFileDestroy(hf);
}


int main(int argc, char **argv) {
    options_t opt;
    FILE *fp;
    HashFile *hf;

    /* process command line arguments of the form -arg */
    opt.directories  = 0;
    opt.verbose      = 0;
    opt.append_mode  = 0;
    opt.prepend_mode = 0;
    opt.basename     = 0;
    opt.header       = NULL;
    opt.footer       = NULL;
    opt.archive      = NULL;
    opt.map          = NULL;

    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-a") == 0 && argc > 1) {
	    opt.archive = argv[1];
	    argv++;
	    argc--;
	}

	if (strcmp(*argv, "-A") == 0)
	    opt.append_mode = 1;

	if (strcmp(*argv, "-O") == 0)
	    opt.prepend_mode = 1;

	if (strcmp(*argv, "-d") == 0)
	    opt.directories = 1;

	if (strcmp(*argv, "-v") == 0)
	    opt.verbose = 1;

	if (strcmp(*argv, "-b") == 0)
	    opt.basename = 1;

	if (strcmp(*argv, "-m") == 0 && argc > 1) {
	    /* Name mapping */
	    opt.map = load_map(argv[1]);
	    if (!opt.map) {
		fprintf(stderr, "Failed to load map '%s'\n", argv[1]);
		return 1;
	    }

	    argv++;
	    argc--;
	}

	if (strcmp(*argv, "-h") == 0 && argc > 1) {
	    /* Common header */
	    hf->headers = (HashFileSection *)
		realloc(hf->headers, (hf->nheaders+1) *
			sizeof(HashFileSection));
	    opt.header = argv[1];
	    hf->nheaders++;
	    argv++;
	    argc--;
	}

	if (strcmp(*argv, "-f") == 0 && argc > 1) {
	    /* Common footer */
	    hf->footers = (HashFileSection *)
		realloc(hf->footers, (hf->nfooters+1) *
			sizeof(HashFileSection));
	    opt.footer = argv[1];
	    hf->nfooters++;
	    argv++;
	    argc--;
	}
    }

    if (argc < 1 && !opt.archive) {
	fprintf(stderr, "Usage: hash_tar [options] [tarfile] > tarfile.hash\n");
	fprintf(stderr, "    -a fname  Tar archive filename: use if reading from stdin\n");
	fprintf(stderr, "    -A        Force no archive name (eg will concat to archive itself)\n");
	fprintf(stderr, "    -O        Set arc. offset to size of hash (use when prepending)\n");
	fprintf(stderr, "    -v        Verbose mode\n");
	fprintf(stderr, "    -d        Index directory names (useless?)\n");
	fprintf(stderr, "    -h name   Set tar entry 'name' to be a file header\n");
	fprintf(stderr, "    -f name   Set tar entry 'name' to be a file footer\n");
	fprintf(stderr, "    -b        Use only the filename portion of a pathname\n");
	fprintf(stderr, "    -m fname  Reads lines of 'old new' and renames entries before indexing.");
	return 1;
    }



    /* Load the tar file index into memory */
    hf = HashFileCreate(0, HASH_DYNAMIC_SIZE);
    if (argc < 1) {
	if (!opt.archive) {
	    fprintf(stderr, "If reading from stdin you must use the "
		    "\"-a archivename\" option\n");
	    return 1;
	}
	
	accumulate(hf, stdin, opt.archive, &opt);
    } else {
	/* Single file mode */
	if (opt.append_mode) {
	    if (argc >= 2) {
		fprintf(stderr, "Can only use append_mode with a single "
			"tar file\n");
		return 1;
	    }
	}

	/* Iterate over all tar files */
	while (argc >= 1) {
	    FILE *fp = fopen(argv[0], "rb");
	    if (fp == NULL) {
		perror(argv[0]);
		return 1;
	    }
	    accumulate(hf, fp, argv[0], &opt);
	    fclose(fp);

	    argc--;
	    argv++;
	}
    }

   
    /*
     * Find the header/footer if specified. For now we only support one of
     * each.
     */
    link_footers(hf, &opt);


    /* Construct the hash */
    construct_hash(hf);


    /* Save hash */
    save_hash(hf, &opt);


    /* Tidy up */
    free(files);

    return 0;
}
