#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <io_lib/tar_format.h>

int main(int argc, char **argv) {
    int directories = 0;
    FILE *fp;
    tar_block blk;
    char member[256];
    size_t size, extra;
    int LongLink = 0;
    size_t offset = 0;
    
    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-d") == 0)
	    directories = 1;
    }

    if (argc != 1) {
	fprintf(stderr, "Usage: index_tar [-d] tarfile > tarfile.index\n");
	return 1;
    }

    /* open the tarfile */
    if (NULL == (fp = fopen(argv[0], "rb"))) {
	perror(argv[0]);
	return 1;
    }

    while(fread(&blk, sizeof(blk), 1, fp) == 1) {
	/*
	 * If a directory is too large to fit in the name (>100) but short
	 * enough to fit in the prefix the name field will be empty, this is
	 * not the case for ordinary files where the name field is always
	 * non-empty
	 */
	if (!blk.header.name[0] && !blk.header.prefix[0])
	    break;

        /* get size of member, rounded to a multiple of TBLOCK */
	size = strtoul(blk.header.size, NULL, 8);
        extra = TBLOCK*((size+TBLOCK-1)/TBLOCK) - size;

        /* skip directories unless requested */
        if (directories || blk.header.typeflag != DIRTYPE || LongLink) {

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
                printf("%lu %.256s\n", (long)offset, member);
                LongLink = 0;
            }
        }

        /* increment offset */
        size += extra;
        fseek(fp, size, SEEK_CUR);
        offset += sizeof(blk) + size;
    }

    fclose(fp);
    return 0;
}
