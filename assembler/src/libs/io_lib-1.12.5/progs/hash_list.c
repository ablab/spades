#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <io_lib/hash_table.h>
#include <io_lib/os.h>

/*
 * Dumps a textual represenation of the hash table to stdout.
 */
void HashTableLongDump(HashFile *hf, FILE *fp, int long_format) {
    HashTable *h = hf->h;
    int i;

    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    HashFileItem *hfi;
	    hfi = (HashFileItem *)hi->data.p;
	    if (long_format) {
		char *aname;
		
		if (hf->archives &&
		    hfi->archive >= 0 &&
		    hfi->archive < hf->narchives) {
		    aname = hf->archives[hfi->archive];
		} else {
		    aname = "?";
		}
		fprintf(fp, "%10"PRId64" %6"PRId32" %.*s %s\n",
			hfi->pos, hfi->size, hi->key_len, hi->key,
			aname);
		/*
		fprintf(fp, "%10ld %6d %.*s\n",
			hfi->pos, hfi->size, hi->key_len, hi->key);
		*/
	    } else {
		fprintf(fp, "%.*s\n", hi->key_len, hi->key);
	    }
	}
    }
}


/*
 * Lists the contents of a .hash file
 */
int main(int argc, char **argv) {
    FILE *fp;
    HashFile *hf;
    int long_format = 0;

    /* process command line arguments of the form -arg */
    if (argc >= 2 && strcmp(argv[1], "-l") == 0) {
	long_format = 1;
	argc--;
	argv++;
    }
    if (argc >= 2) {
	fp = fopen(argv[1], "rb");
	if (NULL == fp) {
	    perror(argv[1]);
	    return 1;
	}
    } else {
	fp = stdin;
    }

    hf = HashFileLoad(fp);
    if (hf) {
	HashTableLongDump(hf, stdout, long_format);
	HashFileDestroy(hf);
    }

    return 0;
}
