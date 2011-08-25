#ifdef HAVE_CONFIG_H
#  include "io_lib_config.h"
#endif

#if !(defined(_MSC_VER) || defined(__MINGW32__))
#  define TRACE_ARCHIVE
#  ifndef HAVE_LIBCURL
#    define USE_WGET
#  endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <limits.h>
#include "io_lib/os.h"
#ifdef TRACE_ARCHIVE
#  include <sys/socket.h>
#  include <netinet/in.h>
#  include <netdb.h>
#  include <sys/time.h>
#  include <errno.h>
#endif
#ifdef USE_WGET
#  include <sys/wait.h>
#endif
#ifndef PATH_MAX
#  define PATH_MAX 1024
#endif
#ifdef HAVE_LIBCURL
#  include <curl/curl.h>
#endif

#include "io_lib/open_trace_file.h"
#include "io_lib/misc.h"
#include "io_lib/tar_format.h"
#include "io_lib/compress.h"
#include "io_lib/hash_table.h"
#include "io_lib/sff.h"
#include "io_lib/srf.h"

/*
 * Supported compression extensions. See the magics array in compress.c for
 * the full structure.
 */
static char *magics[] = {"", ".bz", ".gz", ".Z", ".z", ".bz2", ".sz"};

/*
 * Tokenises the search path splitting on colons (unix) or semicolons (windows).
 * We also  explicitly add a "./" to the end of the search path
 *
 * Returns: A new search path with items separated by nul chars. Two nul
 *          chars in a row represent the end of the tokenised path.
 * Returns NULL for a failure.
 *
 * The returned data has been malloced. It is up to the caller to free this
 * memory.
 */
static char *tokenise_search_path(char *searchpath) {
    char *newsearch;
    unsigned int i, j;
    size_t len;
#ifdef _WIN32
    char path_sep = ';';
#else
    char path_sep = ':';
#endif

    if (!searchpath)
	searchpath="";

    newsearch = (char *)malloc((len = strlen(searchpath))+5);
    if (!newsearch)
	return NULL;

    for (i = 0, j = 0; i < len; i++) {
	/* "::" => ":". Used for escaping colons in http://foo */
	if (i < len-1 && searchpath[i] == ':' && searchpath[i+1] == ':') {
	    newsearch[j++] = ':';
	    i++;
	    continue;
	}

	if (searchpath[i] == path_sep) {
	    /* Skip blank path components */
	    if (j && newsearch[j-1] != 0)
		newsearch[j++] = 0;
	} else {
	    newsearch[j++] = searchpath[i];
	}
    }

    if (j)
	newsearch[j++] = 0;
    newsearch[j++] = '.';
    newsearch[j++] = '/';
    newsearch[j++] = 0;
    newsearch[j++] = 0;
    
    return newsearch;
}

/*
 * Searches for file in the tar pointed to by tarname. If it finds it, it
 * copies it out and returns a file pointer to the temporary file,
 * otherwise we return NULL.
 *
 * If 'tarname'.index exists we will use this as a fast lookup method,
 * otherwise we just do a sequential search through the tar.
 *
 * Offset specifies a starting search position. Set this to zero if you want
 * to search through the entire tar file, otherwise set it to the byte offset
 * into the file of the tar header block for the desired file to extract.
 * (Note that the tar index file overrides this value.)
 *
 * Returns mFILE pointer if found
 *         NULL if not.
 */
static mFILE *find_file_tar(char *file, char *tarname, size_t offset) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    char path[PATH_MAX+101];
    FILE *fp;
    tar_block blk;
    int size;
    int name_len = strlen(file);

    /* Maximum name length for a tar file */
    if (name_len > 100)
	return NULL;

    /* Search the .index file */
    sprintf(path, "%s.index", tarname);
    if (file_exists(path)) {
	FILE *fpind = fopen(path, "r");
	char *cp;
	int tmp_off;
	int found = 0;
	
	if (fpind) {
	    while (fgets(path, PATH_MAX+100, fpind)) {
		if (cp = strchr(path, '\n'))
		    *cp = 0;
		tmp_off = strtol(path, &cp, 10);
		while (isspace(*cp))
		    cp++;
		if (strncmp(cp, file, name_len) == 0) {
		    int i;
		    for (i = 0; i < num_magics; i++) {
			if (strcmp(&cp[name_len], magics[i]) == 0) {
			    offset = tmp_off;
			    found = 1;
			    break;
			}
		    }
		    if (found)
			break;
		}
	    }
	    fclose(fpind);

	    /* Not in index */
	    if (!found)
		return NULL;
	}
    }

    if (NULL == (fp = fopen(tarname, "rb")))
	return NULL;

    /*
     * Search through the tar file (starting from index position) looking
     * for our filename. If there was no index then we start from position 0.
     */
    fseek(fp, offset, SEEK_SET);
    while(fread(&blk, sizeof(blk), 1, fp) == 1) {
	if (!blk.header.name[0])
	    break;

	size = strtol(blk.header.size, NULL, 8);

	/* start with the same name... */
	if (strncmp(blk.header.name, file, name_len) == 0) {
	    char *data;
	    int i;

	    /* ... but does it end with a known compression extension? */
	    for (i = 0; i < num_magics; i++) {
		if (strcmp(&blk.header.name[name_len], magics[i]) == 0) {
		    break;
		}
	    }
	    /* ... apparently not? continue then */
	    if (i == num_magics)
		continue;

	    /* Found it - copy out the data to an mFILE */
	    if (NULL == (data = (char *)malloc(size)))
		return NULL;
	    if (size != fread(data, 1, size, fp)) {
		free(data);
		return NULL;
	    }
	    return mfcreate(data, size);
	}

	fseek(fp, TBLOCK*((size+TBLOCK-1)/TBLOCK), SEEK_CUR);
    }

    fclose(fp);
    return NULL;
}

/*
 * Reads a hash file to look for a filename. The hash file contains the
 * (relative) pathname for the file it is an index for along with the
 * positions and sizes of each file contained within it. The file format
 * of the archive itself is irrelevant provided that the data is not
 * internally compressed in some manner specific to that archive.
 *
 * Return mFILE pointer if found
 *        NULL if not
 */
static mFILE *find_file_hash(char *file, char *hashfile) {
    size_t size;
    static HashFile *hf = NULL;
    static char hf_name[1024];
    char *data;

    /* Cache an open HashFile for fast accesing */
    if (strcmp(hashfile, hf_name) != 0) {
	if (hf)
	    HashFileDestroy(hf);
	hf = HashFileOpen(hashfile);

	if (!hf)
	    return NULL;
	strcpy(hf_name, hashfile);
    }

    /* Search */
    if (NULL == (data = HashFileExtract(hf, file, &size)))
	return NULL;

    /* Found, so copy the contents to a fake FILE pointer */
    return mfcreate(data, size);
}

/*
 * Extracts a single trace from an SRF file.
 *
 * Return mFILE pointer if found
 *        NULL if not
 */
static mFILE *find_file_srf(char *tname, char *srffile) {
    srf_t *srf;
    uint64_t cpos, hpos, dpos;
    mFILE *mf = NULL;
    char *cp;

    if (NULL == (srf = srf_open(srffile, "r")))
	return NULL;

    if (NULL != (cp = strrchr(tname, '/')))
    	tname = cp+1;

    if (0 == srf_find_trace(srf, tname, &cpos, &hpos, &dpos)) {
	char *data = malloc(srf->th.trace_hdr_size + srf->tb.trace_size);
	if (!data) {
	    srf_destroy(srf, 1);
	    return NULL;
	}
	memcpy(data, srf->th.trace_hdr, srf->th.trace_hdr_size);
	memcpy(data + srf->th.trace_hdr_size,
	       srf->tb.trace, srf->tb.trace_size);
	mf = mfcreate(data, srf->th.trace_hdr_size + srf->tb.trace_size);
    }

    srf_destroy(srf, 1);
    return mf;
}

#ifdef TRACE_ARCHIVE
/*
 * Searches for file in the ensembl trace archive pointed to by arcname.
 * If it finds it, it copies it out and returns a file pointer to the
 * temporary file, otherwise we return NULL.
 *
 * Arcname has the form address:port, eg "titan/22100"
 *
 * Returns mFILE pointer if found
 *         NULL if not.
 */
#define RDBUFSZ 8192
static mFILE *find_file_archive(char *file, char *arcname) {
    char server[1024], *cp;
    int port;
    struct hostent *host;
    struct sockaddr_in saddr;
    int s = 0;
    char msg[1024];
    ssize_t msg_len;
    char buf[RDBUFSZ];
    mFILE *fpout;
    int block_count;

    /* Split arc name into server and port */
    if (!(cp = strchr(arcname, '/')))
	return NULL;
    strncpy(server, arcname, 1023);
    server[MIN(1023,cp-arcname)] = 0;
    port = atoi(cp+1);

    /* Make and connect socket */
    if (NULL == (host = gethostbyname(server))) {
	perror("gethostbyname()");
	return NULL;
    }
    saddr.sin_port = htons(port);
    saddr.sin_family = host->h_addrtype;
    memcpy(&saddr.sin_addr,host->h_addr_list[0], host->h_length);
    if ((s = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == -1) {
	perror("socket()");
	return NULL;
    }
    if (connect(s, (struct sockaddr *)&saddr, sizeof(saddr)) == -1) {
	perror("connect()");
	return NULL;
    }

    /* The minimal message to send down is "--scf tracename" */
    sprintf(msg, "--scf %.*s\n", 1000, file);
    msg_len = strlen(msg);
    if (send(s, msg, msg_len, 0) != msg_len) {
	/*
	 * partial request sent, but requests are short so if this
	 * happens it's unlikely we'll cure it by sending multiple
	 * fragments.
	 */
	/* close(s); */
	return NULL;
    }

    /*
     * Create a fake FILE (mFILE) and write to it.
     */
    fpout = mfcreate(NULL, 0);

    /*
     * Read the data back, in multiple blocks if necessary and write it
     * to our temporary file. We use a blocking read with a low timeout to
     * prevent locking up the application indefinitely.
     */
    {
	struct timeval tv = {0, 10000};
	setsockopt(s, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv, sizeof(tv));
    }
    errno = 0;
    block_count = 200;
    while ((msg_len = read(s, buf, RDBUFSZ)) > 0 ||
	   (errno == EWOULDBLOCK && --block_count)) {
	errno = 0;
	if (msg_len > 0)
	    mfwrite(buf, 1, msg_len, fpout);
    }
    close(s);

    if (!block_count) {
	mfclose(fpout);
	return NULL;
    }

    mrewind(fpout);

    return fpout;
}
#endif

#ifdef USE_WGET
/* NB: non-reentrant due to reuse of handle */
static mFILE *find_file_url(char *file, char *url) {
    char buf[8192], *cp;
    mFILE *fp;
    int pid;
    int maxlen = 8190 - strlen(file);
    char *fname = tempnam(NULL, NULL);
    int status;

    /* Expand %s for the trace name */
    for (cp = buf; *url && cp - buf < maxlen; url++) {
	if (*url == '%' && *(url+1) == 's') {
	    url++;
	    cp += strlen(strcpy(cp, file));
	} else {
	    *cp++ = *url;
	}
    }
    *cp++ = 0;

    /* Execute wget */
    if ((pid = fork())) {
	waitpid(pid, &status, 0);
    } else {
	execlp("wget", "wget", "-q", "-O", fname, buf, NULL);
    }

    /* Return a filepointer to the result (if it exists) */
    fp = (!status && file_size(fname) != 0) ? mfopen(fname, "rb") : NULL;
    remove(fname);
    free(fname);

    return fp;
}
#endif

#ifdef HAVE_LIBCURL
static mFILE *find_file_url(char *file, char *url) {
    char buf[8192], *cp;
    mFILE *mf = NULL, *headers = NULL;
    int maxlen = 8190 - strlen(file);
    static CURL *handle = NULL;
    static int curl_init = 0;
    char errbuf[CURL_ERROR_SIZE];

    *errbuf = 0;

    if (!curl_init) {
	if (curl_global_init(CURL_GLOBAL_ALL))
	    return NULL;

	if (NULL == (handle = curl_easy_init()))
	    goto error;

	curl_init = 1;
    }

    /* Expand %s for the trace name */
    for (cp = buf; *url && cp - buf < maxlen; url++) {
	if (*url == '%' && *(url+1) == 's') {
	    url++;
	    cp += strlen(strcpy(cp, file));
	} else {
	    *cp++ = *url;
	}
    }
    *cp++ = 0;

    /* Setup the curl */
    if (NULL == (mf = mfcreate(NULL, 0)) ||
	NULL == (headers = mfcreate(NULL, 0)))
	return NULL;

    if (0 != curl_easy_setopt(handle, CURLOPT_URL, buf))
	goto error;
    if (0 != curl_easy_setopt(handle, CURLOPT_TIMEOUT, 10L))
	goto error;
    if (0 != curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, mfwrite))
	goto error;
    if (0 != curl_easy_setopt(handle, CURLOPT_WRITEDATA, mf))
	goto error;
    if (0 != curl_easy_setopt(handle, CURLOPT_HEADERFUNCTION, mfwrite))
	goto error;
    if (0 != curl_easy_setopt(handle, CURLOPT_WRITEHEADER, headers))
	goto error;
    if (0 != curl_easy_setopt(handle, CURLOPT_ERRORBUFFER, errbuf))
	goto error;
    
    /* Fetch! */
    if (0 != curl_easy_perform(handle))
	goto error;
    
    /* Report errors is approproate. 404 is silent as it may have just been
     * a search via RAWDATA path, everything else is worth reporting.
     */
    {
	float version;
	int response;
	char nul = 0;
	mfwrite(&nul, 1, 1, headers);
	if (2 == sscanf(headers->data, "HTTP/%f %d", &version, &response)) {
	    if (response != 200) {
		if (response != 404)
		    fprintf(stderr, "%.*s\n",
			    (int)headers->size, headers->data);
		goto error;
	    }
	}
    }

    if (mftell(mf) == 0)
	goto error;

    mfdestroy(headers);

    mrewind(mf);
    return mf;

 error:
    if (mf)
	mfdestroy(mf);
    if (headers)
	mfdestroy(headers);
    if (*errbuf)
	fprintf(stderr, "%s\n", errbuf);
    return NULL;
}
#endif

/*
 * Takes an SFF file in 'data' and edits the header to ensure
 * that it has no index listed and only claims to contain a single entry.
 * This isn't strictly necessary for the sff/sff.c reading code, but it is
 * the 'Right Thing' to do.
 *
 * Returns an mFILE on success or NULL on failure.
 */
static mFILE *sff_single(char *data, size_t size) {
    *(uint64_t *)(data+8)  = be_int8(0); /* index offset */
    *(uint32_t *)(data+16) = be_int4(0); /* index size */
    *(uint32_t *)(data+20) = be_int4(1); /* number of reads */

    return mfcreate(data, size);
}

/* Hash (.hsh) format index searching for SFF files */
static mFILE *sff_hash_query(char *sff, char *entry, FILE *fp) {
    static HashFile *hf = NULL;
    static char sff_copy[1024];
    static FILE *fp_copy = NULL;
    char *data;
    size_t size;

    /* Cache an open HashFile for fast accessing */
    if (strcmp(sff, sff_copy) != 0) {
	if (hf) {
	    hf->afp = hf->hfp = NULL; /* will be closed by our parent */
	    HashFileDestroy(hf);
	}
	fseek(fp, -4, SEEK_CUR);
	if (NULL == (hf = HashFileFopen(fp)))
	    return NULL;

	strcpy(sff_copy, sff);
	fp_copy = fp;
    }

    data = HashFileExtract(hf, entry, &size);
    
    return data ? sff_single(data, size) : NULL;
}


/*
 * getuint4_255
 *
 * A function to convert a 4-byte TVF/SFF value into an integer, where
 * the bytes are base 255 numbers.  This is used to store the index offsets.
 */
static unsigned int getuint4_255(unsigned char *b)
{
    return
        ((unsigned int) b[0]) * 255 * 255 * 255 +
        ((unsigned int) b[1]) * 255 * 255 +
        ((unsigned int) b[2]) * 255 +
        ((unsigned int) b[3]);
}

/*
 * 454 sorted format (.srt) index searching for SFF files.
 * Uses a binary search.
 * This function and getuint4_255 above are taken with permission
 * from 454's getsff.c with the following licence:
 *
 * Copyright (c)[2001-2005] 454 Life Sciences Corporation. All Rights Reserved.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 *
 * IN NO EVENT SHALL LICENSOR BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided
 * that this copyright and notice appears in all copies.
 */
static mFILE *sff_sorted_query(char *sff, char *accno, FILE *fp,
			       uint32_t index_length) {
    static unsigned char *index;
    static char sff_copy[1024];
    unsigned char *us;
    uint32_t start, end;
    uint32_t offset;
    char *data = NULL;
    static char chdr[1024];
    static int chdrlen = 0, nflows = 0;
    char rhdr[1024];
    int rhdrlen;
    int nbases, dlen;
    int bytes_per_flow = 2;

    /* Cache index if we're querying the same SFF file */
    if (strcmp(sff_copy, sff) != 0) {
	if (index)
	    xfree(index);
	if (NULL == (index = (unsigned char *)xmalloc(index_length)))
	    return NULL;

	if (index_length != fread(index, 1, index_length, fp)) {
	    xfree(index);
	    return NULL;
	}
	strcpy(sff_copy, sff);

	/* Read the common header too - minimal decoding necessary */
	fseek(fp, 0, SEEK_SET);
	if (31 != fread(chdr, 1, 31, fp))
	    return NULL;
	chdrlen = be_int2(*(uint16_t *)(chdr+24));
	nflows  = be_int2(*(uint16_t *)(chdr+28));
	if (chdrlen-31 != fread(chdr+31, 1, chdrlen-31, fp))
	    return NULL;
    }

    /*
     * Perform a binary search of the index, stopping when the search
     * region becomes relatively small.  This assumes that no accession
     * number is near 200 characters.
     */
    start = 0;
    end = index_length;
    while (end - start > 200) {
	uint32_t mid;
	int val;
        mid = (start + end) / 2;

        /*
         * From the byte midpoint, scan backwards to the beginning of the
         * index record that covers that byte midpoint.
         */
        while (mid > start && index[mid-1] != 255) {
            mid--;
        }
        val = strcmp(accno, (char *)(index+mid));

        if (val == 0) {
            break;
        } else if (val < 0) {
            end = mid;
        } else {
            start = mid;
        }
    }

    /*
     * Scan through the small search region, looking for the accno.
     */
    while (start < end) {
        if (strcmp(accno, (char *)(index+start)) == 0) {
            /*
             * If the accno is found, skip the accno characters,
             * then get the record offset.
             */
            for (us=index+start; *us; us++,start++) ;
            us++;
            start++;

            offset = getuint4_255(us);
            if (us[4] != 255) {
		return NULL;
            }

	    /*
	     * The original getsff.c here computed the record size by
	     * looking at the next index item and comparing it's offset to
	     * this one, or the end of file position if this is the last
	     * item. This has two problems:
	     * 1: It means the index itself cannot be added to the end of
	     *    the file.
	     * 2: It means that we cannot simply add an index to a SFF
	     *    file without also reordering all of the items within it.
	     *
	     * We solve this by reading the read header to work out the
	     * object size instead.
	     */
	    break;
        }

        /*
         * Skip to the beginning of the next index element.
         */
        while (start < end && index[start] != 255) {
            start++;
        }
        start++;
    }

    /*
     * Now offset indicates the position of the SFF entry. Read and decode
     * header to get data length. Then read this too.
     */
    fseek(fp, offset, SEEK_SET);
    if (16 != fread(rhdr, 1, 16, fp))
	return NULL;

    rhdrlen   = be_int2(*(uint16_t *)rhdr);
    nbases   = be_int4(*(uint32_t *)(rhdr+4));
    
    if (rhdrlen-16 != fread(rhdr+16, 1, rhdrlen-16, fp))
	return NULL;
    dlen = (nflows * bytes_per_flow + nbases * 3 + 7) & ~7;

    /* Built up the fake SFF entry */
    if (NULL == (data = (char *)xmalloc(chdrlen + rhdrlen + dlen)))
	return NULL;

    memcpy(data, chdr, chdrlen);
    memcpy(data + chdrlen, rhdr, rhdrlen);
    if (dlen != fread(data + chdrlen + rhdrlen, 1, dlen, fp)) {
	xfree(data);
	return NULL;
    }

    /* Convert to mFILE */
    return sff_single(data, chdrlen + rhdrlen + dlen);
}


/*
 * This returns an mFILE containing an SFF entry.
 *
 * This does the minimal decoding necessary to skip through the SFF
 * container to find an entry. In this respect it is a semi-duplication
 * of sff/sff.[ch], but implemented for efficiency.
 *
 * Having found an entry it packs the common header, the read specific
 * header and the read data into a single block of memory and returns this
 * as an mFILE. In essence it produces a single-read SFF archive. This
 * is then decoded by the normal sff parsing code representing a small
 * amount of redundancy, but one which is swamped by the I/O time.
 */
static mFILE *find_file_sff(char *entry, char *sff) {
    static FILE *fp = NULL;
    static char sff_copy[1024];
    char chdr[65536], rhdr[65536]; /* generous, but worst case */
    uint32_t nkey, nflows, chdrlen, rhdrlen, dlen, magic;
    uint64_t file_pos;
    static uint64_t index_offset = 0;
    static uint32_t index_length = 0;
    static char index_format[8];
    uint32_t nreads, i;
    size_t entry_len = strlen(entry);
    int bytes_per_flow = 2;
    char *fake_file;

    /*
     * Check cached information so rapid queries to the same archive are
     * fast.
     * ASSUMPTION: we won't externally replace the sff file with another of
     * the same name.
     */
    if (strcmp(sff, sff_copy) == 0) {
	if (memcmp(index_format, ".hsh1.00", 8) == 0) {
	    return sff_hash_query(sff, entry, fp);
	} else if (memcmp(index_format, ".srt1.00", 8) == 0) {
	    return sff_sorted_query(sff, entry, fp, index_length-8);
	}
    }

    if (fp)
	fclose(fp);

    strcpy(sff_copy, sff);
    *index_format = 0;


    /* Read the common header */
    if (NULL == (fp = fopen(sff, "rb")))
	return NULL;
    if (31 != fread(chdr, 1, 31, fp))
	return NULL;

    /* Check magic & vers: TODO */
    magic = be_int4(*(uint32_t *)chdr);
    if (magic != SFF_MAGIC)
	return NULL;
    if (memcmp(chdr+4, SFF_VERSION, 4) != 0)
	return NULL;

    /* If we have an index, use it, otherwise search linearly */
    index_offset = be_int8(*(uint64_t *)(chdr+8));
    index_length = be_int4(*(uint32_t *)(chdr+16));
    if (index_length != 0) {
	long orig_pos = ftell(fp);
	fseek(fp, index_offset, SEEK_SET);
	fread(index_format, 1, 8, fp);

	if (memcmp(index_format, ".hsh1.00", 8) == 0) {
	    /* HASH index v1.00 */
	    return sff_hash_query(sff, entry, fp);

	} else if (memcmp(index_format, ".srt1.00", 8) == 0) {
	    /* 454 sorted v1.00 */
	    return sff_sorted_query(sff, entry, fp, index_length-8);
	} else {
	    /* Unknown index: revert back to a slow linear scan */
	    fseek(fp, orig_pos, SEEK_SET);
	}
    }

    nreads  = be_int4(*(uint32_t *)(chdr+20));
    chdrlen = be_int2(*(uint16_t *)(chdr+24));
    nkey    = be_int2(*(uint16_t *)(chdr+26));
    nflows  = be_int2(*(uint16_t *)(chdr+28));

    /* Read the remainder of the header */
    if (chdrlen-31 != fread(chdr+31, 1, chdrlen-31, fp))
	return NULL;

    file_pos = chdrlen;

    /* Loop until we find the correct entry */
    for (i = 0; i < nreads; i++) {
	uint16_t name_len;
	uint32_t nbases;

	/* Index could be between common header and first read - skip */
	if (file_pos == index_offset) {
	    fseek(fp, index_length, SEEK_CUR);
	    file_pos += index_length;
	}

	/* Read 16 bytes to get name length */
	if (16 != fread(rhdr, 1, 16, fp))
	    return NULL;
	rhdrlen   = be_int2(*(uint16_t *)rhdr);
	name_len = be_int2(*(uint16_t *)(rhdr+2));
	nbases   = be_int4(*(uint32_t *)(rhdr+4));

	/* Read the rest of the header */
	if (rhdrlen-16 != fread(rhdr+16, 1, rhdrlen-16, fp))
	    return NULL;

	file_pos += rhdrlen;

	dlen = (nflows * bytes_per_flow + nbases * 3 + 7) & ~7;

	if (name_len == entry_len  && 0 == memcmp(rhdr+16, entry, entry_len))
	    break;

	/* This is not the read you are looking for... */
	fseek(fp, dlen, SEEK_CUR);
    }

    if (i == nreads) {
	/* Not found */
	return NULL;
    }

    /*
     * Although we've decoded some bits already, we take the more modular
     * approach of packing the sections together and passing the entire
     * data structure off as a single-read SFF file to be decoded fully
     * by the sff reading code.
     */
    if (NULL == (fake_file = (char *)xmalloc(chdrlen + rhdrlen + dlen)))
	return NULL;

    memcpy(fake_file, chdr, chdrlen);
    memcpy(fake_file+chdrlen, rhdr, rhdrlen);
    if (dlen != fread(fake_file+chdrlen+rhdrlen, 1, dlen, fp)) {
	xfree(fake_file);
	return NULL;
    }

    /* Convert to an mFILE and return */
    return sff_single(fake_file, chdrlen+rhdrlen+dlen);
}

/*
 * Searches for file in the directory 'dirname'. If it finds it, it opens
 * it. This also searches for compressed versions of the file in dirname
 * too.
 *
 * Returns mFILE pointer if found
 *         NULL if not
 */
static mFILE *find_file_dir(char *file, char *dirname) {
    char path[PATH_MAX+1], path2[PATH_MAX+1];
    size_t len = strlen(dirname);
    char *cp;

    if (dirname[len-1] == '/')
	len--;

    /* Special case for "./" or absolute filenames */
    if (*file == '/' || (len==1 && *dirname == '.'))
	sprintf(path, "%s", file);
    else
	sprintf(path, "%.*s/%s", (int)len, dirname, file);

    if (is_file(path)) {
	return mfopen(path, "rb");
    }

    /*
     * Given a pathname /a/b/c if a/b is a file and not a directory then
     * we'd get an ENOTDIR error. Instead we assume that a/b is an archive
     * and we attempt to work out what type by reading the first and last
     * bits of the file.
     */
    if (cp = strrchr(file, '/')) {
	strcpy(path2, path); /* path contains / too as it's from file */
	*strrchr(path2, '/') = 0;

	if (is_file(path2)) {
	    /* Open the archive to test for magic numbers */
	    char magic[8];
	    FILE *fp;
	    enum archive_type_t {
		NONE, HASH, TAR, SFF, SRF
	    } type = NONE;

	    if (NULL == (fp = fopen(path2, "rb")))
		return NULL;
	    memcpy(magic, "\0\0\0\0\0\0", 4);
	    fread(magic, 1, 4, fp);

	    /* .hsh or .sff at start */
	    if (memcmp(magic, ".hsh", 4) == 0)
		type = HASH;
	    else if (memcmp(magic, ".sff", 4) == 0)
		type = SFF;

	    /* Or .hsh or Ihsh at the end */
	    if (NONE == type) {
		fseek(fp, -16, SEEK_END);
		fread(magic, 1, 8, fp);
		if (memcmp(magic+4, ".hsh", 4) == 0)
		    type = HASH;
		else if (memcmp(magic, "Ihsh", 4) == 0)
		    type = SRF;
	    }

	    /* or ustar 257 bytes in to indicate un-hashed tar */
	    if (NONE == type) {
		fseek(fp, 257, SEEK_SET);
		fread(magic, 1, 5, fp);
		if (memcmp(magic, "ustar", 5) == 0)
		    type = TAR;
	    }
	    fclose(fp);

	    switch (type) {
	    case HASH:
		return find_file_hash(cp+1, path2);
	    case TAR:
		return find_file_tar(cp+1, path2, 0);
	    case SFF:
		return find_file_sff(cp+1, path2);
	    case SRF:
		return find_file_srf(cp+1, path2);
	    case NONE:
		break;
	    }

	    return NULL;
	}
    }

    return NULL;
}

/*
 * ------------------------------------------------------------------------
 * Public functions below.
 */

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in 'path' (which is a colon separated list).
 * If 'path' is NULL it uses the RAWDATA environment variable instead.
 *
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_path_mfile(char *file, char *path, char *relative_to) {
    char *newsearch;
    char *ele;
    mFILE *fp;

    /* Use path first */
    if (!path)
	path = getenv("RAWDATA");
    if (NULL == (newsearch = tokenise_search_path(path)))
	return NULL;
    
    /*
     * Step through the search path testing out each component.
     * We now look through each path element treating some prefixes as
     * special, otherwise we treat the element as a directory.
     */
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
	int i;
	char *suffix[6] = {"", ".gz", ".bz2", ".sz", ".Z", ".bz2"};
	for (i = 0; i < 6; i++) {
	    char file2[1024];
	    char *ele2;
	    int valid = 1;

	    /*
	     * '|' prefixing a path component indicates that we do not
	     * wish to perform the compression extension searching in that
	     * location.
	     */
	    if (*ele == '|') {
		ele2 = ele+1;
		valid = (i == 0);
	    } else {
		ele2 = ele;
	    }

	    sprintf(file2, "%s%s", file, suffix[i]);

	    if (0 == strncmp(ele2, "TAR=", 4)) {
		if (valid && (fp = find_file_tar(file2, ele2+4, 0))) {
		    free(newsearch);
		    return fp;
		}

	    } else if (0 == strncmp(ele2, "HASH=", 5)) {
		if (valid && (fp = find_file_hash(file2, ele2+5))) {
		    free(newsearch);
		    return fp;
		}
#ifdef TRACE_ARCHIVE
	    } else if (0 == strncmp(ele2, "ARC=", 4)) {
		if (valid && (fp = find_file_archive(file2, ele2+4))) {
		    free(newsearch);
		    return fp;
		}
#endif
#if defined(USE_WGET) || defined(HAVE_LIBCURL)
	    } else if (0 == strncmp(ele2, "URL=", 4)) {
		if (valid && (fp = find_file_url(file2, ele2+4))) {
		    free(newsearch);
		    return fp;
		}
#endif
	    } else if (0 == strncmp(ele2, "SFF=", 4)) {
		if (valid && (fp = find_file_sff(file2, ele2+4))) {
		    free(newsearch);
		    return fp;
		}

	    } else if (0 == strncmp(ele2, "SRF=", 4)) {
		if (valid && (fp = find_file_srf(file2, ele2+4))) {
		    free(newsearch);
		    return fp;
		}

	    } else {
		if (valid && (fp = find_file_dir(file2, ele2))) {
		    free(newsearch);
		    return fp;
		}
	    }
	}
    }

    free(newsearch);

    /* Look in the same location as the incoming 'relative_to' filename */
    if (relative_to) {
	char *cp;
	char relative_path[PATH_MAX+1];
	strcpy(relative_path, relative_to);
	if (cp = strrchr(relative_path, '/'))
	    *cp = 0;
	if (fp = find_file_dir(file, relative_path))
	    return fp;
    }

    return NULL;
}

FILE *open_path_file(char *file, char *path, char *relative_to) {
    mFILE *mf = open_path_mfile(file, path, relative_to);
    FILE *fp;

    if (!mf)
	return NULL;

    if (mf->fp)
	return mf->fp;

    /* Secure temporary file generation */
    if (NULL == (fp = tmpfile()))
	return NULL;

    /* Copy the data */
    fwrite(mf->data, 1, mf->size, fp);
    rewind(fp);
    mfclose(mf);

    return fp;
}

static char *exp_path = NULL;
static char *trace_path = NULL;

void  iolib_set_trace_path(char *path) { trace_path = path; }
char *iolib_get_trace_path(void)       { return trace_path; }
void  iolib_set_exp_path  (char *path) { exp_path = path; }
char *iolib_get_exp_path  (void)       { return exp_path; }

/*
 * Trace file functions: uses TRACE_PATH environment variable.
 */
mFILE *open_trace_mfile(char *file, char *rel_to) {
    return open_path_mfile(file, trace_path ? trace_path
			                    : getenv("TRACE_PATH"), rel_to);
}

FILE *open_trace_file(char *file, char *rel_to) {
    return open_path_file(file, trace_path ? trace_path
			                   : getenv("TRACE_PATH"), rel_to);
}

/*
 * Trace file functions: uses EXP_PATH environment variable.
 */
mFILE *open_exp_mfile(char *file, char *relative_to) {
    return open_path_mfile(file, exp_path ? exp_path
			                  : getenv("EXP_PATH"), relative_to);
}

FILE *open_exp_file(char *file, char *relative_to) {
    return open_path_file(file, exp_path ? exp_path
			                 : getenv("EXP_PATH"), relative_to);
}

