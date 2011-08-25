#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "io_lib/Read.h"
#include "io_lib/misc.h"
#include "io_lib/ztr.h"
#include "io_lib/hash_table.h"

#include "io_lib/srf.h"

/*
 * ---------------------------------------------------------------------------
 * Object creation / destruction.
 */

/*
 * Allocates and returns an srf_t structure.
 * fp is the file point associated with this SRF file. Pass as NULL if unknown
 * at this stage.
 *
 * Returns malloced memory on success
 *         NULL on failure
 */
srf_t *srf_create(FILE *fp) {
    srf_t *srf = (srf_t *)calloc(1, sizeof(*srf));

    if (srf)
	srf->fp = fp;
    return srf;
}

/*
 * Opens an SRF archive for reading or writing as defined by 'mode'. Mode is
 * passed directly on to fopen() and so uses the same flags.
 *
 * Returns a srf_t struct pointer on success.
 *         NULL on failure
 */
srf_t *srf_open(char *fn, char *mode) {
    FILE *fp;
    char bmode[11];
    size_t l, i;

    /* Enforce binary mode for windows */
    if ((l = strlen(mode)) < 9) {
	int binary = 0;
	for (i = 0; i < l; i++) {
	    if ('b' == (bmode[i] = mode[i]))
		binary=1;
	}
	if (!binary)
	    bmode[i++] = 'b';
	bmode[i] = 0;
	mode = bmode;
    }
    
    return (fp = fopen(fn, mode)) ? srf_create(fp) : NULL;
}

/*
 * Deallocates an srf_t struct. If auto_close is true then it also closes
 * any associated FILE pointer.
 */
void srf_destroy(srf_t *srf, int auto_close) {
    if (!srf)
	return;
    
    if (auto_close && srf->fp) {
	if (-1 == fclose(srf->fp))
	    perror("fclose(srf->fp)");
    }

    if(srf->th.trace_hdr)
	free(srf->th.trace_hdr);

    if (srf->mf)
	mfdestroy(srf->mf);

    if (srf->ztr)
	delete_ztr(srf->ztr);

    free(srf);
}

/*
 * ---------------------------------------------------------------------------
 * Base data type I/O.
 */

/* 
 * Writes a null-terminated C string in pascal-string form.
 * Returns the number of bytes written or
 *         -1 for failure.
 */
int srf_write_pstring(srf_t *srf, char *str) {
    size_t l = str ? strlen(str) : 0;
    if (l > 255)
	return -1;

    if (l)
	return fprintf(srf->fp, "%c%s", (int)l, str);
    else
	return fprintf(srf->fp, "%c", (int)l);
}

/*
 * As per srf_write_pstring but 'str' may hold binary data including nul
 * chars, hence we pass over length as a separate paremeter.
 */
int srf_write_pstringb(srf_t *srf, char *str, int length) {
    if (length > 255 || length < 0)
	return -1;

    if (length)
	return fprintf(srf->fp, "%c%s", (int)length, str);
    else
	return fprintf(srf->fp, "%c", (int)length);
}

/*
 * Reads a pascal-style string from the srf file.
 * 'str' passed in needs to be at least 256 bytes long. The string read
 * will be stored there and nul terminated.
 *
 * Returns the length of the string read (minus nul char),
 *         -1 for failure
 */
int srf_read_pstring(srf_t *srf, char *str) {
    int len;

    if (EOF == (len = fgetc(srf->fp)))
	return -1;
    if (len != fread(str, 1, len, srf->fp))
	return -1;
    str[len] = '\0';

    return len;
}

/*
 * Read/write unsigned 32-bit and 64-bit values in big-endian format
 * (ie the same as ZTR and SCF endian uses).
 * Functions all return 0 for success, -1 for failure.
 */
int srf_read_uint32(srf_t *srf, uint32_t *val) {
    unsigned char d[4];
    if (1 != fread(d, 4, 1, srf->fp))
	return -1;

    *val = (d[0] << 24) | (d[1] << 16) | (d[2] << 8) | (d[3] << 0);
    return 0;
}

int srf_write_uint32(srf_t *srf, uint32_t val) {
    unsigned char d[4];
    d[0] = (val >> 24) & 0xff;
    d[1] = (val >> 16) & 0xff;
    d[2] = (val >>  8) & 0xff;
    d[3] = (val >>  0) & 0xff;

    return fwrite(d, 4, 1, srf->fp) ? 0 : -1;
}

int srf_read_uint64(srf_t *srf, uint64_t *val) {
    unsigned char d[8];
    if (1 != fread(d, 8, 1, srf->fp))
	return -1;

    *val = ((uint64_t)d[0] << 56)
	 | ((uint64_t)d[1] << 48)
	 | ((uint64_t)d[2] << 40)
	 | ((uint64_t)d[3] << 32)
	 | ((uint64_t)d[4] << 24)
	 | ((uint64_t)d[5] << 16)
	 | ((uint64_t)d[6] <<  8)
	 | ((uint64_t)d[7] <<  0);
    return 0;
}

int srf_write_uint64(srf_t *srf, uint64_t val) {
    unsigned char d[8];
    d[0] = (val >> 56) & 0xff;
    d[1] = (val >> 48) & 0xff;
    d[2] = (val >> 40) & 0xff;
    d[3] = (val >> 32) & 0xff;
    d[4] = (val >> 24) & 0xff;
    d[5] = (val >> 16) & 0xff;
    d[6] = (val >>  8) & 0xff;
    d[7] = (val >>  0) & 0xff;

    return fwrite(d, 8, 1, srf->fp) ? 0 : -1;
}

/*
 * ---------------------------------------------------------------------------
 * Mid level I/O - srf block type handling
 */


/*
 * Allocates and initialises a container header structure. An existing
 * srf_cont_hdr_t may be passed in to avoid allocation of a new object.
 *
 * Returns: allocated cont_header on success, to be freed using
 *          srf_destroy_cont_hdr() (if allocated here),
 *          NULL on failure.
 */
srf_cont_hdr_t *srf_construct_cont_hdr(srf_cont_hdr_t *ch,
				       char *bc,
				       char *bc_version) {
    if (!ch) {
	if (NULL == (ch = (srf_cont_hdr_t *)calloc(1, sizeof(*ch))))
	    return NULL;
    }

    ch->block_type = SRFB_CONTAINER;
    strcpy(ch->version, SRF_VERSION);
    ch->container_type = 'Z';
    strncpy(ch->base_caller, bc, 255);
    strncpy(ch->base_caller_version, bc_version, 255);

    return ch;
}

/*
 * Deallocates an srf_cont_hdr_t constructed by srf_construct_cont_hdr().
 */
void srf_destroy_cont_hdr(srf_cont_hdr_t *ch) {
    if (ch)
	free(ch);
}

/*
 * Reads a container header and stores the result in 'ch'.
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch) {
    char magic[3];
    uint32_t sz;

    if (!ch)
	return -1;

    /* Check block type */
    if (EOF == (ch->block_type = fgetc(srf->fp)))
	return -1;
    if (ch->block_type != SRFB_CONTAINER)
	return -1;

    /* Check magic number && version */
    if (3 != fread(magic, 1, 3, srf->fp))
	return -1;
    if (0 != srf_read_uint32(srf, &sz))
	return -1;
    if (srf_read_pstring(srf, ch->version) < 0)
	return -1;
    if (strncmp(magic, "SRF", 3) || strcmp(ch->version, SRF_VERSION))
	return -1;
    
    /* Containter type, base caller bits */
    if (EOF == (ch->container_type = fgetc(srf->fp)) ||
	srf_read_pstring(srf, ch->base_caller) < 0||
	srf_read_pstring(srf, ch->base_caller_version) < 0)
	return -1;

    return 0;
}

/*
 * Writes a container header to disk.
 *
 * 4 Block type + magic number ("SSRF")
 * 1+n: pString for version number ("1.0")
 * 1: "Z" => container type is ZTR
 * 1+n: pString for base-caller ("eg Bustard")
 * 1+n: pString for base-caller version (eg "1.8.28")
 * 4: uint32 distance from start of header to first data block.
 * 4: uint32 distance from start of block to index block (0). FIXME
 *
 * Returns 0 for success
 *        -1 for failure
 */
int srf_write_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch) {
    uint32_t sz = 0;

    if (!ch)
	return -1;

    /* Magic number && version */
    if (4 != fwrite(SRF_MAGIC, 1, 4, srf->fp))
	return -1;

    /* Header size */
    sz =  9
	+ (ch->version ? strlen(ch->version) : 0) + 1
	+ (ch->base_caller ? strlen(ch->base_caller) : 0) + 1
	+ (ch->base_caller_version ? strlen(ch->base_caller_version) : 0) + 1;
    if (0 != srf_write_uint32(srf, sz))
	return -1;

    if (srf_write_pstring(srf, ch->version) < 0)
	return -1;
    
    /* Containter type, base caller bits */
    if (EOF == fputc(ch->container_type, srf->fp))
	return -1;
    if (srf_write_pstring(srf, ch->base_caller) < 0)
	return -1;
    if (srf_write_pstring(srf, ch->base_caller_version) < 0)
	return -1;

    return ferror(srf->fp) ? -1 : 0;
}

/*
 * Reads an XML TraceInfo block
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_xml(srf_t *srf, srf_xml_t *xml) {
    int block_type;

    if (EOF == (block_type = fgetc(srf->fp)))
	return -1;
    if (block_type != SRFB_XML)
	return -1;

    if (0 != srf_read_uint32(srf, &xml->xml_len))
	return -1;
    xml->xml_len -= 5;

    if (NULL == (xml->xml = (char *)realloc(xml->xml, xml->xml_len+1)))
	return -1;
    if (xml->xml_len != fread(xml->xml, 1, xml->xml_len, srf->fp))
	return -1;
    xml->xml[xml->xml_len] = 0;

    return 0;
}

/*
 * Writes an XML TraceInfo block
 * Returns 0 for success
 *        -1 for failure
 */
int srf_write_xml(srf_t *srf, srf_xml_t *xml) {
    if (!srf->fp)
	return -1;

    if (EOF == fputc(SRFB_XML, srf->fp))
	return -1;

    if (-1 == srf_write_uint32(srf, xml->xml_len+5))
	return -1;

    if (xml->xml_len != fwrite(xml->xml, 1, xml->xml_len, srf->fp))
	return -1;

    return ferror(srf->fp) ? -1 : 0;
}

/*
 * Initialises a srf_trace_header_t and inserts some passed in values.
 * If the supplied th is NULL then a new structure is allocated.
 *
 * Returns a pointer to the dh passed in (or allocated if NULL) on success
 *         NULL on failure
 */
srf_trace_hdr_t *srf_construct_trace_hdr(srf_trace_hdr_t *th,
					 char *prefix,
					 unsigned char *header,
					 uint32_t header_sz) {
    if (!th) {
	if (NULL == (th = (srf_trace_hdr_t *)calloc(1, sizeof(*th))))
	    return NULL;
    }

    th->block_type = SRFB_TRACE_HEADER;
    strncpy(th->id_prefix, prefix, 255);
    th->trace_hdr_size = header_sz;
    th->trace_hdr = header;
    th->read_prefix_type = 'E';

    return th;
}

/*
 * Deallocates a srf_trace_hdr_t if allocated by
 * srf_construct_trace_hdr().
 * Do not use this if you passed in a static srf_trace_hdr to the construct
 * function.
 */
void srf_destroy_trace_hdr(srf_trace_hdr_t *th) {
    if (th) {
	if (th->trace_hdr)
	    free(th->trace_hdr);
	free(th);
    }
}

/*
 * Reads a data header and stores the result in 'th'.
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_trace_hdr(srf_t *srf, srf_trace_hdr_t *th) {
    int z;

    /* Check block type */
    if (EOF == (th->block_type = fgetc(srf->fp)))
	return -1;

    if (th->block_type != SRFB_TRACE_HEADER)
	return -1;
    if (0 != srf_read_uint32(srf, &th->trace_hdr_size))
	return -1;
    th->trace_hdr_size -= 1 + 4 + 1;

    /* Read-id prefix */
    if (EOF == (th->read_prefix_type = fgetc(srf->fp)))
	return -1;
    if ((z = srf_read_pstring(srf, th->id_prefix)) < 0)
	return -1;
    th->trace_hdr_size -= z+1;

    /* The data header itself */
    if (th->trace_hdr)
	free(th->trace_hdr);
    if (th->trace_hdr_size) {
	if (NULL == (th->trace_hdr = malloc(th->trace_hdr_size)))
	    return -1;
	if (th->trace_hdr_size != fread(th->trace_hdr, 1,
					  th->trace_hdr_size, srf->fp)) {
	    free(th->trace_hdr);
	    th->trace_hdr = NULL;
	    return -1;
	}
    } else {
	th->trace_hdr = NULL;
    }

    return 0;
}

/*
 * Writes a srf_trace_hdr_t structure to disk.
 * Returns 0 for sucess
 *        -1 for failure
 */
int srf_write_trace_hdr(srf_t *srf, srf_trace_hdr_t *th) {
    uint32_t sz;

    if (!srf->fp)
	return -1;

    if (EOF == fputc(th->block_type, srf->fp))
	return -1;

    /* Size */
    sz = 1 + 4 + 1
	+ (th->id_prefix ? strlen(th->id_prefix) : 0) + 1
	+ th->trace_hdr_size;
    if (-1 == srf_write_uint32(srf, sz))
	return -1;

    /* Prefix */
    if (EOF == fputc(th->read_prefix_type, srf->fp))
	return -1;
    if (-1 == srf_write_pstring(srf, th->id_prefix))
	return -1;

    /* The ztr header blob itself... */
    if (th->trace_hdr_size !=
	fwrite(th->trace_hdr, 1, th->trace_hdr_size, srf->fp))
	return -1;

    return ferror(srf->fp) ? -1 : 0;
}


srf_trace_body_t *srf_construct_trace_body(srf_trace_body_t *tb,
					   char *suffix,
					   int suffix_len,
					   unsigned char *body,
					   uint32_t body_size,
					   unsigned char flags) {
    if (!tb) {
	if (NULL == (tb = (srf_trace_body_t *)calloc(1, sizeof(*tb))))
	    return NULL;
    }
    tb->block_type = SRFB_TRACE_BODY;
    if (suffix_len == -1) {
	suffix_len = strlen(suffix);
	if (suffix_len > 255)
	    suffix_len = 255;
    }
    memcpy(tb->read_id, suffix, suffix_len);
    tb->read_id[suffix_len] = 0;
    tb->read_id_length = suffix_len;

    tb->trace = body;
    tb->trace_size = body_size;
    tb->flags = flags;

    return tb;
}

void srf_destroy_trace_body(srf_trace_body_t *tb) {
    if (tb)
	free(tb);
}

/*
 * Writes a new trace body.
 *
 * Returns: 0 on success
 *          -1 on failure
 */
int srf_write_trace_body(srf_t *srf, srf_trace_body_t *tb) {
    uint32_t sz;

    if (!srf->fp)
	return -1;

    if (EOF == fputc(tb->block_type, srf->fp))
	return -1;

    /* Size */
    sz = 6 + tb->read_id_length+1 + tb->trace_size;
    if (0 != srf_write_uint32(srf, sz))
	return -1;

    /* Flags and name */
    if (EOF == (fputc(tb->flags, srf->fp)))
	return -1;
    if (-1 == srf_write_pstringb(srf, tb->read_id, tb->read_id_length))
	return -1;

    /* Tbe ztr footer blob itself... */
    if (tb->trace_size != fwrite(tb->trace, 1, tb->trace_size, srf->fp))
	return -1;

    return ferror(srf->fp) ? -1 : 0;
}

/*
 * Reads a trace header + trace 'blob' and stores the result in 'th'
 * If no_trace is true then it skips loading the trace data itself.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_trace_body(srf_t *srf, srf_trace_body_t *tb, int no_trace) {
    int z;

    /* Check block type */
    if (EOF == (tb->block_type = fgetc(srf->fp)))
	return -1;
    if (tb->block_type != SRFB_TRACE_BODY)
	return -1;

    /* Size */
    if (0 != srf_read_uint32(srf, &tb->trace_size))
	return -1;
    tb->trace_size -= 6;

    /* Flags */
    if (EOF == (z = fgetc(srf->fp)))
	return -1;
    tb->flags = z;

    /* Read-id suffix */
    if ((z = srf_read_pstring(srf, tb->read_id)) < 0)
	return -1;
    tb->read_id_length = z;
    tb->trace_size -= z+1;

    /* The trace data itself */
    if (!no_trace) {
	if (tb->trace_size) {
	    if (NULL == (tb->trace = malloc(tb->trace_size)))
		return -1;
	    if (tb->trace_size != fread(tb->trace, 1, tb->trace_size,
					srf->fp)) {
		free(tb->trace);
		tb->trace = NULL;
		return -1;
	    }
	} else {
	    tb->trace = NULL;
	}
    } else {
	/* Skip */
	fseeko(srf->fp, tb->trace_size, SEEK_CUR);
	tb->trace = NULL;
    }

    return 0;
}


/*
 * Reads a SRF index header. See srf_write_index_hdr for the format.
 * If no_seek is true it reads the header starting at the current file
 * offset, otherwise it seeks to the end of the file and reads that
 * header instead.
 *
 * Returns 0 on success and fills out *hdr
 *         -1 on failure
 */
int srf_read_index_hdr(srf_t *srf, srf_index_hdr_t *hdr, int no_seek) {
    int sz, z;

    /* Load footer */
    if (!no_seek) {
	if (0 != fseeko(srf->fp, -16, SEEK_END))
	    return -1;
	
	if (4 != fread(hdr->magic,   1, 4, srf->fp))
	    return -1;
	if (4 != fread(hdr->version, 1, 4, srf->fp))
	    return -1;
	if (0 != srf_read_uint64(srf, &hdr->size))
	    return -1;

	/* Check for validity */
        if (memcmp(hdr->magic,	 SRF_INDEX_MAGIC,   4) ||
	    memcmp(hdr->version, SRF_INDEX_VERSION, 4))
	    return -1;

	/* Seek to index header and re-read */
	if (0 != fseeko(srf->fp, -hdr->size, SEEK_END))
	    return -1;
    }

    if (4 != fread(hdr->magic,   1, 4, srf->fp))
	return -1;
    if (4 != fread(hdr->version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_read_uint64(srf, &hdr->size))
	return -1;

    /* Check once more */
    if (memcmp(hdr->magic,   SRF_INDEX_MAGIC,   4) ||
	memcmp(hdr->version, SRF_INDEX_VERSION, 4))
	return -1;

    /* And finally the remainder of the index header details */
    if (EOF == (hdr->index_type         = fgetc(srf->fp)))
	return -1;
    if (EOF == (hdr->dbh_pos_stored_sep = fgetc(srf->fp)))
	return -1;

    if (0 != srf_read_uint32(srf, &hdr->n_container))
	return -1;
    if (0 != srf_read_uint32(srf, &hdr->n_data_block_hdr))
	return -1;
    if (0 != srf_read_uint64(srf, &hdr->n_buckets))
	return -1;

    sz = 34; /* fixed size of the above records */

    if ((z = srf_read_pstring(srf, hdr->dbh_file)) < 0)
	return -1;
    sz += z+1;
    if ((z = srf_read_pstring(srf, hdr->cont_file)) < 0)
	return -1;
    sz += z+1;

    hdr->index_hdr_sz = sz;

    return 0;
}

/*
 * Writes a SRF index header.
 *
 * Header:
 *   x4    magic number, starting with 'I'.
 *   x4    version code (eg "1.00")
 *   x8    index size
 *   x1    index type ('E' normally)
 *   x1    dbh_pos_stored_sep (indicates if the item list contains the
 *         "data block header" index number).
 *   x4    number of containers
 *   x4    number of DBHs
 *   x8    number of hash buckets
 *
 *   x*    dbhFile  p-string (NULL if held within the same file)
 *   x*    contFile p-string (NULL if held within the same file)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int srf_write_index_hdr(srf_t *srf, srf_index_hdr_t *hdr) {
    if (4 != fwrite(hdr->magic,   1, 4, srf->fp))
	return -1;
    if (4 != fwrite(hdr->version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_write_uint64(srf, hdr->size))
	return -1;

    if (EOF == fputc(hdr->index_type, srf->fp))
	return -1;
    if (EOF == fputc(hdr->dbh_pos_stored_sep, srf->fp))
	return -1;

    if (0 != srf_write_uint32(srf, hdr->n_container))
	return -1;
    if (0 != srf_write_uint32(srf, hdr->n_data_block_hdr))
	return -1;
    if (0 != srf_write_uint64(srf, hdr->n_buckets))
	return -1;

    if (-1 == srf_write_pstring(srf, hdr->dbh_file))
	return -1;
    if (-1 == srf_write_pstring(srf, hdr->cont_file))
	return -1;

    return ferror(srf->fp) ? -1 : 0;
}

/* Position in index - internal struct used for code below only */
typedef struct {
    uint64_t pos;
    uint32_t dbh;
} pos_dbh;

/*
 * This allocates and initialises an srf_index_t struct filling out the
 * fields to default values or the supplied parameters. It does not
 * actually write anything to disc itself.
 *
 * Note: non-NULL values for ch_file and th_file are not implemented yet.
 *
 * ch_file is the container header file. If specified as non-NULL it is the 
 * name of the file storing the container and DB records (trace bodies).
 * NULL implies all the data is in the same file.
 *
 * th_file is the filename where we store the DBH records (trace headers).
 * NULL implies all the data is in the same file.
 *
 * dbh_sep is a boolean value used to indicate whether we store the
 * location of DBH+DB per trace or just the DB record. The latter uses less
 * space and is most generally used, but the former is required if DBH and DB
 * are split apart into two files (ch_file and th_file).
 *
 * Returns srf_index_t pointer on success.
 *         NULL on failure
 */
srf_index_t *srf_index_create(char *ch_file, char *th_file, int dbh_sep) {
    srf_index_t *idx = (srf_index_t *)malloc(sizeof(srf_index_t));
    if (!idx)
	return NULL;

    if (ch_file) {
	strncpy(idx->ch_file, ch_file, PATH_MAX);
	idx->ch_file[PATH_MAX] = 0;
    } else {
	idx->ch_file[0] = 0;
    }

    if (th_file) {
	strncpy(idx->th_file, th_file, PATH_MAX);
	idx->th_file[PATH_MAX] = 0;
    } else {
	idx->th_file[0] = 0;
    }

    idx->dbh_pos_stored_sep = dbh_sep;

    /* Create the arrays and hash table */
    if (!(idx->ch_pos = ArrayCreate(sizeof(uint64_t), 0)))
	return NULL;

    if (!(idx->th_pos = ArrayCreate(sizeof(uint64_t), 0)))
	return NULL;

    if (!(idx->name_blocks = ArrayCreate(sizeof(srf_name_block_t), 0)))
        return NULL;

    if (!(idx->db_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE |
					    HASH_FUNC_JENKINS3 |
					    HASH_NONVOLATILE_KEYS |
					    HASH_POOL_ITEMS)))
	return NULL;

    return idx;
}


/*
 * Deallocates memory used by an srf_index_t structure.
 */
void srf_index_destroy(srf_index_t *idx) {
    size_t i;

    if (!idx)
	return;

    if (idx->db_hash)
	HashTableDestroy(idx->db_hash, 0);
    if (idx->ch_pos)
	ArrayDestroy(idx->ch_pos);
    if (idx->th_pos)
	ArrayDestroy(idx->th_pos);
    if (idx->name_blocks) {
        for (i = 0; i < ArrayMax(idx->name_blocks); i++) {
	    if (NULL != arr(srf_name_block_t, idx->name_blocks, i).names)
	        free(arr(srf_name_block_t, idx->name_blocks, i).names);
        }
	ArrayDestroy(idx->name_blocks);
    }

    free(idx);
}


/*
 * Dumps out some statistics on the index to fp, or stderr if fp
 * is NULL.
 */
void srf_index_stats(srf_index_t *idx, FILE *fp) {
    HashTableStats(idx->db_hash, fp ? fp : stderr);
}


/*
 * Adds a container header (CH block) to the srf index.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int srf_index_add_cont_hdr(srf_index_t *idx, uint64_t pos) {
    uint64_t *ip = ARRP(uint64_t, idx->ch_pos, ArrayMax(idx->ch_pos));
    return ip ? (*ip = pos, 0) : -1;
}


/*
 * Adds a trace header (DBH block) to the srf index.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int srf_index_add_trace_hdr(srf_index_t *idx, uint64_t pos) {
    uint64_t *ip = ARRP(uint64_t, idx->th_pos, ArrayMax(idx->th_pos));
    return ip ? (*ip = pos, 0) : -1;
}


/*
 * Adds a trace body (DB block) to the srf index.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int srf_index_add_trace_body(srf_index_t *idx, char *name, uint64_t pos) {
    HashData hd;
    pos_dbh *pdbh;
    srf_name_block_t *blockp;
    char *name_copy;
    size_t name_len;
    int new;

    if (idx->dbh_pos_stored_sep) {
	if (NULL == (pdbh = (pos_dbh *)malloc(sizeof(*pdbh))))
	    return -1;
	pdbh->pos = pos;
	pdbh->dbh = ArrayMax(idx->th_pos);
	hd.p = pdbh;
    } else {
	hd.i = pos;
    }

    name_len = strlen(name) + 1; /* Include NULL */

    /* Allocate more space for names if needed */
    if (ArrayMax(idx->name_blocks) == 0
	|| arr(srf_name_block_t, idx->name_blocks,
	       ArrayMax(idx->name_blocks) -  1).space <= name_len) {
        blockp = ARRP(srf_name_block_t, idx->name_blocks,
		      ArrayMax(idx->name_blocks));
	if (NULL == blockp) return -1;

	blockp->used = 0;
	blockp->space = (name_len < SRF_INDEX_NAME_BLOCK_SIZE
			 ? SRF_INDEX_NAME_BLOCK_SIZE
			 : name_len);
	blockp->names = xmalloc(blockp->space);
	if (NULL == blockp->names) {
	    ArrayMax(idx->name_blocks)--;
	    return -1;
	}	
    }

    blockp = ARRP(srf_name_block_t, idx->name_blocks,
		  ArrayMax(idx->name_blocks) - 1);
    name_copy = blockp->names + blockp->used;
    memcpy(name_copy, name, name_len);
    blockp->used  += name_len;
    blockp->space -= name_len;

    if (NULL == HashTableAdd(idx->db_hash, name_copy, name_len - 1, hd, &new)){
        return -1;
    }
    if (0 == new) {
	fprintf(stderr, "duplicate read name %s\n", name);
	return -1;
    }

    return 0;
}


/*
 * Writes the HashTable structures to 'fp'.
 * This is a specialisation of the HashTable where the HashData is a
 * position,size tuple.
 *
 * Header:
 *   x4    magic number, starting with 'I'.
 *   x4    version code (eg "1.00")
 *   x8    index size (should be x4 as we assume bucket locs are x4?)
 *
 *   x1    index type ('E' normally)
 *   x1    dbh_pos_stored_sep (indicates if the item list contains the
 *         "data block header" index number).
 *
 *   x4    number of containers
 *   x4    number of DBHs
 *   x8    number of hash buckets
 *
 *   x*    dbhFile  p-string (NULL if held within the same file)
 *   x*    contFile p-string (NULL if held within the same file)
 *
 * Containers: (1 entry per container)
 *   x8    file position of container header
 *
 * Data Block Headers: (1 entry per DBH)
 *   x8    file position of container header
 *
 * Buckets: (1 entry per bucket)
 *   x8    8-byte offset of linked list pos,  rel. to the start of the hdr
 *
 * Items: (1 per trace)
 *   x1    name disambiguation hash, top-most bit set => last item in list
 *   x8    data position
 *  (x4)  (dbh_index - optional; present if dbh_pos_stored_sep is 1)
 *
 * Footer:
 *   x4    magic number
 *   x4    version
 *   x8    index size
 *
 * Returns: the number of bytes written on success
 *         -1 for error
 */
int srf_index_write(srf_t *srf, srf_index_t *idx) {
    unsigned int i, j;
    srf_index_hdr_t hdr;
    uint64_t *bucket_pos;
    int item_sz;
    HashTable *h = idx->db_hash;

    /* Option: whether to store dbh positions directly in the index */
    hdr.dbh_pos_stored_sep = idx->dbh_pos_stored_sep;

    /* Compute index size and bucket offsets */
    hdr.size = 34 +
	1 + (idx->ch_file ? strlen(idx->ch_file) : 0) +
	1 + (idx->th_file ? strlen(idx->th_file) : 0);
    hdr.size += 8*(ArrayMax(idx->ch_pos) +
		   ArrayMax(idx->th_pos) +
		   h->nbuckets);
    if (NULL == (bucket_pos = (uint64_t *)calloc(h->nbuckets,
						 sizeof(*bucket_pos))))
	return -1;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	if (!(hi = h->bucket[i]))
	    continue;
	bucket_pos[i] = hdr.size;
	item_sz = 1 + 8 + (hdr.dbh_pos_stored_sep ? 4 : 0);
	for (; hi; hi = hi->next)
	    hdr.size += item_sz;
    }
    hdr.size += 16; /* footer */

    /* Construct and write out the index header */
    memcpy(hdr.magic,   SRF_INDEX_MAGIC,   4);
    memcpy(hdr.version, SRF_INDEX_VERSION, 4);
    hdr.index_type = 'E';
    hdr.n_container = ArrayMax(idx->ch_pos);
    hdr.n_data_block_hdr = ArrayMax(idx->th_pos);
    hdr.n_buckets = h->nbuckets;
    if (idx->th_file)
	strncpy(hdr.dbh_file,  idx->th_file, 255);
    else
	hdr.dbh_file[0] = 0;
    if (idx->ch_file)
	strncpy(hdr.cont_file, idx->ch_file, 255);
    else
	hdr.cont_file[0] = 0;
    if (0 != srf_write_index_hdr(srf, &hdr))
	return -1;

    /* Write the container and data block header arrays */
    j = ArrayMax(idx->ch_pos);
    for (i = 0; i < j; i++) {
	if (0 != srf_write_uint64(srf, arr(uint64_t, idx->ch_pos, i)))
	    return -1;
    }

    j = ArrayMax(idx->th_pos);
    for (i = 0; i < j; i++) {
	if (0 != srf_write_uint64(srf, arr(uint64_t, idx->th_pos, i)))
	    return -1;
    }

    /* Write out buckets */
    for (i = 0; i < h->nbuckets; i++) {
	if (0 != srf_write_uint64(srf, bucket_pos[i]))
	    return -1;
    }

    /* Write out the trace locations themselves */
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	/*
	fprintf(stderr, "Bucket %d offset %lld vs %lld\n",
		i, ftell(srf->fp), bucket_pos[i]);
	*/
	if (!(hi = h->bucket[i]))
	    continue;
	for (; hi; hi = hi->next) {
	    uint64_t pos;
	    uint32_t dbh = 0;
	    uint32_t h7;

	    if (hdr.dbh_pos_stored_sep) {
		pos_dbh *pdbh = (pos_dbh *)hi->data.p;
		pos = pdbh->pos;
		dbh = pdbh->dbh;
	    } else {
		pos = hi->data.i;
	    }

	    /* Rehash key in 7 bits;  */
	    h7 = hash64(h->options & HASH_FUNC_MASK,
			(uint8_t *)hi->key, hi->key_len) >> 57;
	    if (!hi->next)
		h7 |= 0x80;
	    /*
	    fprintf(stderr, "\t%.*s => %x @ %lld\n",
		    hi->key_len, hi->key, h7, pos);
	    */
	    if (fputc(h7, srf->fp) < 0)
		return -1;
	    if (0 != srf_write_uint64(srf, pos))
		return -1;

	    if (hdr.dbh_pos_stored_sep)
		if (0 != srf_write_uint32(srf, dbh))
		    return -1;
	}
    }

    /* Footer */
    if (4 != fwrite(hdr.magic,   1, 4, srf->fp))
	return -1;
    if (4 != fwrite(hdr.version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_write_uint64(srf, hdr.size))
	return -1;

    free(bucket_pos);

    return 0;
}

/*
 * ---------------------------------------------------------------------------
 * Trace name codec details.
 */

/*
 * Reads up to 32-bits worth of data and returns. Updates the block
 * byte and bit values to indicate the current 'read' position.
 *
 * Returns unsigned value on success (>=0)
 *         -1 on failure
 */
static uint32_t get_hi_bits(block_t *block, int nbits) {
    unsigned int val, bnum = 0;

    if (block->byte*8 + block->bit + nbits > block->alloc * 8)
        return -1;

    /* Fetch the partial byte of data */
    val = (block->data[block->byte]) & (1<<(8-block->bit))-1;
    bnum = 8 - block->bit;

    if (bnum >= nbits) {
	val >>= bnum-nbits;
	val &= (1<<nbits)-1;
	block->bit += nbits;
	return val;
    }

    /* And additional entire bytes worth as required */
    while (bnum+8 <= nbits && bnum+8 < 32) {
	val <<= 8;
        val |= block->data[++block->byte];
        bnum += 8;
    }

    /* The remaining partial byte */
    val <<= nbits-bnum;
    val |= (block->data[++block->byte] >> (8-(nbits-bnum)))
	& (1<<(nbits-bnum))-1;
    block->bit = nbits-bnum;

    return val;
}

/*
 * Stores up to 32-bits of data in a block
 */
static void set_hi_bits(block_t *block, uint32_t val, int nbits) {
    unsigned int curr = block->data[block->byte];
    
    /* Pack first partial byte */
    if (nbits > 8-block->bit) {
	curr |= val >> (nbits -= 8-block->bit);
	block->data[block->byte] = curr;
	block->data[++block->byte] = 0;
	block->bit = 0;
    } else {
	curr |= val << (8-block->bit - nbits);
	block->data[block->byte] = curr;
	if ((block->bit += nbits) == 8) {
	    block->bit = 0;
	    block->data[++block->byte] = 0;
	}
	return;
    }

    /* Handle whole bytes worth */
    while (nbits > 8) {
	block->data[block->byte++] = (val >> (nbits-=8)) & 0xff;
    }

    /* And finally any remaining bits left */
    block->data[block->byte] = (val & ((1<<nbits) - 1)) << (8-nbits);
    block->bit = nbits;
}


/*
 * Formats are specified embedded in 'fmt' using a percent-rule, much
 * like printf().
 *
 * Both fmt and suffix are C-style nul terminated strings.
 *
 * The format consists of:
 *
 * '%' <field-width> <bits-used> <format-code>
 *
 * Field-width is a numerical value indicating the number of characters we
 * wish to print. It is optional as without specifying this we emit as
 * many characters as are needed to describe the data. If specified the
 * output is padded to be at least field-width in size. The padding character
 * may vary on the otuput format, but will typically be '0'.
 *
 * Bits-used consists of '.' (a full stop) followed by a numerical value
 * indicating the number of bits to read from the suffix, starting from
 * bit 0 or the next free bit following a previous format. If not specified
 * generally all bits are used (8 * suffix_len) unless otherwise indicated
 * below.
 * 
 * Format-code may be one of:
 *    'd' decimal values (0-9)
 *    'o' octal values (0-7)
 *    'x' hexidecimal values, lowercase
 *    'X' hexidecimal values, uppercase
 *    'j' base-36 encoding, lowercase
 *    'J' base-36 encoding, uppercase (454)
 *    'c' a single character (default bits used = 8)
 *    's' string (all bits used, treated as ascii)
 *    '%' a literal percent character (no bits used).
 *
 * Returns the number of bytes written to 'name' on success
 *         -1 on failure
 */
#define emit(c) \
    if (out_pos < name_len-1) { \
        name[out_pos++] = (c); \
    } else { \
        block_destroy(blk, 1);\
        return name_len; \
    }

int construct_trace_name(char *fmt,
			 unsigned char *suffix, int suffix_len,
			 char *name, int name_len) {
    block_t *blk = block_create(suffix, suffix_len);
    int out_pos = 0;
    int percent = 0;

    /* Default nul-terminate for abort cases */
    name[name_len-1] = '\0';

    for(; *fmt; fmt++) {
	switch(*fmt) {

	/* A format specifier */
	case '%': {
	    int width = 0;
	    int bits = 0;
	    uint32_t val;

	    fmt++;
	    percent++;

	    /* Width specifier */
	    if (0 == (width = strtol(fmt, &fmt, 10)))
		width = 1; /* minimum width */

	    /* Bit size specifier */
	    if ('.' == *fmt) {
		fmt++;
		bits = strtol(fmt, &fmt, 10);
	    }

	    /* The format code */
	    switch (*fmt) {
		int i;

	    case '%':
		for (i = 0; i < width; i++) {
		    emit('%');
		}
		break;

	    case 'o':
	    case 'd':
	    case 'x':
	    case 'X':
	    case 'j':
	    case 'J': {
		/* One of the integer encoding formats */
		char *digits = "0123456789abcdef";
		int d = 0, tmp_ind = 0;
		char tmp[1024];

		switch(*fmt) {
		case 'o': d=8;	break;
		case 'd': d=10; break;
		case 'x': d=16; break;

		case 'X':
		    d=16;
		    digits = "0123456789ABCDEF";
		    break;

		case 'j':
		    d=36;
		    digits = "abcdefghijklmnopqrstuvwxyz0123456789";
		    break;
		    
		case 'J': /* Used by 454 */
		    d=36;
		    digits = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
		    break;
		}

		while (bits > 0) {
		    int32_t sv;
		    int nb = bits > 32 ? 32 : bits;
		    if (-1 == (sv = get_hi_bits(blk, nb)))
			return -1;
		    val = sv;

		    do {
			tmp[tmp_ind++] = digits[val % d];
			val /= d;
		    } while (val != 0);

		    bits -= nb;
		}

		/* Pad to requested size */
		for (i = width; i > tmp_ind; i--)
		    emit(*digits);

		/* Output the formatted value itself */
		do {
		    emit(tmp[--tmp_ind]);
		} while (tmp_ind != 0);

		break;
	    }

	    case 'c':
		/* A single n-bit character */
		if (!bits)
		    bits = 8;
		if (-1 == (val = get_hi_bits(blk, bits)))
		    return -1;

		emit(val);
		break;

	    case 's': {
		/* A string, to the end of suffix, of n-bit characters */
		if (!bits)
		    bits = 8;
		/* Reading n-bits at a time to produce a string */
		while (-1 != (val = get_hi_bits(blk, bits)))
		    emit(val);
		break;
	    }

	    default:
		fprintf(stderr, "Unknown arg: %c\n", *fmt);
	    }
	    
	}

	case '\0':
	    break;

	default:
	    emit(*fmt);
	}
    }

    /*
     * No percent rule found implies the name is a simple string
     * concatenation of prefix and suffix
     */
    if (!percent) {
	int i;

	/* A strncpy would be more efficient here */
	for (i = 0; suffix[i]; i++) {
	    emit(suffix[i]);
	}
    }

    emit('\0');

    block_destroy(blk, 1);
    return out_pos;
}

/*
 * The opposite of above.
 * Given a format and a set of arguments packs the values into the
 * supplied 'suffix' string. This should be 256 characters long and the
 * first byte will consist of the real length.
 *
 * Returns 0 on success
 *        -1 on failuare
 */
int pack_trace_suffix(unsigned char *suffix, char *fmt, ...) {
    block_t *blk = block_create(NULL, 256);
    va_list args;

    va_start(args, fmt);

    for(; *fmt; fmt++) {
	switch(*fmt) {

	/* A format specifier */
	case '%': {
	    int width = 0;
	    int bits = 0;
	    signed int val;

	    fmt++;

	    /* Width specifier - ignore this */
	    width = strtol(fmt, &fmt, 10);

	    /* Bit size specifier */
	    if ('.' == *fmt) {
		fmt++;
		bits = strtol(fmt, &fmt, 10);
	    }

	    /* The format code */
	    switch (*fmt) {
	    case '%':
		/* A literal percent */
		break;

	    case 'o':
	    case 'd':
	    case 'x':
	    case 'X':
	    case 'j':
	    case 'J':
		/* A numeric value - it doesn't matter what format it
		 * is specified as, just how many bits.
		 */
		val = (uint32_t)va_arg(args, int);
		set_hi_bits(blk, val, bits);
		break;

	    case 'c':
		if (!bits)
		    bits = 8;
		val = (unsigned char)va_arg(args, int);
		set_hi_bits(blk, val, bits);
		break;

	    case 's': {
		char *s = (char *)va_arg(args, char *);
		if (!bits)
		    bits = 8;

		for (; *s; s++) {
		    set_hi_bits(blk, *s, bits);
		}
		break;
	    }

	    default:
		fprintf(stderr, "Unknown arg: %c\n", *fmt);
	    }
	}
	}
    }

    if (blk->byte >= 256)
	return -1;

    *suffix = blk->byte + (blk->bit > 0);
    memcpy(suffix+1, blk->data, *suffix);

    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * Higher level I/O functions
 */

/*
 * Fetches the next trace from an SRF container as a "memory-FILE".
 * Name, if defined (which should be a buffer of at least 512 bytes long)
 * will be filled out to contain the read name.
 * 
 * Returns mFILE containing trace on success
 *         NULL on failure.
 */
mFILE *srf_next_trace(srf_t *srf, char *name) {
    do {
	int type;

	switch(type = srf_next_block_type(srf)) {
	case -1:
	    /* EOF */
	    return NULL;

	case SRFB_NULL_INDEX: {
	    /*
	     * Maybe the last 8 bytes of a the file (or previously was
	     * last 8 bytes prior to concatenating SRF files together).
	     * If so it's the index length and should always be 8 zeros.
	     */
	    uint64_t ilen;
	    if (1 != fread(&ilen, 8, 1, srf->fp))
		return NULL;
	    if (ilen != 0)
		return NULL;
	    break;
	}

	case SRFB_CONTAINER:
	    if (0 != srf_read_cont_hdr(srf, &srf->ch))
		return NULL;
	    break;

	case SRFB_XML:
	    if (0 != srf_read_xml(srf, &srf->xml))
		return NULL;
	    break;

	case SRFB_TRACE_HEADER:
	    if (0 != srf_read_trace_hdr(srf, &srf->th))
		return NULL;

	    break;

	case SRFB_TRACE_BODY: {
	    mFILE *mf = mfcreate(NULL, 0);
	    srf_trace_body_t tb;
	    tb.trace = NULL;

	    if (!mf || 0 != srf_read_trace_body(srf, &tb, 0))
		return NULL;

	    if (name) {
		if (-1 == construct_trace_name(srf->th.id_prefix,
					       (unsigned char *)tb.read_id,
					       tb.read_id_length,
					       name, 512)) {
		    return NULL;
		}
	    }

	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, mf);
	    if (tb.trace_size)
		mfwrite(tb.trace, 1, tb.trace_size, mf);

	    if (tb.trace)
		free(tb.trace);

	    mrewind(mf);
	    return mf;
	}

	case SRFB_INDEX: {
	    off_t pos = ftell(srf->fp);
	    srf_read_index_hdr(srf, &srf->hdr, 1);

	    /* Skip the index body */
	    fseeko(srf->fp, pos + srf->hdr.size, SEEK_SET);
	    break;
	}

	default:
	    fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	    return NULL;
	}
    } while (1);

    return NULL;
}

/*
 * Decodes a partial ZTR file consisting of data in 'mf'.
 * Note that mf may contain a partial chunk, so we need to be careful on
 * error checking.
 *
 * If a ztr object is passed in (in 'z') then we assume we've already
 * loaded the ZTR header and get straight down to decoding the remaining
 * chunks. Otherwise we also decode the header.
 *
 * If no chunk is visible at all then we'll return NULL and rewind mf.
 * Otherwise we'll leave the file pointer at the start of the next 
 * partial chunk (or EOF if none) and return the ztr_t pointer.
 */
ztr_t *partial_decode_ztr(srf_t *srf, mFILE *mf, ztr_t *z) {
    ztr_t *ztr;
    ztr_chunk_t *chunk;
    long pos = 0;

    if (z) {
	/* Use existing ZTR object => already loaded header */
	ztr = z;

    } else {
	/* Allocate or use existing ztr */
	if (NULL == (ztr = new_ztr()))
	    return NULL;

	/* Read the header */
	if (-1 == ztr_read_header(mf, &ztr->header)) {
	    if (!z)
		delete_ztr(ztr);
	    mrewind(mf);
	    return NULL;
	}

	/* Check magic number and version */
	if (memcmp(ztr->header.magic, ZTR_MAGIC, 8) != 0) {
	    if (!z)
		delete_ztr(ztr);
	    mrewind(mf);
	    return NULL;
	}

	if (ztr->header.version_major != ZTR_VERSION_MAJOR) {
	    if (!z)
		delete_ztr(ztr);
	    mrewind(mf);
	    return NULL;
	}
    }

    /* Load chunks */
    pos = mftell(mf);
    while (chunk = ztr_read_chunk_hdr(mf)) {
	chunk->data = (char *)xmalloc(chunk->dlength);
	if (chunk->dlength != mfread(chunk->data, 1, chunk->dlength, mf))
	    break;

	ztr->nchunks++;
	ztr->chunk = (ztr_chunk_t *)xrealloc(ztr->chunk, ztr->nchunks *
					     sizeof(ztr_chunk_t));
	memcpy(&ztr->chunk[ztr->nchunks-1], chunk, sizeof(*chunk));
	xfree(chunk);
	pos = mftell(mf);
    }

    /*
     * At this stage we're 'pos' into the mFILE mf with any remainder being
     * a partial block.
     */
    if (0 == ztr->nchunks) {
	if (!z)
	    delete_ztr(ztr);
	mrewind(mf);
	return NULL;
    }

    /* Ensure we exit at the start of a ztr CHUNK */
    mfseek(mf, pos, SEEK_SET);

    /* If this is the header part, ensure we uncompress and init. data */
    if (!z) {
	/* Force caching of huffman code_sets */
	ztr_find_hcode(ztr, CODE_USER);

	/* And uncompress the rest */
	uncompress_ztr(ztr);
    }

    return ztr;
}

/*
 * Creates a copy of ztr_t 'src' and returns it. The newly returned ztr_t
 * will consist of shared components where src and dest overlap, but freeing
 * dest will know what's appropriate to free and what is not.
 */
ztr_t *ztr_dup(ztr_t *src) {
    ztr_t *dest = new_ztr();
    int i;

    if (!dest)
	return NULL;

    /* Basics */
    *dest = *src;

    /* Mirror chunks */
    dest->chunk = (ztr_chunk_t *)malloc(src->nchunks * sizeof(ztr_chunk_t));
    for (i = 0; i < src->nchunks; i++) {
	dest->chunk[i] = src->chunk[i];
	dest->chunk[i].ztr_owns = 0; /* src owns the data/meta_data */
    }

    /* Mirror text_segments; no overlap here */
    dest->text_segments = (ztr_text_t *)malloc(src->ntext_segments *
					       sizeof(ztr_text_t));
    for (i = 0; i < src->ntext_segments; i++) {
	dest->text_segments[i] = src->text_segments[i];
    }

    /* huffman hcodes */
    dest->hcodes = (ztr_hcode_t *)malloc(src->nhcodes * sizeof(ztr_hcode_t));
    for (i = 0; i < src->nhcodes; i++) {
	dest->hcodes[i] = src->hcodes[i];
	dest->hcodes[i].ztr_owns = 0;
    }

    return dest;
}

/*
 * Fetches the next trace from an SRF container as a ZTR object.
 * This is more efficient than srf_next_trace() if we are serially
 * reading through many traces as we decode ZTR data less often and can
 * cache data from one trace to the next.
 *
 * Name, if defined (which should be a buffer of at least 512 bytes long)
 * will be filled out to contain the read name.
 *
 * filter_mask should consist of zero or more SRF_READ_FLAG_* bits.
 * Reads with one or more flags matching these bits will be skipped over.
 *
 * flags, if non-NULL, is filled out on exit to contain the SRF flags
 * from the Data Block structure.
 *
 * Returns ztr_t * on success
 *         NULL on failure.
 */
ztr_t *srf_next_ztr_flags(srf_t *srf, char *name, int filter_mask,
			  int *flags) {
    do {
	int type;

	switch(type = srf_next_block_type(srf)) {
	case -1:
	    /* EOF */
	    return NULL;

	case SRFB_NULL_INDEX: {
	    /*
	     * Maybe the last 8 bytes of a the file (or previously was
	     * last 8 bytes prior to concatenating SRF files together).
	     * If so it's the index length and should always be 8 zeros.
	     */
	    uint64_t ilen;
	    if (1 != fread(&ilen, 8, 1, srf->fp))
		return NULL;
	    if (ilen != 0)
		return NULL;
	    break;
	}

	case SRFB_CONTAINER:
	    if (0 != srf_read_cont_hdr(srf, &srf->ch))
		return NULL;
	    break;

	case SRFB_XML:
	    if (0 != srf_read_xml(srf, &srf->xml))
		return NULL;
	    break;

	case SRFB_TRACE_HEADER:
	    if (0 != srf_read_trace_hdr(srf, &srf->th))
		return NULL;

	    /* Decode ZTR chunks in the header */
	    if (srf->mf)
		mfdestroy(srf->mf);

	    srf->mf = mfcreate(NULL, 0);
	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, srf->mf);
	    if (srf->ztr)
		delete_ztr(srf->ztr);
	    mrewind(srf->mf);

	    if (NULL != (srf->ztr = partial_decode_ztr(srf, srf->mf, NULL))) {
		srf->mf_pos = mftell(srf->mf);
	    } else {
		/* Maybe not enough to decode or no headerBlob. */
		/* So delay until decoding the body. */
		srf->mf_pos = 0;
	    }
	    mfseek(srf->mf, 0, SEEK_END);
	    srf->mf_end = mftell(srf->mf);

	    break;

	case SRFB_TRACE_BODY: {
	    srf_trace_body_t tb;
	    ztr_t *ztr_tmp;

	    if (!srf->mf || 0 != srf_read_trace_body(srf, &tb, 0))
		return NULL;

	    if (name) {
		if (-1 == construct_trace_name(srf->th.id_prefix,
					       (unsigned char *)tb.read_id,
					       tb.read_id_length,
					       name, 512)) {
		    return NULL;
		}
	    }

	    mfseek(srf->mf, srf->mf_end, SEEK_SET);
	    if (tb.trace_size) {
		mfwrite(tb.trace, 1, tb.trace_size, srf->mf);
		free(tb.trace);
		tb.trace = NULL;
	    }
	    mftruncate(srf->mf, mftell(srf->mf));
	    mfseek(srf->mf, srf->mf_pos, SEEK_SET);

	    if (tb.flags & filter_mask) {
		break; /* Filtered, so skip it */
	    } else {
		if (flags)
		    *flags = tb.flags;
		if (srf->ztr)
		    ztr_tmp = ztr_dup(srf->ztr); /* inefficient, but simple */
		else
		    ztr_tmp = NULL;

		return partial_decode_ztr(srf, srf->mf, ztr_tmp);
	    }
	}

	case SRFB_INDEX: {
	    off_t pos = ftell(srf->fp);
	    srf_read_index_hdr(srf, &srf->hdr, 1);

	    /* Skip the index body */
	    fseeko(srf->fp, pos + srf->hdr.size, SEEK_SET);
	    break;
	}

	default:
	    fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	    return NULL;
	}
    } while (1);

    return NULL;
}

ztr_t *srf_next_ztr(srf_t *srf, char *name, int filter_mask) {
    return srf_next_ztr_flags(srf, name, filter_mask, NULL);
}

/*
 * Returns the type of the next block.
 * -1 for none (EOF)
 */
int srf_next_block_type(srf_t *srf) {
    int c = fgetc(srf->fp);
    if (c == EOF)
	return -1;

    ungetc(c, srf->fp);

    return c;
}

/*
 * Reads the next SRF block from an archive and returns the block type.
 * If the block is a trace it'll return the full trace name too (maximum
 * 512 bytes).
 *
 * Returns block type on success, writing to pos and name as appropriate
 *         -1 on EOF
 *         -2 on failure
 */
int srf_next_block_details(srf_t *srf, uint64_t *pos, char *name) {
    int type;
    *pos = ftell(srf->fp);

    switch(type = srf_next_block_type(srf)) {
    case -1:
	/* EOF */
	return -1;

    case SRFB_NULL_INDEX: {
	/*
	 * Maybe the last 8 bytes of a the file (or previously was
	 * last 8 bytes prior to concatenating SRF files together).
	 * If so it's the index length and should always be 8 zeros.
	 */
	uint64_t ilen;
	if (1 != fread(&ilen, 8, 1, srf->fp))
	    return -2;
	if (ilen != 0)
	    return -2;
	break;
    }

    case SRFB_CONTAINER:
	if (0 != srf_read_cont_hdr(srf, &srf->ch))
	    return -2;
	break;

    case SRFB_TRACE_HEADER:
	if (0 != srf_read_trace_hdr(srf, &srf->th))
	    return -2;
	
	break;

    case SRFB_TRACE_BODY:
	/* Inefficient, but it'll do for testing purposes */
	if (0 != srf_read_trace_body(srf, &srf->tb, 1))
	    return -2;

	if (name) {
	     if (-1 == construct_trace_name(srf->th.id_prefix,
					    (unsigned char *)srf->tb.read_id,
					    srf->tb.read_id_length,
					    name, 512)) {
		 return -2;
	    }
	}
	
	break;

    case SRFB_INDEX:
	srf_read_index_hdr(srf, &srf->hdr, 1);

	/* Skip the index body */
	fseeko(srf->fp, *pos + srf->hdr.size, SEEK_SET);
	break;


    default:
	fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	return -2;
    }

    return type;
}

/*
 * Searches through 'nitems' 8-byte values stored in 'srf' at file offset
 * 'start' onwards for the closest value <= 'query'.
 *
 * Returns 0 on success, setting *res
 *        -1 on failure
 */
static int binary_scan(srf_t *srf, int nitems, uint64_t start, uint64_t query,
		       uint64_t *res) {
    int min = 0;
    int max = nitems;
    int guess, i;
    uint64_t pos = 0, best = 0;

    if (nitems <= 0)
	return -1;

    /* Binary search on disk for approx location */
    while (max - min > 100) {
	guess = (max - min) / 2 + min;

	if (guess == max)
	    guess = max-1;

	if (-1 == fseeko(srf->fp, guess * 8 + start, SEEK_SET))
	    return -1;
	if (0 != srf_read_uint64(srf, &pos))
	    return -1;
	if (pos > query) {
	    max = guess;
	} else {
	    min = guess;
	}
    }

    /* Within a small distance => linear scan now to avoid needless disk IO */
    if (-1 == fseeko(srf->fp, min * 8 + start, SEEK_SET))
	return -1;
    for (i = min; i < max; i++) {
	if (0 != srf_read_uint64(srf, &pos))
	    return -1;
	if (pos > query) {
	    break;
	} else {
	    best = pos;
	}
    }

    assert(best <= query);
    *res = best;

    return 0;
}

/*
 * Searches in an SRF index for a trace of a given name.
 * If found it sets the file offsets for the container (cpos), data block
 * header (hpos) and data block (dpos).
 *
 * On a test with 2 containers and 12 headers this averaged at 6.1 reads per
 * trace fetch and 8.0 seeks.
 *
 * Returns 0 on success
 *        -1 on failure (eg no index)
 *        -2 on trace not found in index.
 */
int srf_find_trace(srf_t *srf, char *tname,
		   uint64_t *cpos, uint64_t *hpos, uint64_t *dpos) {
    srf_index_hdr_t hdr;
    uint64_t hval, bnum;
    uint64_t bucket_pos;
    off_t ipos, skip;
    int item_sz = 8;

    /* Check for valid index */
    if (0 != srf_read_index_hdr(srf, &hdr, 0)) {
	return -1;
    }
    ipos = ftello(srf->fp);
    skip = hdr.n_container * 8 + hdr.n_data_block_hdr * 8;
    if (hdr.dbh_pos_stored_sep)
	item_sz += 4;

    /* Hash and load the bucket */
    hval = hash64(HASH_FUNC_JENKINS3, (unsigned char *)tname, strlen(tname));
    bnum = hval & (hdr.n_buckets - 1);
    if (-1 == fseeko(srf->fp, ipos + skip + bnum * 8, SEEK_SET))
	return -1;

    if (0 != srf_read_uint64(srf, &bucket_pos))
	return -1;
    if (!bucket_pos)
	return -2;

    /* Secondary hash is the top 7-bits */
    hval >>= 57;

    /* Jump to the item list */
    if (-1 == fseeko(srf->fp, ipos-hdr.index_hdr_sz + bucket_pos, SEEK_SET))
	return -1;
    for (;;) {
	char name[1024];
	int h = fgetc(srf->fp);
	off_t saved_pos;
	uint64_t dbh_ind = 0;
	
	if ((h & 0x7f) != hval) {
	    if (h & 0x80)
		return -2; /* end of list and not found */
	    /*
	     * fseeko(srf->fp, 8, SEEK_CUR);
	     * Use fread instead as it's likely already cached and linux
	     * fseeko involves a real system call (lseek).
	     */
	    fread(dpos, 1, item_sz, srf->fp);
	    continue;
	}

	/* Potential hit - investigate to see if it's the real one: */
	/* Seek to dpos and get trace id suffix. Compare to see if valid */
	if (0 != srf_read_uint64(srf, dpos))
	    return -1;
	if (hdr.dbh_pos_stored_sep) {
	    if (0 != srf_read_uint64(srf, &dbh_ind))
		return -1;
	}
	saved_pos = ftello(srf->fp);
	if (-1 == fseeko(srf->fp, (off_t)*dpos, SEEK_SET))
	    return -1;
	if (0 != srf_read_trace_body(srf, &srf->tb, 0))
	    return -1;

	/* Identify the matching hpos (trace header) for this trace body */
	if (hdr.dbh_pos_stored_sep) {
	    /* Hack for now - binary scan through 1 object */
	    if (0 != binary_scan(srf, 1,
				 ipos + hdr.n_container * 8 + dbh_ind * 8,
				 *dpos, hpos))
		return -1;
	} else {
	    if (0 != binary_scan(srf, hdr.n_data_block_hdr,
				 ipos + hdr.n_container * 8,
				 *dpos, hpos))
		return -1;
	}

	/* Check the trace name matches */
	if (-1 == fseeko(srf->fp, *hpos, SEEK_SET))
	    return -1;
	if (0 != srf_read_trace_hdr(srf, &srf->th))
	    return -1;

	if (-1 == construct_trace_name(srf->th.id_prefix,
				       (unsigned char *)srf->tb.read_id,
				       srf->tb.read_id_length,
				       name, 1024))
	    return -1;

	if (strcmp(name, tname)) {
	    /* Not found, continue with next item in list */
	    if (h & 0x80)
		return -2;
	    if (-1 == fseeko(srf->fp, saved_pos, SEEK_SET))
		return -1;
	    continue;
	}
	
	/* Matches, so fetch the container data and return out trace */
	if (0 != binary_scan(srf, hdr.n_container,
			     ipos, *dpos, cpos))
	    return -1;

	/* FIXME: what to do with base-caller and cpos */

	break;
    }

    return 0;
}
