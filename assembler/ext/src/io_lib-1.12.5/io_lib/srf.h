#ifndef _SRF_H_
#define _SRF_H_

#include "io_lib/hash_table.h"
#include "io_lib/ztr.h"
#include "io_lib/mFILE.h"

#define SRF_MAGIC		"SSRF"
#define SRF_VERSION             "1.3"

#define SRFB_CONTAINER 		'S'
#define SRFB_XML		'X'
#define SRFB_TRACE_HEADER	'H'
#define SRFB_TRACE_BODY		'R'
#define SRFB_INDEX		'I'

/* Lack of index => 8 zero bytes at end of file to indicate zero length */
#define SRFB_NULL_INDEX		'\0'

/*--- Public structures */

/* Container header - several per file */
typedef struct {
    int block_type;
    char version[256];
    char container_type;
    char base_caller[256];
    char base_caller_version[256];
} srf_cont_hdr_t;

/* Trace header - several per container */
typedef struct {
    int block_type; 
    char read_prefix_type;
    char id_prefix[256];
    uint32_t trace_hdr_size;
    unsigned char *trace_hdr;
} srf_trace_hdr_t;

/* Trace body - several per trace header */
typedef struct {
    int block_type;
    int read_id_length;
    char read_id[256];
    unsigned char flags;
    uint32_t trace_size;
    unsigned char *trace;
} srf_trace_body_t;

/* XML - NCBI TraceInfo data block */
typedef struct {
    uint32_t xml_len;
    char *xml;
} srf_xml_t;

#define SRF_READ_FLAG_BAD_MASK       (1<<0)
#define SRF_READ_FLAG_WITHDRAWN_MASK (1<<1)
#define SRF_READ_FLAG_USER_MASK      (7<<5)

/* Indexing */
typedef struct {
    char     magic[4];
    char     version[4];
    uint64_t size;
    uint32_t n_container;
    uint32_t n_data_block_hdr;
    uint64_t n_buckets;
    int8_t   index_type;
    int8_t   dbh_pos_stored_sep;
    char     dbh_file[256];
    char     cont_file[256];
    int      index_hdr_sz; /* size of the above data on disk */
} srf_index_hdr_t;

/* In-memory index itself */
#define SRF_INDEX_NAME_BLOCK_SIZE 10000000

typedef struct {
  size_t  used;
  size_t  space;
  char   *names;
} srf_name_block_t;

typedef struct {
    char ch_file[PATH_MAX+1];
    char th_file[PATH_MAX+1];
    Array ch_pos;
    Array th_pos;
    Array name_blocks;
    int dbh_pos_stored_sep;
    HashTable *db_hash;
} srf_index_t;

/* Master SRF object */
typedef struct {
    FILE *fp;

    /* Cached copies of each of the most recent chunk types loaded */
    srf_cont_hdr_t    ch;
    srf_trace_hdr_t   th;
    srf_trace_body_t  tb;
    srf_xml_t         xml;
    srf_index_hdr_t   hdr;

    /* Private: cached data for use by srf_next_ztr */
    ztr_t *ztr;
    mFILE *mf;
    long mf_pos, mf_end;
} srf_t;

#define SRF_INDEX_MAGIC    "Ihsh"
#define SRF_INDEX_VERSION  "1.01"


/*--- Initialisation */
srf_t *srf_create(FILE *fp);
srf_t *srf_open(char *fn, char *mode);
void srf_destroy(srf_t *srf, int auto_close);

/*--- Base type I/O methods */

int srf_write_pstring(srf_t *srf, char *str);
int srf_write_pstringb(srf_t *srf, char *str, int length);
int srf_read_pstring(srf_t *srf, char *str);

int srf_read_uint32(srf_t *srf, uint32_t *val);
int srf_write_uint32(srf_t *srf, uint32_t val);

int srf_read_uint64(srf_t *srf, uint64_t *val);
int srf_write_uint64(srf_t *srf, uint64_t val);


/*--- Mid level I/O - srf block */
srf_cont_hdr_t *srf_construct_cont_hdr(srf_cont_hdr_t *ch,
				       char *bc,
				       char *bc_version);
void srf_destroy_cont_hdr(srf_cont_hdr_t *ch);
int srf_read_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch);
int srf_write_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch);

int srf_read_xml(srf_t *srf, srf_xml_t *xml);
int srf_write_xml(srf_t *srf, srf_xml_t *xml);

srf_trace_hdr_t *srf_construct_trace_hdr(srf_trace_hdr_t *th,
					 char *prefix,
					 unsigned char *header,
					 uint32_t header_sz);
void srf_destroy_trace_hdr(srf_trace_hdr_t *th);
int srf_read_trace_hdr(srf_t *srf, srf_trace_hdr_t *th);
int srf_write_trace_hdr(srf_t *srf, srf_trace_hdr_t *th);

srf_trace_body_t *srf_construct_trace_body(srf_trace_body_t *th,
					   char *suffix,
					   int suffix_len,
					   unsigned char *body,
					   uint32_t body_size,
					   unsigned char flags);
void srf_destroy_trace_body(srf_trace_body_t *th);
int srf_write_trace_body(srf_t *srf, srf_trace_body_t *th);
int srf_read_trace_body(srf_t *srf, srf_trace_body_t *th, int no_trace);

int srf_read_index_hdr(srf_t *srf, srf_index_hdr_t *hdr, int no_seek);
int srf_write_index_hdr(srf_t *srf, srf_index_hdr_t *hdr);
srf_index_t *srf_index_create(char *ch_file, char *th_file, int dbh_sep);
void srf_index_destroy(srf_index_t *idx);
void srf_index_stats(srf_index_t *idx, FILE *fp);
int srf_index_add_cont_hdr(srf_index_t *idx, uint64_t pos);
int srf_index_add_trace_hdr(srf_index_t *idx, uint64_t pos);
int srf_index_add_trace_body(srf_index_t *idx, char *name, uint64_t pos);
int srf_index_write(srf_t *srf, srf_index_t *idx);

/*--- Higher level I/O functions */
mFILE *srf_next_trace(srf_t *srf, char *name);
ztr_t *srf_next_ztr_flags(srf_t *srf, char *name, int filter_mask, int *flags);
ztr_t *srf_next_ztr(srf_t *srf, char *name, int filter_mask);

ztr_t *partial_decode_ztr(srf_t *srf, mFILE *mf, ztr_t *z);
ztr_t *ztr_dup(ztr_t *src);

int srf_next_block_type(srf_t *srf); /* peek ahead */
int srf_next_block_details(srf_t *srf, uint64_t *pos, char *name);

int srf_find_trace(srf_t *srf, char *trace,
		   uint64_t *cpos, uint64_t *hpos, uint64_t *dpos);

int construct_trace_name(char *fmt,
			 unsigned char *suffix, int suffix_len,
			 char *name, int name_len);

#endif /* _SRF_H_ */
