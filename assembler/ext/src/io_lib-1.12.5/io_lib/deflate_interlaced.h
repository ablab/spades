#ifndef _DEFLATE_SIMPLE_H_
#define _DEFLATE_SIMPLE_H_

#ifdef __cplusplus
extern "C" {
#endif

/* inlined codes */
#define CODE_INLINE	  0   

/* predefined codes */
#define CODE_DNA	  1   /* DNA, uppercase only */
#define CODE_DNA_AMBIG	  2   /* DNA, uc with ambiguity codes */
#define CODE_ENGLISH      3   /* English text */
#define NCODES_STATIC     4   /* Cheat, but we count from 0 for ease */

/* predefined elsewhere in HUFF chunks, 128 onwards */
#define CODE_USER	  128 

#define MAX_CODE_LEN	  15 /* maximum allowed by RFC 1951 */

#ifndef ZTR_FORM_STHUFF
#  define ZTR_FORM_STHUFF  77
#endif

/* A single symbol and it's encoding */
typedef struct {
    signed int symbol; /* 0-255 character, 256 = exception code, 257 = EOF */
    int nbits;
    unsigned int code;
    int freq;
} huffman_code_t;

/* A collection of huffman_code_t along with decoding optimisations */
typedef struct {
    huffman_code_t *codes;
    int ncodes;
    int codes_static;
    huffman_code_t lookup[258]; /* Mapping of symbol character to code */
    int max_code_len;
} huffman_codes_t;

/* Use for store_bits() and get_bits() */
typedef struct block {
    unsigned char *data;
    size_t alloc;
    size_t byte;
    int bit;
} block_t;

/* Tree and jump-table data structures used for fast decoding. */
typedef struct {
    /* Graph construction */
    unsigned short c[2]; /* child node */
      signed short l[2]; /* symbol to emit on transition. -1 => none */
} htree_t;

typedef struct {
    /* Byte-wise jumping table */
    unsigned short jump;
    unsigned char symbol[4];
    unsigned char nsymbols;
    unsigned char top_bit;   /* bit 9 of symbol[] */
} h_jump4_t;


/* A collection of huffman_codes_t, for use with the multi-code codec */
typedef struct {
    huffman_codes_t **codes;
    int ncodes;
    int code_set; /* (128-255) The user specified number for this encoding */

    /* Cached binary version of codeset, assumes last block */
    block_t *blk;
    int      bit_num; /* if 1st block, which bit will stored codes end on */

    /* Cache huffman_multi_decode parameters */
    h_jump4_t (*decode_J4)[16];
    htree_t *decode_t;
} huffman_codeset_t;

block_t *block_create(unsigned char *data, size_t size);
void block_destroy(block_t *blk, int keep_data);
int block_resize(block_t *blk, size_t size);
void store_bytes(block_t *block, unsigned char *val, int nbytes);


int huffman_encode(block_t *blk, huffman_codes_t *c, int code_set,
		   unsigned char *data, int len);

block_t *huffman_decode(block_t *in, huffman_codes_t *c);

int huffman_multi_encode(block_t *blk, huffman_codeset_t *cs,
			 int code_set, unsigned char *data, int len);

block_t *huffman_multi_decode(block_t *in, huffman_codeset_t *cs);

huffman_codeset_t *codes2codeset(huffman_code_t *codes, int ncodes,
				 int code_num);
huffman_codeset_t *generate_code_set(int code_set, int ncodes,
				     unsigned char *data, int len,
				     int eof, int max_code_len,
				     int all_codes);

int store_codes(block_t *out,
		huffman_codeset_t *c,
		int last_block);
huffman_codeset_t *restore_codes(block_t *block, int *bfinal);
void huffman_codes_destroy(huffman_codes_t *c);
void huffman_codeset_destroy(huffman_codeset_t *cs);

#ifdef __cplusplus
}
#endif

#endif /* _DEFLATE_SIMPLE_H_ */
