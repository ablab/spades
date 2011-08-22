/*
 * This code implements "interlaced deflate", and is *based* on a
 * simplistic implementation of the Deflate algorithm.
 * See http://www.ietf.org/rfc/rfc1951.txt for details on the original
 * Deflate algorithm.
 *
 * It differs from RFC1951 in two important ways:
 *
 * 1) It only supports the huffman encoding step and does not attempt
 *    to do any LZ-style string matching to generate distance codes.
 *    (These generally do not improve data compression for our desired
 *    use.)
 *
 * 2) It optionally allows interleaving of multiple huffman trees for
 *    a single data stream. NB: when multiple codes are used this is
 *    incompatible with RFC1951.
 *
 * It has been written here, instead of using zlib, so that we can separate
 * out the encoding of the huffman tree from the compression of the data
 * stream into separate memory sections with the intent to optimise
 * compression of very small blocks of data by sharing one set of frequency
 * tables (ie huffman tree) with multiple sets of compressed data blocks.
 *
 * James Bonfield, 2007
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#define NDEBUG /* disable asserts for production use */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <unistd.h>

#include "io_lib/deflate_interlaced.h"

#ifndef MIN
#    define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#    define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

/* #define TEST_MAIN */

/*
 * ---------------------------------------------------------------------------
 * Local structs & defines
 */

/* Used in tree construction only */
typedef struct node {
    int count;
    int sym; /* char or SYM_EOF */
    struct node *parent;
    struct node *next;
} node_t;

#define SYM_EOF 256

static void output_code_set(FILE *fp, huffman_codes_t *codes);
static void output_code_set2(FILE *fp, huffman_codes_t *codes);
int next_symbol(block_t *in, int *htab);

/*
 * ---------------------------------------------------------------------------
 * Our standard precomputed tables, for DNA, English text, etc.
 */

/* DNA */
/*
 * A   00
 * C   01
 * G   110
 * T   10
 * N   1110
 * EOF 11110
 * ?   11111*
 */
static huffman_code_t codes_dna[] = {
    {'A',  2}, {'C',  2}, {'T',  2}, {'G',  3}, {'N',  4}, {  0,  5},
    {SYM_EOF,  6},
    {  1, 13}, {  2, 13}, {  3, 13}, {  4, 13}, {  5, 13}, {  6, 13},
    {  7, 14}, {  8, 14}, {  9, 14}, { 10, 14}, { 11, 14}, { 12, 14},
    { 13, 14}, { 14, 14}, { 15, 14}, { 16, 14}, { 17, 14}, { 18, 14},
    { 19, 14}, { 20, 14}, { 21, 14}, { 22, 14}, { 23, 14}, { 24, 14},
    { 25, 14}, { 26, 14}, { 27, 14}, { 28, 14}, { 29, 14}, { 30, 14},
    { 31, 14}, { 32, 14}, { 33, 14}, { 34, 14}, { 35, 14}, { 36, 14},
    { 37, 14}, { 38, 14}, { 39, 14}, { 40, 14}, { 41, 14}, { 42, 14},
    { 43, 14}, { 44, 14}, { 45, 14}, { 46, 14}, { 47, 14}, {'0', 14},
    {'1', 14}, {'2', 14}, {'3', 14}, {'4', 14}, {'5', 14}, {'6', 14},
    {'7', 14}, {'8', 14}, {'9', 14}, { 58, 14}, { 59, 14}, { 60, 14},
    { 61, 14}, { 62, 14}, { 63, 14}, { 64, 14}, {'B', 14}, {'D', 14},
    {'E', 14}, {'F', 14}, {'H', 14}, {'I', 14}, {'J', 14}, {'K', 14},
    {'L', 14}, {'M', 14}, {'O', 14}, {'P', 14}, {'Q', 14}, {'R', 14},
    {'S', 14}, {'U', 14}, {'V', 14}, {'W', 14}, {'X', 14}, {'Y', 14},
    {'Z', 14}, { 91, 14}, { 92, 14}, { 93, 14}, { 94, 14}, { 95, 14},
    { 96, 14}, {'a', 14}, {'b', 14}, {'c', 14}, {'d', 14}, {'e', 14},
    {'f', 14}, {'g', 14}, {'h', 14}, {'i', 14}, {'j', 14}, {'k', 14},
    {'l', 14}, {'m', 14}, {'n', 14}, {'o', 14}, {'p', 14}, {'q', 14},
    {'r', 14}, {'s', 14}, {'t', 14}, {'u', 14}, {'v', 14}, {'w', 14},
    {'x', 14}, {'y', 14}, {'z', 14}, {123, 14}, {124, 14}, {125, 14},
    {126, 14}, {127, 14}, {128, 14}, {129, 14}, {130, 14}, {131, 14},
    {132, 14}, {133, 14}, {134, 14}, {135, 14}, {136, 14}, {137, 14},
    {138, 14}, {139, 14}, {140, 14}, {141, 14}, {142, 14}, {143, 14},
    {144, 14}, {145, 14}, {146, 14}, {147, 14}, {148, 14}, {149, 14},
    {150, 14}, {151, 14}, {152, 14}, {153, 14}, {154, 14}, {155, 14},
    {156, 14}, {157, 14}, {158, 14}, {159, 14}, {160, 14}, {161, 14},
    {162, 14}, {163, 14}, {164, 14}, {165, 14}, {166, 14}, {167, 14},
    {168, 14}, {169, 14}, {170, 14}, {171, 14}, {172, 14}, {173, 14},
    {174, 14}, {175, 14}, {176, 14}, {177, 14}, {178, 14}, {179, 14},
    {180, 14}, {181, 14}, {182, 14}, {183, 14}, {184, 14}, {185, 14},
    {186, 14}, {187, 14}, {188, 14}, {189, 14}, {190, 14}, {191, 14},
    {192, 14}, {193, 14}, {194, 14}, {195, 14}, {196, 14}, {197, 14},
    {198, 14}, {199, 14}, {200, 14}, {201, 14}, {202, 14}, {203, 14},
    {204, 14}, {205, 14}, {206, 14}, {207, 14}, {208, 14}, {209, 14},
    {210, 14}, {211, 14}, {212, 14}, {213, 14}, {214, 14}, {215, 14},
    {216, 14}, {217, 14}, {218, 14}, {219, 14}, {220, 14}, {221, 14},
    {222, 14}, {223, 14}, {224, 14}, {225, 14}, {226, 14}, {227, 14},
    {228, 14}, {229, 14}, {230, 14}, {231, 14}, {232, 14}, {233, 14},
    {234, 14}, {235, 14}, {236, 14}, {237, 14}, {238, 14}, {239, 14},
    {240, 14}, {241, 14}, {242, 14}, {243, 14}, {244, 14}, {245, 14},
    {246, 14}, {247, 14}, {248, 14}, {249, 14}, {250, 14}, {251, 14},
    {252, 14}, {253, 14}, {254, 14}, {255, 14},
};

/* DNA with a few ambiguity codes */
static huffman_code_t codes_dna_ambig[] = {
    {'A',  2}, {'C',  2}, {'T',  2}, {'G',  3}, {'N',  4}, {  0,  7},
    { 45,  7}, {'B',  8}, {'D',  8}, {'H',  8}, {'K',  8}, {'M',  8},
    {'R',  8}, {'S',  8}, {'V',  8}, {'W',  8}, {'Y',  8}, {SYM_EOF, 11},
    {226, 14}, {  1, 15}, {  2, 15}, {  3, 15}, {  4, 15}, {  5, 15},
    {  6, 15}, {  7, 15}, {  8, 15}, {  9, 15}, { 10, 15}, { 11, 15},
    { 12, 15}, { 13, 15}, { 14, 15}, { 15, 15}, { 16, 15}, { 17, 15},
    { 18, 15}, { 19, 15}, { 20, 15}, { 21, 15}, { 22, 15}, { 23, 15},
    { 24, 15}, { 25, 15}, { 26, 15}, { 27, 15}, { 28, 15}, { 29, 15},
    { 30, 15}, { 31, 15}, { 32, 15}, { 33, 15}, { 34, 15}, { 35, 15},
    { 36, 15}, { 37, 15}, { 38, 15}, { 39, 15}, { 40, 15}, { 41, 15},
    { 42, 15}, { 43, 15}, { 44, 15}, { 46, 15}, { 47, 15}, {'0', 15},
    {'1', 15}, {'2', 15}, {'3', 15}, {'4', 15}, {'5', 15}, {'6', 15},
    {'7', 15}, {'8', 15}, {'9', 15}, { 58, 15}, { 59, 15}, { 60, 15},
    { 61, 15}, { 62, 15}, { 63, 15}, { 64, 15}, {'E', 15}, {'F', 15},
    {'I', 15}, {'J', 15}, {'L', 15}, {'O', 15}, {'P', 15}, {'Q', 15},
    {'U', 15}, {'X', 15}, {'Z', 15}, { 91, 15}, { 92, 15}, { 93, 15},
    { 94, 15}, { 95, 15}, { 96, 15}, {'a', 15}, {'b', 15}, {'c', 15},
    {'d', 15}, {'e', 15}, {'f', 15}, {'g', 15}, {'h', 15}, {'i', 15},
    {'j', 15}, {'k', 15}, {'l', 15}, {'m', 15}, {'n', 15}, {'o', 15},
    {'p', 15}, {'q', 15}, {'r', 15}, {'s', 15}, {'t', 15}, {'u', 15},
    {'v', 15}, {'w', 15}, {'x', 15}, {'y', 15}, {'z', 15}, {123, 15},
    {124, 15}, {125, 15}, {126, 15}, {127, 15}, {128, 15}, {129, 15},
    {130, 15}, {131, 15}, {132, 15}, {133, 15}, {134, 15}, {135, 15},
    {136, 15}, {137, 15}, {138, 15}, {139, 15}, {140, 15}, {141, 15},
    {142, 15}, {143, 15}, {144, 15}, {145, 15}, {146, 15}, {147, 15},
    {148, 15}, {149, 15}, {150, 15}, {151, 15}, {152, 15}, {153, 15},
    {154, 15}, {155, 15}, {156, 15}, {157, 15}, {158, 15}, {159, 15},
    {160, 15}, {161, 15}, {162, 15}, {163, 15}, {164, 15}, {165, 15},
    {166, 15}, {167, 15}, {168, 15}, {169, 15}, {170, 15}, {171, 15},
    {172, 15}, {173, 15}, {174, 15}, {175, 15}, {176, 15}, {177, 15},
    {178, 15}, {179, 15}, {180, 15}, {181, 15}, {182, 15}, {183, 15},
    {184, 15}, {185, 15}, {186, 15}, {187, 15}, {188, 15}, {189, 15},
    {190, 15}, {191, 15}, {192, 15}, {193, 15}, {194, 15}, {195, 15},
    {196, 15}, {197, 15}, {198, 15}, {199, 15}, {200, 15}, {201, 15},
    {202, 15}, {203, 15}, {204, 15}, {205, 15}, {206, 15}, {207, 15},
    {208, 15}, {209, 15}, {210, 15}, {211, 15}, {212, 15}, {213, 15},
    {214, 15}, {215, 15}, {216, 15}, {217, 15}, {218, 15}, {219, 15},
    {220, 15}, {221, 15}, {222, 15}, {223, 15}, {224, 15}, {225, 15},
    {227, 15}, {228, 15}, {229, 15}, {230, 15}, {231, 15}, {232, 15},
    {233, 15}, {234, 15}, {235, 15}, {236, 15}, {237, 15}, {238, 15},
    {239, 15}, {240, 15}, {241, 15}, {242, 15}, {243, 15}, {244, 15},
    {245, 15}, {246, 15}, {247, 15}, {248, 15}, {249, 15}, {250, 15},
    {251, 15}, {252, 15}, {253, 15}, {254, 15}, {255, 15},
};

/* English text */
static huffman_code_t codes_english[] = {
    { 32,  3}, {'e',  3}, {'a',  4}, {'i',  4}, {'n',  4}, {'o',  4},
    {'s',  4}, {'t',  4}, {'d',  5}, {'h',  5}, {'l',  5}, {'r',  5},
    {'u',  5}, { 10,  6}, { 13,  6}, { 44,  6}, {'c',  6}, {'f',  6},
    {'g',  6}, {'m',  6}, {'p',  6}, {'w',  6}, {'y',  6}, { 46,  7},
    {'b',  7}, {'v',  7}, { 34,  8}, {'I',  8}, {'k',  8}, { 45,  9},
    {'A',  9}, {'N',  9}, {'T',  9}, { 39, 10}, { 59, 10}, { 63, 10},
    {'B', 10}, {'C', 10}, {'E', 10}, {'H', 10}, {'M', 10}, {'S', 10},
    {'W', 10}, {'x', 10}, { 33, 11}, {'0', 11}, {'1', 11}, {'F', 11},
    {'G', 11}, {  0, 15}, {  1, 15}, {  2, 15}, {  3, 15}, {  4, 15},
    {  5, 15}, {  6, 15}, {  7, 15}, {  8, 15}, {  9, 15}, { 11, 15},
    { 12, 15}, { 14, 15}, { 15, 15}, { 16, 15}, { 17, 15}, { 18, 15},
    { 19, 15}, { 20, 15}, { 21, 15}, { 22, 15}, { 23, 15}, { 24, 15},
    { 25, 15}, { 26, 15}, { 27, 15}, { 28, 15}, { 29, 15}, { 30, 15},
    { 31, 15}, { 35, 15}, { 36, 15}, { 37, 15}, { 38, 15}, { 40, 15},
    { 41, 15}, { 42, 15}, { 43, 15}, { 47, 15}, {'2', 15}, {'3', 15},
    {'4', 15}, {'5', 15}, {'6', 15}, {'7', 15}, {'8', 15}, {'9', 15},
    { 58, 15}, { 60, 15}, { 61, 15}, { 62, 15}, { 64, 15}, {'D', 15},
    {'J', 15}, {'K', 15}, {'L', 15}, {'O', 15}, {'P', 15}, {'Q', 15},
    {'R', 15}, {'U', 15}, {'V', 15}, {'X', 15}, {'Y', 15}, {'Z', 15},
    { 91, 15}, { 92, 15}, { 93, 15}, { 94, 15}, { 95, 15}, { 96, 15},
    {'j', 15}, {'q', 15}, {'z', 15}, {123, 15}, {124, 15}, {125, 15},
    {126, 15}, {127, 15}, {128, 15}, {129, 15}, {130, 15}, {131, 15},
    {132, 15}, {133, 15}, {134, 15}, {135, 15}, {136, 15}, {137, 15},
    {138, 15}, {139, 15}, {140, 15}, {141, 15}, {142, 15}, {143, 15},
    {144, 15}, {145, 15}, {146, 15}, {147, 15}, {148, 15}, {149, 15},
    {150, 15}, {151, 15}, {152, 15}, {153, 15}, {154, 15}, {155, 15},
    {156, 15}, {157, 15}, {158, 15}, {159, 15}, {160, 15}, {161, 15},
    {162, 15}, {163, 15}, {164, 15}, {165, 15}, {166, 15}, {167, 15},
    {168, 15}, {169, 15}, {170, 15}, {171, 15}, {172, 15}, {173, 15},
    {174, 15}, {175, 15}, {176, 15}, {177, 15}, {178, 15}, {179, 15},
    {180, 15}, {181, 15}, {182, 15}, {183, 15}, {184, 15}, {185, 15},
    {186, 15}, {187, 15}, {188, 15}, {189, 15}, {190, 15}, {191, 15},
    {192, 15}, {193, 15}, {194, 15}, {195, 15}, {196, 15}, {197, 15},
    {198, 15}, {199, 15}, {200, 15}, {201, 15}, {202, 15}, {203, 15},
    {204, 15}, {205, 15}, {206, 15}, {207, 15}, {208, 15}, {209, 15},
    {210, 15}, {211, 15}, {212, 15}, {213, 15}, {214, 15}, {215, 15},
    {216, 15}, {217, 15}, {218, 15}, {219, 15}, {220, 15}, {221, 15},
    {222, 15}, {223, 15}, {224, 15}, {225, 15}, {226, 15}, {227, 15},
    {228, 15}, {229, 15}, {230, 15}, {231, 15}, {232, 15}, {233, 15},
    {234, 15}, {235, 15}, {236, 15}, {237, 15}, {238, 15}, {239, 15},
    {240, 15}, {241, 15}, {242, 15}, {243, 15}, {244, 15}, {245, 15},
    {246, 15}, {247, 15}, {248, 15}, {249, 15}, {250, 15}, {251, 15},
    {252, 15}, {253, 15}, {254, 15}, {255, 15}, {SYM_EOF, 15},
};

static huffman_codeset_t *static_codeset[NCODES_STATIC];

/*
 * ---------------------------------------------------------------------------
 * Block_t structure support
 */

/*
 * Allocates and returns a new block_t struct of a specified default size.
 * A default 'data' pointer may be passed in, in which it must have
 * been created using malloc(size). Otherwise if data is NULL then
 * size indicates the amount of memory to allocate. Size maybe zero to
 * defer allocation.
 *
 * Returns newly created block_t* on success
 *         NULL on failure
 */
block_t *block_create(unsigned char *data, size_t size) {
    block_t *b = (block_t *)malloc(sizeof(*b));
    if (!b)
	return NULL;

    b->data = data;
    b->alloc = size;
    b->byte = 0;
    b->bit = 0;

    if (size && !data && NULL == (b->data = calloc(size, 1))) {
	free(b);
	return NULL;
    }

    return b;
}

/*
 * Deallocates memory created by block_create().
 * keep_data is a boolean which if true requests that the data held within
 * the block should not be deallocated as it is in use elsewhere.
 */
void block_destroy(block_t *blk, int keep_data) {
    if (!blk)
	return;

    if (!keep_data && blk->data)
	free(blk->data);
    
    free(blk);
}

/*
 * Ensures a block_t holds at least 'size' bytes.
 * Newly allocated data is initialised to zero.
 *
 * Returns 0 on success
 *        -1 on failure, leaving block pointing to the existing data
 */
int block_resize(block_t *blk, size_t size) {
    unsigned char *newp = NULL;

    if (!blk)
	return -1;

    /* Grow size to next power of 2, if we're growing */
    if (size > blk->alloc) {
	size--;
	size |= size >> 1;
	size |= size >> 2;
	size |= size >> 4;
	size |= size >> 8;
	size |= size >> 16;
	size++;
    }

    if (NULL == (newp = realloc(blk->data, size)))
	return -1;
    else
	blk->data = newp;

    if (size > blk->alloc)
	memset(&blk->data[blk->alloc], 0, size - blk->alloc);
    blk->alloc = size;

    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * Tree building and code generation functions
 */

/*
 * Reverses the order of bits in the bottom nbits of val.
 * Returns the bit-reverse value.
 */
unsigned int bit_reverse(unsigned int val, int nbits) {
    unsigned int new = 0, i;

    for (i = 0; i < nbits; i++) {
	new = (new << 1) | (val & 1);
	val >>= 1;
    }

    return new;
}


/*
 * Generates canonical huffman codes given a set of symbol bit lengths.
 * The results are stored within the supplied huffman_codes_t struct.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int canonical_codes(huffman_codes_t *c) {
    int i, j;
    unsigned int code, last_len;
    int clens[33];
    int offs[33];
    huffman_code_t ctmp[258];
    signed int symtab[258];

    /* Sort by bit-length, subfield symbol - much faster than qsort() */
    for (i = 0; i < 258; i++)
	symtab[i] = -1;
    for (i = 0; i < c->ncodes; i++)
	symtab[c->codes[i].symbol] = i;
    for (i = 0; i <= 32; i++)
	offs[i] = clens[i] = 0;
    for (i = 0; i < c->ncodes; i++)
	clens[c->codes[i].nbits]++;
    for (i = 1; i <= 32; i++)
	offs[i] = offs[i-1] + clens[i-1];
    for (i = 0; i < 258; i++) {
	if (symtab[i] != -1)
	    ctmp[offs[c->codes[symtab[i]].nbits]++] = c->codes[symtab[i]];
    }
    memcpy(c->codes, ctmp, c->ncodes * sizeof(huffman_code_t));

    /*
     * Force all codes to be <= max_code_len. This is needed due to the
     * 15-bit length limitation of Deflate literal codes and the 7-bit 
     * limit of the code bit-length table.
     */
    /* Find first point of failure */
    for (i = 0; i < c->ncodes; i++) {
	if (c->codes[i].nbits > c->max_code_len)
	    break;
    }

    /*
     * From here on we shrink the length of the current code by increasing
     * the length of an earlier symbol, at last_code.
     */
    if (i != c->ncodes) {
	int delta = 0;

	/*
	fprintf(stderr, "=== REORDERING %d ===\n", c->code_set);
	output_code_set(stderr, c);
	output_code_set2(stderr, c);
	*/

	for (; i < c->ncodes; i++) {
	    int k, cur_len;

	    c->codes[i].nbits -= delta;
	    if (c->codes[i].nbits <= c->max_code_len)
		continue;

	    for (j = i; j >= 0 && c->codes[j].nbits >= c->max_code_len; j--)
		;
	    if (j < 0) {
		fprintf(stderr,
			"Too many symbols to fit in bit-length requirements\n");
		fprintf(stderr, "=== FAILING ===\n");
		output_code_set(stderr, c);
		output_code_set2(stderr, c);
		abort();
	    }

	    /*
	    fprintf(stderr, "Changing code %d/%d to len %d\n",
		    c->codes[i].symbol, c->codes[j].symbol,
		    c->codes[j].nbits+1);
	    */
	    cur_len = c->codes[i].nbits;
	    c->codes[i].nbits = ++c->codes[j].nbits;

	    /*
	     * Shrink the next code by one, or if none at that bit-length
	     * the next 2, and so on
	     */
	    delta = 1;
	    for (k = i+1; delta && k < c->ncodes; k++) {
		while (c->codes[k].nbits > cur_len) {
		    delta *= 2;
		    cur_len++;
		}
		c->codes[k].nbits--;
		delta--;
	    }
	    assert(delta == 0);
	}

	/*
	fprintf(stderr, "=== REORDERED TO %d ===\n", c->code_set);
	output_code_set(stderr, c);
	output_code_set2(stderr, c);
	*/

	/* Ordering is shot - regenerate via brute force way */
	return canonical_codes(c);
    }


    /* Generate codes */
    code = last_len = 0; /* stop warning */
    for (i = 0; i < c->ncodes; i++) {
	int nbits = c->codes[i].nbits;

	if (i == 0) {
	    code = 0;
	    last_len = nbits;
	} else {
	    code++;
	}
	if (nbits > last_len) {
	    code <<= (nbits - last_len);
	    last_len = nbits;
	}
	c->codes[i].code = bit_reverse(code, nbits);
    }

    /* Reindex so the symbol is the primary index into codes */
    for (i = 0; i <= 257; i++) {
	c->lookup[i].nbits = 0;
    }
    for (i = 0; i < c->ncodes; i++) {
        c->lookup[c->codes[i].symbol] = c->codes[i];
    }

    return 0;
}

static int node_compar2(const void *vp1, const void *vp2) {
    const node_t *n1 = *(const node_t **)vp1;
    const node_t *n2 = *(const node_t **)vp2;

    /*
     * The sort order is vital here. This needs to return the same collating
     * order on all systems so that differing qsort() functions will not
     * swap around symbols with the same bit lengths, hence we sort by both
     * fields to force a unique stable ordering.
     */
    if (n1->count != n2->count)
	return n1->count - n2->count;
    else
	return n2->sym - n1->sym;
}

/*
 * Computes the huffman bit-lengths for a data set. We don't care
 * about the actual tree, just how deep the symbols end up.
 *
 * Huffman trees are constructed by constructing a set of nodes
 * initially containing the symbol and it's frequency. We then merge
 * the two least used nodes to produce a new node with a combined
 * frequency. Repeat until one root node is left.
 *
 * data/len is the input data to analyse.
 *
 * 'eof' is a boolean to indicate whether the EOF symbol should be included
 * in the symbols produced.
 *
 * all_codes is a boolean to indicate whether we should include symbols not
 * found in the input data set. (This was used to create the static lookup
 * tables.)
 *
 * Returns huffman_codes_t* on success
 *         NULL on failure
 */
huffman_codes_t *calc_bit_lengths(unsigned char *data, int len,
				  int eof, int max_code_len, int all_codes,
				  int start, int skip) {
    int i, ncodes;
    node_t nodes[258+257], *head, *new = &nodes[258];
    node_t *n2[258+257];
    int map[258];
    huffman_codes_t *c;
    int hist[256];

    if (NULL == (c = (huffman_codes_t *)malloc(sizeof(*c))))
	return NULL;
    c->codes_static = 0;
    c->max_code_len = max_code_len;

    /* Count frequencies of symbols */
    memset(hist, 0, 256*sizeof(*hist));
    /* Calc freqs */
    for (i = start; i < len; i+=skip) {
	hist[data[i]]++;
    }


    /*
     * Initialise nodes. We build a map of ASCII character code to node
     * number. (By default it's a simple 1:1 mapping unless legal_chars is
     * defined.)
     */
    ncodes = 0;
    if (eof) {
	nodes[ncodes].sym = SYM_EOF;
	nodes[ncodes].count = eof;
	nodes[ncodes].parent = NULL;
	n2[ncodes] = &nodes[ncodes];
	map[SYM_EOF] = ncodes++;
    }

    /* All 256 chars existing at a minimal level */
    if (all_codes) {
	for (i = 0; i < 256; i++) {
	    nodes[ncodes].sym = i;
	    nodes[ncodes].count = hist[i];
	    nodes[ncodes].parent = NULL;
	    n2[ncodes] = &nodes[ncodes];
	    map[i] = ncodes++;
	}
    } else {
	/* Only include non-zero symbols */
	for (i = 0; i < 256; i++) {
	    if (hist[i] == 0)
		continue;
	    nodes[ncodes].sym = i;
	    nodes[ncodes].count = hist[i];
	    nodes[ncodes].parent = NULL;
	    n2[ncodes] = &nodes[ncodes];
	    map[i] = ncodes++;
	}
    }

    /* Sort by counts, smallest first and form a sorted linked list */
    qsort(n2, ncodes, sizeof(*n2), node_compar2);

    /* Skip symbols that do not occur, unless all_codes is true */
    for (i = 0; i < ncodes; i++) {
	n2[i]->next = i+1 < ncodes ? n2[i+1] : NULL;
    }

    /* Repeatedly merge two smallest values */
    head = n2[0];
    while (head && head->next) {
	node_t *after = head->next, *n;
	int sum = head->count + head->next->count;
	
	for (n = head->next->next; n; after = n, n = n->next) {
	    if (sum < n->count)
		break;
	}

	/* Produce a new summation node and link it in place */
	after->next = new;
	new->next = n;
	new->sym = '?';
	new->count = sum;
	new->parent = NULL;
	head->parent = new;
	head->next->parent = new;
	head = head->next->next;

	new++;
    }

    /* Walk up tree computing the bit-lengths for our symbols */
    c->ncodes = ncodes;
    c->codes = (huffman_code_t *)malloc(c->ncodes * sizeof(*c->codes));
    if (NULL == c->codes) {
	free(c);
	return NULL;
    }

    for (i = 0; i < ncodes; i++) {
	int len = 0;
	node_t *n;
	for (n = n2[i]->parent; n; n = n->parent) {
	    len++;
	}

	c->codes[i].symbol = n2[i]->sym;
	c->codes[i].freq   = n2[i]->count;
	c->codes[i].nbits  = len ? len : 1; /* special case, nul input */
    }

    if (0 != canonical_codes(c)) {
	free(c);
	return NULL;
    }

    return c;
}

/*
 * A special case of the generate_code_set() function below, but for
 * creating predefined code sets from bit-length arrays. Useful for
 * code that wants to use a predetermined huffman tree.
 *
 *
 * Returns huffman_codes_t* on success; free using huffman_codes_destroy().
 *         NULL on failure.
 */
huffman_codeset_t *codes2codeset(huffman_code_t *codes, int ncodes,
				 int code_num) {
    huffman_codeset_t *cs;
    huffman_codes_t *c;

    if (NULL == (cs = (huffman_codeset_t *)malloc(sizeof(*cs))))
	return NULL;

    if (NULL == (c = (huffman_codes_t *)malloc(sizeof(*c))))
	return NULL;

    cs->codes = (huffman_codes_t **)malloc(sizeof(*cs->codes));
    cs->codes[0] = c;
    cs->ncodes = 1;
    cs->code_set = code_num;
    cs->blk = NULL;
    cs->bit_num = 0;
    cs->decode_t = NULL;
    cs->decode_J4 = NULL;
    
    c->codes_static = 1;
    c->max_code_len = MAX_CODE_LEN;

    c->codes = codes;
    c->ncodes = ncodes;

    cs->bit_num = 0; /* FIXME: need to know this */

    canonical_codes(c);

    return cs;
}

/*
 * Initialises and returns a huffman_codes_t struct from a specified code_set.
 * If code_set is not one of the standard predefined values then the
 * input data is analysed using calc_bit_lengths() above to produce the
 * optimal set of huffman codes, otherwise we return predefined values.
 *
 * 'eof' is a boolean to indicate whether the EOF symbol should be included
 * in the symbols produced.
 *
 * all_codes is a boolean to indicate whether we should include symbols not
 * found in the input data set. (This was used to create the static lookup
 * tables.)
 *
 * Returns huffman_codes_t* on success; free using huffman_codes_destroy().
 *         NULL on failure.
 */
huffman_codeset_t *generate_code_set(int code_set, int ncodes,
				     unsigned char *data, int len,
				     int eof, int max_code_len,
				     int all_codes) {
    huffman_codeset_t *cs;

    /*
     * Either INLINE or a CODE_USER+ set of codes.
     * => analyse the data and compute a new set of bit-lengths & codes.
     */
    if (code_set >= 128 || code_set == CODE_INLINE) { 
	int i;

	if (NULL == (cs = (huffman_codeset_t *)malloc(sizeof(*cs))))
	    return NULL;

	cs->code_set = code_set;
	cs->ncodes = ncodes;
	cs->codes = (huffman_codes_t **)malloc(cs->ncodes*sizeof(*cs->codes));
	cs->blk = NULL;
	cs->bit_num = 0;
	cs->decode_t = NULL;
	cs->decode_J4 = NULL;

	for (i = 0; i < ncodes; i++) {
	    /*
	     * If requested, include EOF all code sets, but at a
	     * frequency of only '1' occurrance where we predict it
	     * not to be needed.
	     */
	    if (eof && (len+i)%ncodes)
		eof = 1;
	    cs->codes[i] = calc_bit_lengths(data, len, eof, max_code_len,
					    all_codes, i, ncodes);
	    cs->codes[i]->codes_static = 0;
	    if (NULL == cs->codes[i]) {
		/* FIXME: tidy up */
		return NULL;
	    }

	    canonical_codes(cs->codes[i]);
	}

    /*
     * Otherwise we use the determined codes at the top of this file, such
     * as codes_dna and codes_english.
     */
    } else {
	if (code_set < 1 || code_set >= NCODES_STATIC) {
	    fprintf(stderr, "Unknown huffman code set '%d'\n", code_set);
	    return NULL;
	}

	/* If our global codeset hasn't been initialised yet, do so */
	if (!static_codeset[code_set]) {
	    huffman_codes_t *c = (huffman_codes_t *)malloc(sizeof(*c));

	    if (NULL == (cs = (huffman_codeset_t *)malloc(sizeof(*cs))))
		return NULL;

	    cs->codes = (huffman_codes_t **)malloc(sizeof(*cs->codes));
	    cs->codes[0] = c;
	    cs->ncodes = 1;
	    cs->code_set = code_set;
	    cs->blk = NULL;
	    cs->bit_num = 0;
	    cs->decode_t = NULL;
	    cs->decode_J4 = NULL;

	    c->codes_static = 1;
	    c->max_code_len = MAX_CODE_LEN;

	    switch(code_set) {
	    case CODE_DNA:
		c->codes = codes_dna;
		c->ncodes = sizeof(codes_dna)/sizeof(*c->codes);
		cs->bit_num = 5;
		break;

	    case CODE_DNA_AMBIG:
		c->codes = codes_dna_ambig;
		c->ncodes = sizeof(codes_dna_ambig)/sizeof(*c->codes);
		cs->bit_num = 1;
		break;

	    case CODE_ENGLISH:
		c->codes = codes_english;
		c->ncodes = sizeof(codes_english)/sizeof(*c->codes);
		cs->bit_num = 1;
		break;

	    default:
		fprintf(stderr, "Unknown huffman code set '%d'\n", code_set);
		return NULL;
	    }

	    canonical_codes(c);

	    static_codeset[code_set] = cs;
	}

	cs = static_codeset[code_set];
    }

    return cs;
}

void huffman_codes_destroy(huffman_codes_t *c) {
    if (!c)
	return;

    if (!c->codes_static && c->codes)
	free(c->codes);

    free(c);
}

void huffman_codeset_destroy(huffman_codeset_t *cs) {
    int i;

    if (!cs)
	return;

    /* If this codeset is one of the predefined global ones we do nothing */
    if (cs->ncodes == 1 && cs->codes[0]->codes_static)
	return;

    for (i = 0; i < cs->ncodes; i++)
	huffman_codes_destroy(cs->codes[i]);
    if (cs->codes)
	free(cs->codes);
    if (cs->blk)
	block_destroy(cs->blk, 0);
    if (cs->decode_t)
	free(cs->decode_t);
    if (cs->decode_J4)
	free(cs->decode_J4);

    free(cs);
}

/*
 * ---------------------------------------------------------------------------
 * Encoding and decoding related functions
 */

/*
 * Can store up to 24-bits worth of data encoded in an integer value
 * Possibly we'd want to have a less optimal store_bits function when dealing
 * with nbits > 24, but for now we assume the codes generated are never
 * that big. (Given this is only possible with 121392 or more
 * characters with exactly the correct frequency distribution we check
 * for it elsewhere.)
 */
static void store_bits(block_t *block, unsigned int val, int nbits) {
    /* fprintf(stderr, " store_bits: %02x %d\n", val, nbits); */

#if 1
    {
	unsigned int curr = block->data[block->byte];
	curr |= (val & ((1 << nbits)-1)) << block->bit;
	block->bit += nbits;
	while (block->bit >= 8) {
	    block->data[block->byte++] = curr & 0xff;
	    curr >>= 8;
	    block->bit -= 8;
	}
	block->data[block->byte] = curr & 0xff;
    }
    return;

#else

    {
      /* Slow, but simple */
      unsigned int mask = 1;
      int bit = 1 << block->bit;
      do {
	if (val & mask)
	    block->data[block->byte] |= bit;
	/*
	 * Data should be zeroed anyway, so this is not needed.
	 *
	 * else
	 *    block->data[block->byte] &= ~bit;
	*/

	if (++block->bit == 8) {
	    block->bit = 0;
	    block->byte++;
	    block->data[block->byte] = 0;
	    bit = 1;
	} else {
	    bit <<= 1;
	}
	mask <<= 1;
      } while(--nbits);
    }

#endif
}

/*
 * Reads up to 24-bits worth of data and returns. Updates the block
 * byte and bit values to indicate the current 'read' position.
 *
 * Returns unsigned value on success (>=0)
 *         -1 on failure
 */
static signed int get_bits(block_t *block, int nbits) {
    unsigned int val, bnum = 0;

    if (block->byte*8 + block->bit + nbits > block->alloc * 8)
	return -1;

    /* Fetch the partial byte of data */
    val = (block->data[block->byte]) >> block->bit;
    bnum = 8 - block->bit;

    /* And additional entire bytes worth as required */
    while (bnum <= nbits) {
	val |= block->data[++block->byte] << bnum;
	bnum += 8;
    }

    block->bit = (block->bit + nbits) % 8;
    return val & ((1 << nbits) - 1);
}

/* stores nbytes bytes, padding to align on the next byte boundary */
void store_bytes(block_t *block, unsigned char *val, int nbytes) {
    /* Align */
    if (block->bit != 0) {
	block->byte++;
	block->bit = 0;
    }

    /* Resize */
    block_resize(block, block->byte + nbytes + 1);

    /* Store */
    memcpy(&block->data[block->byte], val, nbytes);
    block->byte += nbytes;
}


/*
 * Encodes the huffman symbol bit-lengths as a serialised block of data
 * suitable for storing in a ZTR "ZLBH" chunk. This uses the Deflate
 * storage format defined in RFC1951.
 *
 * Returns: 0 on success
 *         -1 on failure
 */
int store_codes_single(block_t *out, huffman_codes_t *c) {
    int i;
    unsigned char bl_code[257]; /* bit-length codes and for codes 16-18 */
    unsigned char bl_opt[257];  /*     the operand to the blcode */
    unsigned char sorted_codes[258];
    int bl_freq[19]; /* frequency of bit-length codes produced */
    int bl_count;
    huffman_codes_t *bl_cds = NULL;
    int hclen_order[] = {
	16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };
    int hlit, hdist, hclen, hcmap[19];

    /* output_code_set (stderr, c); */

    if (out->alloc < out->byte + 1000) {
	out->alloc = out->byte + 1000;
	if (NULL == (out->data = realloc(out->data, out->alloc)))
	    return -1;
    }

    /*
     *-----------------------------------------------------------------
     * Reformat the dynamic code bit-lengths into an alphabet of 19 
     * "code length" symbols as defined in RFC1951.
     */
    memset(sorted_codes, 0, 258);
    for (i = 0; i < c->ncodes; i++) {
	sorted_codes[c->codes[i].symbol] = c->codes[i].nbits;
    }
    for (i = 0; i < 19; i++)
	bl_freq[i] = 0;

    bl_count = 0;
    for (i = 0; i < 257; ) {
	int j = i+1, n;
	int v = sorted_codes[i];
	while (j < 257 && sorted_codes[j] == v)
	    j++;

	n = j-i; /* n = run-length */
	/* fprintf(stderr, "value=%d, run_len=%d\n", v, n); */
	if (v == 0) {
	    /* bit-len zero => no code and uses code 17/18 for run len */
	    while (n > 0) {
		while (n >= 11) {
		    bl_freq[18]++;
		    bl_code[bl_count] = 18;
		    bl_opt[bl_count] = MIN(n, 138)-11;
		    n -= bl_opt[bl_count]+11;
		    bl_count++;
		}

		while (n >= 3) {
		    bl_freq[17]++;
		    bl_code[bl_count] = 17;
		    bl_opt[bl_count] = MIN(n, 10)-3;
		    n -= bl_opt[bl_count]+3;
		    bl_count++;
		}

		switch (n) {
		case 2: bl_code[bl_count++] = 0; bl_freq[0]++; n--;
		case 1: bl_code[bl_count++] = 0; bl_freq[0]++; n--;
		}
	    }
	} else if (v <= 15) {
	    /* non-zero code, uses code 16 for run-len */
	    if (n >= 4) {
		bl_freq[v]++;
		bl_code[bl_count++] = v;
		n--;
		while (n >= 3) {
		    bl_freq[16]++;
		    bl_code[bl_count] = 16;
		    bl_opt[bl_count] = MIN(n, 6)-3;
		    n -= bl_opt[bl_count]+3;
		    bl_count++;
		}
	    }

	    switch(n) {
	    case 3: bl_code[bl_count++] = v; bl_freq[v]++; n--;
	    case 2: bl_code[bl_count++] = v; bl_freq[v]++; n--;
	    case 1: bl_code[bl_count++] = v; bl_freq[v]++; n--;
	    }
	} else {
	    fprintf(stderr, "Unsupported code length: %d\n", v);
	}

	i = j;
    }
    hlit = 257;

    /* Add a single distance code of zero bits. This means that there
     * are no distance codes used at all.
     */
    bl_code[bl_count++] = 0;
    bl_freq[0]++;
    hdist = 1;

    /* Produce new huffman codes for our code-length symbols. */
    bl_cds = calc_bit_lengths(bl_code, bl_count, 0, 7, 0, 0, 1);
    /* output_code_set (stderr, bl_cds); */

    /*
     *-----------------------------------------------------------------
     * Output the "code length" bit-lengths, 3 bits at a time.
     *
     * Compute how many HCLEN code length values we need, using the
     * predefined order in the RFC.
     */
    for (hclen = 19; hclen > 0; hclen--) {
	if (bl_freq[hclen_order[hclen-1]])
	    break;
    }

    store_bits(out, hlit-257,  5);
    store_bits(out, hdist-1, 5);
    store_bits(out, hclen-4, 4);

    for (i = 0; i < 19; i++)
	hcmap[i] = -1;
    for (i = 0; i < bl_cds->ncodes; i++)
	hcmap[bl_cds->codes[i].symbol] = i;
    for (i = 0; i < hclen; i++) {
	if (hcmap[hclen_order[i]] >= 0) {
	    store_bits(out, bl_cds->codes[hcmap[hclen_order[i]]].nbits, 3);
	} else {
	    store_bits(out, 0, 3);
	}
    }


    /*
     *----------------------------------------------------------------
     * Finally output the original bit-lengths using the code-len codes.
     */
    for (i = 0; i < bl_count; i++) {
	huffman_code_t *c = &bl_cds->codes[hcmap[bl_code[i]]];
	store_bits(out, c->code, c->nbits); 
	/*
	fprintf(stderr, "bl_code %d (opt %d), code %d/%d\n",
		bl_code[i], bl_opt[i], c->code, c->nbits);
	*/
	switch(bl_code[i]) {
	case 18:
	    store_bits(out, bl_opt[i], 7);
	    break;
	case 17:
	    store_bits(out, bl_opt[i], 3);
	    break;
	case 16:
	    store_bits(out, bl_opt[i], 2);
	    break;
	}
    }

    huffman_codes_destroy(bl_cds);

    return 0;
}

/*
 * A wrapper around store_codes_single to output either a single or multiple
 * huffman codes to a block.
 *
 * This also creates a new block and fills out the block header appropriately.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int store_codes(block_t *out, huffman_codeset_t *cs, int last_block) {
    int i;

    if (out->alloc < out->byte + 1000) {
	out->alloc = out->byte + 1000;
	if (NULL == (out->data = realloc(out->data, out->alloc)))
	    return -1;
    }

    /* Header details */
    store_bits(out, last_block != 0,  1); /* last block */
    if (cs->ncodes == 1) {
	store_bits(out, 2,  2); /* dynamic huffman */
    } else {
	int nbits = 0;
	store_bits(out, 3,  2); /* multiple tree dynamic huffman */
	while (1<<nbits <= cs->ncodes-1)
	    nbits++;
	store_bits(out, nbits-1, 4);
	store_bits(out, cs->ncodes-1, nbits);
    }

    for (i = 0; i < cs->ncodes; i++) {
	if (-1 == store_codes_single(out, cs->codes[i]))
	    return -1;
    }

    return 0;
}

/*
 * This is the opposite of the store_codes() function. It loads generates
 * huffman_codes_t structs from the a serialised data stream as presented
 * in the above format.
 *
 * The input data is the data-string. On return the number of bytes
 * consumed will be returned in *len_used (if non NULL).
 * This is to allow stripping off of the huffman codes from a longer
 * array of data (ie probably followed by the STHUFF encoded chunk
 * itself).
 *
 * Returns: malloced huffman_codes_t structure on success.
 *          NULL on failure.
 */
huffman_codes_t *restore_codes_single(block_t *block) {
    int hlit, hdist, hclen;
    int hclen_order[19] = {
	16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };
    int hc_bitlen[19], i;
    huffman_codes_t *bl_cds, *cds;
    int sym, sym_val, last_len;
    int htab[256];

    hlit   = get_bits(block, 5)+257;
    hdist  = get_bits(block, 5)+1;
    hclen  = get_bits(block, 4)+4;

    /*
    fprintf(stderr, "bfinal = %d, btype=%d\n", *bfinal, btype);
    fprintf(stderr, "hlit=0x%x, hdist=0x%x, hclen=0x%x\n",
	    hlit, hdist, hclen);
    */

    /* Read HCLEN code-lengths and construct huffman codes from them */
    for (i = 0; i < hclen; i++)
	hc_bitlen[hclen_order[i]] = get_bits(block, 3);
    for (; i < 19; i++)
	hc_bitlen[hclen_order[i]] = 0;

    bl_cds = (huffman_codes_t *)malloc(sizeof(*bl_cds));
    bl_cds->codes_static = 0;
    bl_cds->ncodes = 0;
    bl_cds->codes = (huffman_code_t *)malloc(19 * sizeof(*bl_cds->codes));
    bl_cds->max_code_len = 7;
    for (i = 0; i < 19; i++) {
	if (hc_bitlen[i]) {
	    bl_cds->codes[bl_cds->ncodes].symbol = i;
	    bl_cds->codes[bl_cds->ncodes].nbits = hc_bitlen[i];
	    bl_cds->ncodes++;
	}
    }

    canonical_codes(bl_cds);
    /* output_code_set (stderr, bl_cds); */

    /* Build a lookup table of possible codes to symbols */
    for (i = 0; i < 256; i++)
	htab[i] = -1;
    for (i = 0; i < bl_cds->ncodes; i++) {
	htab[bit_reverse(bl_cds->codes[i].code, bl_cds->codes[i].nbits)
	     | (1<<bl_cds->codes[i].nbits)]
	    = bl_cds->codes[i].symbol;
    }

    /* Now decode the next HLIT literal codes using bl_cds */
    cds = (huffman_codes_t *)malloc(sizeof(*cds));
    cds->codes_static = 0;
    cds->ncodes = 0;
    cds->codes = (huffman_code_t *)malloc(257 * sizeof(*cds->codes));
    cds->max_code_len = 15;
    sym_val = last_len = 0;
    while ((sym = next_symbol(block, htab)) != -1) {
	int count;
	/* fprintf(stderr, "LIT Sym=%d\n", sym); */

	switch(sym) {
	case 16:
	    count = get_bits(block, 2) + 3;
	    /* fprintf(stderr, "   +%d\n", count); */
	    for (i = 0; i < count; i++) {
		cds->codes[cds->ncodes].symbol = sym_val++;
		cds->codes[cds->ncodes++].nbits = last_len;
	    }
	    break;
	    
	case 17:
	    count = get_bits(block, 3) + 3;
	    /* fprintf(stderr, "   +%d\n", count); */
	    sym_val += count;
	    last_len = 0;
	    break;
	    
	case 18:
	    count = get_bits(block, 7) + 11;
	    /* fprintf(stderr, "   +%d\n", count); */
	    sym_val += count;
	    last_len = 0;
	    break;

	case 0:
	    sym_val++;
	    last_len = 0;
	    break;

	default:
	    cds->codes[cds->ncodes].symbol = sym_val++;
	    last_len = cds->codes[cds->ncodes++].nbits = sym;
	}

	if (sym_val >= hlit)
	    break;
    }
    assert(sym != -1);
    assert(cds->ncodes <= 257);

    /* Skip HDIST codes. Hopefully only 1 of zero length */
    sym_val = 0;
    while ((sym = next_symbol(block, htab)) != -1) {
	/* fprintf(stderr, "DIST Sym=%d\n", sym); */

	switch(sym) {
	case 16:
	    sym_val += get_bits(block, 2) + 3;
	    break;
	    
	case 17:
	    sym_val += get_bits(block, 3) + 3;
	    break;
	    
	case 18:
	    sym_val += get_bits(block, 7) + 11;
	    break;

	default:
	    sym_val++;
	}

	if (sym_val >= hdist)
	    break;
    }
    assert(sym != -1);

    huffman_codes_destroy(bl_cds);
    canonical_codes(cds);
    /* output_code_set(stderr, cds); */

    return cds;
}

huffman_codeset_t *restore_codes(block_t *block, int *bfinal) {
    int btype;
    huffman_codeset_t *cs;
    
    /* Header details */
    if (bfinal)
	*bfinal = get_bits(block, 1);
    else
	get_bits(block, 1);
    btype  = get_bits(block, 2);

    cs = (huffman_codeset_t *)malloc(sizeof(*cs));
    cs->code_set = 0;
    cs->blk = NULL;
    cs->bit_num = 0;
    cs->decode_t = NULL;
    cs->decode_J4 = NULL;

    if (btype == 2) {
	/* Standard Deflate algorithm */
	cs->ncodes = 1;
	cs->codes = (huffman_codes_t **)malloc(cs->ncodes*sizeof(*cs->codes));
	cs->codes[0] = restore_codes_single(block);
    } else if (btype == 3) {
	/* Deflate extension - multiple codes */
	int nbits, i;
	nbits  = get_bits(block, 4) + 1;
	cs->ncodes = get_bits(block, nbits) + 1;
	cs->codes = (huffman_codes_t **)malloc(cs->ncodes*sizeof(*cs->codes));
	for (i = 0; i < cs->ncodes; i++) {
	    cs->codes[i] = restore_codes_single(block);
	}
    } else {
	fprintf(stderr, "restore_codes() only implemented for "
		"BTYPE == DYNAMIC HUFFMAN and INTERLACED HUFFMAN\n");
	return NULL;
    }

    cs->bit_num = block->bit;

    return cs;
}

/*
 * Given a multiple sets of huffman codes and a block of data this
 * compresses and returns the data block. It iterates around each set
 * of huffman codes in a cyclic fashion encoding each byte with the
 * appropriate huffman codes.
 *
 * Returns: 0 on success
 *          -1 on failure
 */
int huffman_multi_encode(block_t *blk, huffman_codeset_t *cs,
			 int code_set, unsigned char *data, int len) {
    int i, nc;
    huffman_code_t *lookup;
    huffman_codes_t **c;

    if (!cs) {
	/* No codes known, so derive our own */
	fprintf(stderr, "FIXME: use generate_code_set() to build our own codes here\n");
	return -1;
    } else {
	c  = cs->codes;
	nc = cs->ncodes;
    }

    /*
     * The maximum size to encode len bytes is < 9 bits per symbol
     * (not quite 8 due to an EOF symbol) plus the overhead of the bit-length
     * tree. That in turn, with alternating 8/9 bit-lengths would max out
     * as 258*8 + 5+5+4 + 19*3 + 258*5 bits (429 bytes), but in practice
     * I'm not even sure if it's possible to construct such a set of code
     * lengths that would compress that poor.
     *
     * This of course assumes we're using appropriate compression codes.
     * Given a user may give a completely inappropriate code we have to
     * assume every symbol is actually 15 bits instead of < 9 on average.
     *
     * We ensure blk here is large enough for the worst case scenario so we
     * don't incur overheads in store_bits().
     */
    if (blk->alloc <= 429 + 2*len + 2 + blk->byte) {
	blk->alloc  = 429 + 2*len + 2 + blk->byte;
	blk->data = realloc(blk->data, blk->alloc);
	if (!blk->data)
	    return -1;
    }

    /*
     * Splitting this special case out is worth it as it's approx 30% faster.
     * Also note that the nc > 1 case is faster with a separate counter and
     * test than using modulus (by a factor of 2). It could be sped up further
     * for powers of 2 using bitwise AND, but the difference is not huge.
     */
    if (nc == 1) {
	lookup = c[0]->lookup;
	for (i = 0; i < len; i++) {
	    assert(lookup[data[i]].nbits > 0);
	    store_bits(blk, lookup[data[i]].code,
		       lookup[data[i]].nbits);
	}
    } else {
	int count = 0;
	for (i = 0; i < len; i++) {
	    lookup = c[count]->lookup;
	    assert(lookup[data[i]].nbits > 0);
	    store_bits(blk, lookup[data[i]].code,
		       lookup[data[i]].nbits);
	    if (++count == nc)
		count = 0;
	}
    }

    lookup = c[i%nc]->lookup;
    store_bits(blk, lookup[SYM_EOF].code, lookup[SYM_EOF].nbits);

    assert(blk->alloc > blk->byte);

    /* We probably massively overallocated, so return some of it back */
    blk->data = realloc(blk->data, blk->byte+1);
    blk->alloc = blk->byte+1;

    return 0;
}

/*
 * The opposite of huffman_encode().
 * Decode a huffman stream from 'block' using huffman codes 'c'.
 *
 * Returns: allocated block_t pointer on success
 *          NULL on failure.
 *
 * Method 1
 * --------
 *
 * At any node in our tree we can precompute a lookup table so that upon
 * reading the next 'k' bits we know the new node we'd end up in and what
 * symbols to export.
 * Then decoding simply works in fixed sets of k bits at a time.
 *
 * We use k=4 for efficient table space (they fit neatly in cache) and ease
 * of decoding 4-bits at a time. k=8 is about 20% faster as reading the input
 * byte by byte is easy, but the setup time is substantially longer
 * (16x at a guess) and the lookup tables no longer fit in the L1 cache.
 */
block_t *huffman_decode(block_t *in, huffman_codes_t *c) {
    block_t *out;
    htree_t t[513]; /* number of internal nodes */
    int i, j, n;
    int new_node, node_num;
    h_jump4_t J4[513][16];
    unsigned char *cp;

    if (NULL == (out = block_create(NULL, 8*in->alloc+8))) {
	block_destroy(in, 0);
	return NULL;
    }

    /* Construct the tree from the codes */
    new_node = 1;
    t[0].l[0] = t[0].l[1] = -1;
    t[0].c[0] = t[0].c[1] = 0;
    for (i = 0; i < c->ncodes; i++) {
	int n = 0;
	unsigned int v = c->codes[i].code;

	for (j = 0; j < c->codes[i].nbits-1; j++) {
	    int b = v & 1;
	    if (t[n].c[b]) {
		n = t[n].c[b];
	    } else {
		n = (t[n].c[b] = new_node++);
		t[n].c[0] = t[n].c[1] = 0;
		t[n].l[0] = t[n].l[1] = -1;
	    }
	    v >>= 1;
	}
	/* last bit */
	t[n].l[v & 1] = c->codes[i].symbol;
    }

    /* Build the 16 wide lookup table per node */
    for (n = 0; n < new_node; n++) {
	for (j = 0; j < 16; j++) {
	    unsigned int v = j;
	    int n2 = n;
	    h_jump4_t *hj = &J4[n][j];
	    hj->nsymbols = 0;
	    hj->top_bit = 0;
	    for (i = 0; i < 4; i++) {
		int b = v & 1;
		if (t[n2].l[b] >= 0) {
		    hj->symbol[hj->nsymbols++] = t[n2].l[b];
		    if (t[n2].l[b] == SYM_EOF)
			if (!hj->top_bit)
			    hj->top_bit |= 1 << (hj->nsymbols-1);
		}
		n2 = t[n2].c[b];
		v >>= 1;
	    }
	    hj->jump = n2;
	}
    }

#if 0
    /* Debug output */
    for (n = 0; n < new_node; n++) {
	printf("Node %d, c[]={%d,%d}, l[]={%d,%d}\n",
	       n, t[n].c[0], t[n].c[1], t[n].l[0], t[n].l[1]);
	for (i = 0; i < 256; i++) {
	    printf("\t%02x %s =>%02d, ", i, print_8rev(i), J4[n][i].jump);
	    for (k = 0; k < J4[n][i].nsymbols; k++) {
		if (isprint(J4[n][i].symbol[k]))
		    printf(" '%c'", J4[n][i].symbol[k]);
		else
		    printf(" %03d", J4[n][i].symbol[k]);
	    }
	    printf("\n");
	}
    }
#endif


    /*
     * Decoding - part 1
     * We're part way through a byte, so decode bit by bit up to the next
     * whole byte and then we start the fast decoding section.
     */
    cp = out->data;
    node_num = 0;
    while (in->bit != 0) {
	int b = get_bits(in, 1);
	if (t[node_num].l[b] != -1) {
	    if (t[node_num].l[b] != SYM_EOF) {
		*cp++ = t[node_num].l[b];
	    } else {
		out->byte = cp - out->data;
		return out;
	    }
	}
	node_num = t[node_num].c[b];
    }

    /*
     * Decoding - part 2
     *
     * We handle data nibble by nibble using the nibble to get an
     * h_jump4_t lookup from the J4[] table.
     * If top_bit is clear then we know we have no funny business (SYM_EOF)
     * so we use a fast decoding technique, otherwise we have to do a slower
     * loop with a check.
     */
    {
	int last_node = 0;
	unsigned char *last_cp = cp;
	h_jump4_t *x = &J4[node_num][in->data[in->byte] & 0x0f];
	int l = x->nsymbols;
	int b;

	/*
	 * This is the tight loop, so we over-optimise here by ignoring EOF
	 * and relying on knowing the length of the input data stream.
	 * This allows us to ignore the 9-bit data and only operate on
	 * the basic 0-255 symbols, glossing over the minor issue that EOF
	 * will look like an ordinary symbol.
	 */
	for (i = in->byte; i < in->alloc; i++) {
	    last_cp = cp;
	    last_node = node_num;

	    x = &J4[node_num][in->data[i] & 0x0f];
	    l = x->nsymbols;

	    for (j = 0; j < l; j++) {	
		*cp++ = x->symbol[j];
	    }
	    node_num = x->jump;

	    if (x->top_bit)
		break;

	    x = &J4[node_num][(in->data[i] >> 4) & 0x0f];
	    l = x->nsymbols;
	     
	    for (j = 0; j < l; j++) {	
		*cp++ = x->symbol[j];
	    }
	    node_num = x->jump;

	    if (x->top_bit)
		break;
	}

    
	/*
	 * Decoding - part 3
	 *
	 * The above optimisation has unfortunately added EOF to our data
	 * along with whatever else was packed in the last byte after the
	 * EOF symbol. So we rewind one byte and finish off decoding
	 * the slow way - walking the tree.
	 */
	cp = last_cp;
	node_num = last_node;
	in->byte = i;
	in->bit = 0;
	while (-1 != (b = get_bits(in, 1))) {
	    if (t[node_num].l[b] != -1) {
		if (t[node_num].l[b] != SYM_EOF) {
		    *cp++ = t[node_num].l[b];
		} else {
		    out->byte = cp - out->data;
		    return out;
		}
	    }
	    node_num = t[node_num].c[b];
	}
     }


    /* We shouldn't reach here */
    return NULL;
}

int init_decode_tables(huffman_codeset_t *cs) {
    int nnodes, i, j, n, nc;
    huffman_codes_t **c;
    int new_node, rec;
    h_jump4_t (*J4)[16] = NULL;
    htree_t *t;
    
    c = cs->codes;
    nc = cs->ncodes;

    /* Allocate memory for internal nodes (nsyms-1 for each code set) */
    for (nnodes = i = 0; i < nc; i++) {
	nnodes += c[i]->ncodes-1;
    }

    if (NULL == (t = (htree_t *)malloc(nnodes * sizeof(*t))))
	goto error;

    if (NULL == (J4 = (h_jump4_t (*)[16])malloc(nnodes * sizeof(*J4))))
	goto error;

    /*
     * Construct the tree from the codes.
     * We have one tree for all 'nc' huffman codes with each tree pointing
     * to the root of the next one (or first) tree whenever we emit a
     * symbol.
     *
     * This then effectively means the decoding step is identical to the
     * single huffman code function.
     */
    new_node = 0;
    for (rec = 0; rec < nc; rec++) {
	int root = new_node++;
	int next_root = rec == nc-1
	    ? 0
	    : root + c[rec]->ncodes-1;
	    
	t[root].l[0] = t[root].l[1] = -1;
	t[root].c[0] = t[root].c[1] = next_root;
	for (i = 0; i < c[rec]->ncodes; i++) {
	    int n = root;
	    unsigned int v = c[rec]->codes[i].code;

	    for (j = 0; j < c[rec]->codes[i].nbits-1; j++) {
		int b = v & 1;
		if (t[n].c[b] != next_root) {
		    n = t[n].c[b];
		} else {
		    n = (t[n].c[b] = new_node++);
		    t[n].c[0] = t[n].c[1] = next_root;
		    t[n].l[0] = t[n].l[1] = -1;
		}
		v >>= 1;
	    }
	    /* last bit */
	    t[n].l[v & 1] = c[rec]->codes[i].symbol;
	}
    }

    /*
    for (i = 0; i < new_node; i++) {
	printf("t[%d] = {(%d,%d), (%d,%d)}\n",
	       i,
	       t[i].l[0], t[i].l[1],
	       t[i].c[0], t[i].c[1]);
    }
    */

    /* Build the 16 wide lookup table per node */
    for (n = 0; n < new_node; n++) {
	for (j = 0; j < 16; j++) {
	    unsigned int v = j;
	    int n2 = n;
	    h_jump4_t *hj = &J4[n][j];
	    hj->nsymbols = 0;
	    hj->top_bit = 0;
	    for (i = 0; i < 4; i++) {
		int b = v & 1;
		if (t[n2].l[b] >= 0) {
		    hj->symbol[hj->nsymbols++] = t[n2].l[b];
		    if (t[n2].l[b] == SYM_EOF)
			if (!hj->top_bit)
			    hj->top_bit |= 1 << (hj->nsymbols-1);
		}
		n2 = t[n2].c[b];
		v >>= 1;
	    }
	    hj->jump = n2;
	    /*
	    printf("J4[%d][%d] = {'%.*s', %d}\n",
		   n, j, hj->nsymbols, hj->symbol, n2);
	    */
	}
    }

    cs->decode_t = t;
    cs->decode_J4 = J4;

    return 0;

 error:
    if (t)
	free(t);

    if (J4)
	free(J4);

    cs->decode_t = NULL;
    cs->decode_J4 = NULL;

    return -1;
}

/*
 * The opposite of huffman_encode().
 * Decode a huffman stream from 'block' using huffman codes 'c'.
 *
 * Returns: allocated block_t pointer on success
 *          NULL on failure.
 *
 * Method 1
 * --------
 *
 * At any node in our tree we can precompute a lookup table so that upon
 * reading the next 'k' bits we know the new node we'd end up in and what
 * symbols to export.
 * Then decoding simply works in fixed sets of k bits at a time.
 *
 * We use k=4 for efficient table space (they fit neatly in cache) and ease
 * of decoding 4-bits at a time. k=8 is about 20% faster as reading the input
 * byte by byte is easy, but the setup time is substantially longer
 * (16x at a guess) and the lookup tables no longer fit in the L1 cache.
 *
 * NB: This version also handles multiple interleaved huffman codes as
 * this support doesn't really slow down the decoding process.
 */
block_t *huffman_multi_decode(block_t *in, huffman_codeset_t *cs) {
    block_t *out = NULL;
    int i, j;
    int node_num;
    unsigned char *cp;
    huffman_codes_t **c;
    int nc;
    h_jump4_t (*J4)[16];
    htree_t *t;

    if (!cs)
	return NULL;
    c = cs->codes;
    nc = cs->ncodes;

    /* Ensure precomputed lookup tables exist */
    if (!cs->decode_t || !cs->decode_J4)
	if (-1 == init_decode_tables(cs))
	    return NULL;

    t  = cs->decode_t;
    J4 = cs->decode_J4;

    if (NULL == (out = block_create(NULL, 9*(in->alloc+1)))) {
	goto error;
    }

    /*
     * Decoding - part 1
     * We're part way through a byte, so decode bit by bit up to the next
     * whole byte and then we start the fast decoding section.
     */
    cp = out->data;
    node_num = 0;
    while (in->bit != 0) {
	int b = get_bits(in, 1);
	htree_t *t2 = &t[node_num];
	if (t2->l[b] != -1) {
	    if (t2->l[b] != SYM_EOF) {
		*cp++ = t2->l[b];
	    } else {
		out->byte = cp - out->data;
		goto success;
	    }
	}
	node_num = t2->c[b];
    }

    /*
     * Decoding - part 2
     *
     * We now handle data nibble by nibble using the nibble to get an
     * h_jump4_t lookup from the J4[] table.
     * If top_bit is clear then we know we have no funny business (SYM_EOF)
     * so we use a fast decoding technique, otherwise we have to do a slower
     * loop with a check.
     */
    {
	int last_node = node_num;
	unsigned char *last_cp = cp;
	h_jump4_t *x = &J4[node_num][in->data[in->byte] & 0x0f];
	int l = x->nsymbols;
	int b;

	/*
	 * This is the tight loop, so we over-optimise here by ignoring EOF
	 * and relying on knowing the length of the input data stream.
	 * This allows us to ignore the 9-bit data and only operate on
	 * the basic 0-255 symbols, glossing over the minor issue that EOF
	 * will look like an ordinary symbol.
	 */
	for (i = in->byte; i < in->alloc; i++) {
	    last_cp = cp;
	    last_node = node_num;

	    x = &J4[node_num][in->data[i] & 0x0f];
	    l = x->nsymbols;

	    /* printf("val=%d\n", in->data[i] & 0x0f); */
	    for (j = 0; j < l; j++) {
		*cp++ = x->symbol[j];
	    }
	    node_num = x->jump;

	    if (x->top_bit)
		break;

	    x = &J4[node_num][(in->data[i] >> 4) & 0x0f];
	    l = x->nsymbols;
	     
	    for (j = 0; j < l; j++) {	
		*cp++ = x->symbol[j];
	    }
	    node_num = x->jump;

	    if (x->top_bit)
		break;
	}

    
	/*
	 * Decoding - part 3
	 *
	 * The above optimisation has unfortunately added EOF to our data
	 * along with whatever else was packed in the last byte after the
	 * EOF symbol. So we rewind one byte and finish off decoding
	 * the slow way - walking the tree.
	 */
	cp = last_cp;
	node_num = last_node;
	in->byte = i;
	in->bit = 0;
	while (-1 != (b = get_bits(in, 1))) {
	    htree_t *t2 = &t[node_num];
	    if (t2->l[b] != -1) {
		if (t2->l[b] != SYM_EOF) {
		    *cp++ = t2->l[b];
		} else {
		    out->byte = cp - out->data;
		    goto success;
		}
	    }
	    node_num = t2->c[b];
	}
     }

 success:
    return out;

 error:
    if (out)
	block_destroy(out, 0);

    return NULL;
}

#if 0
/* A simple to understand (but slow) version of the above function */
block_t *huffman_multi_decode(block_t *in, huffman_codes_t **c, int nc) {
    block_t *out;
    htree_t (*t)[513]; /* number of internal nodes */
    int i, j, n, rec;
    int new_node, node_num;
    unsigned char *cp;
    int bC; /* byte count */

    if (NULL == (t = (htree_t (*)[513])malloc(nc * sizeof(*t))))
	return NULL;

    if (NULL == (out = block_create(NULL, 8*in->alloc+8))) {
	block_destroy(in, 0);
	free(t);
	return NULL;
    }

    /* Construct the tree from the codes */
    for (rec = 0; rec < nc; rec++) {
	new_node = 1;
	t[rec][0].l[0] = t[rec][0].l[1] = -1;
	t[rec][0].c[0] = t[rec][0].c[1] = 0;
	for (i = 0; i < c[rec]->ncodes; i++) {
	    int n = 0;
	    unsigned int v = c[rec]->codes[i].code;

	    for (j = 0; j < c[rec]->codes[i].nbits-1; j++) {
		int b = v & 1;
		if (t[rec][n].c[b]) {
		    n = t[rec][n].c[b];
		} else {
		    n = (t[rec][n].c[b] = new_node++);
		    t[rec][n].c[0] = t[rec][n].c[1] = 0;
		    t[rec][n].l[0] = t[rec][n].l[1] = -1;
		}
		v >>= 1;
	    }
	    /* last bit */
	    t[rec][n].l[v & 1] = c[rec]->codes[i].symbol;
	}
    }

    /*
     * Decoding - the slow way. How to speed up multi-code decoding?
     */
    cp = out->data;
    node_num = 0;
    bC = 0;
    while (in->byte < in->alloc) {
	htree_t *t2;
	int b = get_bits(in, 1);
	t2 = &t[bC%nc][node_num];
	if (t2->l[b] != -1) {
	    if (t2->l[b] != SYM_EOF) {
		*cp++ = t2->l[b];
		bC++;
		node_num = 0;
	    } else {
		out->byte = cp - out->data;
		free(t);
		return out;
	    }
	} else {
	    node_num = t2->c[b];
	}
    }

    /* We shouldn't reach here */
    free(t);
    return NULL;
}
#endif

/*
 * A slow version of the above huffman_decode function. This is designed as
 * a piecemeal decoder for purposes of restoring the huffman codes themselves.
 * NB: this only works for code lengths small enough to keep inside the
 * htab[] dimensions - IT DOES NOT CHECK THIS.
 *
 * Returns the next symbol
 *         -1 for failure
 */
int next_symbol(block_t *in, int *htab) {
    int b, v = 0, c = 1;
    while ((b = get_bits(in, 1)) != -1) {
	v = (v<<1) | b | (c <<= 1);

	if (htab[v] != -1)
	    return htab[v];
    }

    return -1;
}


/*
 * ---------------------------------------------------------------------------
 * Debug code. This turns the library into a stand-alone program for
 * easy debugging.x
 */
static void output_code_set(FILE *fp, huffman_codes_t *cds) {
    int i, j;
    int nbits_in = 0, nbits_out = 0;
    huffman_code_t *codes = cds->codes;
    int ncodes = cds->ncodes;

    fprintf(fp, "static huffman_code_t codes_FIXME[] = {\n");
    for (i = j = 0; i < ncodes; i++) {
	nbits_out += codes[i].nbits * codes[i].freq;
	nbits_in  += 8*codes[i].freq;
	if (j == 0)
	    fprintf(fp, "    ");
	if (codes[i].symbol == SYM_EOF) {
	    fprintf(fp, "{SYM_EOF,%3d}, ", codes[i].nbits);
	    j = 10;
	} else {
	    if (isalnum(codes[i].symbol)) {
		fprintf(fp, "{'%c',%3d}, ", codes[i].symbol, codes[i].nbits);
	    } else {
		fprintf(fp, "{%3d,%3d}, ", codes[i].symbol, codes[i].nbits);
	    }
	}
	j++;
	
	if (j >= 6) {
	    fputc('\n', fp);
	    j = 0;
	}
    }
    if (j)
	fputc('\n', fp);
    fprintf(fp, "};\n");
    fprintf(fp, "/* Expected compression to %f of input */\n",
	    (double)nbits_out/nbits_in);
}

static void output_code_set2(FILE *fp, huffman_codes_t *cds) {
    int i;
    huffman_code_t *codes = cds->codes;
    int ncodes = cds->ncodes;

    fprintf(fp, "huffman_code_t = {\n");
    for (i = 0; i < ncodes; i++) {
	fprintf(fp, "\t%d:\t%3d %c %2d %04x %d\n",
		i,codes[i].symbol, codes[i].symbol,
		codes[i].nbits, codes[i].code,
		codes[i].freq);
    }
    fprintf(fp, "};\n");
}


/*
 * --------------------------------------------------------------------------
 * A test main() to create an application capable of compressing and
 * uncompressing stdin.
 */

#ifdef TEST_MAIN
#include <fcntl.h>
/* #include <zlib.h> */

/*
 * Slurps the entirety of stdin into a malloced buffer and returns a pointer
 * to it.
 *
 * Returns: malloced buffer on success, *lenp equal to length
 *          NULL on failure
 */
static unsigned char *load(int *lenp, char *fn) {
    unsigned char *data = NULL;
    int dsize = 0;
    int dcurr = 0, len;
    int fd = 0;

    if (fn) {
	if (-1 == (fd = open(fn, O_RDONLY, 0))) {
	    perror(fn);
	    return NULL;
	}
    }

    do {
	if (dsize - dcurr < 8192) {
	    dsize = dsize ? dsize * 2 : 8192;
	    if (NULL == (data = realloc(data, dsize))) {
		if (fd)
		    close(fd);
		return NULL;
	    }
	}

	len = read(fd, data + dcurr, 8192);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
	if (fd)
	    close(fd);
	return NULL;
    }

    if (fd)
	close(fd);

    *lenp = dcurr;
    return data;
}

/*
 * Returns 0 for success
 *        -1 for failure.
 */
int decode_main(unsigned char *data, int len, int code_set) {
    huffman_codeset_t *cs = NULL;
    block_t *blk_in;
    unsigned char *out;
    int err, out_len, i;
    int bfinal;

    blk_in = block_create(NULL, 1000 + len);

    if (code_set != 0) {
	cs = generate_code_set(code_set, 1,  /* no. codes */
			       NULL, 0,      /* data + size */
			       1,	     /* eof */
			       MAX_CODE_LEN,
			       0);	     /* all_codes */
	store_codes(blk_in, cs, 1);
    }

    if (blk_in->bit != 0) {
	blk_in->data[blk_in->byte] |= data[0];
	memcpy(&blk_in->data[blk_in->byte+1], data+1, len-1);
    } else {
	memcpy(&blk_in->data[blk_in->byte], data, len);
    }

    /* Do the decoding */
    do {
	block_t *out;
	if (!cs)
	    cs = restore_codes(blk_in, &bfinal);
	out = huffman_multi_decode(blk_in, cs);
	write(1, out->data, out->byte);
	block_destroy(out, 0);
	huffman_codeset_destroy(cs);
	cs = NULL;
    } while (!bfinal);

    block_destroy(blk_in, 0);

    return 0;
}

/*
 * Returns 0 for success
 *        -1 for failure.
 */
int encode_main(unsigned char *data, int len, int code_set, int rec_size,
		int blk_size, int dump_tree, int exit_after_tree) {
    /* Encoding */
    unsigned char *d2 = data;
    block_t *blk;
    huffman_codeset_t *cs;

    blk = block_create(NULL, 8192);
    fprintf(stderr, "Input %d bytes\n", len);

    do {
	int l2 = len > blk_size ? blk_size : len;
	int i;

	if (code_set != 0)
	    l2 = len; /* predefined code-sets have final-block bit set */
	cs = generate_code_set(code_set, rec_size,
			       d2, l2, /* Data and length */
			       1,      /* eof */
			       MAX_CODE_LEN,
			       0);     /* all codes */
	if (!cs)
	    return -1;

	if (dump_tree) {
	    for (i = 0; i < rec_size; i++) {
		printf("==Sub-set %d==\n", i);
		output_code_set(stdout, cs->codes[i]);
		/* output_code_set2(stdout, cs->codes[i]); */
	    }
	    if (exit_after_tree)
		return 0;
	}

	store_codes(blk, cs, l2 == len);
	if (code_set != 0) {
	    blk->data[blk->byte = 0] = 0;  /* starting bit no. preseved */
	} else {
	    /*
	      fprintf(stderr, "codes finished at %d bytes, %d bits\n",
		    blk->byte, blk->bit);
	    */
	}
	if (exit_after_tree) {
	    write(1, blk->data, blk->byte + (blk->bit != 0));
	    return 0;
	}

	huffman_multi_encode(blk, cs, code_set, d2, l2);

	huffman_codeset_destroy(cs);
	len -= l2;
	d2  += l2;

    } while (len > 0);
    write(1, blk->data, blk->byte + (blk->bit != 0));

    fprintf(stderr, "Output %ld bytes\n", blk->byte + (blk->bit != 0));
    block_destroy(blk, 0);

    return 0;
}

int main(int argc, char **argv) {
    unsigned char *data;
    int decode = 0;
    int dump_tree = 0;
    int exit_after_tree = 0;
    int code_set = CODE_INLINE;
    int blk_size = 0x7fff;
    int rec_size = 1;
    int c, r, len;
    char *fn = NULL;

    while ((c = getopt(argc, argv, "c:detxl:b:hr:i:")) != -1) {
	switch (c) {
	case 'b':
	    blk_size = atoi(optarg);
	    break;

	case 'c':
	    code_set = atoi(optarg);
	    break;

	case 'r':
	    rec_size = atoi(optarg);
	    break;

	case 'd':
	    decode = 1;
	    break;

	case 'e':
	    decode = 0;
	    break;

	case 't':
	    dump_tree = 1;
	    break;

	case 'x':
	    exit_after_tree = 1;
	    break;
	    
	case 'i':
	    fn = optarg;
	    break;

	default:
	    fprintf(stderr, "Usage: huffman_static [options] < stdin > stdout\n");
	    fprintf(stderr, "    Decoding options\n");
	    fprintf(stderr, "        -d\tdecode flag\n");
	    fprintf(stderr, "    Encoding options\n");
	    fprintf(stderr, "        -e\tencode flag\n");
	    fprintf(stderr, "        -b size\tspecify the block-size\n");
	    fprintf(stderr, "        -c code\tspecify code-set. 0 => inline\n");
	    fprintf(stderr, "        -t\tpretty-print the code-set used\n");
	    fprintf(stderr, "        -x\texit after outputting code-set\n");
	    exit(c != 'h');
	}
    }

    data = load(&len, fn);

    r = decode
	? decode_main(data, len, code_set)
	: encode_main(data, len, code_set, rec_size,
		      blk_size, dump_tree, exit_after_tree);

    free(data);

    return r == 0 ? 0 : 1;
}
#endif
