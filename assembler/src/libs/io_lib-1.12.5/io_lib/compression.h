#ifndef _COMPRESSION_H_
#define _COMPRESSION_H_

#include "io_lib/os.h"
#include <zlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * zlib_huff()
 *
 * Compresses data using huffman encoding, as implemented by zlib.
 *
 * Arguments:
 *	uncomp		Uncompressed input data
 *	uncomp_len	Length of uncomp data
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *zlib_huff(char *uncomp, int uncomp_len, int strategy, int *comp_len);

/*
 * zlib_dehuff()
 *
 * Uncompresses data using huffman encoding, as implemented by zlib.
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *zlib_dehuff(char *comp, int comp_len, int *uncomp_len);

/*
 * zlib_dehuff2()
 *
 * Uncompresses data using huffman encoding, as implemented by zlib.
 * Similar to zlib_dehuff above, but with the following differences:
 *
 * 1) It pastes together the zlib stream from two components; comp1+comp2
 *    with the last byte of comp1 overlapping (ORed) with the first byte
 *    of comp2. This allows for separation of the huffman codes from
 *    the compressed data itself.
 * 2) It uses the raw Deflate format rather than Zlib's wrapping of it.
 * 3) It uses an EOF symbol to mark the end rather than encoding the
 *    uncompressed size in the header
 * 
 *
 * Arguments:
 *	comp1		Compressed input data part 1
 *	comp1_len	Length of comp1 data
 *	comp2		Compressed input data part 2
 *	comp2_len	Length of comp2 data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *zlib_dehuff2(char *comp1, int comp1_len,
		   char *comp2, int comp2_len,
		   int *uncomp_len);

/*
 * Run length encoding.
 *
 * Any run of 3 or more identical characters (up to 255 in a row) are replaced
 * by a 'guard' byte followed by the number of characters followed by
 * the character value itself.
 * Any single guard value in the input is escaped using 'guard 0'.
 *
 * Specifying guard as -1 will automatically pick one of the least used
 * characters in the input as the guard.
 *
 * Arguments:
 *	uncomp		Input data
 *	uncomp_len	Length of input data 'uncomp'
 *	guard		Guard byte - used to encode "N" copies of data
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *rle(char *uncomp, int uncomp_len, int guard, int *comp_len);

/*
 * Reverses run length encoding.
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *unrle(char *comp, int comp_len, int *uncomp_len);

/*
 * Mutli-byte run length encoding.
 *
 * Any run of 3 or more identical characters (up to 255 in a row) are replaced
 * by a 'guard' byte followed by the number of characters followed by
 * the character value itself.
 * Any single guard value in the input is escaped using 'guard 0'.
 *
 * Specifying guard as -1 will automatically pick one of the least used
 * characters in the input as the guard.
 *
 * Arguments:
 *	uncomp		Input data
 *	uncomp_len	Length of input data 'uncomp'
 *	guard		Guard byte - used to encode "N" copies of data
 *      rsz             Size of blocks to compare for run checking.
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *xrle(char *uncomp, int uncomp_len, int guard, int rsz, int *comp_len);

/*
 * Reverses multi-byte run length encoding.
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *unxrle(char *comp, int comp_len, int *uncomp_len);

/*
 * Mutli-byte run length encoding.
 *
 * Steps along in words of size 'rsz'. Unlike XRLE above this does run-length
 * encoding by writing out an additional "length" word every time 2 or more
 * words in a row are spotted. This removes the need for a guard byte.
 *
 * Additionally this method ensures that both input and output formats remain
 * aligned on words of size 'rsz'.
 *
 * Arguments:
 *	uncomp		Input data
 *	uncomp_len	Length of input data 'uncomp'
 *      rsz             Size of blocks to compare for run checking.
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data if successful
 *	NULL if not successful
 */
char *xrle2(char *uncomp, int uncomp_len, int rsz, int *comp_len);

/*
 * Reverses multi-byte run length encoding (xrle_new).
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Uncompressed data if successful
 *	NULL if not successful
 */
char *unxrle2(char *comp, int comp_len, int *uncomp_len);

/*
 * decorrelate1()
 *
 * Produce successive deltas from a 1-byte array.
 *
 * Arguments:
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data
 *	level		Differencing level (must be 1, 2 or 3)
 *	comp_len	Return: where to store new compressed length
 *
 * Returns:
 *	Success: A decorrelated buffer (malloced)
 *	Failure: NULL
 */
char *decorrelate1(char *uncomp, int uncomp_len, int level, int *comp_len);
char *decorrelate1dyn(char *s_uncomp, int uncomp_len, int *comp_len);

/*
 * recorrelate1()
 *
 * The reverse of decorrelate1()
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Success: uncompressed data
 *	Failure: NULL
 */
char *recorrelate1(char *comp, int comp_len, int *uncomp_len);

/*
 * decorrelate2()
 *
 * Produce successive deltas from a 2-byte array (big endian)
 *
 * Arguments:
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data
 *	level		Differencing level (must be 1, 2 or 3)
 *	comp_len	Return: where to store new compressed length
 *
 * Returns:
 *	Success: A decorrelated buffer (malloced)
 *	Failure: NULL
 */
char *decorrelate2(char *uncomp, int uncomp_len, int level, int *comp_len);
char *decorrelate2dyn(char *s_uncomp, int uncomp_len, int *comp_len);

/*
 * recorrelate2()
 *
 * The reverse of decorrelate2()
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Success: uncompressed data
 *	Failure: NULL
 */
char *recorrelate2(char *comp, int comp_len, int *uncomp_len);

/*
 * decorrelate4()
 *
 * Produce successive deltas from a 4-byte array (big endian)
 *
 * Arguments:
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data
 *	level		Differencing level (must be 1, 2 or 3)
 *	comp_len	Return: where to store new compressed length
 *
 * Returns:
 *	Success: A decorrelated buffer (malloced)
 *	Failure: NULL
 */
char *decorrelate4(char *uncomp, int uncomp_len, int level, int *comp_len);

/*
 * recorrelate4()
 *
 * The reverse of decorrelate4()
 *
 * Arguments:
 *	comp		Compressed input data
 *	comp_len	Length of comp data
 *	uncomp_len	Output: length of uncompressed data
 *
 * Returns:
 *	Success: uncompressed data
 *	Failure: NULL
 */
char *recorrelate4(char *comp, int comp_len, int *uncomp_len);

/*
 * shrink_16to8()
 *
 * Stores an array of 16-bit (big endian) array elements in an 8-bit array.
 * We assume that most 16-bit elements encode numbers that fit in an 8-bit
 * value. When not possible, we store a marker followed by the 16-bit value
 * stored as multiple 8-bit values.
 *
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data (in bytes)
 *	comp_len	Return: where to store new compressed length
 *	
 * Returns:
 *	Success: An 8-bit array (malloced)
 *	Failure: NULL
 */
char *shrink_16to8(char *uncomp, int uncomp_len, int *comp_len);

/*
 * expand_8to16()
 *
 * The opposite of the shrink_16to8() function.
 *
 *	comp		Compressed input data
 *	comp_len	Length of comp data (in bytes)
 *	uncomp_len	Output: length of uncompressed data (in bytes)
 *	
 * Returns:
 *	Success: Uncompressed data (char *)
 *	Failure: NULL
 */
char *expand_8to16(char *comp, int comp_len, int *uncomp_len);

/*
 * shrink_32to8()
 *
 * Stores an array of 32-bit (big endian) array elements in an 8-bit array.
 * We assume that most 32-bit elements encode numbers that fit in an 8-bit
 * value. When not possible, we store a marker followed by the 32-bit value
 * stored as multiple 8-bit values.
 *
 *	uncomp		Uncompressed data
 *	uncomp_len	Length of uncompressed data (in bytes)
 *	comp_len	Return: where to store new compressed length
 *	
 * Returns:
 *	Success: An 8-bit array (malloced)
 *	Failure: NULL
 */
char *shrink_32to8(char *uncomp, int uncomp_len, int *comp_len);

/*
 * expand_8to32()
 *
 * The opposite of the shrink_32to8() function.
 *
 *	comp		Compressed input data
 *	comp_len	Length of comp data (in bytes)
 *	uncomp_len	Output: length of uncompressed data (in bytes)
 *	
 * Returns:
 *	Success: Uncompressed data (char *)
 *	Failure: NULL
 */
char *expand_8to32(char *comp, int comp_len, int *uncomp_len);

char *follow1(char *s_uncomp,
	      int uncomp_len,
	      int *comp_len);

char *unfollow1(char *s_comp,
		int comp_len,
		int *uncomp_len);

char *ichebcomp(char *uncomp,
		int uncomp_len,
		int *data_len);

char *ichebuncomp(char *comp,
		  int comp_len,
		  int *uncomp_len);

/*
 * This is a LOSSY compression. It replaces N with 10 * log2(N).
 */
char *log2_data(char *x_uncomp,
		int uncomp_len,
		int *comp_len);

char *unlog2_data(char *x_comp,
		  int comp_len,
		  int *uncomp_len);

/*
 * Implements compression using a set of static huffman codes stored using
 * the Deflate algorithm (and so in this respect it's similar to zlib).
 *
 * The huffman codes though can be previously stored in the ztr object
 * using ztr_add_hcode(). "cset" indicates which numbered stored huffman
 * code set is to be used, or passing zero will use inline codes (ie they
 * are stored in the data stream itself, just as in standard deflate).
 *
 * Arguments:
 *	ztr		ztr_t pointer; used to find stored code-sets
 *	uncomp		The uncompressed input data
 *	uncomp_len	Length of uncomp
 *	cset		Stored code-set number, zero for inline
 *	recsz		Record size - only used when cset == 0.
 *	comp_len	Output: length of compressed data
 *
 * Returns:
 *	Compressed data stream if successful + comp_len
 *      NULL on failure
 */
char *sthuff(ztr_t *ztr, char *uncomp, int uncomp_len, 
	     int cset, int recsz, int *comp_len);
char *unsthuff(ztr_t *ztr, char *comp, int comp_len, int *uncomp_len);

/*
 * Reorders quality data from its RAW format to an interleaved 4-byte
 * aligned format.
 *
 * Starting with sequence A1 C2 G3 the raw format is quality of called
 * bases followed by quality of remaining bases:
 * 0 (RAW format)
 * Q(A1) Q(C2) Q(G3)
 * Q(C2) Q(A2) Q(A3) 
 * Q(G2) Q(G2) Q(C3) 
 * Q(T2) Q(T2) Q(T3) 
 *
 * We reorder it to:
 * ZTR_FORM_QSHIFT <any> <any> 0(raw)
 * Q(A1) Q(C1) Q(G1) Q(T1)
 * Q(C2) Q(A2) Q(G2) Q(T2)
 * Q(G3) Q(A3) Q(C3) Q(T3)
 * 
 * Returns shifted data on success
 *         NULL on failure
 */
char *qshift(char *qold, int qlen, int *new_len);
char *unqshift(char *qold, int qlen, int *new_len);

/*
 * Given a sequence ACTG this shifts trace data from the order:
 *
 *     A1A2A3A4 C1C2C3C4 G1G2G3G4 T1T2T3T4
 *
 * to
 *
 *     A1C1G1T1 C2A2G2T2 T3A3C3G3 G4C4C4T4
 *
 * Ie for each base it ouputs the signal for the called base first
 * followed by the remaining 3 signals in A,C,G,T order (minus the
 * called signal already output).
 */
char *tshift(ztr_t *ztr, char *told_c, int tlen, int *new_len);
char *untshift(ztr_t *ztr, char *told_c, int tlen, int *new_len);

#ifdef __cplusplus
}
#endif

#endif /* _COMPRESSION_H_ */
