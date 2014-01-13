/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/types.h>

#define BGZF_BLOCK_SIZE     0xff00 // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE 0x10000

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

#define BGZF_ERR_ZLIB   1
#define BGZF_ERR_HEADER 2
#define BGZF_ERR_IO     4
#define BGZF_ERR_MISUSE 8

#ifdef _USE_KNETFILE
#include "knetfile.h"
typedef knetFile *_bgzf_file_t;
#define _bgzf_open(fn, mode) knet_open(fn, mode)
#define _bgzf_dopen(fp, mode) knet_dopen(fp, mode)
#define _bgzf_close(fp) knet_close((_bgzf_file_t)fp)
#define _bgzf_fileno(fp) (((_bgzf_file_t)fp)->fd)
#define _bgzf_tell(fp) knet_tell((_bgzf_file_t)fp)
#define _bgzf_seek(fp, offset, whence) knet_seek((_bgzf_file_t)fp, offset, whence)
#define _bgzf_read(fp, buf, len) knet_read((_bgzf_file_t)fp, buf, len)
#define _bgzf_write(fp, buf, len) knet_write((_bgzf_file_t)fp, buf, len)
#else // ~defined(_USE_KNETFILE)
#if defined(_WIN32) || defined(_MSC_VER)
#define ftello(fp) ftell((_bgzf_file_t)fp)
#define fseeko(fp, offset, whence) fseek((_bgzf_file_t)fp, offset, whence)
#else // ~defined(_WIN32)
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif // ~defined(_WIN32)
typedef FILE *_bgzf_file_t;
#define _bgzf_open(fn, mode) fopen(fn, mode)
#define _bgzf_dopen(fp, mode) fdopen(fp, mode)
#define _bgzf_close(fp) fclose((_bgzf_file_t)fp)
#define _bgzf_fileno(fp) fileno((_bgzf_file_t)fp)
#define _bgzf_tell(fp) ftello((_bgzf_file_t)fp)
#define _bgzf_seek(fp, offset, whence) fseeko((_bgzf_file_t)fp, offset, whence)
#define _bgzf_read(fp, buf, len) fread(buf, 1, len, (_bgzf_file_t)fp)
#define _bgzf_write(fp, buf, len) fwrite(buf, 1, len, (_bgzf_file_t)fp)
#endif // ~define(_USE_KNETFILE)


typedef struct {
    int errcode:16, is_write:2, compress_level:14;
    int cache_size;
    int block_length, block_offset;
    int64_t block_address;
    void *uncompressed_block, *compressed_block;
    void *cache; // a pointer to a hash table
    void *fp; // actual file handler; FILE* on writing; FILE* or knetFile* on reading
    void *mt; // only used for multi-threading
} BGZF;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

	/******************
	 * Basic routines *
	 ******************/

	/**
	 * Open an existing file descriptor for reading or writing.
	 *
	 * @param fd    file descriptor
	 * @param mode  mode matching /[rwu0-9]+/: 'r' for reading, 'w' for writing and a digit specifies
	 *              the zlib compression level; if both 'r' and 'w' are present, 'w' is ignored.
	 * @return      BGZF file handler; 0 on error
	 */
	BGZF* bgzf_dopen(int fd, const char *mode);

	#define bgzf_fdopen(fd, mode) bgzf_dopen((fd), (mode)) // for backward compatibility

	/**
	 * Open the specified file for reading or writing.
	 */
	BGZF* bgzf_open(const char* path, const char *mode);

	/**
	 * Close the BGZF and free all associated resources.
	 *
	 * @param fp    BGZF file handler
	 * @return      0 on success and -1 on error
	 */
	int bgzf_close(BGZF *fp);

	/**
	 * Read up to _length_ bytes from the file storing into _data_.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to read into
	 * @param length size of data to read
	 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
	 */
	ssize_t bgzf_read(BGZF *fp, void *data, ssize_t length);

	/**
	 * Write _length_ bytes from _data_ to the file.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to write
	 * @param length size of data to write
	 * @return       number of bytes actually written; -1 on error
	 */
	ssize_t bgzf_write(BGZF *fp, const void *data, ssize_t length);

	/**
	 * Write the data in the buffer to the file.
	 */
	int bgzf_flush(BGZF *fp);

	/**
	 * Return a virtual file pointer to the current location in the file.
	 * No interpetation of the value should be made, other than a subsequent
	 * call to bgzf_seek can be used to position the file at the same point.
	 * Return value is non-negative on success.
	 */
	#define bgzf_tell(fp) ((fp->block_address << 16) | (fp->block_offset & 0xFFFF))

	/**
	 * Set the file to read from the location specified by _pos_.
	 *
	 * @param fp     BGZF file handler
	 * @param pos    virtual file offset returned by bgzf_tell()
	 * @param whence must be SEEK_SET
	 * @return       0 on success and -1 on error
	 */
	int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence);

	/**
	 * Check if the BGZF end-of-file (EOF) marker is present
	 *
	 * @param fp    BGZF file handler opened for reading
	 * @return      1 if EOF is present; 0 if not or on I/O error
	 */
	int bgzf_check_EOF(BGZF *fp);

	/**
	 * Check if a file is in the BGZF format
	 *
	 * @param fn    file name
	 * @return      1 if _fn_ is BGZF; 0 if not or on I/O error
	 */
	 int bgzf_is_bgzf(const char *fn);

	/*********************
	 * Advanced routines *
	 *********************/

	/**
	 * Set the cache size. Only effective when compiled with -DBGZF_CACHE.
	 *
	 * @param fp    BGZF file handler
	 * @param size  size of cache in bytes; 0 to disable caching (default)
	 */
	void bgzf_set_cache_size(BGZF *fp, int size);

	/**
	 * Flush the file if the remaining buffer size is smaller than _size_ 
	 */
	int bgzf_flush_try(BGZF *fp, ssize_t size);

	/**
	 * Read one byte from a BGZF file. It is faster than bgzf_read()
	 * @param fp     BGZF file handler
	 * @return       byte read; -1 on end-of-file or error
	 */
	int bgzf_getc(BGZF *fp);

	/**
	 * Read one line from a BGZF file. It is faster than bgzf_getc()
	 *
	 * @param fp     BGZF file handler
	 * @param delim  delimitor
	 * @param str    string to write to; must be initialized
	 * @return       length of the string; 0 on end-of-file; negative on error
	 */
	int bgzf_getline(BGZF *fp, int delim, kstring_t *str);

	/**
	 * Read the next BGZF block.
	 */
	int bgzf_read_block(BGZF *fp);

	/**
	 * Enable multi-threading (only effective on writing)
	 *
	 * @param fp          BGZF file handler; must be opened for writing
	 * @param n_threads   #threads used for writing
	 * @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
	 */
	int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);

        inline void
          packInt16(uint8_t* buffer, uint16_t value);
        inline int
          unpackInt16(const uint8_t* buffer);
        inline void
          packInt32(uint8_t* buffer, uint32_t value);
        int
          bgzf_check_header(const uint8_t* header);
        int bgzf_compress(void *_dst, int *dlen, void *src, int slen, int level);

#ifdef __cplusplus
}
#endif

#endif
