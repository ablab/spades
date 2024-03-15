#ifndef PBGZF_H_
#define PBGZF_H_

#include <pthread.h>
#include "bgzf.h"
#include "block.h"
#include "queue.h"
#include "reader.h"
#include "consumer.h"
#include "writer.h"
#include "util.h"

#define DISABLE_BZ2

#ifndef DISABLE_BZ2
#define BZ2_DEFAULT_LEVEL 9
#endif

#ifdef __cplusplus 
extern "C" {
#endif

typedef struct {
    pthread_attr_t attr;
    pthread_t *threads;
    int32_t n;
    consumer_t **c;
} consumers_t;

typedef struct {
    pthread_attr_t attr;
    pthread_t thread;
    reader_t *r;
} producer_t;

typedef struct {
    pthread_attr_t attr;
    pthread_t thread;
    writer_t *w;
} outputter_t;

typedef struct {
    block_t *block; // buffer block
    block_pool_t *pool;

    char open_mode;
    int32_t is_write;
    int32_t queue_size;
    int32_t num_threads;
    int32_t block_offset; // for pbgzf_flush_tr
    int64_t block_address; // for pbgzf_flush_tr
    int32_t eof_ok; // for pbgzf_check_EOF
    int32_t eof;
    int64_t n_blocks;

    queue_t *input;
    queue_t *output;
    reader_t *r;
    writer_t *w;
    consumers_t *c;
    producer_t *p;
    outputter_t *o;
} PBGZF;

#ifdef __cplusplus
extern "C" {
#endif

#define PBGZF_QUEUE_SIZE 100
#define PBGZF_BLOCKS_POOL_NUM 100
//#define PBGZF_USE_LOCAL_POOLS 

/*
 * Sets the number of consumer threads per file handle.  If set to
 * X, the numebr of threads will be X+1, due to the need to have a 
 * single reader or writer thread in addition to the consumer 
 * thread(s).
 */
void
pbgzf_set_num_threads_per(int32_t n);

/*
 * Open an existing file descriptor for reading or writing.
 * Mode must be either "r" or "w".
 * A subsequent pbgzf_close will not close the file descriptor.
 * Returns null on error.
 */
PBGZF* pbgzf_fdopen(int fd, const char* __restrict mode);

/*
 * Open the specified file for reading or writing.
 * Mode must be either "r" or "w".
 * Returns null on error.
 */
PBGZF* pbgzf_open(const char* path, const char* __restrict mode);

/*
 * Close the BGZ file and free all associated resources.
 * Does not close the underlying file descriptor if created with pbgzf_fdopen.
 * Returns zero on success, -1 on error.
 */
int pbgzf_close(PBGZF* fp);

/*
 * Read up to length bytes from the file storing into data.
 * Returns the number of bytes actually read.
 * Returns zero on end of file.
 * Returns -1 on error.
 */
int pbgzf_read(PBGZF* fp, void* data, int length);

/*
 * Write length bytes from data to the file.
 * Returns the number of bytes written.
 * Returns -1 on error.
 */
int pbgzf_write(PBGZF* fp, const void* data, int length);

/*
 * Return a virtual file pointer to the current location in the file.
 * No interpetation of the value should be made, other than a subsequent
 * call to pbgzf_seek can be used to position the file at the same point.
 * Return value is non-negative on success.
 * Returns -1 on error.
 */
int64_t pbgzf_tell(PBGZF *fp);

/*
 * Set the file to read from the location specified by pos, which must
 * be a value previously returned by pbgzf_tell for this file (but not
 * necessarily one returned by this file handle).
 * The where argument must be SEEK_SET.
 * Seeking on a file opened for write is not supported.
 * Returns zero on success, -1 on error.
 */
int64_t pbgzf_seek(PBGZF* fp, int64_t pos, int where);

int pbgzf_check_EOF(PBGZF *fp);

int pbgzf_flush(PBGZF* fp);

int pbgzf_flush_try(PBGZF *fp, int size);

void pbgzf_set_cache_size(PBGZF *fp, int cache_size);

#ifdef __cplusplus
}
#endif

void
pbgzf_main(int f_src, int f_dst, int compress, int compress_level, int compress_type, int queue_size, int num_threads);

#ifdef __cplusplus 
}
#endif

#endif
