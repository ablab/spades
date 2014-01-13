#ifndef WRITER_H_
#define WRITER_H_

typedef struct {
    BGZF *fp_bgzf;
    FILE *fp_file;
    queue_t *output;
    uint8_t compress;
    int32_t compress_type;
    uint8_t is_done;
    uint8_t is_closed;
    block_pool_t *pool_fp;
    block_pool_t *pool_local;
} writer_t;

writer_t*
writer_init(int fd, queue_t *output, uint8_t compress, int32_t compress_level, int32_t compress_type, block_pool_t *pool);

void*
writer_run(void *arg);

void
writer_destroy(writer_t *w);

void
writer_reset(writer_t *w);

#endif
