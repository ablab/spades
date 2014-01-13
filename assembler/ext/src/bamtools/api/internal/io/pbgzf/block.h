#ifndef BLOCK_H_
#define BLOCK_H_

typedef struct {
    int8_t *buffer;
    int32_t block_length;
    int32_t block_offset; // used by bgzf_write and bgzf_flush
    int64_t block_address; // used by bgzf_write and bgzf_flush
    int64_t id; // Used by the queue
    int32_t mem;
} block_t;

block_t*
block_init();

void
block_destroy(block_t *block);

typedef struct {
    int32_t head;
    int32_t n;
    int32_t m;
    block_t **blocks;
    pthread_mutex_t *mut;
} block_pool_t;

// pool is full, with mutex
block_pool_t*
block_pool_init(int32_t m);

// no blocks in the pool, no mutex
block_pool_t*
block_pool_init2(int32_t m);

// add to the end of the queue
int32_t
block_pool_add(block_pool_t *pool, block_t *block);

// get from the front of the queue
block_t*
block_pool_get(block_pool_t *pool);

block_t*
block_pool_peek(block_pool_t *pool);

void
block_pool_destroy(block_pool_t *pool);

void
block_pool_reset(block_pool_t *pool);

#endif
