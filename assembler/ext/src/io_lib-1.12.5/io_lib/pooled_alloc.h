#ifndef _POOLED_ALLOC_H_
#define _POOLED_ALLOC_H_

/*
 * Implements a pooled block allocator where all items are the same size,
 * but we need many of them.
 */
typedef struct {
    void   *pool;
    size_t  used;
} pool_t;

typedef struct {
    size_t dsize;
    size_t npools;
    pool_t *pools;
    void *free;
} pool_alloc_t;

pool_alloc_t *pool_create(size_t dsize);
void pool_destroy(pool_alloc_t *p);
void *pool_alloc(pool_alloc_t *p);
void pool_free(pool_alloc_t *p, void *ptr);


#endif /*_POOLED_ALLOC_H_*/
