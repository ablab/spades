#ifndef QUEUE_H_
#define QUEUE_H_

#define QUEUE_DEBUG

enum {
    QUEUE_STATE_OK = 0,
    QUEUE_STATE_EOF = 1,
    QUEUE_STATE_FLUSH = 2
};

typedef struct {
    block_t **queue;
    int32_t mem;
    int32_t head;
    int32_t tail;
    int32_t n; 
    int32_t length;
    int64_t id;
    int8_t ordered;
    int32_t num_adders;
    int32_t num_getters;
    pthread_mutex_t *mut;
    pthread_cond_t *not_full;
    pthread_cond_t *not_empty;
    pthread_cond_t *is_empty;
    pthread_cond_t *not_flush;
    int8_t state;
#ifdef QUEUE_DEBUG
    int32_t num_waiting[4];
#endif
} queue_t;

queue_t*
queue_init(int32_t capacity, int8_t ordered, int32_t num_adders, int32_t num_getters);

int8_t
queue_add(queue_t *q, block_t *b, int8_t wait);

block_t*
queue_get(queue_t *q, int8_t wait);

void
queue_signal(queue_t *q);

// TODO
/*
int32_t
queue_add_batch(queue_t *q, block_pool_t *pool, int8_t wait);

int32_t
queue_get_batch(queue_t *q, block_pool_t *pool, int8_t wait);
*/

void
queue_wait_until_empty(queue_t *q);

void
queue_wait_until_not_flush(queue_t *q);

#define queue_set_flush(_q) (_q->state = QUEUE_STATE_FLUSH)

void
queue_remove_flush(queue_t *q);

void
queue_close(queue_t *q);

void
queue_destroy(queue_t *q);

void
queue_reset(queue_t *q, int32_t num_adders, int32_t num_getters);

void
queue_wake_all(queue_t *q);

void 
queue_remove_adder(queue_t *q);

void 
queue_remove_getter(queue_t *q);

// DEBUG
void
queue_print_status(queue_t *q, FILE *fp);

#endif
