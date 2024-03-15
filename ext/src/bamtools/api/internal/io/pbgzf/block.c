#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "bgzf.h"
#include "util.h"
#include "block.h"

block_t*
block_init()
{
  block_t *b = calloc(1, sizeof(block_t));
  b->buffer = malloc(sizeof(int8_t) * BGZF_MAX_BLOCK_SIZE);
  b->id = -1;
  return b;
}

void
block_destroy(block_t *block)
{
  if(NULL == block) return;
  free(block->buffer);
  free(block);
}

block_pool_t*
block_pool_init2(int32_t m)
{
  int32_t i;
  block_pool_t *pool = calloc(1, sizeof(block_pool_t));

  pool->m = m;
  pool->n = 0;
  pool->head = 0;
  pool->blocks = calloc(pool->m, sizeof(block_pool_t*));
  for(i=0;i<pool->m;i++) {
      pool->blocks[i] = NULL;
  }

  return pool;
}

block_pool_t*
block_pool_init(int32_t m)
{
  int32_t i;
  block_pool_t *pool = NULL;

  pool = block_pool_init2(m);
  for(i=0;i<pool->m;i++) {
      pool->blocks[i] = block_init();
  }
  pool->n = m;
  
  pool->mut = calloc(1, sizeof(pthread_mutex_t));
  if(0 != pthread_mutex_init(pool->mut, NULL)) {
      fprintf(stderr, "Could not create mutex\n");
      exit(1);
  }

  return pool;
}

static void
block_pool_shift(block_pool_t *pool)
{
  int32_t i;
  /*
  fprintf(stderr, "SHIFTING head=%d n=%d m=%d\n",
          pool->head, pool->n, pool->m);
          */
  if(0 == pool->head) return;
  for(i=0;i<pool->n;i++) {
      pool->blocks[i] = pool->blocks[i+pool->head];
      pool->blocks[i+pool->head] = NULL;
  }
  pool->head = 0;
  /*
  fprintf(stderr, "SHIFTED head=%d n=%d m=%d\n",
          pool->head, pool->n, pool->m);
          */
}

int32_t
block_pool_add(block_pool_t *pool, block_t *block)
{
  if(NULL == block) return 0;
  if(NULL != pool->mut) safe_mutex_lock(pool->mut);
  if(pool->n < pool->m) {
      if(pool->head + pool->n == pool->m) {
          block_pool_shift(pool);
      }
      pool->blocks[pool->head+pool->n] = block;
      pool->n++;
  }
  else { // ignore
      block_destroy(block);
  }
  if(NULL != pool->mut) safe_mutex_unlock(pool->mut);
  return 1;
}

block_t*
block_pool_get(block_pool_t *pool)
{
  block_t *b = NULL;
  if(NULL != pool->mut) safe_mutex_lock(pool->mut);
  if(0 < pool->n) {
      b = pool->blocks[pool->head];
      pool->blocks[pool->head] = NULL;
      pool->n--;
      pool->head++;
      if(pool->head == pool->m) pool->head = 0; // NB: pool->n == 0
  }
  if(NULL != pool->mut) safe_mutex_unlock(pool->mut);
  return b;
}

block_t*
block_pool_peek(block_pool_t *pool)
{
  block_t *b = NULL;
  if(NULL != pool->mut) safe_mutex_lock(pool->mut);
  if(0 < pool->n) {
      b = pool->blocks[pool->head];
  }
  if(NULL != pool->mut) safe_mutex_unlock(pool->mut);
  return b;
}

void
block_pool_destroy(block_pool_t *pool)
{
  int32_t i;
  if(NULL == pool) return;
  for(i=0;i<pool->m;i++) {
      block_destroy(pool->blocks[i]);
  }
  free(pool->blocks);
  free(pool->mut);
  free(pool);
}

void
block_pool_reset(block_pool_t *pool)
{
  int32_t i;
  for(i=0;i<pool->m;i++) {
      block_destroy(pool->blocks[i]);
      pool->blocks[i] = NULL;
  }
  pool->n = 0;
}
