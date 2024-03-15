#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>

#include "util.h"
#include "block.h"
#include "queue.h"

static void
queue_close_nolock(queue_t *q)
{
  if(QUEUE_STATE_EOF == q->state) return;
  q->state = QUEUE_STATE_EOF;
  pthread_cond_broadcast(q->not_full);
  pthread_cond_broadcast(q->not_empty);
}

void
queue_signal(queue_t *q)
{
  if(q->n < q->mem) pthread_cond_signal(q->not_full);
  if(0 == q->n) pthread_cond_signal(q->is_empty);
  if(0 < q->n) pthread_cond_signal(q->not_empty);
}

queue_t*
queue_init(int32_t capacity, int8_t ordered, int32_t num_adders, int32_t num_getters)
{
  queue_t *q = calloc(1, sizeof(queue_t));

  q->mem = capacity;
  q->queue = calloc(q->mem, sizeof(block_t*));
  q->ordered = ordered;

  q->mut = calloc(1, sizeof(pthread_mutex_t));
  q->not_full = calloc(1, sizeof(pthread_cond_t));
  q->not_empty = calloc(1, sizeof(pthread_cond_t));
  q->is_empty = calloc(1, sizeof(pthread_cond_t));
  q->not_flush = calloc(1, sizeof(pthread_cond_t));
  q->state = QUEUE_STATE_OK;
  q->num_adders = num_adders;
  q->num_getters = num_getters;
#ifdef QUEUE_DEBUG
  q->num_waiting[0] = 0;
  q->num_waiting[1] = 0;
  q->num_waiting[2] = 0;
  q->num_waiting[3] = 0;
#endif
    
  if(0 != pthread_mutex_init(q->mut, NULL)) {
      fprintf(stderr, "Could not create mutex\n");
      exit(1);
  }
  if(0 != pthread_cond_init(q->not_full, NULL)) {
      fprintf(stderr, "Could not create condition\n");
      exit(1);
  }
  if(0 != pthread_cond_init(q->not_empty, NULL)) {
      fprintf(stderr, "Could not create condition\n");
      exit(1);
  }
  if(0 != pthread_cond_init(q->is_empty, NULL)) {
      fprintf(stderr, "Could not create condition\n");
      exit(1);
  }
  if(0 != pthread_cond_init(q->not_flush, NULL)) {
      fprintf(stderr, "Could not create condition\n");
      exit(1);
  }

  return q;
}

int8_t
queue_add(queue_t *q, block_t *b, int8_t wait)
{
  safe_mutex_lock(q->mut);
  queue_signal(q);
  if(0 == q->num_getters) { // no more getters
      queue_close_nolock(q); // close
      safe_mutex_unlock(q->mut);
      return 0;
  }
  else if(0 == q->num_adders) { // then why are you adding?
      queue_close_nolock(q);
      safe_mutex_unlock(q->mut);
      return 0;
  }
  while(q->n == q->mem) {
      if(wait && QUEUE_STATE_OK == q->state) {
#ifdef QUEUE_DEBUG
          q->num_waiting[0]++;
#endif
          if(0 != pthread_cond_wait(q->not_full, q->mut)) {
              fprintf(stderr, "Could not condition wait\n");
              exit(1);
          }
#ifdef QUEUE_DEBUG
          q->num_waiting[0]--;
#endif
      }
      else {
          if(0 == q->num_getters) queue_close_nolock(q);
          queue_signal(q);
          safe_mutex_unlock(q->mut);
          return 0;
      }
  }
  if(1 == q->ordered) {
      while(q->id - q->n + q->mem <= b->id) {
          if(wait && QUEUE_STATE_OK == q->state) {
#ifdef QUEUE_DEBUG
              q->num_waiting[1]++;
#endif
              queue_signal(q); // NB: important to signal in case a getter may open this slot 
              if(0 != pthread_cond_wait(q->not_full, q->mut)) {
                  fprintf(stderr, "Could not condition wait\n");
                  exit(1);
              }
#ifdef QUEUE_DEBUG
              q->num_waiting[1]--;
#endif
          }
          else {
              if(0 == q->num_getters) queue_close_nolock(q);
              queue_signal(q);
              safe_mutex_unlock(q->mut);
              return 0;
          }
      }
      if(NULL != q->queue[b->id % q->mem]) {
          fprintf(stderr, "Overwritting an existing block\n");
          exit(1);
      }
      q->queue[b->id % q->mem] = b;
  }
  else {
      b->id = q->id;
      q->id++;
      q->queue[q->tail++] = b;
      if(q->tail == q->mem) q->tail = 0;
  }
  q->n++;
  queue_signal(q);
  safe_mutex_unlock(q->mut);
  return 1;
}

block_t*
queue_get(queue_t *q, int8_t wait)
{
  block_t *b = NULL;
  safe_mutex_lock(q->mut);
  queue_signal(q);
  if(0 == q->num_getters) { // then why are you getting
      queue_close_nolock(q); // close the queue
      safe_mutex_unlock(q->mut);
      return NULL;
  }
  else if(0 == q->n && 0 == q->num_adders) { 
      queue_close_nolock(q);
      safe_mutex_unlock(q->mut);
      return NULL;
  }
  while(0 == q->n) {
      if(1 == wait && QUEUE_STATE_OK == q->state) {
#ifdef QUEUE_DEBUG
          q->num_waiting[2]++;
#endif
          if(0 != pthread_cond_wait(q->not_empty, q->mut)) {
              fprintf(stderr, "Could not condition wait\n");
              exit(1);
          }
#ifdef QUEUE_DEBUG
          q->num_waiting[2]--;
#endif
      }
      else {
          if(0 == q->num_adders && 0 == q->n) queue_close_nolock(q); // close the queue
          queue_signal(q);
          safe_mutex_unlock(q->mut);
          return NULL;
      }
  }
  b = q->queue[q->head];
  if(q->ordered) {
      while(NULL == b) {
          if(1 == wait && QUEUE_STATE_OK == q->state) {
#ifdef QUEUE_DEBUG
              q->num_waiting[3]++;
#endif
              queue_signal(q); // NB: important to signal as an adder may add to this slot
              if(0 != pthread_cond_wait(q->not_empty, q->mut)) {
                  fprintf(stderr, "Could not condition wait\n");
                  exit(1);
              }
#ifdef QUEUE_DEBUG
              q->num_waiting[3]--;
#endif
          }
          else {
              if(0 == q->num_adders && 0 == q->n) queue_close_nolock(q); // close the queue
              queue_signal(q);
              safe_mutex_unlock(q->mut);
              return NULL;
          }
          b = q->queue[q->head];
      }
  }
  q->queue[q->head++] = NULL;
  if(q->head == q->mem) q->head = 0;
  if(q->ordered) q->id++;
  q->n--;
  queue_signal(q);
  safe_mutex_unlock(q->mut);
  return b;
}

void
queue_wait_until_empty(queue_t *q)
{
  safe_mutex_lock(q->mut);
  if(0 < q->n) { // wait
      if(0 != pthread_cond_wait(q->is_empty, q->mut)) {
          fprintf(stderr, "Could not condition wait\n");
          exit(1);
      }
  }
  safe_mutex_unlock(q->mut);
}

void
queue_wait_until_not_flush(queue_t *q)
{
  safe_mutex_lock(q->mut);
  if(QUEUE_STATE_FLUSH == q->state) { // wait
      if(0 != pthread_cond_wait(q->not_flush, q->mut)) {
          fprintf(stderr, "Could not condition wait\n");
          exit(1);
      }
  }
  safe_mutex_unlock(q->mut);
}

void
queue_remove_flush(queue_t *q)
{
  if(QUEUE_STATE_FLUSH != q->state) return;
  safe_mutex_lock(q->mut);
  q->state = QUEUE_STATE_OK;
  pthread_cond_signal(q->not_flush);
  safe_mutex_unlock(q->mut);
}

void
queue_close(queue_t *q)
{
  if(QUEUE_STATE_EOF == q->state) return;
  safe_mutex_lock(q->mut);
  queue_close_nolock(q);
  safe_mutex_unlock(q->mut);
}

void
queue_reset(queue_t *q, int32_t num_adders, int32_t num_getters)
{
  int32_t i;
  safe_mutex_lock(q->mut);
  for(i=0;i<q->mem;i++) {
      if(NULL != q->queue[i]) {
          block_destroy(q->queue[i]);
          q->queue[i] = NULL;
      }
  }
  q->head = q->tail = q->n = 0;
  q->id = 0;
  q->state = QUEUE_STATE_OK;
  q->num_adders = num_adders;
  q->num_getters = num_getters;
  safe_mutex_unlock(q->mut);
}

void
queue_destroy(queue_t *q)
{
  int32_t i;
  if(NULL == q) return;
  queue_close(q);
  for(i=0;i<q->mem;i++) {
      block_destroy(q->queue[i]);
  }
  free(q->queue);
  free(q->mut);
  free(q->not_full);
  free(q->not_empty);
  free(q->is_empty);
  free(q->not_flush);
  free(q);
}

static void
queue_wake_all_no_lock(queue_t *q)
{
  if(0 == q->num_getters || (0 == q->num_adders && 0 == q->n)) {
      q->state = QUEUE_STATE_EOF;
  }
  pthread_cond_signal(q->not_full);
  pthread_cond_signal(q->not_empty);
  pthread_cond_signal(q->is_empty);
  pthread_cond_signal(q->not_flush);
}

void
queue_wake_all(queue_t *q)
{
  safe_mutex_lock(q->mut);
  queue_wake_all_no_lock(q);
  safe_mutex_unlock(q->mut);
}

void 
queue_remove_adder(queue_t *q)
{
  safe_mutex_lock(q->mut);
  q->num_adders--;
  if(0 == q->num_adders && 0 == q->n) queue_wake_all_no_lock(q);
  safe_mutex_unlock(q->mut);
}

void 
queue_remove_getter(queue_t *q)
{
  safe_mutex_lock(q->mut);
  q->num_getters--;
  if(0 == q->num_getters) queue_wake_all_no_lock(q);
  safe_mutex_unlock(q->mut);
}

void
queue_print_status(queue_t *q, FILE *fp)
{
    fprintf(fp, "QUEUE STATUS\n");
    fprintf(fp, "mem=%d head=%d tail=%d n=%d length=%d id=%lld ordered=%d num_adders=%d num_getters=%d\n",
            q->mem, q->head, q->tail, q->n, q->length, q->id, q->ordered, q->num_adders, q->num_getters);
#ifdef QUEUE_DEBUG
    fprintf(fp, "num_waiting=[%d,%d,%d,%d]\n", 
            q->num_waiting[0], q->num_waiting[1],
            q->num_waiting[2], q->num_waiting[3]);
#endif
}
