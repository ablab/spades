#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#include "bgzf.h"
#include "util.h"
#include "block.h"
#include "queue.h"
#include "pbgzf.h"
#include "reader.h"

static const int WINDOW_SIZE = BGZF_MAX_BLOCK_SIZE;

reader_t*
reader_init(int fd, queue_t *input, uint8_t compress, block_pool_t *pool)
{
  reader_t *r = calloc(1, sizeof(reader_t));

  if(0 == compress) {
      r->fp_bgzf = bgzf_fdopen(fd, "r");
  }
  else {
      r->fd_file = fd;
  }
  r->input = input;
  r->compress = compress;
  r->pool = pool;

  return r;
}

static int
reader_read_block(BGZF* fp, block_t *b)
{
  uint8_t header[BLOCK_HEADER_LENGTH];
  int count, size = 0, remaining;
  int64_t block_address = _bgzf_tell(fp->fp);
  count = _bgzf_read(fp->fp, header, sizeof(header));
  if (count == 0) {
      fp->block_length = b->block_length = 0;
      return 0;
  }
  size = count;
  if (count != sizeof(header)) {
      fprintf(stderr, "read failed\n");
      return -1;
  }
  if (!bgzf_check_header(header)) {
      fprintf(stderr, "invalid block header\n");
      return -1;
  }
  b->block_length = unpackInt16((uint8_t*)&header[16]) + 1;
  int8_t* compressed_block = (int8_t*) b->buffer;
  memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
  remaining = b->block_length - BLOCK_HEADER_LENGTH;
  count = _bgzf_read(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
  if (count != remaining) {
      fprintf(stderr, "read failed\n");
      return -1;
  }
  size += count;
  /*
  count = inflate_block(fp, block_length);
  if (count < 0) return -1;
  */
  if (fp->block_length != 0) {
      // Do not reset offset if this read follows a seek.
      fp->block_offset = 0;
  }
  fp->block_address = block_address;
  b->block_address = block_address;
  /*
   // TODO: in the consumer, with a mutex
  fp->block_length = count;
  cache_block(fp, size);
  */
  return 0;
}

void*
reader_run(void *arg)
{
  reader_t *r = (reader_t*)arg;
  block_t *b = NULL;
  int32_t wait;
  uint64_t n = 0;
  block_pool_t *pool;
  
  //fprintf(stderr, "reader staring\n");

  pool = block_pool_init2(PBGZF_BLOCKS_POOL_NUM);

  while(!r->is_done) {
#ifdef PBGZF_USE_LOCAL_POOLS
      // read block
      while(pool->n < pool->m) {
          if(NULL == r->pool || NULL == (b = block_pool_get(r->pool))) {
              b = block_init(); 
          }
          if(0 == r->compress) {
              if(reader_read_block(r->fp_bgzf, b) < 0) {
                  fprintf(stderr, "reader reader_read_block: bug encountered\n");
                  exit(1);
              }
          }
          else { 
              if((b->block_length = read(r->fd_file, b->buffer, WINDOW_SIZE)) < 0) {
                  fprintf(stderr, "reader read: bug encountered\n");
                  exit(1);
              }
          }
          if(NULL == b || 0 == b->block_length) {
              block_pool_add(r->pool, b);
              b = NULL;
              break;
          }
          if(0 == block_pool_add(pool, b)) {
              fprintf(stderr, "reader block_pool_add: bug encountered\n");
              exit(1);
          }
      }
      //fprintf(stderr, "reader: read in pool->n=%d\n", pool->n);

      if(0 == pool->n) {
          break;
      }

      // add to the queue
      while(0 < pool->n) {
          b = block_pool_peek(pool);
          if(NULL == b) {
              fprintf(stderr, "reader block_pool_get: bug encountered\n");
              exit(1);
          }
          wait = (pool->n == pool->m) ? 1 : 0; // NB: only wait if we cannot read in any more...
          if(0 == queue_add(r->input, b, wait)) {
              if(1 == wait) {
                  if(QUEUE_STATE_OK == r->input->state) {
                      fprintf(stderr, "reader queue_add: bug encountered\n");
                      exit(1);
                  }
                  else if(QUEUE_STATE_EOF == r->input->state) { // EOF, quit
                      break;
                  }
                  else {
                      // NB: if the reader has blocks, it does not make sense to
                      // flush
                      fprintf(stderr, "reader queue_add: bug encountered\n");
                      exit(1);
                  }
              }
              else {
                  break;
              }
          }
          block_pool_get(pool); // ignore return
          b = NULL;
          n++;
      }
      //fprintf(stderr, "reader: add to pool->n=%d\n", pool->n);
      //fprintf(stderr, "r->output->n=%d\n", r->input->n);
#else
      // read block
      //fprintf(stderr, "Reader #%d read block\n", 0);
      b = block_init(); 
      if(0 == r->compress) {
          if(reader_read_block(r->fp_bgzf, b) < 0) {
              fprintf(stderr, "reader reader_read_block: bug encountered\n");
              exit(1);
          }
      }
      else { 
          if((b->block_length = read(r->fd_file, b->buffer, WINDOW_SIZE)) < 0) {
              fprintf(stderr, "reader read: bug encountered\n");
              exit(1);
          }
      }
      if(NULL == b || 0 == b->block_length) {
          block_destroy(b);
          b = NULL;
          break;
      }

      // add to the queue
      //fprintf(stderr, "Reader #%d add to queue\n", 0);
      wait = 1;
      if(0 == queue_add(r->input, b, wait)) {
          if(1 == wait) {
              if(QUEUE_STATE_OK == r->input->state) {
                  fprintf(stderr, "reader queue_add: bug encountered\n");
                  exit(1);
              }
              else if(QUEUE_STATE_EOF == r->input->state) { // EOF, quit
                  block_destroy(b);
                  b = NULL;
                  break;
              }
              else {
                  queue_wait_until_not_flush(r->input);
                  continue;
              }
          }
          else {
              block_destroy(b);
              b = NULL;
              break;
          }
      }
      b = NULL;
      n++;
      //fprintf(stderr, "reader read %llu blocks\n", n);
#endif
  }
  block_destroy(b);
  b = NULL;
  
  r->is_done = 1;
  
  // NB: EOF should be handled when the adder is removed
  queue_remove_adder(r->input);
  
  //fprintf(stderr, "reader read %llu blocks\n", n);
  //fprintf(stderr, "reader r->input->n=%d\n", r->input->n);
  //fprintf(stderr, "reader r->input->state=%d QUEUE_STATE_EOF=%d\n", r->input->state, QUEUE_STATE_EOF);

  block_pool_destroy(pool);

  return arg;
}

void
reader_destroy(reader_t *r)
{
  if(NULL == r) return;
  if(0 == r->compress) {
      if(bgzf_close(r->fp_bgzf) < 0) {
          fprintf(stderr, "reader bgzf_close: bug encountered\n");
          exit(1);
      }
  }
  free(r);
}

void
reader_reset(reader_t *r)
{
    r->is_done = r->is_closed = 0;
}
