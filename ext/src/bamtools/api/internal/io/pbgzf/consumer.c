#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#ifndef DISABLE_BZ2
#include <bzlib.h>
#endif
#include <pthread.h>

#include "bgzf.h"
#include "util.h"
#include "block.h"
#include "queue.h"
#include "reader.h"
#include "pbgzf.h"
#include "consumer.h"

#define GZIP_WINDOW_BITS -15 // no zlib header
#define Z_DEFAULT_MEM_LEVEL 8

consumer_t*
consumer_init(queue_t *input,
              queue_t *output,
              reader_t *reader,
              int8_t compress,
              int32_t compress_level,
              int32_t compress_type,
              int32_t cid)
{
  consumer_t *c = calloc(1, sizeof(consumer_t));

  c->input = input;
  c->output = output;
  c->reader = reader;
  c->compress = compress;
#ifndef DISABLE_BZ2
  if (compress_type == 0) {
#endif
      c->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
#ifndef DISABLE_BZ2
  } else {
      c->compress_level = compress_level < 1? BZ2_DEFAULT_LEVEL : compress_level; 
  }
#endif
  c->compress_type = compress_type;
  c->cid = cid;

  c->buffer = malloc(sizeof(uint8_t)*BGZF_MAX_BLOCK_SIZE);

  return c;
}

static int
consumer_inflate_block_gz(consumer_t *c, block_t *block)
{
  z_stream zs;
  int status;
  
  // copy compressed buffer into consumer buffer
  memcpy(c->buffer, block->buffer, block->block_length);

  zs.zalloc = NULL;
  zs.zfree = NULL;
  zs.next_in = c->buffer + 18;
  zs.avail_in = block->block_length - 16;
  zs.next_out = (void*)block->buffer;
  zs.avail_out = BGZF_MAX_BLOCK_SIZE;

  status = inflateInit2(&zs, GZIP_WINDOW_BITS);
  if (status != Z_OK) {
      fprintf(stderr, "inflate init failed\n");
      return -1;
  }

  status = inflate(&zs, Z_FINISH);
  if (status != Z_STREAM_END) {
      inflateEnd(&zs);
      fprintf(stderr, "inflate failed\n");
      return -1;
  }

  status = inflateEnd(&zs);
  if (status != Z_OK) {
      fprintf(stderr, "inflate end failed\n");
      return -1;
  }

  return zs.total_out;
}

#ifndef DISABLE_BZ2
static int
consumer_inflate_block_bz2(consumer_t *c, block_t *b)
{
  int32_t block_length;
  int status;
  unsigned int destLen;

  memcpy(c->buffer, b->buffer, b->block_length); // copy uncompressed to compressed 
  block_length = b->block_length;
  
  destLen = BGZF_MAX_BLOCK_SIZE;
  status = BZ2_bzBuffToBuffDecompress((char*)b->buffer, &destLen, (char*)(c->buffer + 18), block_length-16, 0, 0); 
  if (BZ_OK != status) {
      // TODO: deal with BZ_OUTBUFF_FULL 
      fprintf(stderr, "inflate failed\n");
      return -1;
  }
  b->block_length = destLen;

  return b->block_length;
}
#endif

static int
consumer_inflate_block(consumer_t *c, block_t *block)
{
#ifndef DISABLE_BZ2
  // from bgzf.c
  static uint8_t magic[28] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
  static int32_t magic_l = 28;
  if (0 == c->compress_type) return consumer_inflate_block_gz(c, block);
  else if (magic_l == block->block_length) { // check EOF magic #...
      if (0 != memcmp(magic, block->buffer, magic_l)) return consumer_inflate_block_bz2(c, block);
      else return consumer_inflate_block_gz(c, block);
  } else {
      return consumer_inflate_block_bz2(c, block);
  }
#else
  return consumer_inflate_block_gz(c, block);
#endif
}

static int
consumer_deflate_block_gz(consumer_t *c, block_t *b)
{
  // Deflate the block in fp->uncompressed_block into fp->compressed_block.
  // Also adds an extra field that stores the compressed block length.
  int32_t block_length;

  memcpy(c->buffer, b->buffer, b->block_length); // copy uncompressed to compressed 

  // Notes:
  // fp->compressed_block is now b->buffer
  // fp->uncompressed_block is now c->buffer
  // block_length is now b->block_length
  
  block_length = b->block_length;
  if(0 != bgzf_compress(b->buffer, &b->block_length, c->buffer, block_length, c->compress_level)) {
      fprintf(stderr, "deflate failed\n");
      return -1;
  }
  b->block_offset = 0;

  return b->block_length;
}

#ifndef DISABLE_BZ2
static int
consumer_deflate_block_bz2(consumer_t *c, block_t *b)
{
  // Deflate the block in fp->uncompressed_block into fp->compressed_block.
  // Also adds an extra field that stores the compressed block length.
  int32_t block_length, compressed_length = 0;
  int status, input_length;
  unsigned int destLen;

  memcpy(c->buffer, b->buffer, b->block_length); // copy uncompressed to compressed 
  
  // NB: use gzip header, for now...
  consumer_init_header(b->buffer);

  input_length = block_length = b->block_length;
  destLen = BGZF_MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
  while (1) {
      status = BZ2_bzBuffToBuffCompress((char*)(b->buffer + BLOCK_HEADER_LENGTH), &destLen, (char*)c->buffer, input_length, 9, 0, 0); 
      if (BZ_OK != status) {
          if (status == BZ_OUTBUFF_FULL) {
              // Not enough space in buffer.
              // Can happen in the rare case the input doesn't compress enough.
              // Reduce the amount of input until it fits.
              input_length -= 1024;
              if (input_length <= 0) {
                  // should never happen
                  fprintf(stderr, "input reduction failed\n");
                  return -1;
              }
              continue;
          }
          fprintf(stderr, "deflate failed\n");
          return -1;
      }
      break;
  }
  compressed_length = destLen;
  compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
  if (compressed_length > BGZF_MAX_BLOCK_SIZE) {
      // should never happen
      fprintf(stderr, "deflate overflow\n");
      return -1;
  }
  b->block_length = compressed_length;
  
  consumer_update_header(c, b->buffer, input_length, compressed_length);

  int remaining = block_length - input_length;
  b->block_offset = remaining;
  b->block_offset = 0;
  
  return compressed_length;
}
#endif

static int
consumer_deflate_block(consumer_t *c, block_t *b)
{
#ifndef DISABLE_BZ2
  if(0 == c->compress_type) return consumer_deflate_block_gz(c, b);
  else return consumer_deflate_block_bz2(c, b);
#else
  return consumer_deflate_block_gz(c, b);
#endif
}

void*
consumer_run(void *arg)
{
  consumer_t *c = (consumer_t*)arg;
  block_t *b = NULL;
  int32_t wait;
  block_pool_t *pool_in = NULL, *pool_out = NULL;

  pool_in = block_pool_init2(PBGZF_BLOCKS_POOL_NUM);
  pool_out = block_pool_init2(PBGZF_BLOCKS_POOL_NUM);
  c->n = 0;

  //fprintf(stderr, "consumer #%d starting\n", c->cid);
  //fprintf(stderr, "consumer start pool_in->n=%d\n", pool_in->n);
  //fprintf(stderr, "consumer start pool_out->n=%d\n", pool_out->n);

  while(1) {
#ifdef PBGZF_USE_LOCAL_POOLS
      // get block(s)
      while(pool_in->n < pool_in->m) { // more to read in
          b = queue_get(c->input, (0 == pool_in->n && 0 == pool_out->n) ? 1 : 0); // NB: only wait if the pools are empty
          if(NULL == b) {
              break;
          }
          //fprintf(stderr, "ADDING to pool_in b->id=%d\n", b->id);
          if(0 == block_pool_add(pool_in, b)) {
              fprintf(stderr, "consumer block_pool_add: bug encountered\n");
              exit(1);
          }
          b = NULL;
      }
      //fprintf(stderr, "consumer get blocks pool_in->n=%d\n", pool_in->n);
      //fprintf(stderr, "consumer get blocks pool_out->n=%d\n", pool_out->n);
      if(0 == pool_in->n && 0 == pool_out->n) { // no more data
          if(QUEUE_STATE_FLUSH == c->input->state) {
              queue_wait_until_not_flush(c->input);
              continue;
          }
          else if((NULL == c->reader && QUEUE_STATE_EOF == c->input->state) 
             || (NULL != c->reader && 1 == c->reader->is_done)) { // TODO: does this need to be synced?
              break;
          }
          else {
              fprintf(stderr, "consumer queue_get: bug encountered\n");
              exit(1);
          }
      }

      // inflate/deflate
      while(0 < pool_in->n && pool_out->n < pool_in->m) { // consume while the in has more and the out has room
          b = block_pool_peek(pool_in);
          //fprintf(stderr, "PEEK from pool_in b->id=%d\n", b->id);
          if(NULL == b) {
              fprintf(stderr, "consumer block_pool_get: bug encountered\n");
              exit(1);
          }
          if(0 == c->compress) {
              if((b->block_length = consumer_inflate_block(c, b)) < 0) {
                  fprintf(stderr, "Error decompressing\n");
                  exit(1);
              }
          }
          else if(1 == c->compress) {
              if((b->block_length = consumer_deflate_block(c, b)) < 0) {
                  fprintf(stderr, "Error decompressing\n");
                  exit(1);
              }
          }
          if(0 == block_pool_add(pool_out, b)) {
              fprintf(stderr, "consumer block_pool_add: bug encountered\n");
              exit(1);
          }
          block_pool_get(pool_in); // ignore return
          b = NULL;
      }
      //fprintf(stderr, "consumer inflate/deflate pool_in->n=%d\n", pool_in->n);
      //fprintf(stderr, "consumer inflate/deflate pool_out->n=%d\n", pool_out->n);

      // put back a block
      while(0 < pool_out->n) {
          b = block_pool_peek(pool_out);
          // NB: only wait if the pools are full
          wait = (pool_in->m == pool_in->n && pool_out->m == pool_out->n) ? 1 : 0;
          if(!queue_add(c->output, b, wait)) {
              if(1 == wait && QUEUE_STATE_EOF != c->output->state) {
                  fprintf(stderr, "consumer queue_add: bug encountered\n");
                  exit(1);
              }
              else {
                  break;
              }
          }
          block_pool_get(pool_out); // ignore return
          b = NULL;
          c->n++;
      }
      //fprintf(stderr, "consumer output pool_in->n=%d\n", pool_in->n);
      //fprintf(stderr, "consumer output pool_out->n=%d\n", pool_out->n);
      //fprintf(stderr, "consumer output c->output->n=%d\n", c->output->n);
#else
      // get block
      //fprintf(stderr, "Consumer #%d get block\n", c->cid);
      b = queue_get(c->input, 1);
      if(NULL == b) {
          if(QUEUE_STATE_FLUSH == c->input->state) {
              queue_wait_until_not_flush(c->input);
              continue;
          }
          else if((NULL == c->reader && QUEUE_STATE_EOF == c->input->state) 
             || (NULL != c->reader && 1 == c->reader->is_done)) { // TODO: does this need to be synced?
              break;
          }
          else {
              fprintf(stderr, "consumer queue_get: bug encountered\n");
              exit(1);
          }
      }

      // inflate/deflate
      //fprintf(stderr, "Consumer #%d inflate/deflate\n", c->cid);
      if(0 == c->compress) {
          if((b->block_length = consumer_inflate_block(c, b)) < 0) {
              fprintf(stderr, "Error decompressing\n");
              exit(1);
          }
      }
      else if(1 == c->compress) {
          if((b->block_length = consumer_deflate_block(c, b)) < 0) {
              fprintf(stderr, "Error decompressing\n");
              exit(1);
          }
      }

      // put back a block
      //fprintf(stderr, "Consumer #%d add block\n", c->cid);
      wait = 1;
      if(!queue_add(c->output, b, wait)) {
          if(1 == wait && QUEUE_STATE_EOF != c->output->state) {
              fprintf(stderr, "consumer queue_add: bug encountered\n");
              exit(1);
          }
          else {
              break;
          }
      }
      b = NULL;
      c->n++;
#endif
      
      /*
      fprintf(stderr, "consumer #%d c->input=[%d/%d,%d] c->output=[%d/%d,%d]\n",
              c->cid,
              c->input->n, c->input->mem, c->input->state,
              c->output->n, c->output->mem, c->output->state);
      */
  }

  c->is_done = 1;
  // NB: queue handles waking...
  queue_remove_getter(c->input);
  queue_remove_adder(c->output);

  //fprintf(stderr, "consumer #%d done processed %llu blocks\n", c->cid, c->n);
  //queue_print_status(c->input, stderr);
  //queue_print_status(c->output, stderr);

  // TODO: why do you need to wake all here?
  queue_wake_all(c->input);
  queue_wake_all(c->output);

  // destroy the pool
  block_pool_destroy(pool_in);
  block_pool_destroy(pool_out);

  //fprintf(stderr, "consumer #%d exiting\n", c->cid);
  return arg;
}

void
consumer_destroy(consumer_t *c)
{
  if(NULL == c) return;
  free(c->buffer);
  free(c);
}

void
consumer_reset(consumer_t *c)
{
  c->is_done = 0;
}
