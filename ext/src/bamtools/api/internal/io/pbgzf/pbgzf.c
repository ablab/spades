#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "bgzf.h"
#include "util.h"
#include "block.h"
#include "queue.h"
#include "reader.h"
#include "writer.h"
#include "consumer.h"
#include "pbgzf.h"

static int32_t num_threads_per_pbgzf = -1;

void
pbgzf_set_num_threads_per(int32_t n)
{
  fprintf(stderr, "Setting the number of threads per PBGZF file handle to %d\n", n);
  num_threads_per_pbgzf = n;
}

static consumers_t*
consumers_init(int32_t n, queue_t *input, queue_t *output, reader_t *reader, 
               int32_t compress, int32_t compress_level, int32_t compress_type)
{
  consumers_t *c = NULL;
  int32_t i;

  c = calloc(1, sizeof(consumers_t));
  c->threads = calloc(n, sizeof(pthread_t));
  c->n = n;
  c->c = calloc(n, sizeof(consumer_t*));

  for(i=0;i<n;i++) {
      c->c[i] = consumer_init(input, output, reader, compress, compress_level, compress_type, i);
  }

  pthread_attr_init(&c->attr);
  pthread_attr_setdetachstate(&c->attr, PTHREAD_CREATE_JOINABLE);

  return c;
}

static void
consumers_destroy(consumers_t *c)
{
  int32_t i;
  free(c->threads);
  for(i=0;i<c->n;i++) {
      consumer_destroy(c->c[i]);
  }
  free(c->c);
  free(c);
}

static void
consumers_run(consumers_t *c)
{
  int32_t i;

  for(i=0;i<c->n;i++) {
      if(0 != pthread_create(&c->threads[i], &c->attr, consumer_run, c->c[i])) {
          fprintf(stderr, "failed to create threads");
          exit(1);
      }
  }
}

static void
consumers_join(consumers_t *c)
{
  int32_t i;

  for(i=0;i<c->n;i++) {
      if(0 != pthread_join(c->threads[i], NULL)) {
          fprintf(stderr, "failed to join threads");
          exit(1);
      }
  }
  // flush the output queue
  if(QUEUE_STATE_OK == c->c[0]->output->state) c->c[0]->output->state = QUEUE_STATE_FLUSH;
}

static void
consumers_reset(consumers_t *c)
{
  int32_t i;
  for(i=0;i<c->n;i++) {
      consumer_reset(c->c[i]);
  }
}

static producer_t*
producer_init(reader_t *r)
{
  producer_t *p;

  p = calloc(1, sizeof(producer_t));
  p->r = r;

  pthread_attr_init(&p->attr);
  pthread_attr_setdetachstate(&p->attr, PTHREAD_CREATE_JOINABLE);

  return p;
}

static void
producer_destroy(producer_t *p)
{
  free(p);
}

static void
producer_run(producer_t *p)
{
  if(0 != pthread_create(&p->thread, &p->attr, reader_run, p->r)) {
      fprintf(stderr, "failed to create threads");
      exit(1);
  }
}

static void
producer_join(producer_t *p)
{
  if(0 != pthread_join(p->thread, NULL)) {
      fprintf(stderr, "failed to join threads");
      exit(1);
  }
  // close the input queue
  // TODO: sync?
  if(p->r->input->state == QUEUE_STATE_OK) p->r->input->state = QUEUE_STATE_FLUSH;
}

static void
producer_reset(producer_t *p)
{
  reader_reset(p->r);
}

static outputter_t*
outputter_init(writer_t *w)
{
  outputter_t *o;

  o = calloc(1, sizeof(outputter_t));
  o->w = w;

  pthread_attr_init(&o->attr);
  pthread_attr_setdetachstate(&o->attr, PTHREAD_CREATE_JOINABLE);

  return o;
}

static void
outputter_destroy(outputter_t *o)
{
  free(o);
}

static void
outputter_run(outputter_t *o)
{
  if(0 != pthread_create(&o->thread, &o->attr, writer_run, o->w)) {
      fprintf(stderr, "failed to create threads");
      exit(1);
  }
}

static void
outputter_join(outputter_t *o)
{
  if(0 != pthread_join(o->thread, NULL)) {
      fprintf(stderr, "failed to join threads");
      exit(1);
  }
}

static void
outputter_reset(outputter_t *o)
{
  writer_reset(o->w);
}

static inline
int
pbgzf_min(int x, int y)
{
  return (x < y) ? x : y;
}

static void
pbgzf_run(PBGZF *fp)
{
  if(NULL != fp->p) producer_run(fp->p);
  if(NULL != fp->c) consumers_run(fp->c);
  if(NULL != fp->o) outputter_run(fp->o);
}

static void
pbgzf_join(PBGZF *fp)
{
  // join producer
  if(NULL != fp->p) producer_join(fp->p);
  else if(QUEUE_STATE_OK == fp->input->state) fp->input->state = QUEUE_STATE_FLUSH;

  // wake, just in case
  queue_wake_all(fp->input);
  if('r' == fp->open_mode) {
      queue_wake_all(fp->output);
  }

  // join the consumers
  if(NULL != fp->c) consumers_join(fp->c);
  else if(QUEUE_STATE_OK == fp->output->state) fp->output->state = QUEUE_STATE_FLUSH;

  // wake just in case
  queue_wake_all(fp->output);

  // join the outputter
  if(NULL != fp->o) outputter_join(fp->o);
}

static PBGZF*
pbgzf_init(int fd, const char* __restrict mode)
{
  int i, compress_level = -1, compress_type = 0, is_write;
  char open_mode;
  PBGZF *fp = NULL;

  // set compress_level
  for (i = 0; mode[i]; ++i)
    if (mode[i] >= '0' && mode[i] <= '9') break;
  if (mode[i]) compress_level = (int)mode[i] - '0';
  if (strchr(mode, 'u')) compress_level = 0;

  // set read/write
  if (strchr(mode, 'r') || strchr(mode, 'R')) { /* The reading mode is preferred. */
      open_mode = 'r';
      is_write = 0;
  } else if (strchr(mode, 'w') || strchr(mode, 'W')) {
      open_mode = 'w';
      is_write = 1;
  }
  else {
      return NULL;
  }
  
  // set type
#ifndef DISABLE_BZ2
  if (strchr(mode, 'Z')) {
      compress_type = 0;
  } else if (strchr(mode, 'B')) {
      compress_type = 1;
      if (compress_level < 1) {
          compress_type = 0; // switch to gz no compress mode
          compress_level = 0;
      }
  }
  else {
      compress_type = 0;
  }
#endif

  fp = calloc(1, sizeof(PBGZF));

  // queues
  fp->open_mode = open_mode;
  fp->is_write = is_write;
  if(num_threads_per_pbgzf <= 0) {
      fp->num_threads = detect_cpus(); 
  }
  else {
      fp->num_threads = num_threads_per_pbgzf;
  }
  fprintf(stderr, "%s with %d threads.\n", ('r' == open_mode) ? "Reading" : "Writing", fp->num_threads);
  fp->queue_size = PBGZF_QUEUE_SIZE;
  fp->input = queue_init(fp->queue_size, 0, 1, fp->num_threads);
  fp->output = queue_init(fp->queue_size, 1, fp->num_threads, 1);

  fp->pool = block_pool_init(PBGZF_BLOCKS_POOL_NUM);
  fp->block = NULL;
  fp->n_blocks = 0;

  if('w' == open_mode) { // write to a compressed file
      fp->r = NULL; // do not read
      fp->p = NULL; // do not produce data
      fp->c = consumers_init(fp->num_threads, fp->input, fp->output, fp->r, 1, compress_level, compress_type); // deflate/compress
      fp->w = writer_init(fd, fp->output, 1, compress_level, compress_type, fp->pool); // write data
      fp->o = outputter_init(fp->w);
  }
  else { // read from a compressed file
      if(strchr(mode, 'u')) {// hidden functionality
          fp->r = reader_init(fd, fp->input, 1, fp->pool); // read the uncompressed file
          fp->p = producer_init(fp->r);
          fp->c = consumers_init(fp->num_threads, fp->input, fp->output, fp->r, 2, compress_level, compress_type); // do nothing
      }
      else {
          fp->r = reader_init(fd, fp->input, 0, fp->pool); // read the compressed file
          fp->p = producer_init(fp->r);
          fp->c = consumers_init(fp->num_threads, fp->input, fp->output, fp->r, 0, compress_level, compress_type); // inflate
          fp->eof_ok = bgzf_check_EOF(fp->r->fp_bgzf);
      }
      fp->w = NULL;
      fp->o = NULL; // do not write
  }

  pbgzf_run(fp);

  return fp;
}

PBGZF* pbgzf_fdopen(int fd, const char* __restrict mode)
{
  return pbgzf_init(fd, mode);
}

PBGZF* pbgzf_open(const char* path, const char* __restrict mode)
{
  int fd;
  if (strchr(mode, 'r') || strchr(mode, 'R')) { /* The reading mode is preferred. */
#ifdef _WIN32
      fd = open(path, O_RDONLY | O_BINARY); 
#else 
      fd = open(path, O_RDONLY);
#endif
  }
  else { // write to a compressed file
#ifdef _WIN32
      fd = open(path, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0666); 
#else
      fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0666); 
#endif
  }
  if(-1 == fd) return NULL;
  return pbgzf_fdopen(fd, mode);
}

int 
pbgzf_read(PBGZF* fp, void* data, int length)
{
  if(length <= 0) {
      return 0;
  }
  if(fp->open_mode != 'r') {
      fprintf(stderr, "file not open for reading\n");
      return -1;
  }
  if(fp->eof == 1) return 0;

  int bytes_read = 0, available = 0;
  int8_t* output = data;
  while(bytes_read < length) {
      int copy_length;

      available = (NULL == fp->block) ? 0 : (fp->block->block_length - fp->block->block_offset);
      if(0 == available) {
          if(NULL != fp->block) {
              block_destroy(fp->block);
              fp->n_blocks++;
          }
          fp->block = queue_get(fp->output, 1);
      }
      available = (NULL == fp->block) ? 0 : (fp->block->block_length - fp->block->block_offset);
      if(available <= 0) {
          break;
      }

      int8_t *buffer;
      copy_length = pbgzf_min(length-bytes_read, available);
      buffer = (int8_t*)fp->block->buffer;
      memcpy(output, buffer + fp->block->block_offset, copy_length);
      fp->block->block_offset += copy_length;
      output += copy_length;
      bytes_read += copy_length;
  }
  // Try to get a new block and reset address/offset
  if(NULL == fp->block || fp->block->block_offset == fp->block->block_length) {
      if(NULL != fp->block) {
          block_destroy(fp->block);
          fp->block = NULL;
          fp->n_blocks++;
      }
      if(0 < available) {
          //fp->block = queue_get(fp->output, 1-fp->r->is_done);
          // NB: do not wait if there will be no more data
          fp->block = queue_get(fp->output, (QUEUE_STATE_EOF == fp->output->state) ? 0 : 1);
      } // TODO: otherwise EOF?
      if(NULL == fp->block) {
          fp->block_offset = 0;
          fp->block_address = bgzf_tell(fp->r->fp_bgzf);
      }
      else {
          fp->block_offset = fp->block->block_offset;
          fp->block_address = fp->block->block_address;
      }
  }
  else {
      fp->block_offset = fp->block->block_offset;
      fp->block_address = fp->block->block_address;
  }

  if(0 == bytes_read) fp->eof = 1;

  return bytes_read;
}

int 
pbgzf_write(PBGZF* fp, const void* data, int length)
{
  const int8_t *input = data;
  int block_length, bytes_written;

  if(fp->open_mode != 'w') {
      fprintf(stderr, "file not open for writing\n");
      return -1;
  }

  if(NULL == fp->block) {
      fp->block = block_init();
      fp->block->block_length = BGZF_MAX_BLOCK_SIZE;
  }

  input = data;
  block_length = fp->block->block_length;
  bytes_written = 0;
  while (bytes_written < length) {
      int copy_length = pbgzf_min(block_length - fp->block->block_offset, length - bytes_written);
      int8_t* buffer = fp->block->buffer;
      memcpy(buffer + fp->block->block_offset, input, copy_length);
      fp->block->block_offset += copy_length;
      fp->block_offset += copy_length;
      input += copy_length;
      bytes_written += copy_length;
      /*
         fprintf(stderr, "fp->block_offset=%d copy_length=%d bytes_written=%d\n", 
         fp->block_offset, copy_length, bytes_written);
         */
      if (fp->block->block_offset == block_length) {
          // add to the queue
          if(!queue_add(fp->input, fp->block, 1)) {
              fprintf(stderr, "pbgzf_write queue_add: bug encountered\n");
              exit(1);
          }
          fp->block = NULL;
          fp->block = block_init();
          fp->block_offset = 0;
          fp->block->block_length = BGZF_MAX_BLOCK_SIZE;
          fp->n_blocks++;
      }
  }
  return bytes_written;
}

int64_t 
pbgzf_tell(PBGZF *fp)
{
  if('w' == fp->open_mode) {
      // wait until the input queue and output queue are empty
      queue_wait_until_empty(fp->input);
      queue_wait_until_empty(fp->output);
      return bgzf_tell(fp->w->fp_bgzf);
  }
  else { // reading
      return ((fp->block_address << 16) | (fp->block_offset & 0xFFFF));
  }
}

int64_t 
pbgzf_seek(PBGZF* fp, int64_t pos, int where)
{
  int32_t i;

  if(fp->open_mode != 'r') {
      fprintf(stderr, "file not open for read\n");
      return -1;
  }
  if (where != SEEK_SET) {
      fprintf(stderr, "unimplemented seek option\n");
      return -1;
  }

  // signal and join
  // signal other threads to finish
  pthread_cond_signal(fp->input->not_full);
  pthread_cond_signal(fp->input->not_empty);
  if(NULL != fp->w) fp->w->is_done = 1;
  if(NULL != fp->c) {
      for(i=0;i<fp->c->n;i++) {
          fp->c->c[i]->is_done = 1;
      }
  }
  if(NULL != fp->r) fp->r->is_done = 1;

  // join
  pbgzf_join(fp);

  // seek
  if(bgzf_seek(fp->r->fp_bgzf, pos, where) < 0) {
      return -1;
  };

  // reset the producer/consumer/outputter
  if(NULL != fp->p) producer_reset(fp->p);
  if(NULL != fp->c) consumers_reset(fp->c);
  if(NULL != fp->o) outputter_reset(fp->o);

  // reset the queues
  queue_reset(fp->input, 1, fp->num_threads);
  queue_reset(fp->output, fp->num_threads, 1);
  
  // restart threads
  pbgzf_run(fp);

  // get a block
  if(NULL != fp->block) block_destroy(fp->block);
  fp->block = queue_get(fp->output, 1);
  if(NULL == fp->block) fp->eof = 1; // must be EOF
  else fp->eof = 0;

  // reset block offset/address
  fp->block_offset = pos & 0xFFFF;
  fp->block_address = (pos >> 16) & 0xFFFFFFFFFFFFLL;
  if(NULL != fp->block) fp->block->block_offset = fp->block_offset;

  return 0;
}

int 
pbgzf_check_EOF(PBGZF *fp)
{
  if('r' != fp->open_mode) {
      fprintf(stderr, "file not open for reading\n");
      exit(1);
  }
  return fp->eof_ok;
}

static int 
pbgzf_flush_aux(PBGZF* fp, int32_t restart)
{
  int ret;

  if('w' != fp->open_mode) {
      fprintf(stderr, "file not open for writing\n");
      exit(1);
  }

  // flush
  if(0 < fp->block->block_offset) {
      fp->block->block_length = fp->block->block_offset;
      if(!queue_add(fp->input, fp->block, 1)) {
          fprintf(stderr, "pbgzf_flush_aux queue_add: bug encountered\n");
          exit(1);
      }
      fp->block = NULL;

      // wait until the input is empty
      queue_wait_until_empty(fp->input);

      // reset block
      fp->block = block_init();
      fp->block_offset = 0;
      fp->block->block_length = BGZF_MAX_BLOCK_SIZE;
  }
  else {
      // wait until the input is empty
      queue_wait_until_empty(fp->input);
  }

  // close the input queue 
  queue_close(fp->input);

  // wake all
  queue_wake_all(fp->input);

  // join
  pbgzf_join(fp);

  // flush the underlying stream
  ret = bgzf_flush(fp->w->fp_bgzf);
  if(0 != ret) return ret;

  if(1 == restart) {
      // reset the producer/consumer/outputter
      if(NULL != fp->p) producer_reset(fp->p);
      if(NULL != fp->c) consumers_reset(fp->c);
      if(NULL != fp->o) outputter_reset(fp->o);

      // reset the queues
      queue_reset(fp->input, 1, fp->num_threads);
      queue_reset(fp->output, fp->num_threads, 1);

      // restart threads
      pbgzf_run(fp);
  }

  return 0;
}

int 
pbgzf_flush(PBGZF* fp)
{
  return pbgzf_flush_aux(fp, 1);
}

int 
pbgzf_flush_try(PBGZF *fp, int size)
{
  if (fp->block->block_offset + size > fp->block->block_length) { 
      //NB: no need to restart the threads, just flush the current block

      if('w' != fp->open_mode) {
          fprintf(stderr, "file not open for writing\n");
          exit(1);
      }

      // flush
      if(0 < fp->block->block_offset) {
          fp->block->block_length = fp->block->block_offset;
          if(!queue_add(fp->input, fp->block, 1)) {
              fprintf(stderr, "pbgzf_flush_try queue_add: bug encountered\n");
              exit(1);
          }
          fp->block = NULL;

          // reset block
          fp->block = block_init();
          fp->block_offset = 0;
          fp->block->block_length = BGZF_MAX_BLOCK_SIZE;
      }
  }
  return -1;
}

int 
pbgzf_close(PBGZF* fp)
{
  int32_t i;
  if(NULL == fp) return 0;
  if('w' == fp->open_mode) {
      // flush the data to the output file
      pbgzf_flush_aux(fp, 0); // NB: no need to restart the threads
  }
  else {
      // shut down the reader, regardless of where it is reading
      fp->r->is_done = 1;

      // shut down the consumers, regardless of what they are consuming
      for(i=0;i<fp->c->n;i++) {
          fp->c->c[i]->is_done = 1;
      }
  
      // close the input queue
      queue_close(fp->input); // NB: consumers should shut down when the input queue has EOF set
      queue_close(fp->output); // NB: consumers should shut down when the input queue has EOF set

      // wake all of the threads
      queue_wake_all(fp->input);
      queue_wake_all(fp->output);

      // TODO: faster signal to the consumers that computation is complete
      // join
      pbgzf_join(fp);
  }

  // destroy
  if(NULL != fp->c) consumers_destroy(fp->c);
  if(NULL != fp->p) producer_destroy(fp->p);
  if(NULL != fp->o) outputter_destroy(fp->o);
  if(NULL != fp->input) queue_destroy(fp->input);
  if(NULL != fp->output) queue_destroy(fp->output);
  if(NULL != fp->r) reader_destroy(fp->r);
  if(NULL != fp->w) writer_destroy(fp->w);

  if(NULL != fp->block) block_destroy(fp->block);
  fp->block = NULL;
  block_pool_destroy(fp->pool);
  free(fp);

  return 0;
}

void pbgzf_set_cache_size(PBGZF *fp, int cache_size)
{
  if(fp && 'r' == fp->open_mode) bgzf_set_cache_size(fp->r->fp_bgzf, cache_size);
}

void
pbgzf_main(int f_src, int f_dst, int compress, int compress_level, int compress_type, int queue_size, int num_threads)
{
  // NB: this gives us greater control over queue size and the like
  queue_t *input = NULL;
  queue_t *output = NULL;
  reader_t *r = NULL;
  writer_t *w = NULL;
  consumers_t *c = NULL;
  producer_t *p = NULL;
  outputter_t *o = NULL;
  block_pool_t *pool = NULL;

  pool = block_pool_init(PBGZF_BLOCKS_POOL_NUM);
  input = queue_init(queue_size, 0, 1, num_threads);
  output = queue_init(queue_size, 1, num_threads, 1);

  r = reader_init(f_src, input, compress, pool);
  w = writer_init(f_dst, output, compress, compress_level, compress_type, pool);
  c = consumers_init(num_threads, input, output, r, compress, compress_level, compress_type);
  p = producer_init(r);
  o = outputter_init(w);

  producer_run(p);
  consumers_run(c);
  outputter_run(o);

  producer_join(p);
  consumers_join(c);
  outputter_join(o);

  consumers_destroy(c);
  producer_destroy(p);
  outputter_destroy(o);

  queue_destroy(input);
  queue_destroy(output);
  reader_destroy(r);
  writer_destroy(w);
  block_pool_destroy(pool);
}
