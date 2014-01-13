#ifndef CONSUMER_H_
#define CONSUMER_H_

typedef struct {
    queue_t *input;
    queue_t *output;
    reader_t *reader;
    uint8_t *buffer;
    int8_t is_done;
    int8_t compress;
    int32_t compress_level;
    int32_t compress_type;
    int16_t cid;
    int64_t n;
} consumer_t;

consumer_t*
consumer_init(queue_t *input,
              queue_t *output,
              reader_t *reader,
              int8_t compress,
              int32_t compress_level,
              int32_t compress_type,
              int32_t cid);

void*
consumer_run(void *arg);

void
consumer_destroy(consumer_t *c);

void
consumer_reset(consumer_t *c);

#endif
