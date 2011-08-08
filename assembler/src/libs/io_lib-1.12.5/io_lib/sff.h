#ifndef _SFF_H_
#define _SFF_H_

#include "io_lib/Read.h"
#include "io_lib/os.h"
#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This mirrors the order that the SFF header has on disc. We have one of
 * these only.
 *
 * The flow and key are variable length fields as defined by flow_len
 * and key_len.
 * The on-disc structure is then also padded out with zeros to the 8-byte
 * word boundary.
 */
typedef struct {
    uint32_t magic;
    char      version[4];
    uint64_t  index_offset;
    uint32_t  index_len;
    uint32_t  nreads;
    uint16_t  header_len;
    uint16_t  key_len;
    uint16_t  flow_len;
    uint8_t   flowgram_format;
    char     *flow;
    char     *key;
} sff_common_header;

#define SFF_MAGIC   0x2e736666 /* ".sff" */
#define SFF_VERSION "\0\0\0\1"

/*
 * We have one read_header per "reading" in the SFF archive.
 * It too is padded to an 8-byte boundary.
 */
typedef struct {
    uint16_t  header_len;
    uint16_t  name_len;
    uint32_t  nbases;
    uint16_t  clip_qual_left;
    uint16_t  clip_qual_right;
    uint16_t  clip_adapter_left;
    uint16_t  clip_adapter_right;
    char     *name;
} sff_read_header;

/*
 * We have one read_data section per reading, following the read_header.
 * It is padded to an 8-byte boundary.
 */
typedef struct {
    uint16_t *flowgram;   /* x 100.0 */
    uint8_t  *flow_index; /* relative to last */
    char     *bases;
    uint8_t  *quality;
} sff_read_data;


/*
 * Low level functions to decode SFF internals
 */
sff_common_header *decode_sff_common_header(unsigned char *buf);
sff_common_header *read_sff_common_header(mFILE *mf);
void free_sff_common_header(sff_common_header *h);
sff_read_header *decode_sff_read_header(unsigned char *buf);
void free_sff_read_header(sff_read_header *h);
sff_read_header *read_sff_read_header(mFILE *mf);
void free_sff_read_data(sff_read_data *d);
sff_read_data *read_sff_read_data(mFILE *mf, int nflows, int nbases);


/*
 * Loads the first SFF sequence from an SFF container and returns as a Read.
 */
Read *mfread_sff(mFILE *mf);

#ifdef __cplusplus
}
#endif

#endif /* _SFF_H_ */
