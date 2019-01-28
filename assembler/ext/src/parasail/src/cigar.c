#include "config.h"

#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"

#if defined(_MSC_VER)
#define snprintf _snprintf
#endif

/*                     0123456789 */
#define BAM_CIGAR_STR "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4u

/* array index is an ASCII character value from a CIGAR,
   element value is the corresponding integer opcode between 0 and 9 */
const uint8_t parasail_cigar_encoded_ops[] = {
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0 /*   */, 0 /* ! */, 0 /* " */, 0 /* # */,
    0 /* $ */, 0 /* % */, 0 /* & */, 0 /* ' */,
    0 /* ( */, 0 /* ) */, 0 /* * */, 0 /* + */,
    0 /* , */, 0 /* - */, 0 /* . */, 0 /* / */,
    0 /* 0 */, 0 /* 1 */, 0 /* 2 */, 0 /* 3 */,
    0 /* 4 */, 0 /* 5 */, 0 /* 6 */, 0 /* 7 */,
    0 /* 8 */, 0 /* 9 */, 0 /* : */, 0 /* ; */,
    0 /* < */, 7 /* = */, 0 /* > */, 0 /* ? */,
    0 /* @ */, 0 /* A */, 0 /* B */, 0 /* C */,
    2 /* D */, 0 /* E */, 0 /* F */, 0 /* G */,
    5 /* H */, 1 /* I */, 0 /* J */, 0 /* K */,
    0 /* L */, 0 /* M */, 3 /* N */, 0 /* O */,
    6 /* P */, 0 /* Q */, 0 /* R */, 4 /* S */,
    0 /* T */, 0 /* U */, 0 /* V */, 0 /* W */,
    8 /* X */, 0 /* Y */, 0 /* Z */, 0 /* [ */,
    0 /* \ */, 0 /* ] */, 0 /* ^ */, 0 /* _ */,
    0 /* ` */, 0 /* a */, 0 /* b */, 0 /* c */,
    0 /* d */, 0 /* e */, 0 /* f */, 0 /* g */,
    0 /* h */, 0 /* i */, 0 /* j */, 0 /* k */,
    0 /* l */, 0 /* m */, 0 /* n */, 0 /* o */,
    0 /* p */, 0 /* q */, 0 /* r */, 0 /* s */,
    0 /* t */, 0 /* u */, 0 /* v */, 0 /* w */,
    0 /* x */, 0 /* y */, 0 /* z */, 0 /* { */,
    0 /* | */, 0 /* } */, 0 /* ~ */, 0 /*  */
};

uint32_t parasail_cigar_encode(uint32_t length, char op_letter)
{
    return (length << BAM_CIGAR_SHIFT) | (parasail_cigar_encoded_ops[(int)op_letter]);
}

parasail_cigar_t* parasail_cigar_encode_string(const char *cigar)
{
    int sscanf_retcode = 0;
    size_t offset = 0;
    int chars_read = 0;
    unsigned int len = 0;
    char op = 'M';
    int done = 0;
    size_t string_length = 0;
    size_t size = 0;
    parasail_cigar_t *ret = NULL;

    string_length = strlen(cigar);
    size = sizeof(uint32_t)*string_length;
    ret = malloc(sizeof(parasail_cigar_t));
    ret->seq = malloc(size);
    ret->len = 0;

    while (!done) {
        sscanf_retcode = sscanf(
                &cigar[offset], "%u%c%n", &len, &op, &chars_read);
        if (2 != sscanf_retcode) {
            fprintf(stderr, "invalid CIGAR string\n");
            parasail_cigar_free(ret);
            return NULL;
        }
        offset += chars_read;
        ret->len += 1;
        if ((size_t)ret->len >= size) {
            size *= 2;
            ret->seq = realloc(ret->seq, size);
        }
        ret->seq[ret->len-1] = parasail_cigar_encode(len, op);
        if (offset >= string_length) {
            done = 1;
        }
    }

    return ret;
}

char parasail_cigar_decode_op(uint32_t cigar_int) {
    return (cigar_int & 0xfU) > 9 ? 'M': BAM_CIGAR_STR[cigar_int & 0xfU];
}

uint32_t parasail_cigar_decode_len(uint32_t cigar_int) {
    return cigar_int >> BAM_CIGAR_SHIFT;
}

char* parasail_cigar_decode(parasail_cigar_t *cigar)
{
#define SIZE 40
    char *ret = NULL;
    size_t retlen = 0;
    size_t size = 0;
    int i = 0;

    /* initial allocation for 1 op and 3 number characters per cigar int */
    size = sizeof(char)*cigar->len*4;
    ret = malloc(size+1);
    ret[0] = '\0';

    for (i=0; i<cigar->len; ++i) {
        char tmp[SIZE];
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        int snprintf_retcode = snprintf(tmp, SIZE, "%u%c", len, op);
        if (snprintf_retcode >= SIZE) {
            fprintf(stderr, "invalid CIGAR\n");
            free(ret);
            return NULL;
        }
        retlen += snprintf_retcode;
        if (retlen >= size) {
            size *= 2;
            ret = realloc(ret, size+1);
        }
        strcat(ret, tmp);
    }

    return ret;
}

void parasail_cigar_free(parasail_cigar_t *cigar)
{
    free(cigar->seq);
    free(cigar);
}

#define CONCAT_(X, Y) X##Y
#define CONCAT(X, Y) CONCAT_(X, Y)
#define CONCAT3_(X, Y, Z) X##Y##Z
#define CONCAT3(X, Y, Z) CONCAT3_(X, Y, Z)
#define LOC_NOVEC int64_t loc = i*lenb + j;
#define LOC_STRIPED int64_t loc = j*segLen*segWidth + (i%segLen)*segWidth + (i/segLen);

#define T 8
#include "cigar_template.c"
#undef T

#define T 8
#define STRIPED
#include "cigar_template.c"
#undef T
#undef STRIPED

#define T 16
#include "cigar_template.c"
#undef T

#define T 16
#define STRIPED
#include "cigar_template.c"
#undef T
#undef STRIPED

#define T 32
#include "cigar_template.c"
#undef T

#define T 32
#define STRIPED
#include "cigar_template.c"
#undef T
#undef STRIPED

#define T 64
#include "cigar_template.c"
#undef T

#define T 64
#define STRIPED
#include "cigar_template.c"
#undef T
#undef STRIPED

parasail_cigar_t* parasail_result_get_cigar(
        parasail_result_t *result,
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix)
{
    assert(parasail_result_is_trace(result));

    if (result->flag & PARASAIL_FLAG_STRIPED || result->flag & PARASAIL_FLAG_SCAN) {
        if (result->flag & PARASAIL_FLAG_BITS_8) {
            return parasail_cigar_striped_8(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_16) {
            return parasail_cigar_striped_16(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_32) {
            return parasail_cigar_striped_32(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_64) {
            return parasail_cigar_striped_64(seqA, lena, seqB, lenb, matrix, result);
        }
    }
    else {
        return parasail_cigar_8(seqA, lena, seqB, lenb, matrix, result);
    }

    /* should not get here, but to silence warnings */
    return NULL;
}

