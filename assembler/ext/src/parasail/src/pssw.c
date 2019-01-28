/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"

parasail_profile_t* parasail_ssw_init(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix, const int8_t score_size)
{
    parasail_profile_t *profile8 = NULL;
    parasail_profile_t *profile16 = NULL;

    if (score_size == 0 || score_size == 2) {
        profile8 = parasail_profile_create_8(s1, s1Len, matrix);
    }
    if (score_size == 1 || score_size == 2) {
        profile16 = parasail_profile_create_16(s1, s1Len, matrix);
    }

    if (profile8 && profile16) {
        profile8->profile16 = profile16->profile16;
        free(profile16);
        return profile8;
    }
    else if (profile8) {
        return profile8;
    }
    else if (profile16) {
        return profile16;
    }

    return NULL;
}

parasail_result_ssw_t* parasail_ssw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t* matrix)
{
    parasail_profile_t *profile = parasail_ssw_init(s1, s1Len, matrix, 2);
    parasail_result_ssw_t *result = parasail_ssw_profile(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_ssw_t* parasail_ssw_profile(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len_,
        const int open, const int gap)
{
    const char *s1 = profile->s1;
    const parasail_matrix_t *matrix = profile->matrix;
    parasail_result_t *result_forward = NULL;
    parasail_result_t *result_reverse = NULL;
    parasail_result_t *result_final = NULL;
    parasail_result_ssw_t *result_ssw = NULL;
    int has8 = 0;
    int has16 = 0;
    char *s1_reverse = NULL;
    char *s2_reverse = NULL;
    int s1Len = 0;
    int s2Len = 0;
    int s1Off = 0;
    int s2Off = 0;
    parasail_cigar_t *cigar = NULL;

    /* find the end loc first with the faster implementation */
    has8 = (NULL != profile->profile8.score);
    has16 = (NULL != profile->profile16.score);

    /* find the end loc first with the faster implementation */
    if (has8) {
        result_forward = parasail_sw_striped_profile_8(profile, s2, s2Len_, open, gap);
        if (parasail_result_is_saturated(result_forward)) {
            has8 = 0;
            parasail_result_free(result_forward);
            result_forward = NULL;
        }
    }
    if (NULL == result_forward && has16) {
        result_forward = parasail_sw_striped_profile_16(profile, s2, s2Len_, open, gap);
        if (parasail_result_is_saturated(result_forward)) {
            parasail_result_free(result_forward);
            result_forward = NULL;
        }
    }
    /* 8- and/or 16-bit options could fail */
    if (NULL == result_forward) {
        return NULL;
    }

    /* find beginning loc by going in reverse */
    s1Len = result_forward->end_query+1;
    s2Len = result_forward->end_ref+1;
    s1_reverse = parasail_reverse(s1, s1Len);
    s2_reverse = parasail_reverse(s2, s2Len);
    if (has8) {
        result_reverse = parasail_sw_striped_8(
                s1_reverse, s1Len,
                s2_reverse, s2Len,
                open, gap, matrix);
        if (parasail_result_is_saturated(result_reverse)) {
            has8 = 0;
            parasail_result_free(result_reverse);
            result_reverse = NULL;
        }
    }
    if (NULL == result_reverse && has16) {
        result_reverse = parasail_sw_striped_16(
                s1_reverse, s1Len,
                s2_reverse, s2Len,
                open, gap, matrix);
        if (parasail_result_is_saturated(result_reverse)) {
            parasail_result_free(result_reverse);
            result_reverse = NULL;
        }
    }
    free(s2_reverse);
    free(s1_reverse);

    /* 8- and/or 16-bit options could fail */
    if (NULL == result_reverse) {
        parasail_result_free(result_forward); /* free the forward result */
        return NULL;
    }

    /* run trace version of sw on just the aligned portion */
    s1Off = result_forward->end_query - result_reverse->end_query;
    s2Off = result_forward->end_ref - result_reverse->end_ref;
    s1Len = result_reverse->end_query+1;
    s2Len = result_reverse->end_ref+1;
    if (has8) {
        result_final = parasail_sw_trace_striped_8(
                &s1[s1Off], s1Len,
                &s2[s2Off], s2Len,
                open, gap, matrix);
        if (parasail_result_is_saturated(result_final)) {
            has8 = 0;
            parasail_result_free(result_final);
            result_final = NULL;
        }
    }
    if (NULL == result_final && has16) {
        result_final = parasail_sw_trace_striped_16(
                &s1[s1Off], s1Len,
                &s2[s2Off], s2Len,
                open, gap, matrix);
        if (parasail_result_is_saturated(result_final)) {
            parasail_result_free(result_final);
            result_final = NULL;
        }
    }

    /* 8- and/or 16-bit options could fail */
    if (NULL == result_final) {
        parasail_result_free(result_reverse); /* free reverse result */
        parasail_result_free(result_forward); /* free forward result */
        return NULL;
    }

    /* get cigar from traceback */
    cigar = parasail_result_get_cigar(result_final,
            &s1[s1Off], s1Len,
            &s2[s2Off], s2Len,
            matrix);

    result_ssw = (parasail_result_ssw_t*)malloc(sizeof(parasail_result_ssw_t));
    result_ssw->score1 = result_forward->score;
    result_ssw->ref_begin1 = s2Off;
    result_ssw->ref_end1 = result_forward->end_ref;
    result_ssw->read_begin1 = s1Off;
    result_ssw->read_end1 = result_forward->end_query;
    result_ssw->cigar = cigar->seq;
    result_ssw->cigarLen = cigar->len;

    free(cigar);
    parasail_result_free(result_final);
    parasail_result_free(result_reverse);
    parasail_result_free(result_forward);

    return result_ssw;;
}

void parasail_result_ssw_free(parasail_result_ssw_t *result)
{
    free(result->cigar);
    free(result);
}

