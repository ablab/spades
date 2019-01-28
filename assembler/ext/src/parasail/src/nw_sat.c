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
#include <string.h>

#include "parasail/parasail.h"


parasail_result_t* parasail_nw_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

parasail_result_t* parasail_nw_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_scan_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = parasail_nw_scan_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_striped_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = parasail_nw_striped_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = parasail_nw_table_scan_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = parasail_nw_table_striped_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = parasail_nw_rowcol_scan_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = parasail_nw_rowcol_striped_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_scan_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = parasail_nw_stats_scan_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_striped_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = parasail_nw_stats_striped_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = parasail_nw_stats_table_scan_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = parasail_nw_stats_table_striped_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = parasail_nw_stats_rowcol_scan_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = parasail_nw_stats_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = parasail_nw_stats_rowcol_striped_profile_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}

