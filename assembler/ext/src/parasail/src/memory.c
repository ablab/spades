/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"

void* parasail_memalign(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(HAVE__ALIGNED_MALLOC)
    ptr = _aligned_malloc(size, alignment);
#elif defined(HAVE_POSIX_MEMALIGN)
    int retcode = posix_memalign(&ptr, alignment, size);
    assert(0 == retcode);
#elif defined(HAVE_ALIGNED_ALLOC)
    ptr = aligned_alloc(alignment, size);
#elif defined(HAVE_MEMALIGN)
    ptr = memalign(alignment, size);
#else
#error "No suitable memory alignment routine found."
#endif
    assert(NULL != ptr);
    return ptr;
}

void parasail_free(void *ptr)
{
#if defined(HAVE__ALIGNED_MALLOC)
     _aligned_free(ptr);
#else
    free(ptr);
#endif
}

void parasail_free_unaligned(void *ptr)
{
    free(ptr);
}

int * parasail_memalign_int(size_t alignment, size_t size)
{
    return (int *) parasail_memalign(alignment, size*sizeof(int));
}

int8_t * parasail_memalign_int8_t(size_t alignment, size_t size)
{
    return (int8_t *) parasail_memalign(alignment, size*sizeof(int8_t));
}

int16_t * parasail_memalign_int16_t(size_t alignment, size_t size)
{
    return (int16_t *) parasail_memalign(alignment, size*sizeof(int16_t));
}

int32_t * parasail_memalign_int32_t(size_t alignment, size_t size)
{
    return (int32_t *) parasail_memalign(alignment, size*sizeof(int32_t));
}

int64_t * parasail_memalign_int64_t(size_t alignment, size_t size)
{
    return (int64_t *) parasail_memalign(alignment, size*sizeof(int64_t));
}

void parasail_memset(void *b, int c, size_t len)
{
    (void)memset(b, c, len);
}

void parasail_memset_int(int *b, int c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int8_t(int8_t *b, int8_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int16_t(int16_t *b, int16_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int32_t(int32_t *b, int32_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int64_t(int64_t *b, int64_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

parasail_result_t* parasail_result_new()
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    result = (parasail_result_t*)malloc(sizeof(parasail_result_t));
    assert(result);

    result->score = 0;
    result->end_query = 0;
    result->end_ref = 0;
    result->flag = 0;
    result->extra = NULL;

    return result;
}

parasail_result_t* parasail_result_new_stats()
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* allocate struct to hold memory */
    result = parasail_result_new();

    /* allocate only tables */
    result->stats = (parasail_result_extra_stats_t*)
        malloc(sizeof(parasail_result_extra_stats_t));
    assert(result->stats);

    return result;
}

parasail_result_t* parasail_result_new_table1(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new();

    /* allocate only score table */
    result->tables = (parasail_result_extra_tables_t*)malloc(sizeof(parasail_result_extra_tables_t));
    assert(result->tables);
    result->tables->score_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->tables->score_table);

    return result;
}

parasail_result_t* parasail_result_new_rowcol1(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new();

    /* allocate only score col and row */
    result->rowcols = (parasail_result_extra_rowcols_t*)malloc(sizeof(parasail_result_extra_rowcols_t));
    assert(result->rowcols);
    result->rowcols->score_row = (int *)malloc(sizeof(int)*b);
    assert(result->rowcols->score_row);
    result->rowcols->score_col = (int *)malloc(sizeof(int)*a);
    assert(result->rowcols->score_col);

    return result;
}

parasail_result_t* parasail_result_new_table3(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new();
    
    /* allocate only tables */
    result->stats = (parasail_result_extra_stats_t*)
        malloc(sizeof(parasail_result_extra_stats_t));
    assert(result->stats);

    result->stats->tables = (parasail_result_extra_stats_tables_t*)
        malloc(sizeof(parasail_result_extra_stats_tables_t));
    assert(result->stats->tables);

    result->stats->tables->score_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->stats->tables->score_table);
    result->stats->tables->matches_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->stats->tables->matches_table);
    result->stats->tables->similar_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->stats->tables->similar_table);
    result->stats->tables->length_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->stats->tables->length_table);

    return result;
}

parasail_result_t* parasail_result_new_trace(const int a, const int b, const size_t alignment, const size_t size)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);

    /* allocate struct to hold memory */
    result = parasail_result_new();

    result->trace = (parasail_result_extra_trace_t*)malloc(sizeof(parasail_result_extra_trace_t));
    assert(result->trace);
    result->trace->trace_table = parasail_memalign(alignment, size*a*b);
    assert(result->trace->trace_table);
    result->trace->trace_ins_table = NULL;
    result->trace->trace_del_table = NULL;

    return result;
}

parasail_result_t* parasail_result_new_rowcol3(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new();
    
    result->stats = (parasail_result_extra_stats_t*)
        malloc(sizeof(parasail_result_extra_stats_t));
    assert(result->stats);

    result->stats->rowcols = (parasail_result_extra_stats_rowcols_t*)
        malloc(sizeof(parasail_result_extra_stats_rowcols_t));
    assert(result->stats->rowcols);

    result->stats->rowcols->score_row = (int *)malloc(sizeof(int)*b);
    assert(result->stats->rowcols->score_row);
    result->stats->rowcols->matches_row = (int *)malloc(sizeof(int)*b);
    assert(result->stats->rowcols->matches_row);
    result->stats->rowcols->similar_row = (int *)malloc(sizeof(int)*b);
    assert(result->stats->rowcols->similar_row);
    result->stats->rowcols->length_row = (int *)malloc(sizeof(int)*b);
    assert(result->stats->rowcols->length_row);

    result->stats->rowcols->score_col = (int *)malloc(sizeof(int)*a);
    assert(result->stats->rowcols->score_col);
    result->stats->rowcols->matches_col = (int *)malloc(sizeof(int)*a);
    assert(result->stats->rowcols->matches_col);
    result->stats->rowcols->similar_col = (int *)malloc(sizeof(int)*a);
    assert(result->stats->rowcols->similar_col);
    result->stats->rowcols->length_col = (int *)malloc(sizeof(int)*a);
    assert(result->stats->rowcols->length_col);

    return result;
}

void parasail_result_free(parasail_result_t *result)
{
    /* validate inputs */
    assert(NULL != result);
    
    if (result->flag & PARASAIL_FLAG_STATS) {
        if (result->flag & PARASAIL_FLAG_TABLE) {
            free(result->stats->tables->score_table);
            free(result->stats->tables->matches_table);
            free(result->stats->tables->similar_table);
            free(result->stats->tables->length_table);
            free(result->stats->tables);
        }
        if (result->flag & PARASAIL_FLAG_ROWCOL) {
            free(result->stats->rowcols->score_row);
            free(result->stats->rowcols->matches_row);
            free(result->stats->rowcols->similar_row);
            free(result->stats->rowcols->length_row);
            free(result->stats->rowcols->score_col);
            free(result->stats->rowcols->matches_col);
            free(result->stats->rowcols->similar_col);
            free(result->stats->rowcols->length_col);
            free(result->stats->rowcols);
        }
        free(result->stats);
    }
    else {
        if (result->flag & PARASAIL_FLAG_TABLE) {
            free(result->tables->score_table);
            free(result->tables);
        }
        if (result->flag & PARASAIL_FLAG_ROWCOL) {
            free(result->rowcols->score_row);
            free(result->rowcols->score_col);
            free(result->rowcols);
        }
    }
    if (result->flag & PARASAIL_FLAG_TRACE) {
        parasail_free(result->trace->trace_table);
        if (NULL != result->trace->trace_ins_table)
            parasail_free(result->trace->trace_ins_table);
        if (NULL != result->trace->trace_del_table)
            parasail_free(result->trace->trace_del_table);
        free(result->trace);
    }

    free(result);
}

void parasail_version(int *major, int *minor, int *patch)
{
    *major = PARASAIL_VERSION_MAJOR;
    *minor = PARASAIL_VERSION_MINOR;
    *patch = PARASAIL_VERSION_PATCH;
}

parasail_matrix_t* parasail_matrix_create(
        const char *alphabet, const int match, const int mismatch)
{
    parasail_matrix_t *retval = NULL;
    int *matrix = NULL;
    int *mapper = NULL;
    size_t size = 0;
    size_t size1 = 0;
    size_t i = 0;
    size_t j = 0;
    size_t c = 0;

    size = strlen(alphabet);
    assert(size < INT_MAX);
    size1 = size + 1;

    matrix = (int*)malloc(sizeof(int)*size1*size1);
    assert(matrix);
    for (i=0; i<size; ++i) {
        for (j=0; j<size; ++j) {
            if (i == j) {
                matrix[c++] = match;
            }
            else {
                matrix[c++] = mismatch;
            }
        }
        matrix[c++] = 0;
    }
    for (j=0; j<size1; ++j) {
        matrix[c++] = 0;
    }

    mapper = (int*)malloc(sizeof(int)*256);
    assert(mapper);
    parasail_memset_int(mapper, (int)size, 256);
    for (i=0; i<size; ++i) {
        mapper[toupper((unsigned char)alphabet[i])] = (int)i;
        mapper[tolower((unsigned char)alphabet[i])] = (int)i;
    }

    retval = (parasail_matrix_t*)malloc(sizeof(parasail_matrix_t));
    assert(retval);
    retval->name = "";
    retval->matrix = matrix;
    retval->mapper = mapper;
    retval->size = (int)size1;
    retval->max = match > mismatch ? match : mismatch;
    retval->min = match > mismatch ? mismatch : match;
    retval->user_matrix = matrix;
    return retval;
}

parasail_matrix_t* parasail_matrix_copy(const parasail_matrix_t *original)
{
    parasail_matrix_t *retval = NULL;

    retval = (parasail_matrix_t*)malloc(sizeof(parasail_matrix_t));
    assert(retval);
    retval->name = original->name;
    retval->size = original->size;
    retval->max = original->max;
    retval->min = original->min;

    {
        size_t matrix_size = sizeof(int)*original->size*original->size;
        size_t mapper_size = sizeof(int)*256;
        int *new_mapper = NULL;
        int *new_matrix = NULL;

        new_mapper = (int*)malloc(mapper_size);
        assert(new_mapper);
        (void)memcpy(new_mapper, original->mapper, mapper_size);

        new_matrix = (int*)malloc(matrix_size);
        assert(new_matrix);
        (void)memcpy(new_matrix, original->matrix, matrix_size);

        retval->mapper = new_mapper;
        retval->matrix = new_matrix;
        retval->user_matrix = new_matrix;
    }

    return retval;
}

void parasail_matrix_set_value(parasail_matrix_t *matrix, int row, int col, int value)
{
    assert(matrix);

    if (NULL == matrix->user_matrix) {
        fprintf(stderr, "attempted to set value of built-in matrix '%s'\n",
                matrix->name);
        return;
    }

    matrix->user_matrix[row*matrix->size + col] = value;
    if (value > matrix->max) {
        matrix->max = value;
    }
    if (value < matrix->min) {
        matrix->min = value;
    }
}

void parasail_matrix_free(parasail_matrix_t *matrix)
{
    /* validate inputs */
    assert(NULL != matrix);
    if (NULL != matrix->user_matrix) {
        free((void*)matrix->matrix);
        free((void*)matrix->mapper);
        free(matrix);
    }
    else {
        fprintf(stderr, "attempted to free built-in matrix '%s'\n",
                matrix->name);
    }
}

parasail_profile_t* parasail_profile_new(
        const char * s1, const int s1Len, const parasail_matrix_t *matrix)
{
    /* declare all variables */
    parasail_profile_t *profile = NULL;

    profile = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
    assert(profile);

    profile->s1 = s1;
    profile->s1Len = s1Len;
    profile->matrix = matrix;
    profile->profile8.score = NULL;
    profile->profile8.matches = NULL;
    profile->profile8.similar = NULL;
    profile->profile16.score = NULL;
    profile->profile16.matches = NULL;
    profile->profile16.similar = NULL;
    profile->profile32.score = NULL;
    profile->profile32.matches = NULL;
    profile->profile32.similar = NULL;
    profile->profile64.score = NULL;
    profile->profile64.matches = NULL;
    profile->profile64.similar = NULL;
    profile->free = NULL;
    profile->stop = INT32_MAX;

    return profile;
}

void parasail_profile_free(parasail_profile_t *profile)
{
    if (NULL != profile->profile8.score) {
        profile->free(profile->profile8.score);
    }
    if (NULL != profile->profile8.matches) {
        profile->free(profile->profile8.matches);
    }
    if (NULL != profile->profile8.similar) {
        profile->free(profile->profile8.similar);
    }

    if (NULL != profile->profile16.score) {
        profile->free(profile->profile16.score);
    }
    if (NULL != profile->profile16.matches) {
        profile->free(profile->profile16.matches);
    }
    if (NULL != profile->profile16.similar) {
        profile->free(profile->profile16.similar);
    }

    if (NULL != profile->profile32.score) {
        profile->free(profile->profile32.score);
    }
    if (NULL != profile->profile32.matches) {
        profile->free(profile->profile32.matches);
    }
    if (NULL != profile->profile32.similar) {
        profile->free(profile->profile32.similar);
    }

    if (NULL != profile->profile64.score) {
        profile->free(profile->profile64.score);
    }
    if (NULL != profile->profile64.matches) {
        profile->free(profile->profile64.matches);
    }
    if (NULL != profile->profile64.similar) {
        profile->free(profile->profile64.similar);
    }

    free(profile);
}

char* parasail_reverse(const char *s, size_t length)
{
    char *r = NULL;
    size_t i = 0;
    size_t j = 0;

    r = (char*)malloc(sizeof(char)*(length + 1));
    r[length] = '\0';
    for (i=0,j=length-1; i<length; ++i,--j) {
        r[i] = s[j];
    }

    return r;
}

uint32_t* parasail_reverse_uint32_t(const uint32_t *s, size_t length)
{
    uint32_t *r = NULL;
    size_t i = 0;
    size_t j = 0;

    r = (uint32_t*)malloc(sizeof(uint32_t)*(length));
    for (i=0,j=length-1; i<length; ++i,--j) {
        r[i] = s[j];
    }

    return r;
}

int parasail_result_is_nw(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_NW;
}

int parasail_result_is_sg(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_SG;
}

int parasail_result_is_sw(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_SW;
}

int parasail_result_is_saturated(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_SATURATED;
}

int parasail_result_is_banded(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_BANDED;
}

int parasail_result_is_scan(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_SCAN;
}

int parasail_result_is_striped(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_STRIPED;
}

int parasail_result_is_diag(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_DIAG;
}

int parasail_result_is_blocked(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_BLOCKED;
}

int parasail_result_is_stats(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_STATS;
}

int parasail_result_is_stats_table(const parasail_result_t * const restrict result)
{
    return (result->flag & PARASAIL_FLAG_STATS) && (result->flag & PARASAIL_FLAG_TABLE);
}

int parasail_result_is_stats_rowcol(const parasail_result_t * const restrict result)
{
    return (result->flag & PARASAIL_FLAG_STATS) && (result->flag & PARASAIL_FLAG_ROWCOL);
}

int parasail_result_is_table(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_TABLE;
}

int parasail_result_is_rowcol(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_ROWCOL;
}

int parasail_result_is_trace(const parasail_result_t * const restrict result)
{
    return result->flag & PARASAIL_FLAG_TRACE;
}

int parasail_result_get_score(const parasail_result_t * const restrict result)
{
    return result->score;
}

int parasail_result_get_end_query(const parasail_result_t * const restrict result)
{
    return result->end_query;
}

int parasail_result_get_end_ref(const parasail_result_t * const restrict result)
{
    return result->end_ref;
}

int parasail_result_get_matches(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats(result));
    return result->stats->matches;
}

int parasail_result_get_similar(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats(result));
    return result->stats->similar;
}

int parasail_result_get_length(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats(result));
    return result->stats->length;
}

int* parasail_result_get_score_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_table(result) || parasail_result_is_stats_table(result));
    if (parasail_result_is_stats_table(result)) {
        return result->stats->tables->score_table;
    }
    if (parasail_result_is_table(result)) {
        return result->tables->score_table;
    }
    return NULL; /* should not reach */
}

int* parasail_result_get_matches_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_table(result));
    return result->stats->tables->matches_table;
}

int* parasail_result_get_similar_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_table(result));
    return result->stats->tables->similar_table;
}

int* parasail_result_get_length_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_table(result));
    return result->stats->tables->length_table;
}

int* parasail_result_get_score_row(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result) || parasail_result_is_rowcol(result));
    if (parasail_result_is_stats_rowcol(result)) {
        return result->stats->rowcols->score_row;
    }
    if (parasail_result_is_rowcol(result)) {
        return result->rowcols->score_row;
    }
    return NULL; /* should not reach */
}

int* parasail_result_get_matches_row(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->matches_row;
}

int* parasail_result_get_similar_row(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->similar_row;
}

int* parasail_result_get_length_row(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->length_row;
}

int* parasail_result_get_score_col(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result) || parasail_result_is_rowcol(result));
    if (parasail_result_is_stats_rowcol(result)) {
        return result->stats->rowcols->score_col;
    }
    if (parasail_result_is_rowcol(result)) {
        return result->rowcols->score_col;
    }
    return NULL; /* should not reach */
}

int* parasail_result_get_matches_col(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->matches_col;
}

int* parasail_result_get_similar_col(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->similar_col;
}

int* parasail_result_get_length_col(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->length_col;
}

int* parasail_result_get_trace_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_trace(result));
    return result->trace->trace_table;
}

int* parasail_result_get_trace_ins_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_trace(result));
    return result->trace->trace_ins_table;
}

int* parasail_result_get_trace_del_table(const parasail_result_t * const restrict result)
{
    assert(parasail_result_is_trace(result));
    return result->trace->trace_del_table;
}

