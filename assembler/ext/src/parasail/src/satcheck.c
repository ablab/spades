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
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_trace_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_trace_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_trace_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_trace_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_trace_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_trace_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_trace_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_trace_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_trace_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_trace_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_trace_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_trace_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_trace_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_trace_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_trace_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_trace_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_trace_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_trace_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_trace_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_trace_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_trace_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_trace_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_trace_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_trace_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_trace_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_trace_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_trace_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_trace_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_trace_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_trace_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_trace_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_trace_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_trace_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_trace_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_trace_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_trace_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_trace_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_trace_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_trace_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_trace_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_trace_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_diag_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_trace_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_diag_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_trace_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_diag_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_trace_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_diag_altivec_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_altivec_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_altivec_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_trace_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_diag_neon_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_neon_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_diag_neon_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_table_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_table_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_table_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_table_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_rowcol_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_rowcol_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_stats_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_stats_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_stats_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_stats_rowcol_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_stats_rowcol_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_trace_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_trace_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_trace_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_trace_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_trace_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_trace_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_nw_trace_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_nw_trace_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_nw_trace_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_nw_trace_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_nw_trace_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_nw_trace_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_nw_trace_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_nw_trace_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_table_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_table_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_table_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_table_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_rowcol_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_rowcol_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_stats_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_stats_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_stats_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_stats_rowcol_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_stats_rowcol_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_trace_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_trace_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_trace_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_trace_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_trace_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_trace_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sg_trace_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sg_trace_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sg_trace_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sg_trace_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sg_trace_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sg_trace_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sg_trace_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sg_trace_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_table_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_table_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_table_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_table_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_rowcol_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_rowcol_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_stats_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_stats_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_stats_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_stats_rowcol_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_stats_rowcol_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_trace_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_trace_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_trace_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_trace_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_trace_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_trace_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_scan_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_scan_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif


parasail_result_t* parasail_sw_trace_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_profile_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_32(profile, s2, s2Len, open, gap);
    }

    return result;
}


#if HAVE_SSE2
parasail_result_t* parasail_sw_trace_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_sse2_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_SSE41
parasail_result_t* parasail_sw_trace_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_profile_sse41_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_sse41_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_sse41_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_AVX2
parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_avx2_256_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_ALTIVEC
parasail_result_t* parasail_sw_trace_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_profile_altivec_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_altivec_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_altivec_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

#if HAVE_NEON
parasail_result_t* parasail_sw_trace_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;

    result = parasail_sw_trace_striped_profile_neon_128_8(profile, s2, s2Len, open, gap);
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_neon_128_16(profile, s2, s2Len, open, gap);
    }
    if (parasail_result_is_saturated(result)) {
        parasail_result_free(result);
        result = parasail_sw_trace_striped_profile_neon_128_32(profile, s2, s2Len, open, gap);
    }

    return result;
}
#endif

