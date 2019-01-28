#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <parasail.h>
#include <parasail/matrices/blosum62.h>

void compare(const char *s1, const char *s2)
{
    parasail_profile_t *profile = NULL;
    parasail_result_t *result = NULL;
    parasail_cigar_t *cigar = NULL;
    char *cigar_str = NULL;
    int s1Len = (int)strlen(s1);
    int s2Len = (int)strlen(s2);
    char *sat = "OKAY";
    
    printf("comparing '%s'\n"
           "          '%s'\n", s1, s2);

    profile = parasail_profile_create_sat(s1, s1Len, &parasail_blosum62);
    result = parasail_nw_trace_scan_profile_sat(profile, s2, s2Len, 9, 1);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "prof sat", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
    parasail_profile_free(profile);

    result = parasail_nw_trace_striped_sse41_128_8(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "striped 8", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_striped_sse41_128_16(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "striped 16", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_striped_sse41_128_32(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "striped 32", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_striped_sse41_128_64(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "striped 64", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_scan_sse41_128_8(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "scan 8", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_scan_sse41_128_16(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "scan 16", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_scan_sse41_128_32(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "scan 32", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_scan_sse41_128_64(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "scan 64", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    result = parasail_nw_trace_scan(s1, s1Len, s2, s2Len, 9, 1, &parasail_blosum62);
    cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, &parasail_blosum62);
    cigar_str = parasail_cigar_decode(cigar);
    sat = parasail_result_is_saturated(result) ? "FAIL" : "OKAY";
    printf("%s %10s cigar_str='%s'\n", sat, "scan", cigar_str);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
}

int main(int argc, char **argv)
{
    const char *s1 = NULL;
    const char *s2 = NULL;

    s1 = "AAA";
    s2 = "AAA";
    compare(s1, s2);

    s1 = "AAAAAAAAAAAA";
    s2 = "AAAAAAAAAAAA";
    compare(s1, s2);

    s1 = "AAAAAAAAAAAAAAAAA";
    s2 = "AAAAAAAAAAAAAAAAA";
    compare(s1, s2);

    s1 = "AAAAAAAAAAAAAAAAAAAAAAAAA";
    s2 = "AAAAAAAAAAAAAAAAAAAAAAAAA";
    compare(s1, s2);

    s1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    s2 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    compare(s1, s2);

    s1 = "AAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAA";
    s2 = "AAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAA";
    compare(s1, s2);

    s1 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    s2 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    compare(s1, s2);

    return 0;
}

