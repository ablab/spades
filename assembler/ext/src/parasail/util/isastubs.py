#!/usr/bin/env python

# header stuff
print """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * These functions borrow heavily, if not verbatim, from the files in
 * the parsail contrib directory. Please see contrib for details and
 * copyrights of the respective authors/sources.
 */
#include "config.h"

#include "parasail/parasail.h"

#include <errno.h>
#include <stdlib.h>

#define UNUSED(expr) do { (void)(expr); } while (0)"""

# so we don't end up with an empty source file if no ISAs are supported
print """
extern
parasail_result_t* parasail_isastub_dummy() {
    errno = ENOSYS;
    return NULL;
}
"""

def isa_to_guard(isa):
    if '_sse_' in isa:
        print '#if HAVE_SSE2 || HAVE_SSE41'
        print '#else'
    elif '_avx_' in isa:
        print "#if HAVE_AVX2"
        print '#else'
    elif 'sse2' in isa:
        print "#if HAVE_SSE2"
        print '#else'
    elif 'sse41' in isa:
        print "#if HAVE_SSE41"
        print '#else'
    elif 'avx2' in isa:
        print "#if HAVE_AVX2"
        print '#else'
    elif 'altivec' in isa:
        print "#if HAVE_ALTIVEC"
        print '#else'
    elif 'neon' in isa:
        print "#if HAVE_NEON"
        print '#else'
    else:
        print isa
        assert False

def body1():
    print """{
    UNUSED(s1);
    UNUSED(s1Len);
    UNUSED(s2);
    UNUSED(s2Len);
    UNUSED(open);
    UNUSED(gap);
    UNUSED(matrix);
    errno = ENOSYS;
    return NULL;
}
#endif"""

def body2():
    print """{
    UNUSED(profile);
    UNUSED(s2);
    UNUSED(s2Len);
    UNUSED(open);
    UNUSED(gap);
    errno = ENOSYS;
    return NULL;
}
#endif"""

def body3():
    print """{
    UNUSED(s1);
    UNUSED(s1Len);
    UNUSED(matrix);
    errno = ENOSYS;
    return NULL;
}
#endif"""

# vectorized implementations (3x2x3x3x13 = 702 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat"
    ]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for i in isa:
                    print ""
                    isa_to_guard(i)
                    print "extern"
                    print "parasail/parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix)"
                    body1()

# vectorized profile implementations (3x2x3x2x13 = 468 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat"
    ]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for i in isa:
                    print ""
                    isa_to_guard(i)
                    print "extern"
                    print "parasail/parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const parasail_profile_t * const restrict profile,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap)"
                    body2()

# vectorized implementations of blocked (1x1x3x1x2 = 6 impl)
alg = ["sw"]
stats = [""]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_blocked"]
isa = ["_sse41_128_32", "_sse41_128_16"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for i in isa:
                    print ""
                    isa_to_guard(i)
                    print "extern"
                    print "parasail/parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix)"
                    body1()

# profile creation functions (2x13 = 26 impl)
stats = ["", "_stats"]
isa = [
    "_sse_128_64", "_sse_128_32", "_sse_128_16", "_sse_128_8", "_sse_128_sat",
    "_avx_256_64", "_avx_256_32", "_avx_256_16", "_avx_256_8", "_avx_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat",
    ]
for s in stats:
    for i in isa:
        print ""
        isa_to_guard(i)
        print "extern"
        print "parasail/parasail_profile_t* parasail_profile_create"+s+i+'('
        print " "*8+"const char * const restrict s1, const int s1Len,"
        print " "*8+"const parasail_matrix_t* matrix)"
        body3()

print # for newline at end of file
