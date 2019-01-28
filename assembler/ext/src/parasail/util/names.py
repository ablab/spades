#!/usr/bin/env python

# serial reference implementations (3x2x3 = 18 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            print ""
            print "extern parasail_result_t* parasail_"+a+s+t+'('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# serial scan reference implementations (3x2x3 = 18 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            print ""
            print "extern parasail_result_t* parasail_"+a+s+t+'_scan('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

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
                    print "extern parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

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
                    print "extern parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const parasail_profile_t * const restrict profile,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap);"

# vectorized implementations of blocked (1x1x3x1x2 = 6 impl)
alg = ["sw"]
stats = [""]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_blocked"]
isa = ["_sse41_128_32", "_sse41_128_16"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t: continue
            for p in par:
                for i in isa:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# dispatching implementations (3x2x3x3x4 = 216 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for w in width:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# dispatching profile implementations (3x2x3x2x4 = 144 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for w in width:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const parasail_profile_t * const restrict profile,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap);"

# profile creation functions (2x13 = 26 impl)
stats = ["", "_stats"]
isa = [
    "_sse_128_64", "_sse_128_32", "_sse_128_16", "_sse_128_8", "_sse_128_sat",
    "_avx_256_64", "_avx_256_32", "_avx_256_16", "_avx_256_8", "_avx_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat",
    "_64", "_32", "_16", "_8", "_sat"
    ]
for s in stats:
    for i in isa:
        print ""
        print "extern parasail_profile_t* parasail_profile_create"+s+i+'('
        print " "*8+"const char * const restrict s1, const int s1Len,"
        print " "*8+"const parasail_matrix_t* matrix);"

# dispatching saturation check implementations (3x2x3x3 = 54 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                print ""
                print "extern parasail_result_t* parasail_"+a+s+t+p+'_sat('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"const parasail_matrix_t* matrix);"

