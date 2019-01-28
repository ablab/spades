#!/usr/bin/python

# Generate the various function dispatcher implemementations.

import os

def codegen(alg):
    txt = """/**
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
#include "parasail/parasail/cpuid.h"

/* forward declare the dispatcher functions */
"""
    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    txt += "parasail/parasail_function_t %s_dispatcher;\n" % prefix

    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan_profile", "striped_profile"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    txt += "parasail/parasail_pfunction_t %s_dispatcher;\n" % prefix

    txt += """
/* declare and initialize the pointer to the dispatcher function */
"""
    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    txt += "parasail/parasail_function_t * %s_pointer = %s_dispatcher;\n"%(
                            prefix, prefix)

    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan_profile", "striped_profile"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    txt += "parasail/parasail_pfunction_t * %s_pointer = %s_dispatcher;\n"%(
                            prefix, prefix)

    txt += """
/* dispatcher function implementations */
"""
    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    prefix2 = "parasail/parasail_%s%s%s_%s"%(
                        alg, stats, table, par)
                    base = alg
                    if par == "scan":
                        base += "_scan"
                    params = {
                            "ALG": alg,
                            "BASE": base,
                            "PREFIX": prefix,
                            "PREFIX2": prefix2,
                            "TABLE": table,
                            "STATS": stats,
                            "PAR": par,
                            "WIDTH": width
                    }
                    txt += """
parasail_result_t* %(PREFIX)s_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_avx2_256_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_sse41_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_sse2_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_altivec_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_neon_128_%(WIDTH)s;
    }
    else
#endif
    {
        %(PREFIX)s_pointer = parasail_%(BASE)s;
    }
    return %(PREFIX)s_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}
""" % params

    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan_profile", "striped_profile"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    prefix2 = "parasail/parasail_%s%s%s_%s"%(
                        alg, stats, table, par)
                    params = {
                            "ALG": alg,
                            "PREFIX": prefix,
                            "PREFIX2": prefix2,
                            "TABLE": table,
                            "STATS": stats,
                            "PAR": par,
                            "WIDTH": width
                    }
                    txt += """
parasail_result_t* %(PREFIX)s_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_avx2_256_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_sse41_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_sse2_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_altivec_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_neon_128_%(WIDTH)s;
    }
    else
#endif
    {
        %(PREFIX)s_pointer = NULL;
    }
    return %(PREFIX)s_pointer(profile, s2, s2Len, open, gap);
}
""" % params

    txt += """
/* implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct impl */
"""
    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    params = {
                            "ALG": alg,
                            "PREFIX": prefix,
                            "TABLE": table,
                            "STATS": stats,
                            "PAR": par,
                            "WIDTH": width
                    }
                    txt += """
parasail_result_t* %(PREFIX)s(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return %(PREFIX)s_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}
""" % params

    for table in ["", "_table", "_rowcol", "_trace"]:
        for stats in ["", "_stats"]:
            if 'stats' in stats and 'trace' in table: continue
            for par in ["scan_profile", "striped_profile"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail/parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    params = {
                            "ALG": alg,
                            "PREFIX": prefix,
                            "TABLE": table,
                            "STATS": stats,
                            "PAR": par,
                            "WIDTH": width
                    }
                    txt += """
parasail_result_t* %(PREFIX)s(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return %(PREFIX)s_pointer(profile, s2, s2Len, open, gap);
}
""" % params

    return txt

output_dir = "generated/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for alg in ["nw", "sg", "sw"]:
    output_filename = "%s%s_dispatch.c" % (output_dir, alg)
    writer = open(output_filename, "w")
    writer.write(codegen(alg))
    writer.write("\n")
    writer.close()
