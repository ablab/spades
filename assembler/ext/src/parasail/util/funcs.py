#!/usr/bin/env python

print """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_FUNCTION_TYPE_H_
#define _PARASAIL_FUNCTION_TYPE_H_

#include "parasail/parasail.h"

#ifdef __cplusplus
extern "C" {
#endif
"""


def print_fmt(*args):
    fmt = '{%-36s %-38s %5s %10s %-8s %6s %5s %3s %1s %1s %1s %1s %1s},'
    new_args = [arg for arg in args]
    new_args[0] = '%s,'   % new_args[0]
    new_args[1] = '"%s",' % new_args[1]
    new_args[2] = '"%s",' % new_args[2]
    new_args[3] = '"%s",' % new_args[3]
    new_args[4] = '"%s",' % new_args[4]
    new_args[5] = '"%s",' % new_args[5]
    new_args[6] = '"%s",' % new_args[6]
    new_args[7] = '%d,'   % new_args[7]
    new_args[8] = '%d,'   % new_args[8]
    new_args[9] = '%d,'   % new_args[9]
    new_args[10]= '%d,'   % new_args[10]
    new_args[11]= '%d,'   % new_args[11]
    new_args[12]= '%d'    % new_args[12]
    print fmt % tuple(new_args)

def print_pfmt(*args):
    fmt = '{%-36s %-36s %-38s %5s %10s %-8s %6s %5s %3s %1s %1s %1s %1s %1s},'
    new_args = [arg for arg in args]
    new_args[0] = '%s,'   % new_args[0]
    new_args[1] = '%s,'   % new_args[1]
    new_args[2] = '"%s",' % new_args[2]
    new_args[3] = '"%s",' % new_args[3]
    new_args[4] = '"%s",' % new_args[4]
    new_args[5] = '"%s",' % new_args[5]
    new_args[6] = '"%s",' % new_args[6]
    new_args[7] = '"%s",' % new_args[7]
    new_args[8] = '%d,'   % new_args[8]
    new_args[9] = '%d,'   % new_args[9]
    new_args[10]= '%d,'   % new_args[10]
    new_args[11]= '%d,'   % new_args[11]
    new_args[12]= '%d,'   % new_args[12]
    new_args[13]= '%d'    % new_args[13]
    print fmt % tuple(new_args)

def print_null():
    fmt = '{%s, "%s", "%s", "%s", "%s", "%s", "%s", %d, %d, %d, %d, %d, %d},'
    print fmt[:-1] % ("NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0)

def print_pnull():
    fmt = '{%s, %s, "%s", "%s", "%s", "%s", "%s", "%s", %d, %d, %d, %d, %d, %d},'
    print fmt[:-1] % ("NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0)

isa_to_bits = {
    "sse2"    : 128,
    "sse41"   : 128,
    "avx2"    : 256,
    "altivec" : 128,
    "neon"    : 128,
}

print "static const parasail_function_info_t functions[] = {"

for table in ["", "_table", "_rowcol", "_trace"]:
    for stats in ["", "_stats"]:
        if 'stats' in stats and 'trace' in table: continue
        for alg in ["nw", "sg", "sw"]:
            is_table = 0
            is_rowcol = 0
            is_trace = 0
            if table == "_table":
                is_table = 1
            elif table == "_rowcol":
                is_rowcol = 1
            elif table == "_trace":
                is_trace = 1
            is_stats = 0
            if stats:
                is_stats = 1
            pre = "parasail/parasail_"+alg+stats+table
            print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_rowcol, is_trace, is_stats, 1)
            print_fmt(pre+"_scan", pre+"_scan", alg+stats, "scan", "NA", "32", "32", 1, is_table, is_rowcol, is_trace, is_stats, 0)
            for isa in ["sse2", "sse41", "avx2", "altivec", "neon"]:
                print "#if HAVE_%s" % isa.upper()
                bits = isa_to_bits[isa]
                for par in ["scan", "striped", "diag"]:
                    for width in [64, 32, 16, 8, "sat"]:
                        elem = 0
                        if "sat" == width:
                            elem = bits/8
                        else:
                            elem = bits/width
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits,
                                width, elem, is_table, is_rowcol, is_trace, is_stats, 0)
                # blocked implementations only exist for sw sse41 32 and 16 bit
                if (isa == "sse41" and alg == "sw"
                        and not stats and "trace" not in table):
                    par = "blocked"
                    for width in [32, 16]:
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits, width, bits/width, is_table, is_rowcol, is_trace, is_stats, 0)
                print "#endif"
            # also print the dispatcher function
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8, "sat"]:
                    name = "%s_%s_%s" % (pre, par, width)
                    print_fmt(name, name, alg+stats, par, "disp", "NA", width, -1, is_table, is_rowcol, is_trace, is_stats, 0)

print_null()
print "};"

print
print "static const parasail_pfunction_info_t pfunctions[] = {"

for table in ["", "_table", "_rowcol", "_trace"]:
    for stats in ["", "_stats"]:
        if 'stats' in stats and 'trace' in table: continue
        for alg in ["nw", "sg", "sw"]:
            is_table = 0
            is_rowcol = 0
            is_trace = 0
            if table == "_table":
                is_table = 1
            elif table == "_rowcol":
                is_rowcol = 1
            elif table == "_trace":
                is_trace = 1
            is_stats = 0
            if stats:
                is_stats = 1
            pre = "parasail/parasail_"+alg+stats+table
            for isa in ["sse2", "sse41", "avx2", "altivec", "neon"]:
                print "#if HAVE_%s" % isa.upper()
                bits = isa_to_bits[isa]
                for par in ["scan_profile", "striped_profile"]:
                    for width in [64, 32, 16, 8, "sat"]:
                        elem = 0
                        if "sat" == width:
                            elem = bits/8
                        else:
                            elem = bits/width
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        if 'altivec' in isa:
                            creator = "parasail/parasail_profile_create%s_%s_%s_%s" % (stats, isa, bits, width)
                        elif 'neon' in isa:
                            creator = "parasail/parasail_profile_create%s_%s_%s_%s" % (stats, isa, bits, width)
                        else:
                            creator = "parasail/parasail_profile_create%s_%s_%s_%s" % (stats, isa[:3], bits, width)
                        print_pfmt(name, creator, name, alg+stats, par, isa, bits, width, elem, is_table, is_rowcol, is_trace, is_stats, 0)
                print "#endif"
            # also print the dispatcher function
            for par in ["scan_profile", "striped_profile"]:
                for width in [64, 32, 16, 8, "sat"]:
                    name = "%s_%s_%s" % (pre, par, width)
                    creator = "parasail/parasail_profile_create%s_%s" % (stats, width)
                    print_pfmt(name, creator, name, alg+stats, par, "disp", "NA", width, -1, is_table, is_rowcol, is_trace, is_stats, 0)

print_pnull()
print "};"

print """
#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_FUNCTION_TYPE_H_ */
"""
