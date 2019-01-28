#!/usr/bin/env python

print """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_FUNCTION_GROUP_TABLE_H_
#define _PARASAIL_FUNCTION_GROUP_TABLE_H_

#include "parasail/parasail.h"

typedef struct parasail_function_group {
    const char * name;
    parasail_function_info_t *fs;
} parasail_function_group_t;
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

def print_null():
    fmt = '{%s, "%s", "%s", "%s", "%s", "%s", "%s", %d, %d, %d, %d, %d, %d},'
    print fmt[:-1] % ("NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0)

isa_to_bits = {
    "sse2"    : 128,
    "sse41"   : 128,
    "avx2"    : 256,
    "altivec" : 128,
    "neon"    : 128,
}

for table in ["_table"]:
    for stats in ["", "_stats"]:
        for alg in ["nw", "sg", "sw"]:
            is_table = 0
            is_rowcol = 0
            is_trace = 0
            is_stats = 0
            if table:
                is_table = 1
            if stats:
                is_stats = 1
            pre = "parasail/parasail_"+alg+stats+table
            for isa in ["sse2", "sse41", "avx2", "altivec", "neon"]:
                print "#if HAVE_%s" % isa.upper()
                print "static parasail_function_info_t %s_%s_functions[] = {" % (pre, isa)
                print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_rowcol, is_trace, is_stats, 1)
                bits = isa_to_bits[isa]
                for par in ["scan", "striped", "diag"]:
                    widths = [64, 32, 16, 8]
                    for width in widths:
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits, width, bits/width, is_table, is_rowcol, is_trace, is_stats, 0)
                print_null()
                print "};"
                print 'static parasail_function_group_t %s_%s = {"%s_%s", %s_%s_functions};' % ((pre, isa)*3)
                print "#endif"
            # non-isa-specific functions
            isa = "disp"
            print "static parasail_function_info_t %s_%s_functions[] = {" % (pre, isa)
            print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_rowcol, is_trace, is_stats, 1)
            # also print the dispatcher function
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    name = "%s_%s_%s" % (pre, par, width)
                    print_fmt(name, name, alg+stats, par, "disp", "NA", width, -
1, is_table, is_rowcol, is_trace, is_stats, 0)
            # also print the saturation check function
            for par in ["scan", "striped", "diag"]:
                name = "%s_%s_sat" % (pre, par)
                print_fmt(name, name, alg+stats, par, "sat", "NA", 8, -1, is_table, is_rowcol, is_trace, is_stats, 0)
            print_null()
            print "};"
            print 'static parasail_function_group_t %s_%s = {"%s_%s", %s_%s_functions};' % ((pre, isa)*3)
            # non-vectorized functions
            isa = "serial"
            print "static parasail_function_info_t %s_%s_functions[] = {" % (pre, isa)
            print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_rowcol, is_trace, is_stats, 1)
            print_fmt(pre+"_scan", pre+"_scan", alg+stats, "scan", "NA", "32", "32", 1, is_table, is_rowcol, is_trace, is_stats, 0)
            print_null()
            print "};"
            print 'static parasail_function_group_t %s_%s = {"%s_%s", %s_%s_functions};' % ((pre, isa)*3)


print """
#endif /* _PARASAIL_FUNCTION_GROUP_TABLE_H_ */
"""
