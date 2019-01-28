/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_FUNCTION_GROUP_H_
#define _PARASAIL_FUNCTION_GROUP_H_

#include "parasail/parasail.h"

typedef struct parasail_function_group {
    const char * name;
    parasail_function_info_t *fs;
} parasail_function_group_t;

#if HAVE_SSE2
static parasail_function_info_t parasail_nw_sse2_functions[] = {
{parasail_nw,                         "parasail/parasail_nw",                         "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_nw_scan,                    "parasail/parasail_nw_scan",                    "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse2_128_64,        "parasail/parasail_nw_scan_sse2_128_64",        "nw",    "scan", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse2_128_32,        "parasail/parasail_nw_scan_sse2_128_32",        "nw",    "scan", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse2_128_16,        "parasail/parasail_nw_scan_sse2_128_16",        "nw",    "scan", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse2_128_8,         "parasail/parasail_nw_scan_sse2_128_8",         "nw",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse2_128_64,     "parasail/parasail_nw_striped_sse2_128_64",     "nw", "striped", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse2_128_32,     "parasail/parasail_nw_striped_sse2_128_32",     "nw", "striped", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse2_128_16,     "parasail/parasail_nw_striped_sse2_128_16",     "nw", "striped", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse2_128_8,      "parasail/parasail_nw_striped_sse2_128_8",      "nw", "striped", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse2_128_64,        "parasail/parasail_nw_diag_sse2_128_64",        "nw",    "diag", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse2_128_32,        "parasail/parasail_nw_diag_sse2_128_32",        "nw",    "diag", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse2_128_16,        "parasail/parasail_nw_diag_sse2_128_16",        "nw",    "diag", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse2_128_8,         "parasail/parasail_nw_diag_sse2_128_8",         "nw",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_sse2 = {"parasail/parasail_nw_sse2", parasail_nw_sse2_functions};
#endif
#if HAVE_SSE41
static parasail_function_info_t parasail_nw_sse41_functions[] = {
{parasail_nw,                         "parasail/parasail_nw",                         "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_nw_scan,                    "parasail/parasail_nw_scan",                    "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse41_128_64,       "parasail/parasail_nw_scan_sse41_128_64",       "nw",    "scan", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse41_128_32,       "parasail/parasail_nw_scan_sse41_128_32",       "nw",    "scan", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse41_128_16,       "parasail/parasail_nw_scan_sse41_128_16",       "nw",    "scan", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_scan_sse41_128_8,        "parasail/parasail_nw_scan_sse41_128_8",        "nw",    "scan", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse41_128_64,    "parasail/parasail_nw_striped_sse41_128_64",    "nw", "striped", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse41_128_32,    "parasail/parasail_nw_striped_sse41_128_32",    "nw", "striped", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse41_128_16,    "parasail/parasail_nw_striped_sse41_128_16",    "nw", "striped", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_striped_sse41_128_8,     "parasail/parasail_nw_striped_sse41_128_8",     "nw", "striped", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse41_128_64,       "parasail/parasail_nw_diag_sse41_128_64",       "nw",    "diag", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse41_128_32,       "parasail/parasail_nw_diag_sse41_128_32",       "nw",    "diag", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse41_128_16,       "parasail/parasail_nw_diag_sse41_128_16",       "nw",    "diag", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_diag_sse41_128_8,        "parasail/parasail_nw_diag_sse41_128_8",        "nw",    "diag", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_sse41 = {"parasail/parasail_nw_sse41", parasail_nw_sse41_functions};
#endif
#if HAVE_AVX2
static parasail_function_info_t parasail_nw_avx2_functions[] = {
{parasail_nw,                         "parasail/parasail_nw",                         "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_nw_scan,                    "parasail/parasail_nw_scan",                    "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_nw_scan_avx2_256_64,        "parasail/parasail_nw_scan_avx2_256_64",        "nw",    "scan", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_nw_scan_avx2_256_32,        "parasail/parasail_nw_scan_avx2_256_32",        "nw",    "scan", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_nw_scan_avx2_256_16,        "parasail/parasail_nw_scan_avx2_256_16",        "nw",    "scan", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_nw_scan_avx2_256_8,         "parasail/parasail_nw_scan_avx2_256_8",         "nw",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{parasail_nw_striped_avx2_256_64,     "parasail/parasail_nw_striped_avx2_256_64",     "nw", "striped", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_nw_striped_avx2_256_32,     "parasail/parasail_nw_striped_avx2_256_32",     "nw", "striped", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_nw_striped_avx2_256_16,     "parasail/parasail_nw_striped_avx2_256_16",     "nw", "striped", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_nw_striped_avx2_256_8,      "parasail/parasail_nw_striped_avx2_256_8",      "nw", "striped", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{parasail_nw_diag_avx2_256_64,        "parasail/parasail_nw_diag_avx2_256_64",        "nw",    "diag", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_nw_diag_avx2_256_32,        "parasail/parasail_nw_diag_avx2_256_32",        "nw",    "diag", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_nw_diag_avx2_256_16,        "parasail/parasail_nw_diag_avx2_256_16",        "nw",    "diag", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_nw_diag_avx2_256_8,         "parasail/parasail_nw_diag_avx2_256_8",         "nw",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_avx2 = {"parasail/parasail_nw_avx2", parasail_nw_avx2_functions};
#endif
#if HAVE_ALTIVEC
static parasail_function_info_t parasail_nw_altivec_functions[] = {
{parasail_nw,                         "parasail/parasail_nw",                         "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_nw_scan,                    "parasail/parasail_nw_scan",                    "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_nw_scan_altivec_128_64,     "parasail/parasail_nw_scan_altivec_128_64",     "nw",    "scan", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_scan_altivec_128_32,     "parasail/parasail_nw_scan_altivec_128_32",     "nw",    "scan", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_scan_altivec_128_16,     "parasail/parasail_nw_scan_altivec_128_16",     "nw",    "scan", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_scan_altivec_128_8,      "parasail/parasail_nw_scan_altivec_128_8",      "nw",    "scan", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_striped_altivec_128_64,  "parasail/parasail_nw_striped_altivec_128_64",  "nw", "striped", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_striped_altivec_128_32,  "parasail/parasail_nw_striped_altivec_128_32",  "nw", "striped", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_striped_altivec_128_16,  "parasail/parasail_nw_striped_altivec_128_16",  "nw", "striped", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_striped_altivec_128_8,   "parasail/parasail_nw_striped_altivec_128_8",   "nw", "striped", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_diag_altivec_128_64,     "parasail/parasail_nw_diag_altivec_128_64",     "nw",    "diag", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_diag_altivec_128_32,     "parasail/parasail_nw_diag_altivec_128_32",     "nw",    "diag", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_diag_altivec_128_16,     "parasail/parasail_nw_diag_altivec_128_16",     "nw",    "diag", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_diag_altivec_128_8,      "parasail/parasail_nw_diag_altivec_128_8",      "nw",    "diag", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_altivec = {"parasail/parasail_nw_altivec", parasail_nw_altivec_functions};
#endif
#if HAVE_NEON
static parasail_function_info_t parasail_nw_neon_functions[] = {
{parasail_nw,                         "parasail/parasail_nw",                         "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_nw_scan,                    "parasail/parasail_nw_scan",                    "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_nw_scan_neon_128_64,        "parasail/parasail_nw_scan_neon_128_64",        "nw",    "scan", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_scan_neon_128_32,        "parasail/parasail_nw_scan_neon_128_32",        "nw",    "scan", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_scan_neon_128_16,        "parasail/parasail_nw_scan_neon_128_16",        "nw",    "scan", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_scan_neon_128_8,         "parasail/parasail_nw_scan_neon_128_8",         "nw",    "scan", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_striped_neon_128_64,     "parasail/parasail_nw_striped_neon_128_64",     "nw", "striped", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_striped_neon_128_32,     "parasail/parasail_nw_striped_neon_128_32",     "nw", "striped", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_striped_neon_128_16,     "parasail/parasail_nw_striped_neon_128_16",     "nw", "striped", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_striped_neon_128_8,      "parasail/parasail_nw_striped_neon_128_8",      "nw", "striped", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_nw_diag_neon_128_64,        "parasail/parasail_nw_diag_neon_128_64",        "nw",    "diag", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_nw_diag_neon_128_32,        "parasail/parasail_nw_diag_neon_128_32",        "nw",    "diag", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_nw_diag_neon_128_16,        "parasail/parasail_nw_diag_neon_128_16",        "nw",    "diag", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_nw_diag_neon_128_8,         "parasail/parasail_nw_diag_neon_128_8",         "nw",    "diag", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_neon = {"parasail/parasail_nw_neon", parasail_nw_neon_functions};
#endif
static parasail_function_info_t parasail_nw_disp_functions[] = {
{parasail_nw,                         "parasail/parasail_nw",                         "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_nw_scan,                    "parasail/parasail_nw_scan",                    "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_nw_scan_64,                 "parasail/parasail_nw_scan_64",                 "nw",    "scan", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_nw_scan_32,                 "parasail/parasail_nw_scan_32",                 "nw",    "scan", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_nw_scan_16,                 "parasail/parasail_nw_scan_16",                 "nw",    "scan", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_nw_scan_8,                  "parasail/parasail_nw_scan_8",                  "nw",    "scan", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_nw_striped_64,              "parasail/parasail_nw_striped_64",              "nw", "striped", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_nw_striped_32,              "parasail/parasail_nw_striped_32",              "nw", "striped", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_nw_striped_16,              "parasail/parasail_nw_striped_16",              "nw", "striped", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_nw_striped_8,               "parasail/parasail_nw_striped_8",               "nw", "striped", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_nw_diag_64,                 "parasail/parasail_nw_diag_64",                 "nw",    "diag", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_nw_diag_32,                 "parasail/parasail_nw_diag_32",                 "nw",    "diag", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_nw_diag_16,                 "parasail/parasail_nw_diag_16",                 "nw",    "diag", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_nw_diag_8,                  "parasail/parasail_nw_diag_8",                  "nw",    "diag", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_nw_scan_sat,                "parasail/parasail_nw_scan_sat",                "nw",    "scan", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_nw_striped_sat,             "parasail/parasail_nw_striped_sat",             "nw", "striped", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_nw_diag_sat,                "parasail/parasail_nw_diag_sat",                "nw",    "diag", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_disp = {"parasail/parasail_nw_disp", parasail_nw_disp_functions};
#if HAVE_SSE2
static parasail_function_info_t parasail_sg_sse2_functions[] = {
{parasail_sg,                         "parasail/parasail_sg",                         "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sg_scan,                    "parasail/parasail_sg_scan",                    "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse2_128_64,        "parasail/parasail_sg_scan_sse2_128_64",        "sg",    "scan", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse2_128_32,        "parasail/parasail_sg_scan_sse2_128_32",        "sg",    "scan", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse2_128_16,        "parasail/parasail_sg_scan_sse2_128_16",        "sg",    "scan", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse2_128_8,         "parasail/parasail_sg_scan_sse2_128_8",         "sg",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse2_128_64,     "parasail/parasail_sg_striped_sse2_128_64",     "sg", "striped", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse2_128_32,     "parasail/parasail_sg_striped_sse2_128_32",     "sg", "striped", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse2_128_16,     "parasail/parasail_sg_striped_sse2_128_16",     "sg", "striped", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse2_128_8,      "parasail/parasail_sg_striped_sse2_128_8",      "sg", "striped", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse2_128_64,        "parasail/parasail_sg_diag_sse2_128_64",        "sg",    "diag", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse2_128_32,        "parasail/parasail_sg_diag_sse2_128_32",        "sg",    "diag", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse2_128_16,        "parasail/parasail_sg_diag_sse2_128_16",        "sg",    "diag", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse2_128_8,         "parasail/parasail_sg_diag_sse2_128_8",         "sg",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_sse2 = {"parasail/parasail_sg_sse2", parasail_sg_sse2_functions};
#endif
#if HAVE_SSE41
static parasail_function_info_t parasail_sg_sse41_functions[] = {
{parasail_sg,                         "parasail/parasail_sg",                         "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sg_scan,                    "parasail/parasail_sg_scan",                    "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse41_128_64,       "parasail/parasail_sg_scan_sse41_128_64",       "sg",    "scan", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse41_128_32,       "parasail/parasail_sg_scan_sse41_128_32",       "sg",    "scan", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse41_128_16,       "parasail/parasail_sg_scan_sse41_128_16",       "sg",    "scan", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_scan_sse41_128_8,        "parasail/parasail_sg_scan_sse41_128_8",        "sg",    "scan", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse41_128_64,    "parasail/parasail_sg_striped_sse41_128_64",    "sg", "striped", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse41_128_32,    "parasail/parasail_sg_striped_sse41_128_32",    "sg", "striped", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse41_128_16,    "parasail/parasail_sg_striped_sse41_128_16",    "sg", "striped", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_striped_sse41_128_8,     "parasail/parasail_sg_striped_sse41_128_8",     "sg", "striped", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse41_128_64,       "parasail/parasail_sg_diag_sse41_128_64",       "sg",    "diag", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse41_128_32,       "parasail/parasail_sg_diag_sse41_128_32",       "sg",    "diag", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse41_128_16,       "parasail/parasail_sg_diag_sse41_128_16",       "sg",    "diag", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_diag_sse41_128_8,        "parasail/parasail_sg_diag_sse41_128_8",        "sg",    "diag", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_sse41 = {"parasail/parasail_sg_sse41", parasail_sg_sse41_functions};
#endif
#if HAVE_AVX2
static parasail_function_info_t parasail_sg_avx2_functions[] = {
{parasail_sg,                         "parasail/parasail_sg",                         "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sg_scan,                    "parasail/parasail_sg_scan",                    "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sg_scan_avx2_256_64,        "parasail/parasail_sg_scan_avx2_256_64",        "sg",    "scan", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_sg_scan_avx2_256_32,        "parasail/parasail_sg_scan_avx2_256_32",        "sg",    "scan", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_sg_scan_avx2_256_16,        "parasail/parasail_sg_scan_avx2_256_16",        "sg",    "scan", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_sg_scan_avx2_256_8,         "parasail/parasail_sg_scan_avx2_256_8",         "sg",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{parasail_sg_striped_avx2_256_64,     "parasail/parasail_sg_striped_avx2_256_64",     "sg", "striped", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_sg_striped_avx2_256_32,     "parasail/parasail_sg_striped_avx2_256_32",     "sg", "striped", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_sg_striped_avx2_256_16,     "parasail/parasail_sg_striped_avx2_256_16",     "sg", "striped", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_sg_striped_avx2_256_8,      "parasail/parasail_sg_striped_avx2_256_8",      "sg", "striped", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{parasail_sg_diag_avx2_256_64,        "parasail/parasail_sg_diag_avx2_256_64",        "sg",    "diag", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_sg_diag_avx2_256_32,        "parasail/parasail_sg_diag_avx2_256_32",        "sg",    "diag", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_sg_diag_avx2_256_16,        "parasail/parasail_sg_diag_avx2_256_16",        "sg",    "diag", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_sg_diag_avx2_256_8,         "parasail/parasail_sg_diag_avx2_256_8",         "sg",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_avx2 = {"parasail/parasail_sg_avx2", parasail_sg_avx2_functions};
#endif
#if HAVE_ALTIVEC
static parasail_function_info_t parasail_sg_altivec_functions[] = {
{parasail_sg,                         "parasail/parasail_sg",                         "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sg_scan,                    "parasail/parasail_sg_scan",                    "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sg_scan_altivec_128_64,     "parasail/parasail_sg_scan_altivec_128_64",     "sg",    "scan", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_scan_altivec_128_32,     "parasail/parasail_sg_scan_altivec_128_32",     "sg",    "scan", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_scan_altivec_128_16,     "parasail/parasail_sg_scan_altivec_128_16",     "sg",    "scan", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_scan_altivec_128_8,      "parasail/parasail_sg_scan_altivec_128_8",      "sg",    "scan", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_striped_altivec_128_64,  "parasail/parasail_sg_striped_altivec_128_64",  "sg", "striped", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_striped_altivec_128_32,  "parasail/parasail_sg_striped_altivec_128_32",  "sg", "striped", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_striped_altivec_128_16,  "parasail/parasail_sg_striped_altivec_128_16",  "sg", "striped", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_striped_altivec_128_8,   "parasail/parasail_sg_striped_altivec_128_8",   "sg", "striped", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_diag_altivec_128_64,     "parasail/parasail_sg_diag_altivec_128_64",     "sg",    "diag", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_diag_altivec_128_32,     "parasail/parasail_sg_diag_altivec_128_32",     "sg",    "diag", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_diag_altivec_128_16,     "parasail/parasail_sg_diag_altivec_128_16",     "sg",    "diag", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_diag_altivec_128_8,      "parasail/parasail_sg_diag_altivec_128_8",      "sg",    "diag", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_altivec = {"parasail/parasail_sg_altivec", parasail_sg_altivec_functions};
#endif
#if HAVE_NEON
static parasail_function_info_t parasail_sg_neon_functions[] = {
{parasail_sg,                         "parasail/parasail_sg",                         "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sg_scan,                    "parasail/parasail_sg_scan",                    "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sg_scan_neon_128_64,        "parasail/parasail_sg_scan_neon_128_64",        "sg",    "scan", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_scan_neon_128_32,        "parasail/parasail_sg_scan_neon_128_32",        "sg",    "scan", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_scan_neon_128_16,        "parasail/parasail_sg_scan_neon_128_16",        "sg",    "scan", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_scan_neon_128_8,         "parasail/parasail_sg_scan_neon_128_8",         "sg",    "scan", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_striped_neon_128_64,     "parasail/parasail_sg_striped_neon_128_64",     "sg", "striped", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_striped_neon_128_32,     "parasail/parasail_sg_striped_neon_128_32",     "sg", "striped", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_striped_neon_128_16,     "parasail/parasail_sg_striped_neon_128_16",     "sg", "striped", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_striped_neon_128_8,      "parasail/parasail_sg_striped_neon_128_8",      "sg", "striped", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sg_diag_neon_128_64,        "parasail/parasail_sg_diag_neon_128_64",        "sg",    "diag", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sg_diag_neon_128_32,        "parasail/parasail_sg_diag_neon_128_32",        "sg",    "diag", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sg_diag_neon_128_16,        "parasail/parasail_sg_diag_neon_128_16",        "sg",    "diag", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sg_diag_neon_128_8,         "parasail/parasail_sg_diag_neon_128_8",         "sg",    "diag", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_neon = {"parasail/parasail_sg_neon", parasail_sg_neon_functions};
#endif
static parasail_function_info_t parasail_sg_disp_functions[] = {
{parasail_sg,                         "parasail/parasail_sg",                         "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sg_scan,                    "parasail/parasail_sg_scan",                    "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sg_scan_64,                 "parasail/parasail_sg_scan_64",                 "sg",    "scan", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_sg_scan_32,                 "parasail/parasail_sg_scan_32",                 "sg",    "scan", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_sg_scan_16,                 "parasail/parasail_sg_scan_16",                 "sg",    "scan", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_sg_scan_8,                  "parasail/parasail_sg_scan_8",                  "sg",    "scan", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sg_striped_64,              "parasail/parasail_sg_striped_64",              "sg", "striped", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_sg_striped_32,              "parasail/parasail_sg_striped_32",              "sg", "striped", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_sg_striped_16,              "parasail/parasail_sg_striped_16",              "sg", "striped", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_sg_striped_8,               "parasail/parasail_sg_striped_8",               "sg", "striped", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sg_diag_64,                 "parasail/parasail_sg_diag_64",                 "sg",    "diag", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_sg_diag_32,                 "parasail/parasail_sg_diag_32",                 "sg",    "diag", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_sg_diag_16,                 "parasail/parasail_sg_diag_16",                 "sg",    "diag", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_sg_diag_8,                  "parasail/parasail_sg_diag_8",                  "sg",    "diag", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sg_scan_sat,                "parasail/parasail_sg_scan_sat",                "sg",    "scan", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sg_striped_sat,             "parasail/parasail_sg_striped_sat",             "sg", "striped", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sg_diag_sat,                "parasail/parasail_sg_diag_sat",                "sg",    "diag", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_disp = {"parasail/parasail_sg_disp", parasail_sg_disp_functions};
#if HAVE_SSE2
static parasail_function_info_t parasail_sw_sse2_functions[] = {
{parasail_sw,                         "parasail/parasail_sw",                         "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sw_scan,                    "parasail/parasail_sw_scan",                    "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse2_128_64,        "parasail/parasail_sw_scan_sse2_128_64",        "sw",    "scan", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse2_128_32,        "parasail/parasail_sw_scan_sse2_128_32",        "sw",    "scan", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse2_128_16,        "parasail/parasail_sw_scan_sse2_128_16",        "sw",    "scan", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse2_128_8,         "parasail/parasail_sw_scan_sse2_128_8",         "sw",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse2_128_64,     "parasail/parasail_sw_striped_sse2_128_64",     "sw", "striped", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse2_128_32,     "parasail/parasail_sw_striped_sse2_128_32",     "sw", "striped", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse2_128_16,     "parasail/parasail_sw_striped_sse2_128_16",     "sw", "striped", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse2_128_8,      "parasail/parasail_sw_striped_sse2_128_8",      "sw", "striped", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse2_128_64,        "parasail/parasail_sw_diag_sse2_128_64",        "sw",    "diag", "sse2",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse2_128_32,        "parasail/parasail_sw_diag_sse2_128_32",        "sw",    "diag", "sse2",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse2_128_16,        "parasail/parasail_sw_diag_sse2_128_16",        "sw",    "diag", "sse2",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse2_128_8,         "parasail/parasail_sw_diag_sse2_128_8",         "sw",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_sse2 = {"parasail/parasail_sw_sse2", parasail_sw_sse2_functions};
#endif
#if HAVE_SSE41
static parasail_function_info_t parasail_sw_sse41_functions[] = {
{parasail_sw,                         "parasail/parasail_sw",                         "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sw_scan,                    "parasail/parasail_sw_scan",                    "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse41_128_64,       "parasail/parasail_sw_scan_sse41_128_64",       "sw",    "scan", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse41_128_32,       "parasail/parasail_sw_scan_sse41_128_32",       "sw",    "scan", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse41_128_16,       "parasail/parasail_sw_scan_sse41_128_16",       "sw",    "scan", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_scan_sse41_128_8,        "parasail/parasail_sw_scan_sse41_128_8",        "sw",    "scan", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse41_128_64,    "parasail/parasail_sw_striped_sse41_128_64",    "sw", "striped", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse41_128_32,    "parasail/parasail_sw_striped_sse41_128_32",    "sw", "striped", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse41_128_16,    "parasail/parasail_sw_striped_sse41_128_16",    "sw", "striped", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_striped_sse41_128_8,     "parasail/parasail_sw_striped_sse41_128_8",     "sw", "striped", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse41_128_64,       "parasail/parasail_sw_diag_sse41_128_64",       "sw",    "diag", "sse41", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse41_128_32,       "parasail/parasail_sw_diag_sse41_128_32",       "sw",    "diag", "sse41", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse41_128_16,       "parasail/parasail_sw_diag_sse41_128_16",       "sw",    "diag", "sse41", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_diag_sse41_128_8,        "parasail/parasail_sw_diag_sse41_128_8",        "sw",    "diag", "sse41", "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_sse41 = {"parasail/parasail_sw_sse41", parasail_sw_sse41_functions};
#endif
#if HAVE_AVX2
static parasail_function_info_t parasail_sw_avx2_functions[] = {
{parasail_sw,                         "parasail/parasail_sw",                         "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sw_scan,                    "parasail/parasail_sw_scan",                    "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sw_scan_avx2_256_64,        "parasail/parasail_sw_scan_avx2_256_64",        "sw",    "scan", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_sw_scan_avx2_256_32,        "parasail/parasail_sw_scan_avx2_256_32",        "sw",    "scan", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_sw_scan_avx2_256_16,        "parasail/parasail_sw_scan_avx2_256_16",        "sw",    "scan", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_sw_scan_avx2_256_8,         "parasail/parasail_sw_scan_avx2_256_8",         "sw",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{parasail_sw_striped_avx2_256_64,     "parasail/parasail_sw_striped_avx2_256_64",     "sw", "striped", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_sw_striped_avx2_256_32,     "parasail/parasail_sw_striped_avx2_256_32",     "sw", "striped", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_sw_striped_avx2_256_16,     "parasail/parasail_sw_striped_avx2_256_16",     "sw", "striped", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_sw_striped_avx2_256_8,      "parasail/parasail_sw_striped_avx2_256_8",      "sw", "striped", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{parasail_sw_diag_avx2_256_64,        "parasail/parasail_sw_diag_avx2_256_64",        "sw",    "diag", "avx2",  "256", "64",  4, 0, 0, 0, 0, 0},
{parasail_sw_diag_avx2_256_32,        "parasail/parasail_sw_diag_avx2_256_32",        "sw",    "diag", "avx2",  "256", "32",  8, 0, 0, 0, 0, 0},
{parasail_sw_diag_avx2_256_16,        "parasail/parasail_sw_diag_avx2_256_16",        "sw",    "diag", "avx2",  "256", "16", 16, 0, 0, 0, 0, 0},
{parasail_sw_diag_avx2_256_8,         "parasail/parasail_sw_diag_avx2_256_8",         "sw",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_avx2 = {"parasail/parasail_sw_avx2", parasail_sw_avx2_functions};
#endif
#if HAVE_ALTIVEC
static parasail_function_info_t parasail_sw_altivec_functions[] = {
{parasail_sw,                         "parasail/parasail_sw",                         "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sw_scan,                    "parasail/parasail_sw_scan",                    "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sw_scan_altivec_128_64,     "parasail/parasail_sw_scan_altivec_128_64",     "sw",    "scan", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_scan_altivec_128_32,     "parasail/parasail_sw_scan_altivec_128_32",     "sw",    "scan", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_scan_altivec_128_16,     "parasail/parasail_sw_scan_altivec_128_16",     "sw",    "scan", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_scan_altivec_128_8,      "parasail/parasail_sw_scan_altivec_128_8",      "sw",    "scan", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_striped_altivec_128_64,  "parasail/parasail_sw_striped_altivec_128_64",  "sw", "striped", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_striped_altivec_128_32,  "parasail/parasail_sw_striped_altivec_128_32",  "sw", "striped", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_striped_altivec_128_16,  "parasail/parasail_sw_striped_altivec_128_16",  "sw", "striped", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_striped_altivec_128_8,   "parasail/parasail_sw_striped_altivec_128_8",   "sw", "striped", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_diag_altivec_128_64,     "parasail/parasail_sw_diag_altivec_128_64",     "sw",    "diag", "altivec", "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_diag_altivec_128_32,     "parasail/parasail_sw_diag_altivec_128_32",     "sw",    "diag", "altivec", "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_diag_altivec_128_16,     "parasail/parasail_sw_diag_altivec_128_16",     "sw",    "diag", "altivec", "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_diag_altivec_128_8,      "parasail/parasail_sw_diag_altivec_128_8",      "sw",    "diag", "altivec", "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_altivec = {"parasail/parasail_sw_altivec", parasail_sw_altivec_functions};
#endif
#if HAVE_NEON
static parasail_function_info_t parasail_sw_neon_functions[] = {
{parasail_sw,                         "parasail/parasail_sw",                         "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sw_scan,                    "parasail/parasail_sw_scan",                    "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sw_scan_neon_128_64,        "parasail/parasail_sw_scan_neon_128_64",        "sw",    "scan", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_scan_neon_128_32,        "parasail/parasail_sw_scan_neon_128_32",        "sw",    "scan", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_scan_neon_128_16,        "parasail/parasail_sw_scan_neon_128_16",        "sw",    "scan", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_scan_neon_128_8,         "parasail/parasail_sw_scan_neon_128_8",         "sw",    "scan", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_striped_neon_128_64,     "parasail/parasail_sw_striped_neon_128_64",     "sw", "striped", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_striped_neon_128_32,     "parasail/parasail_sw_striped_neon_128_32",     "sw", "striped", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_striped_neon_128_16,     "parasail/parasail_sw_striped_neon_128_16",     "sw", "striped", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_striped_neon_128_8,      "parasail/parasail_sw_striped_neon_128_8",      "sw", "striped", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{parasail_sw_diag_neon_128_64,        "parasail/parasail_sw_diag_neon_128_64",        "sw",    "diag", "neon",  "128", "64",  2, 0, 0, 0, 0, 0},
{parasail_sw_diag_neon_128_32,        "parasail/parasail_sw_diag_neon_128_32",        "sw",    "diag", "neon",  "128", "32",  4, 0, 0, 0, 0, 0},
{parasail_sw_diag_neon_128_16,        "parasail/parasail_sw_diag_neon_128_16",        "sw",    "diag", "neon",  "128", "16",  8, 0, 0, 0, 0, 0},
{parasail_sw_diag_neon_128_8,         "parasail/parasail_sw_diag_neon_128_8",         "sw",    "diag", "neon",  "128",  "8", 16, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_neon = {"parasail/parasail_sw_neon", parasail_sw_neon_functions};
#endif
static parasail_function_info_t parasail_sw_disp_functions[] = {
{parasail_sw,                         "parasail/parasail_sw",                         "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 0, 1},
{parasail_sw_scan,                    "parasail/parasail_sw_scan",                    "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 0, 0},
{parasail_sw_scan_64,                 "parasail/parasail_sw_scan_64",                 "sw",    "scan", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_sw_scan_32,                 "parasail/parasail_sw_scan_32",                 "sw",    "scan", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_sw_scan_16,                 "parasail/parasail_sw_scan_16",                 "sw",    "scan", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_sw_scan_8,                  "parasail/parasail_sw_scan_8",                  "sw",    "scan", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sw_striped_64,              "parasail/parasail_sw_striped_64",              "sw", "striped", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_sw_striped_32,              "parasail/parasail_sw_striped_32",              "sw", "striped", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_sw_striped_16,              "parasail/parasail_sw_striped_16",              "sw", "striped", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_sw_striped_8,               "parasail/parasail_sw_striped_8",               "sw", "striped", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sw_diag_64,                 "parasail/parasail_sw_diag_64",                 "sw",    "diag", "disp",   "NA", "64", -1, 0, 0, 0, 0, 0},
{parasail_sw_diag_32,                 "parasail/parasail_sw_diag_32",                 "sw",    "diag", "disp",   "NA", "32", -1, 0, 0, 0, 0, 0},
{parasail_sw_diag_16,                 "parasail/parasail_sw_diag_16",                 "sw",    "diag", "disp",   "NA", "16", -1, 0, 0, 0, 0, 0},
{parasail_sw_diag_8,                  "parasail/parasail_sw_diag_8",                  "sw",    "diag", "disp",   "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sw_scan_sat,                "parasail/parasail_sw_scan_sat",                "sw",    "scan", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sw_striped_sat,             "parasail/parasail_sw_striped_sat",             "sw", "striped", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{parasail_sw_diag_sat,                "parasail/parasail_sw_diag_sat",                "sw",    "diag", "sat",    "NA",  "8", -1, 0, 0, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_disp = {"parasail/parasail_sw_disp", parasail_sw_disp_functions};
#if HAVE_SSE2
static parasail_function_info_t parasail_nw_stats_sse2_functions[] = {
{parasail_nw_stats,                   "parasail/parasail_nw_stats",                   "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_nw_stats_scan,              "parasail/parasail_nw_stats_scan",              "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse2_128_64,  "parasail/parasail_nw_stats_scan_sse2_128_64",  "nw_stats",    "scan", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse2_128_32,  "parasail/parasail_nw_stats_scan_sse2_128_32",  "nw_stats",    "scan", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse2_128_16,  "parasail/parasail_nw_stats_scan_sse2_128_16",  "nw_stats",    "scan", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse2_128_8,   "parasail/parasail_nw_stats_scan_sse2_128_8",   "nw_stats",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse2_128_64, "parasail/parasail_nw_stats_striped_sse2_128_64", "nw_stats", "striped", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse2_128_32, "parasail/parasail_nw_stats_striped_sse2_128_32", "nw_stats", "striped", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse2_128_16, "parasail/parasail_nw_stats_striped_sse2_128_16", "nw_stats", "striped", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse2_128_8, "parasail/parasail_nw_stats_striped_sse2_128_8", "nw_stats", "striped", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse2_128_64,  "parasail/parasail_nw_stats_diag_sse2_128_64",  "nw_stats",    "diag", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse2_128_32,  "parasail/parasail_nw_stats_diag_sse2_128_32",  "nw_stats",    "diag", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse2_128_16,  "parasail/parasail_nw_stats_diag_sse2_128_16",  "nw_stats",    "diag", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse2_128_8,   "parasail/parasail_nw_stats_diag_sse2_128_8",   "nw_stats",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_stats_sse2 = {"parasail/parasail_nw_stats_sse2", parasail_nw_stats_sse2_functions};
#endif
#if HAVE_SSE41
static parasail_function_info_t parasail_nw_stats_sse41_functions[] = {
{parasail_nw_stats,                   "parasail/parasail_nw_stats",                   "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_nw_stats_scan,              "parasail/parasail_nw_stats_scan",              "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse41_128_64, "parasail/parasail_nw_stats_scan_sse41_128_64", "nw_stats",    "scan", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse41_128_32, "parasail/parasail_nw_stats_scan_sse41_128_32", "nw_stats",    "scan", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse41_128_16, "parasail/parasail_nw_stats_scan_sse41_128_16", "nw_stats",    "scan", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sse41_128_8,  "parasail/parasail_nw_stats_scan_sse41_128_8",  "nw_stats",    "scan", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse41_128_64, "parasail/parasail_nw_stats_striped_sse41_128_64", "nw_stats", "striped", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse41_128_32, "parasail/parasail_nw_stats_striped_sse41_128_32", "nw_stats", "striped", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse41_128_16, "parasail/parasail_nw_stats_striped_sse41_128_16", "nw_stats", "striped", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sse41_128_8, "parasail/parasail_nw_stats_striped_sse41_128_8", "nw_stats", "striped", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse41_128_64, "parasail/parasail_nw_stats_diag_sse41_128_64", "nw_stats",    "diag", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse41_128_32, "parasail/parasail_nw_stats_diag_sse41_128_32", "nw_stats",    "diag", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse41_128_16, "parasail/parasail_nw_stats_diag_sse41_128_16", "nw_stats",    "diag", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sse41_128_8,  "parasail/parasail_nw_stats_diag_sse41_128_8",  "nw_stats",    "diag", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_stats_sse41 = {"parasail/parasail_nw_stats_sse41", parasail_nw_stats_sse41_functions};
#endif
#if HAVE_AVX2
static parasail_function_info_t parasail_nw_stats_avx2_functions[] = {
{parasail_nw_stats,                   "parasail/parasail_nw_stats",                   "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_nw_stats_scan,              "parasail/parasail_nw_stats_scan",              "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_avx2_256_64,  "parasail/parasail_nw_stats_scan_avx2_256_64",  "nw_stats",    "scan", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_avx2_256_32,  "parasail/parasail_nw_stats_scan_avx2_256_32",  "nw_stats",    "scan", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_avx2_256_16,  "parasail/parasail_nw_stats_scan_avx2_256_16",  "nw_stats",    "scan", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_avx2_256_8,   "parasail/parasail_nw_stats_scan_avx2_256_8",   "nw_stats",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_avx2_256_64, "parasail/parasail_nw_stats_striped_avx2_256_64", "nw_stats", "striped", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_avx2_256_32, "parasail/parasail_nw_stats_striped_avx2_256_32", "nw_stats", "striped", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_avx2_256_16, "parasail/parasail_nw_stats_striped_avx2_256_16", "nw_stats", "striped", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_avx2_256_8, "parasail/parasail_nw_stats_striped_avx2_256_8", "nw_stats", "striped", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_avx2_256_64,  "parasail/parasail_nw_stats_diag_avx2_256_64",  "nw_stats",    "diag", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_avx2_256_32,  "parasail/parasail_nw_stats_diag_avx2_256_32",  "nw_stats",    "diag", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_avx2_256_16,  "parasail/parasail_nw_stats_diag_avx2_256_16",  "nw_stats",    "diag", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_avx2_256_8,   "parasail/parasail_nw_stats_diag_avx2_256_8",   "nw_stats",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_stats_avx2 = {"parasail/parasail_nw_stats_avx2", parasail_nw_stats_avx2_functions};
#endif
#if HAVE_ALTIVEC
static parasail_function_info_t parasail_nw_stats_altivec_functions[] = {
{parasail_nw_stats,                   "parasail/parasail_nw_stats",                   "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_nw_stats_scan,              "parasail/parasail_nw_stats_scan",              "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_altivec_128_64, "parasail/parasail_nw_stats_scan_altivec_128_64", "nw_stats",    "scan", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_altivec_128_32, "parasail/parasail_nw_stats_scan_altivec_128_32", "nw_stats",    "scan", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_altivec_128_16, "parasail/parasail_nw_stats_scan_altivec_128_16", "nw_stats",    "scan", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_altivec_128_8, "parasail/parasail_nw_stats_scan_altivec_128_8", "nw_stats",    "scan", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_altivec_128_64, "parasail/parasail_nw_stats_striped_altivec_128_64", "nw_stats", "striped", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_altivec_128_32, "parasail/parasail_nw_stats_striped_altivec_128_32", "nw_stats", "striped", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_altivec_128_16, "parasail/parasail_nw_stats_striped_altivec_128_16", "nw_stats", "striped", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_altivec_128_8, "parasail/parasail_nw_stats_striped_altivec_128_8", "nw_stats", "striped", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_altivec_128_64, "parasail/parasail_nw_stats_diag_altivec_128_64", "nw_stats",    "diag", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_altivec_128_32, "parasail/parasail_nw_stats_diag_altivec_128_32", "nw_stats",    "diag", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_altivec_128_16, "parasail/parasail_nw_stats_diag_altivec_128_16", "nw_stats",    "diag", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_altivec_128_8, "parasail/parasail_nw_stats_diag_altivec_128_8", "nw_stats",    "diag", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_stats_altivec = {"parasail/parasail_nw_stats_altivec", parasail_nw_stats_altivec_functions};
#endif
#if HAVE_NEON
static parasail_function_info_t parasail_nw_stats_neon_functions[] = {
{parasail_nw_stats,                   "parasail/parasail_nw_stats",                   "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_nw_stats_scan,              "parasail/parasail_nw_stats_scan",              "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_neon_128_64,  "parasail/parasail_nw_stats_scan_neon_128_64",  "nw_stats",    "scan", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_neon_128_32,  "parasail/parasail_nw_stats_scan_neon_128_32",  "nw_stats",    "scan", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_neon_128_16,  "parasail/parasail_nw_stats_scan_neon_128_16",  "nw_stats",    "scan", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_neon_128_8,   "parasail/parasail_nw_stats_scan_neon_128_8",   "nw_stats",    "scan", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_neon_128_64, "parasail/parasail_nw_stats_striped_neon_128_64", "nw_stats", "striped", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_neon_128_32, "parasail/parasail_nw_stats_striped_neon_128_32", "nw_stats", "striped", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_neon_128_16, "parasail/parasail_nw_stats_striped_neon_128_16", "nw_stats", "striped", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_neon_128_8, "parasail/parasail_nw_stats_striped_neon_128_8", "nw_stats", "striped", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_neon_128_64,  "parasail/parasail_nw_stats_diag_neon_128_64",  "nw_stats",    "diag", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_neon_128_32,  "parasail/parasail_nw_stats_diag_neon_128_32",  "nw_stats",    "diag", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_neon_128_16,  "parasail/parasail_nw_stats_diag_neon_128_16",  "nw_stats",    "diag", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_neon_128_8,   "parasail/parasail_nw_stats_diag_neon_128_8",   "nw_stats",    "diag", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_stats_neon = {"parasail/parasail_nw_stats_neon", parasail_nw_stats_neon_functions};
#endif
static parasail_function_info_t parasail_nw_stats_disp_functions[] = {
{parasail_nw_stats,                   "parasail/parasail_nw_stats",                   "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_nw_stats_scan,              "parasail/parasail_nw_stats_scan",              "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_64,           "parasail/parasail_nw_stats_scan_64",           "nw_stats",    "scan", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_32,           "parasail/parasail_nw_stats_scan_32",           "nw_stats",    "scan", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_16,           "parasail/parasail_nw_stats_scan_16",           "nw_stats",    "scan", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_8,            "parasail/parasail_nw_stats_scan_8",            "nw_stats",    "scan", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_64,        "parasail/parasail_nw_stats_striped_64",        "nw_stats", "striped", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_32,        "parasail/parasail_nw_stats_striped_32",        "nw_stats", "striped", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_16,        "parasail/parasail_nw_stats_striped_16",        "nw_stats", "striped", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_8,         "parasail/parasail_nw_stats_striped_8",         "nw_stats", "striped", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_64,           "parasail/parasail_nw_stats_diag_64",           "nw_stats",    "diag", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_32,           "parasail/parasail_nw_stats_diag_32",           "nw_stats",    "diag", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_16,           "parasail/parasail_nw_stats_diag_16",           "nw_stats",    "diag", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_8,            "parasail/parasail_nw_stats_diag_8",            "nw_stats",    "diag", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_scan_sat,          "parasail/parasail_nw_stats_scan_sat",          "nw_stats",    "scan", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_striped_sat,       "parasail/parasail_nw_stats_striped_sat",       "nw_stats", "striped", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_nw_stats_diag_sat,          "parasail/parasail_nw_stats_diag_sat",          "nw_stats",    "diag", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_nw_stats_disp = {"parasail/parasail_nw_stats_disp", parasail_nw_stats_disp_functions};
#if HAVE_SSE2
static parasail_function_info_t parasail_sg_stats_sse2_functions[] = {
{parasail_sg_stats,                   "parasail/parasail_sg_stats",                   "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sg_stats_scan,              "parasail/parasail_sg_stats_scan",              "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse2_128_64,  "parasail/parasail_sg_stats_scan_sse2_128_64",  "sg_stats",    "scan", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse2_128_32,  "parasail/parasail_sg_stats_scan_sse2_128_32",  "sg_stats",    "scan", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse2_128_16,  "parasail/parasail_sg_stats_scan_sse2_128_16",  "sg_stats",    "scan", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse2_128_8,   "parasail/parasail_sg_stats_scan_sse2_128_8",   "sg_stats",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse2_128_64, "parasail/parasail_sg_stats_striped_sse2_128_64", "sg_stats", "striped", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse2_128_32, "parasail/parasail_sg_stats_striped_sse2_128_32", "sg_stats", "striped", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse2_128_16, "parasail/parasail_sg_stats_striped_sse2_128_16", "sg_stats", "striped", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse2_128_8, "parasail/parasail_sg_stats_striped_sse2_128_8", "sg_stats", "striped", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse2_128_64,  "parasail/parasail_sg_stats_diag_sse2_128_64",  "sg_stats",    "diag", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse2_128_32,  "parasail/parasail_sg_stats_diag_sse2_128_32",  "sg_stats",    "diag", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse2_128_16,  "parasail/parasail_sg_stats_diag_sse2_128_16",  "sg_stats",    "diag", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse2_128_8,   "parasail/parasail_sg_stats_diag_sse2_128_8",   "sg_stats",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_stats_sse2 = {"parasail/parasail_sg_stats_sse2", parasail_sg_stats_sse2_functions};
#endif
#if HAVE_SSE41
static parasail_function_info_t parasail_sg_stats_sse41_functions[] = {
{parasail_sg_stats,                   "parasail/parasail_sg_stats",                   "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sg_stats_scan,              "parasail/parasail_sg_stats_scan",              "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse41_128_64, "parasail/parasail_sg_stats_scan_sse41_128_64", "sg_stats",    "scan", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse41_128_32, "parasail/parasail_sg_stats_scan_sse41_128_32", "sg_stats",    "scan", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse41_128_16, "parasail/parasail_sg_stats_scan_sse41_128_16", "sg_stats",    "scan", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sse41_128_8,  "parasail/parasail_sg_stats_scan_sse41_128_8",  "sg_stats",    "scan", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse41_128_64, "parasail/parasail_sg_stats_striped_sse41_128_64", "sg_stats", "striped", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse41_128_32, "parasail/parasail_sg_stats_striped_sse41_128_32", "sg_stats", "striped", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse41_128_16, "parasail/parasail_sg_stats_striped_sse41_128_16", "sg_stats", "striped", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sse41_128_8, "parasail/parasail_sg_stats_striped_sse41_128_8", "sg_stats", "striped", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse41_128_64, "parasail/parasail_sg_stats_diag_sse41_128_64", "sg_stats",    "diag", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse41_128_32, "parasail/parasail_sg_stats_diag_sse41_128_32", "sg_stats",    "diag", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse41_128_16, "parasail/parasail_sg_stats_diag_sse41_128_16", "sg_stats",    "diag", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sse41_128_8,  "parasail/parasail_sg_stats_diag_sse41_128_8",  "sg_stats",    "diag", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_stats_sse41 = {"parasail/parasail_sg_stats_sse41", parasail_sg_stats_sse41_functions};
#endif
#if HAVE_AVX2
static parasail_function_info_t parasail_sg_stats_avx2_functions[] = {
{parasail_sg_stats,                   "parasail/parasail_sg_stats",                   "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sg_stats_scan,              "parasail/parasail_sg_stats_scan",              "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_avx2_256_64,  "parasail/parasail_sg_stats_scan_avx2_256_64",  "sg_stats",    "scan", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_avx2_256_32,  "parasail/parasail_sg_stats_scan_avx2_256_32",  "sg_stats",    "scan", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_avx2_256_16,  "parasail/parasail_sg_stats_scan_avx2_256_16",  "sg_stats",    "scan", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_avx2_256_8,   "parasail/parasail_sg_stats_scan_avx2_256_8",   "sg_stats",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_avx2_256_64, "parasail/parasail_sg_stats_striped_avx2_256_64", "sg_stats", "striped", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_avx2_256_32, "parasail/parasail_sg_stats_striped_avx2_256_32", "sg_stats", "striped", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_avx2_256_16, "parasail/parasail_sg_stats_striped_avx2_256_16", "sg_stats", "striped", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_avx2_256_8, "parasail/parasail_sg_stats_striped_avx2_256_8", "sg_stats", "striped", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_avx2_256_64,  "parasail/parasail_sg_stats_diag_avx2_256_64",  "sg_stats",    "diag", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_avx2_256_32,  "parasail/parasail_sg_stats_diag_avx2_256_32",  "sg_stats",    "diag", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_avx2_256_16,  "parasail/parasail_sg_stats_diag_avx2_256_16",  "sg_stats",    "diag", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_avx2_256_8,   "parasail/parasail_sg_stats_diag_avx2_256_8",   "sg_stats",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_stats_avx2 = {"parasail/parasail_sg_stats_avx2", parasail_sg_stats_avx2_functions};
#endif
#if HAVE_ALTIVEC
static parasail_function_info_t parasail_sg_stats_altivec_functions[] = {
{parasail_sg_stats,                   "parasail/parasail_sg_stats",                   "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sg_stats_scan,              "parasail/parasail_sg_stats_scan",              "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_altivec_128_64, "parasail/parasail_sg_stats_scan_altivec_128_64", "sg_stats",    "scan", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_altivec_128_32, "parasail/parasail_sg_stats_scan_altivec_128_32", "sg_stats",    "scan", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_altivec_128_16, "parasail/parasail_sg_stats_scan_altivec_128_16", "sg_stats",    "scan", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_altivec_128_8, "parasail/parasail_sg_stats_scan_altivec_128_8", "sg_stats",    "scan", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_altivec_128_64, "parasail/parasail_sg_stats_striped_altivec_128_64", "sg_stats", "striped", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_altivec_128_32, "parasail/parasail_sg_stats_striped_altivec_128_32", "sg_stats", "striped", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_altivec_128_16, "parasail/parasail_sg_stats_striped_altivec_128_16", "sg_stats", "striped", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_altivec_128_8, "parasail/parasail_sg_stats_striped_altivec_128_8", "sg_stats", "striped", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_altivec_128_64, "parasail/parasail_sg_stats_diag_altivec_128_64", "sg_stats",    "diag", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_altivec_128_32, "parasail/parasail_sg_stats_diag_altivec_128_32", "sg_stats",    "diag", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_altivec_128_16, "parasail/parasail_sg_stats_diag_altivec_128_16", "sg_stats",    "diag", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_altivec_128_8, "parasail/parasail_sg_stats_diag_altivec_128_8", "sg_stats",    "diag", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_stats_altivec = {"parasail/parasail_sg_stats_altivec", parasail_sg_stats_altivec_functions};
#endif
#if HAVE_NEON
static parasail_function_info_t parasail_sg_stats_neon_functions[] = {
{parasail_sg_stats,                   "parasail/parasail_sg_stats",                   "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sg_stats_scan,              "parasail/parasail_sg_stats_scan",              "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_neon_128_64,  "parasail/parasail_sg_stats_scan_neon_128_64",  "sg_stats",    "scan", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_neon_128_32,  "parasail/parasail_sg_stats_scan_neon_128_32",  "sg_stats",    "scan", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_neon_128_16,  "parasail/parasail_sg_stats_scan_neon_128_16",  "sg_stats",    "scan", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_neon_128_8,   "parasail/parasail_sg_stats_scan_neon_128_8",   "sg_stats",    "scan", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_neon_128_64, "parasail/parasail_sg_stats_striped_neon_128_64", "sg_stats", "striped", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_neon_128_32, "parasail/parasail_sg_stats_striped_neon_128_32", "sg_stats", "striped", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_neon_128_16, "parasail/parasail_sg_stats_striped_neon_128_16", "sg_stats", "striped", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_neon_128_8, "parasail/parasail_sg_stats_striped_neon_128_8", "sg_stats", "striped", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_neon_128_64,  "parasail/parasail_sg_stats_diag_neon_128_64",  "sg_stats",    "diag", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_neon_128_32,  "parasail/parasail_sg_stats_diag_neon_128_32",  "sg_stats",    "diag", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_neon_128_16,  "parasail/parasail_sg_stats_diag_neon_128_16",  "sg_stats",    "diag", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_neon_128_8,   "parasail/parasail_sg_stats_diag_neon_128_8",   "sg_stats",    "diag", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_stats_neon = {"parasail/parasail_sg_stats_neon", parasail_sg_stats_neon_functions};
#endif
static parasail_function_info_t parasail_sg_stats_disp_functions[] = {
{parasail_sg_stats,                   "parasail/parasail_sg_stats",                   "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sg_stats_scan,              "parasail/parasail_sg_stats_scan",              "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_64,           "parasail/parasail_sg_stats_scan_64",           "sg_stats",    "scan", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_32,           "parasail/parasail_sg_stats_scan_32",           "sg_stats",    "scan", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_16,           "parasail/parasail_sg_stats_scan_16",           "sg_stats",    "scan", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_8,            "parasail/parasail_sg_stats_scan_8",            "sg_stats",    "scan", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_64,        "parasail/parasail_sg_stats_striped_64",        "sg_stats", "striped", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_32,        "parasail/parasail_sg_stats_striped_32",        "sg_stats", "striped", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_16,        "parasail/parasail_sg_stats_striped_16",        "sg_stats", "striped", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_8,         "parasail/parasail_sg_stats_striped_8",         "sg_stats", "striped", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_64,           "parasail/parasail_sg_stats_diag_64",           "sg_stats",    "diag", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_32,           "parasail/parasail_sg_stats_diag_32",           "sg_stats",    "diag", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_16,           "parasail/parasail_sg_stats_diag_16",           "sg_stats",    "diag", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_8,            "parasail/parasail_sg_stats_diag_8",            "sg_stats",    "diag", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_scan_sat,          "parasail/parasail_sg_stats_scan_sat",          "sg_stats",    "scan", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_striped_sat,       "parasail/parasail_sg_stats_striped_sat",       "sg_stats", "striped", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sg_stats_diag_sat,          "parasail/parasail_sg_stats_diag_sat",          "sg_stats",    "diag", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sg_stats_disp = {"parasail/parasail_sg_stats_disp", parasail_sg_stats_disp_functions};
#if HAVE_SSE2
static parasail_function_info_t parasail_sw_stats_sse2_functions[] = {
{parasail_sw_stats,                   "parasail/parasail_sw_stats",                   "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sw_stats_scan,              "parasail/parasail_sw_stats_scan",              "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse2_128_64,  "parasail/parasail_sw_stats_scan_sse2_128_64",  "sw_stats",    "scan", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse2_128_32,  "parasail/parasail_sw_stats_scan_sse2_128_32",  "sw_stats",    "scan", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse2_128_16,  "parasail/parasail_sw_stats_scan_sse2_128_16",  "sw_stats",    "scan", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse2_128_8,   "parasail/parasail_sw_stats_scan_sse2_128_8",   "sw_stats",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse2_128_64, "parasail/parasail_sw_stats_striped_sse2_128_64", "sw_stats", "striped", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse2_128_32, "parasail/parasail_sw_stats_striped_sse2_128_32", "sw_stats", "striped", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse2_128_16, "parasail/parasail_sw_stats_striped_sse2_128_16", "sw_stats", "striped", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse2_128_8, "parasail/parasail_sw_stats_striped_sse2_128_8", "sw_stats", "striped", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse2_128_64,  "parasail/parasail_sw_stats_diag_sse2_128_64",  "sw_stats",    "diag", "sse2",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse2_128_32,  "parasail/parasail_sw_stats_diag_sse2_128_32",  "sw_stats",    "diag", "sse2",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse2_128_16,  "parasail/parasail_sw_stats_diag_sse2_128_16",  "sw_stats",    "diag", "sse2",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse2_128_8,   "parasail/parasail_sw_stats_diag_sse2_128_8",   "sw_stats",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_stats_sse2 = {"parasail/parasail_sw_stats_sse2", parasail_sw_stats_sse2_functions};
#endif
#if HAVE_SSE41
static parasail_function_info_t parasail_sw_stats_sse41_functions[] = {
{parasail_sw_stats,                   "parasail/parasail_sw_stats",                   "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sw_stats_scan,              "parasail/parasail_sw_stats_scan",              "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse41_128_64, "parasail/parasail_sw_stats_scan_sse41_128_64", "sw_stats",    "scan", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse41_128_32, "parasail/parasail_sw_stats_scan_sse41_128_32", "sw_stats",    "scan", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse41_128_16, "parasail/parasail_sw_stats_scan_sse41_128_16", "sw_stats",    "scan", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sse41_128_8,  "parasail/parasail_sw_stats_scan_sse41_128_8",  "sw_stats",    "scan", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse41_128_64, "parasail/parasail_sw_stats_striped_sse41_128_64", "sw_stats", "striped", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse41_128_32, "parasail/parasail_sw_stats_striped_sse41_128_32", "sw_stats", "striped", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse41_128_16, "parasail/parasail_sw_stats_striped_sse41_128_16", "sw_stats", "striped", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sse41_128_8, "parasail/parasail_sw_stats_striped_sse41_128_8", "sw_stats", "striped", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse41_128_64, "parasail/parasail_sw_stats_diag_sse41_128_64", "sw_stats",    "diag", "sse41", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse41_128_32, "parasail/parasail_sw_stats_diag_sse41_128_32", "sw_stats",    "diag", "sse41", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse41_128_16, "parasail/parasail_sw_stats_diag_sse41_128_16", "sw_stats",    "diag", "sse41", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sse41_128_8,  "parasail/parasail_sw_stats_diag_sse41_128_8",  "sw_stats",    "diag", "sse41", "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_stats_sse41 = {"parasail/parasail_sw_stats_sse41", parasail_sw_stats_sse41_functions};
#endif
#if HAVE_AVX2
static parasail_function_info_t parasail_sw_stats_avx2_functions[] = {
{parasail_sw_stats,                   "parasail/parasail_sw_stats",                   "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sw_stats_scan,              "parasail/parasail_sw_stats_scan",              "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_avx2_256_64,  "parasail/parasail_sw_stats_scan_avx2_256_64",  "sw_stats",    "scan", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_avx2_256_32,  "parasail/parasail_sw_stats_scan_avx2_256_32",  "sw_stats",    "scan", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_avx2_256_16,  "parasail/parasail_sw_stats_scan_avx2_256_16",  "sw_stats",    "scan", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_avx2_256_8,   "parasail/parasail_sw_stats_scan_avx2_256_8",   "sw_stats",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_avx2_256_64, "parasail/parasail_sw_stats_striped_avx2_256_64", "sw_stats", "striped", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_avx2_256_32, "parasail/parasail_sw_stats_striped_avx2_256_32", "sw_stats", "striped", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_avx2_256_16, "parasail/parasail_sw_stats_striped_avx2_256_16", "sw_stats", "striped", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_avx2_256_8, "parasail/parasail_sw_stats_striped_avx2_256_8", "sw_stats", "striped", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_avx2_256_64,  "parasail/parasail_sw_stats_diag_avx2_256_64",  "sw_stats",    "diag", "avx2",  "256", "64",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_avx2_256_32,  "parasail/parasail_sw_stats_diag_avx2_256_32",  "sw_stats",    "diag", "avx2",  "256", "32",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_avx2_256_16,  "parasail/parasail_sw_stats_diag_avx2_256_16",  "sw_stats",    "diag", "avx2",  "256", "16", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_avx2_256_8,   "parasail/parasail_sw_stats_diag_avx2_256_8",   "sw_stats",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_stats_avx2 = {"parasail/parasail_sw_stats_avx2", parasail_sw_stats_avx2_functions};
#endif
#if HAVE_ALTIVEC
static parasail_function_info_t parasail_sw_stats_altivec_functions[] = {
{parasail_sw_stats,                   "parasail/parasail_sw_stats",                   "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sw_stats_scan,              "parasail/parasail_sw_stats_scan",              "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_altivec_128_64, "parasail/parasail_sw_stats_scan_altivec_128_64", "sw_stats",    "scan", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_altivec_128_32, "parasail/parasail_sw_stats_scan_altivec_128_32", "sw_stats",    "scan", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_altivec_128_16, "parasail/parasail_sw_stats_scan_altivec_128_16", "sw_stats",    "scan", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_altivec_128_8, "parasail/parasail_sw_stats_scan_altivec_128_8", "sw_stats",    "scan", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_altivec_128_64, "parasail/parasail_sw_stats_striped_altivec_128_64", "sw_stats", "striped", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_altivec_128_32, "parasail/parasail_sw_stats_striped_altivec_128_32", "sw_stats", "striped", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_altivec_128_16, "parasail/parasail_sw_stats_striped_altivec_128_16", "sw_stats", "striped", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_altivec_128_8, "parasail/parasail_sw_stats_striped_altivec_128_8", "sw_stats", "striped", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_altivec_128_64, "parasail/parasail_sw_stats_diag_altivec_128_64", "sw_stats",    "diag", "altivec", "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_altivec_128_32, "parasail/parasail_sw_stats_diag_altivec_128_32", "sw_stats",    "diag", "altivec", "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_altivec_128_16, "parasail/parasail_sw_stats_diag_altivec_128_16", "sw_stats",    "diag", "altivec", "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_altivec_128_8, "parasail/parasail_sw_stats_diag_altivec_128_8", "sw_stats",    "diag", "altivec", "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_stats_altivec = {"parasail/parasail_sw_stats_altivec", parasail_sw_stats_altivec_functions};
#endif
#if HAVE_NEON
static parasail_function_info_t parasail_sw_stats_neon_functions[] = {
{parasail_sw_stats,                   "parasail/parasail_sw_stats",                   "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sw_stats_scan,              "parasail/parasail_sw_stats_scan",              "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_neon_128_64,  "parasail/parasail_sw_stats_scan_neon_128_64",  "sw_stats",    "scan", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_neon_128_32,  "parasail/parasail_sw_stats_scan_neon_128_32",  "sw_stats",    "scan", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_neon_128_16,  "parasail/parasail_sw_stats_scan_neon_128_16",  "sw_stats",    "scan", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_neon_128_8,   "parasail/parasail_sw_stats_scan_neon_128_8",   "sw_stats",    "scan", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_neon_128_64, "parasail/parasail_sw_stats_striped_neon_128_64", "sw_stats", "striped", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_neon_128_32, "parasail/parasail_sw_stats_striped_neon_128_32", "sw_stats", "striped", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_neon_128_16, "parasail/parasail_sw_stats_striped_neon_128_16", "sw_stats", "striped", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_neon_128_8, "parasail/parasail_sw_stats_striped_neon_128_8", "sw_stats", "striped", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_neon_128_64,  "parasail/parasail_sw_stats_diag_neon_128_64",  "sw_stats",    "diag", "neon",  "128", "64",  2, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_neon_128_32,  "parasail/parasail_sw_stats_diag_neon_128_32",  "sw_stats",    "diag", "neon",  "128", "32",  4, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_neon_128_16,  "parasail/parasail_sw_stats_diag_neon_128_16",  "sw_stats",    "diag", "neon",  "128", "16",  8, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_neon_128_8,   "parasail/parasail_sw_stats_diag_neon_128_8",   "sw_stats",    "diag", "neon",  "128",  "8", 16, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_stats_neon = {"parasail/parasail_sw_stats_neon", parasail_sw_stats_neon_functions};
#endif
static parasail_function_info_t parasail_sw_stats_disp_functions[] = {
{parasail_sw_stats,                   "parasail/parasail_sw_stats",                   "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 0, 0, 1, 1},
{parasail_sw_stats_scan,              "parasail/parasail_sw_stats_scan",              "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_64,           "parasail/parasail_sw_stats_scan_64",           "sw_stats",    "scan", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_32,           "parasail/parasail_sw_stats_scan_32",           "sw_stats",    "scan", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_16,           "parasail/parasail_sw_stats_scan_16",           "sw_stats",    "scan", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_8,            "parasail/parasail_sw_stats_scan_8",            "sw_stats",    "scan", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_64,        "parasail/parasail_sw_stats_striped_64",        "sw_stats", "striped", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_32,        "parasail/parasail_sw_stats_striped_32",        "sw_stats", "striped", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_16,        "parasail/parasail_sw_stats_striped_16",        "sw_stats", "striped", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_8,         "parasail/parasail_sw_stats_striped_8",         "sw_stats", "striped", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_64,           "parasail/parasail_sw_stats_diag_64",           "sw_stats",    "diag", "disp",   "NA", "64", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_32,           "parasail/parasail_sw_stats_diag_32",           "sw_stats",    "diag", "disp",   "NA", "32", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_16,           "parasail/parasail_sw_stats_diag_16",           "sw_stats",    "diag", "disp",   "NA", "16", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_8,            "parasail/parasail_sw_stats_diag_8",            "sw_stats",    "diag", "disp",   "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_scan_sat,          "parasail/parasail_sw_stats_scan_sat",          "sw_stats",    "scan", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_striped_sat,       "parasail/parasail_sw_stats_striped_sat",       "sw_stats", "striped", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{parasail_sw_stats_diag_sat,          "parasail/parasail_sw_stats_diag_sat",          "sw_stats",    "diag", "sat",    "NA",  "8", -1, 0, 0, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0, 0, 0}
};
static parasail_function_group_t parasail_sw_stats_disp = {"parasail/parasail_sw_stats_disp", parasail_sw_stats_disp_functions};

#endif /* _PARASAIL_FUNCTION_GROUP_H_ */

