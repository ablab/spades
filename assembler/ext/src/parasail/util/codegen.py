#!/usr/bin/python
#
# This file uses the C templates located in ./templates and fills them
# in based on the infomation and script contained in this file.
# The template is filled in as appropriate for various instruction sets
# and integer precisions.
#
# author jeff.daily@pnnl.gov
#
# Copyright (c) 2015 Battelle Memorial Institute.

import copy
import os
import re
import string
import sys

from isa import sse2
from isa import sse41
from isa import avx2
from isa import altivec
from isa import neon

keys = sse2.keys()

# gather templates
template_dir = "templates/"

template_filenames = [
"nw_diag.c",
"nw_scan.c",
"nw_striped.c",
"sg_diag.c",
"sg_scan.c",
"sg_striped.c",
"sw_diag.c",
"sw_scan.c",
"sw_striped.c",

"nw_stats_diag.c",
"nw_stats_scan.c",
"nw_stats_striped.c",
"sg_stats_diag.c",
"sg_stats_scan.c",
"sg_stats_striped.c",
"sw_stats_diag.c",
"sw_stats_scan.c",
"sw_stats_striped.c",

"nw_trace_diag.c",
"nw_trace_scan.c",
"nw_trace_striped.c",
"sg_trace_diag.c",
"sg_trace_scan.c",
"sg_trace_striped.c",
"sw_trace_diag.c",
"sw_trace_scan.c",
"sw_trace_striped.c",
]

special_templates = [
"sg_diag_8.c",
"sw_diag_8.c",
"sw_stats_diag_8.c",
"sg_trace_diag_8.c",
"sw_trace_diag_8.c",
]

bias_templates = [
"sw_striped_bias.c",
"sw_stats_striped_bias.c",
]


output_dir = "generated/"
if not os.path.exists(output_dir):
        os.makedirs(output_dir)


def generate_H(params):
    text = ""
    if "striped" in params["NAME"]:
        params["PVH_VAR"] = "pvHStore"
    else:
        params["PVH_VAR"] = "pvH"
    if "neon" in params["ISA"]:
        text = """    /* initialize H */
    {
        %(INDEX)s index = 0;
        for (i=0; i<segLen; ++i) {
            %(INDEX)s segNum = 0;
            %(VTYPE)s h;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.i%(WIDTH)s[segNum] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
            }
            %(VSTORE)s(&%(PVH_VAR)s[index], h);
            ++index;
        }
    }""" % params
    else:
       text = """    /* initialize H */
    {
        %(INDEX)s index = 0;
        for (i=0; i<segLen; ++i) {
            %(INDEX)s segNum = 0;
            %(VTYPE)s_%(WIDTH)s_t h;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
            }
            %(VSTORE)s(&%(PVH_VAR)s[index], h.m);
            ++index;
        }
    }""" % params
    params["INIT_H"] = text
    return params


def generate_H_and_E(params):
    text = ""
    if "striped" in params["NAME"]:
        params["PVH_VAR"] = "pvHStore"
    else:
        params["PVH_VAR"] = "pvH"
    params["PVEA_STORE"] = ""
    if "striped" in params["NAME"] and "trace" in params["NAME"]:
        if "neon" in params["ISA"]:
            params["E_M_VAR"] = "e"
        else:
            params["E_M_VAR"] = "e.m"
        params["PVEA_STORE"] = """
            %(VSTORE)s(&pvEaStore[index], %(E_M_VAR)s);""" % params
    if "neon" in params["ISA"]:
        text = """    /* initialize H and E */
    {
        %(INDEX)s index = 0;
        for (i=0; i<segLen; ++i) {
            %(INDEX)s segNum = 0;
            %(VTYPE)s h;
            %(VTYPE)s e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.i%(WIDTH)s[segNum] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
                tmp = tmp - open;
                e.i%(WIDTH)s[segNum] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
            }
            %(VSTORE)s(&%(PVH_VAR)s[index], h);
            %(VSTORE)s(&pvE[index], e);%(PVEA_STORE)s
            ++index;
        }
    }""" % params
    else:
       text = """    /* initialize H and E */
    {
        %(INDEX)s index = 0;
        for (i=0; i<segLen; ++i) {
            %(INDEX)s segNum = 0;
            %(VTYPE)s_%(WIDTH)s_t h;
            %(VTYPE)s_%(WIDTH)s_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
            }
            %(VSTORE)s(&%(PVH_VAR)s[index], h.m);
            %(VSTORE)s(&pvE[index], e.m);%(PVEA_STORE)s
            ++index;
        }
    }""" % params
    params["INIT_H_AND_E"] = text
    return params


def generate_printer(params):
    text = ""
    trace = ""
    bias = ""
    rowcol = ""
    bias_rowcol = ""
    if "striped" in params["NAME"] or "scan" in params["NAME"]:
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            if params["LANES"] / 10:
                text += "    array[1LL*(%(LANE)2d*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)2d);\n" % params
                trace += "    array[1LL*(%(LANE)2d*seglen+t)*dlen + d] = (int8_t)%(VEXTRACT)s(vH, %(LANE)2d);\n" % params
            else:
                text += "    array[1LL*(%(LANE)s*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)s);\n" % params
                trace += "    array[1LL*(%(LANE)s*seglen+t)*dlen + d] = (int8_t)%(VEXTRACT)s(vH, %(LANE)s);\n" % params
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            if params["LANES"] / 10:
                rowcol += "    col[%(LANE)2d*seglen+t] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)2d);\n" % params
            else:
                rowcol += "    col[%(LANE)s*seglen+t] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)s);\n" % params
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            if params["LANES"] / 10:
                bias += "    array[1LL*(%(LANE)2d*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)2d) - bias;\n" % params
            else:
                bias += "    array[1LL*(%(LANE)s*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)s) - bias;\n" % params
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            if params["LANES"] / 10:
                bias_rowcol += "    col[%(LANE)2d*seglen+t] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)2d) - bias;\n" % params
            else:
                bias_rowcol += "    col[%(LANE)s*seglen+t] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)s) - bias;\n" % params
    elif "diag" in params["NAME"]:
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            params["LANE_END"] = params["LANES"]-lane-1
            text += """
    if (0 <= i+%(LANE)s && i+%(LANE)s < s1Len && 0 <= j-%(LANE)s && j-%(LANE)s < s2Len) {
        array[1LL*(i+%(LANE)s)*s2Len + (j-%(LANE)s)] = (%(INT)s)%(VEXTRACT)s(vWH, %(LANE_END)s);
    }\n"""[1:] % params
            trace += """
    if (0 <= i+%(LANE)s && i+%(LANE)s < s1Len && 0 <= j-%(LANE)s && j-%(LANE)s < s2Len) {
        array[1LL*(i+%(LANE)s)*s2Len + (j-%(LANE)s)] = (int8_t)%(VEXTRACT)s(vWH, %(LANE_END)s);
    }\n"""[1:] % params
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            params["LANE_END"] = params["LANES"]-lane-1
            rowcol += """
    if (i+%(LANE)s == s1Len-1 && 0 <= j-%(LANE)s && j-%(LANE)s < s2Len) {
        row[j-%(LANE)s] = (%(INT)s)%(VEXTRACT)s(vWH, %(LANE_END)s);
    }\n"""[1:] % params
            rowcol += """
    if (j-%(LANE)s == s2Len-1 && 0 <= i+%(LANE)s && i+%(LANE)s < s1Len) {
        col[(i+%(LANE)s)] = (%(INT)s)%(VEXTRACT)s(vWH, %(LANE_END)s);
    }\n"""[1:] % params
    else:
        print "bad printer name"
        sys.exit(1)
    params["PRINTER"] = text[:-1] # remove last newline
    params["PRINTER_TRACE"] = trace[:-1] # remove last newline
    params["PRINTER_BIAS"] = bias[:-1] # remove last newline
    params["PRINTER_ROWCOL"] = rowcol[:-1] # remove last newline
    params["PRINTER_BIAS_ROWCOL"] = bias_rowcol[:-1] # remove last newline
    return params


def generate_saturation_check_old(params):
    width = params["WIDTH"]
    if width == 8:

        params["SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vSaturationCheck = %(VSET0)s();
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);""".strip() % params

        params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheck = %(VOR)s(vSaturationCheck,
                        %(VOR)s(
                            %(VCMPEQ)s(vH, vNegLimit),
                            %(VCMPEQ)s(vH, vPosLimit)));
            }""".strip() % params

        params["SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(vSaturationCheck)) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
    }""".strip() % params

        params["STATS_SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vSaturationCheck = %(VSET0)s();
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);""".strip() % params

        params["STATS_SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheck = %(VOR)s(vSaturationCheck,
                        %(VOR)s(
                            %(VCMPEQ)s(vH, vNegLimit),
                            %(VCMPEQ)s(vH, vPosLimit)));
            }""".strip() % params

        params["STATS_SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(vSaturationCheck)) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
    }""".strip() % params

        params["NEG_INF"] = "INT8_MIN"
        params["VADD"] = params["VADDSx8"]
        params["VSUB"] = params["VSUBSx8"]
    else:
        params["SATURATION_CHECK_INIT"] = ""
        params["SATURATION_CHECK_MID"] = ""
        params["SATURATION_CHECK_FINAL"] = ""
        params["STATS_SATURATION_CHECK_INIT"] = ""
        params["STATS_SATURATION_CHECK_MID"] = ""
        params["STATS_SATURATION_CHECK_FINAL"] = ""
    return params


def generate_saturation_check(params):
    width = params["WIDTH"]
    if width == 8:
        params["SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;""".strip() % params
        if "diag" in params["NAME"]:
            params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWH);
            }""".strip() % params
        elif "scan" in params["NAME"]:
            params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            }""".strip() % params
        elif "striped" in params["NAME"]:
            params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vE);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vF);
            }""".strip() % params
        else:
            params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            }""".strip() % params

        params["SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPEQ)s(vSaturationCheckMin, vNegLimit),
            %(VCMPEQ)s(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }""".strip() % params

        params["STATS_SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;""".strip() % params
        if "diag" in params["NAME"]:
            params["STATS_SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWL);
            }""".strip() % params
        elif "scan" in params["NAME"]:
            params["STATS_SATURATION_CHECK_MID1"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            }""".strip() % params
            params["STATS_SATURATION_CHECK_MID2"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vL);
            }""".strip() % params
        else:
            params["STATS_SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
            }""".strip() % params

        params["STATS_SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPEQ)s(vSaturationCheckMin, vNegLimit),
            %(VCMPEQ)s(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }""".strip() % params

        params["NEG_INF"] = "INT8_MIN"
        params["VADD"] = params["VADDSx8"]
        params["VSUB"] = params["VSUBSx8"]
        if "sw" in params["NAME"] and "striped" in params["NAME"]:
            pass
        else:
            for p in ["VMAX", "VMIN"]:
                if (params[p].endswith("_rpl")
                        and params[p] not in params["FIXES"]):
                    params["FIXES"] += params[params[p]]
    else:
        params["SATURATION_CHECK_INIT"] = ""
        params["SATURATION_CHECK_MID"] = ""
        params["SATURATION_CHECK_FINAL"] = ""
        params["STATS_SATURATION_CHECK_INIT"] = ""
        params["STATS_SATURATION_CHECK_MID"] = ""
        params["STATS_SATURATION_CHECK_MID1"] = ""
        params["STATS_SATURATION_CHECK_MID2"] = ""
        params["STATS_SATURATION_CHECK_FINAL"] = ""
    return params


def generated_params_diag(params):
    lanes = params["LANES"]
    params["DIAG_I"] = ",".join(["%d"%i for i in range(lanes)])
    params["DIAG_ILO"] = ",".join(["%d"%i for i in range(lanes/2,lanes)])
    params["DIAG_IHI"] = ",".join(["%d"%i for i in range(lanes/2)])
    params["DIAG_J"] = ",".join(["%d"%-i for i in range(lanes)])
    params["DIAG_JLO"] = ",".join(["%d"%-i for i in range(lanes/2,lanes)])
    params["DIAG_JHI"] = ",".join(["%d"%-i for i in range(lanes/2)])
    params["DIAG_IBoundary"] = "            ".join(
            ["-open-%d*gap,\n"%(i)
                for i in range(lanes)])[:-2]
    params["DIAG_VS1"] = "                ".join(
            ["s1[i+%d],\n"%i
                for i in range(lanes)])[:-2]
    params["DIAG_MATROW_DECL"] = "        ".join(
            ["const int * const restrict matrow%d = &matrix->matrix[matrix->size*s1[i+%d]];\n"%(i,i)
                for i in range(lanes)])[:-1]
    params["DIAG_MATROW_USE"] = "                    ".join(
            ["matrow%d[s2[j-%d]],\n"%(i,i)
                for i in range(lanes)])[:-2]
    return params


def generated_params_striped(params):
    params["STRIPED_INSERT_MASK"] = "0,"*(params["LANES"]-1)+"1"
    params["POSITION_MASK"] = ",".join([str(i) for i in range(params["LANES"])])
    return params


def generated_params_scan(params):
    lanes = params["LANES"]
    if params["LANES"] / 10:
        params["STATS_SCAN_VFT"] = ("    "*3).join(
                ["tmp.v[%2d] = MAX(tmp.v[%2d]-segLen*gap, tmp.v[%2d]);\n"%(i,i-1,i)
                    for i in range(1,lanes)])[:-1]
        params["STATS_SCAN_UMP"] = ("    "*4).join(
                ["uMp.v[%2d] = uC.v[%2d] ? uMp.v[%2d] : uMp.v[%2d];\n"%(i,i,i-1,i)
                    for i in range(1,lanes)])[:-1]
        params["STATS_SCAN_USP"] = ("    "*4).join(
                ["uSp.v[%2d] = uC.v[%2d] ? uSp.v[%2d] : uSp.v[%2d];\n"%(i,i,i-1,i)
                    for i in range(1,lanes)])[:-1]
        params["STATS_SCAN_ULP"] = ("    "*4).join(
                ["uLp.v[%2d] = uC.v[%2d] ? uLp.v[%2d] + uLp.v[%2d] : uLp.v[%2d];\n"%(
                    i,i,i,i-1,i)
                    for i in range(1,lanes)])[:-1]
    else:
        params["STATS_SCAN_VFT"] = ("    "*3).join(
                ["tmp.v[%d] = MAX(tmp.v[%d]-segLen*gap, tmp.v[%d]);\n"%(i,i-1,i)
                    for i in range(1,lanes)])[:-1]
        params["STATS_SCAN_UMP"] = ("    "*4).join(
                ["uMp.v[%d] = uC.v[%d] ? uMp.v[%d] : uMp.v[%d];\n"%(i,i,i-1,i)
                    for i in range(1,lanes)])[:-1]
        params["STATS_SCAN_USP"] = ("    "*4).join(
                ["uSp.v[%d] = uC.v[%d] ? uSp.v[%d] : uSp.v[%d];\n"%(i,i,i-1,i)
                    for i in range(1,lanes)])[:-1]
        params["STATS_SCAN_ULP"] = ("    "*4).join(
                ["uLp.v[%d] = uC.v[%d] ? uLp.v[%d] + uLp.v[%d] : uLp.v[%d];\n"%(
                    i,i,i,i-1,i)
                    for i in range(1,lanes)])[:-1]
    params["STATS_SCAN_INSERT_MASK"] = "0,"*(params["LANES"]-1)+"1"
    params["SCAN_INSERT_MASK"] = "1"+",0"*(params["LANES"]-1)
    params["SCAN_NEG_INF_FRONT"] = "0,"*(params["LANES"]-1)+"NEG_LIMIT"
    params["POSITION_MASK"] = ",".join([str(i) for i in range(params["LANES"])])
    if "avx" in params["ISA"]:
        params["SCAN_AVX2_BLENDV_FIX"] = """

/* clang optimization broke blendv in this code */
#if defined(__clang__) && defined(__OPTIMIZE__)
#define _mm256_blendv_epi8 _mm256_blendv_epi8_rpl
static inline __m256i _mm256_blendv_epi8_rpl(__m256i a, __m256i b, __m256i mask) {
    a = _mm256_andnot_si256(mask, a);
    a = _mm256_or_si256(a, _mm256_and_si256(mask, b));
    return a;
}
#endif
    """
    else:
        params["SCAN_AVX2_BLENDV_FIX"] = ""

    return params


def generated_params(template, params):
    # some params are generated from given params
    bits = params["BITS"]
    width = params["WIDTH"]
    params["INDEX"] = "int32_t"
    params["ALIGNMENT"] = bits/8
    params["BYTES"] = width/8
    params["LANES"] = bits/width
    params["LAST_POS"] = params["LANES"]-1
    params["INT"] = "int%(WIDTH)s_t" % params
    params["NEG_INF"] = "(INT%(WIDTH)s_MIN/(%(INT)s)(2))" % params
    if "diag" in params["NAME"]:
        params = generated_params_diag(params)
    elif "striped" in params["NAME"]:
        params = generated_params_striped(params)
    elif "scan" in params["NAME"]:
        params = generated_params_scan(params)
    # select appropriate vector functions for given width
    suffix = "x%s" % width
    for key in keys:
        if key.endswith(suffix):
            new_key = key.split('x')[0]
            params[new_key] = params[key]
    fixes = ""
    template_params = re.findall(r'%\([A-Za-z0-9]+\)s', template)
    for param in params:
        wrapped_param = r'%%(%s)s'%param
        if (wrapped_param in template_params
                and str(params[param]).endswith("_rpl")):
            fixes += params[params[param]]
    if ("trace" in params["NAME"]
            and ("scan" in params["NAME"] or "striped" in params["NAME"])
            and "nw" not in params["NAME"]
            and "sg" not in params["NAME"]):
        if params["VEXTRACT"].endswith("_rpl"):
            fixes = params[params["VEXTRACT"]] + fixes
    params["FIXES"] = fixes
    params = generate_printer(params)
    params = generate_saturation_check(params)
    params = generate_H(params)
    params = generate_H_and_E(params)
    return params


for template_filename in template_filenames:
    template = open(template_dir+template_filename).read()
    for width in [64,32,16,8]:
        for isa in [sse2,sse41,avx2,altivec,neon]:
            params = copy.deepcopy(isa)
            params["WIDTH"] = width
            prefix = template_filename[:-2]
            prefix_prof = prefix + "_profile"
            parts = prefix.split('_')
            table_prefix = ""
            rowcol_prefix = ""
            trace_prefix = ""
            if len(parts) == 2:
                table_prefix = "%s_table_%s" % (parts[0], parts[1])
                rowcol_prefix = "%s_rowcol_%s" % (parts[0], parts[1])
                trace_prefix = "%s_trace_%s" % (parts[0], parts[1])
            if len(parts) == 3:
                table_prefix = "%s_%s_table_%s" % (parts[0], parts[1], parts[2])
                rowcol_prefix = "%s_%s_rowcol_%s" % (parts[0], parts[1], parts[2])
                trace_prefix = "%s_trace_%s" % (parts[0], parts[2])
            table_prefix_prof = table_prefix + "_profile"
            rowcol_prefix_prof = rowcol_prefix + "_profile"
            trace_prefix_prof = trace_prefix + "_profile"
            function_name = "%s_%s%s_%s_%s" % (prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_table_name = "%s_%s%s_%s_%s" % (table_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_rowcol_name = "%s_%s%s_%s_%s" % (rowcol_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_trace_name = "%s_%s%s_%s_%s" % (trace_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_pname = "%s_%s%s_%s_%s" % (prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_table_pname = "%s_%s%s_%s_%s" % (table_prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_rowcol_pname = "%s_%s%s_%s_%s" % (rowcol_prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_trace_pname = "%s_%s%s_%s_%s" % (trace_prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            params["NAME"] = "parasail/parasail_"+function_name
            params["NAME_BASE"] = string.replace(params["NAME"], "_stats", "")
            params["NAME_TABLE"] = "parasail/parasail_"+function_table_name
            params["NAME_ROWCOL"] = "parasail/parasail_"+function_rowcol_name
            params["NAME_TRACE"] = "parasail/parasail_"+function_trace_name
            params["PNAME"] = "parasail/parasail_"+function_pname
            params["PNAME_BASE"] = string.replace(params["PNAME"], "_stats", "")
            params["PNAME_TABLE"] = "parasail/parasail_"+function_table_pname
            params["PNAME_ROWCOL"] = "parasail/parasail_"+function_rowcol_pname
            params["PNAME_TRACE"] = "parasail/parasail_"+function_trace_pname
            params = generated_params(template, params)
            output_filename = "%s%s.c" % (output_dir, function_name)
            result = template % params
            writer = open(output_filename, "w")
            writer.write(template % params)
            writer.write("\n")
            writer.close()


# some templates have specializations for certain bit widths, e.g., 8
for template_filename in special_templates:
    template = open(template_dir+template_filename).read()
    prefix = template_filename[:-2]
    parts = prefix.split('_')
    width = int(parts[-1])
    parts = parts[:-1]
    prefix = "_".join(parts)
    table_prefix = ""
    rowcol_prefix = ""
    trace_prefix = ""
    if len(parts) == 2:
        table_prefix = "%s_table_%s" % (parts[0], parts[1])
        rowcol_prefix = "%s_rowcol_%s" % (parts[0], parts[1])
        trace_prefix = "%s_trace_%s" % (parts[0], parts[1])
    if len(parts) == 3:
        table_prefix = "%s_%s_table_%s" % (parts[0], parts[1], parts[2])
        rowcol_prefix = "%s_%s_rowcol_%s" % (parts[0], parts[1], parts[2])
        trace_prefix = "%s_%s_%s" % (parts[0], parts[1], parts[2])
    for isa in [sse2,sse41,avx2,altivec,neon]:
        params = copy.deepcopy(isa)
        params["WIDTH"] = width
        function_name = "%s_%s%s_%s_%s" % (prefix,
                isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
        function_table_name = "%s_%s%s_%s_%s" % (table_prefix,
                isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
        function_rowcol_name = "%s_%s%s_%s_%s" % (rowcol_prefix,
                isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
        function_trace_name = "%s_%s%s_%s_%s" % (trace_prefix,
                isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
        params["NAME"] = "parasail/parasail_"+function_name
        params["NAME_TABLE"] = "parasail/parasail_"+function_table_name
        params["NAME_ROWCOL"] = "parasail/parasail_"+function_rowcol_name
        params["NAME_TRACE"] = "parasail/parasail_"+function_trace_name
        params = generated_params(template, params)
        output_filename = "%s%s.c" % (output_dir, function_name)
        result = template % params
        writer = open(output_filename, "w")
        writer.write(template % params)
        writer.write("\n")
        writer.close()

# some templates have specializations for using bias
for template_filename in bias_templates:
    template = open(template_dir+template_filename).read()
    prefix = template_filename[:-2]
    parts = prefix.split('_')
    parts = parts[:-1]
    prefix = "_".join(parts)
    prefix_prof = prefix + "_profile"
    table_prefix = ""
    rowcol_prefix = ""
    trace_prefix = ""
    if len(parts) == 2:
        table_prefix = "%s_table_%s" % (parts[0], parts[1])
        rowcol_prefix = "%s_rowcol_%s" % (parts[0], parts[1])
        trace_prefix = "%s_trace_%s" % (parts[0], parts[1])
    if len(parts) == 3:
        table_prefix = "%s_%s_table_%s" % (parts[0], parts[1], parts[2])
        rowcol_prefix = "%s_%s_rowcol_%s" % (parts[0], parts[1], parts[2])
        trace_prefix = "%s_%s_trace_%s" % (parts[0], parts[1], parts[2])
    table_prefix_prof = table_prefix + "_profile"
    rowcol_prefix_prof = rowcol_prefix + "_profile"
    trace_prefix_prof = trace_prefix + "_profile"
    for width in [16,8]:
        for isa in [sse2,sse41,avx2,altivec,neon]:
            params = copy.deepcopy(isa)
            params["WIDTH"] = width
            function_name = "%s_%s%s_%s_%s" % (prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_table_name = "%s_%s%s_%s_%s" % (table_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_rowcol_name = "%s_%s%s_%s_%s" % (rowcol_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_trace_name = "%s_%s%s_%s_%s" % (trace_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_pname = "%s_%s%s_%s_%s" % (prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_table_pname = "%s_%s%s_%s_%s" % (table_prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_rowcol_pname = "%s_%s%s_%s_%s" % (rowcol_prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_trace_pname = "%s_%s%s_%s_%s" % (trace_prefix_prof,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            params["NAME"] = "parasail/parasail_"+function_name
            params["NAME_BASE"] = string.replace(params["NAME"], "_stats", "")
            params["NAME_TABLE"] = "parasail/parasail_"+function_table_name
            params["NAME_ROWCOL"] = "parasail/parasail_"+function_rowcol_name
            params["NAME_TRACE"] = "parasail/parasail_"+function_trace_name
            params["PNAME"] = "parasail/parasail_"+function_pname
            params["PNAME_BASE"] = string.replace(params["PNAME"], "_stats", "")
            params["PNAME_TABLE"] = "parasail/parasail_"+function_table_pname
            params["PNAME_ROWCOL"] = "parasail/parasail_"+function_rowcol_pname
            params["PNAME_TRACE"] = "parasail/parasail_"+function_trace_pname
            params = generated_params(template, params)
            params["VADD"] = params["VADDSx%d"%width]
            params["VSUB"] = params["VSUBSx%d"%width]
            output_filename = "%s%s.c" % (output_dir, function_name)
            result = template % params
            writer = open(output_filename, "w")
            writer.write(template % params)
            writer.write("\n")
            writer.close()

