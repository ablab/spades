############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from string_dist_utils import multi_lcs

__author__ = 'anton'

import sys

def CommonPrefix(s1, s2):
    n = 0
    while n < len(s1) and n < len(s2) and s1[n] == s2[n]:
        n += 1
    return n

def CommonSuffix(s1, s2):
    n = 0
    while n < len(s1) and n < len(s2) and s1[-n - 1] == s2[-n - 1]:
        n += 1
    return n

def FindCommon(lines):
    if len(lines) == 0:
        return 0, 0
    left = len(lines[0])
    right = len(lines[0])
    min_len = len(lines[0])
    for line in lines:
        l, r = CommonPrefix(line, lines[0]), CommonSuffix(line, lines[0])
        left = min(left, l)
        right = min(right, r)
        min_len = min(min_len, len(line))
    return left, min(right, min_len - left)

def generate_ids(lines):
    l, r = FindCommon(lines)
    lines = [line[l: len(line) - r] for line in lines]
    id_candidates = generate_id_candidates(lines)
    return select_ids_from_candidates(id_candidates)

def select_ids_from_candidates(id_candidates):
    if len(id_candidates) == 1:
        return [""]
    ids = [""] * len(id_candidates)
    for i in range(len(id_candidates[0])):
        for bcid in range(len(id_candidates)):
            ids[bcid] += id_candidates[bcid][i]
        if len(set(ids)) == len(ids):
            return ids
    return None


def generate_id_candidates(lines):
    all_lcs = multi_lcs(lines)
    id_candidates = []
    for line in lines:
        id_candidates.append(generate_id_candidates_for_barcode(all_lcs, line))
    return id_candidates


def generate_id_candidates_for_barcode(all_lcs, line):
    id_candidate = []
    cur = ""
    pos = 0
    for i in range(len(line)):
        if pos < len(all_lcs) and line[i] == all_lcs[pos]:
            pos += 1
            if cur != "":
                id_candidate.append(cur)
                cur = ""
        else:
            cur = cur + line[i]
    if cur != "":
        id_candidate.append(cur)
    return id_candidate
