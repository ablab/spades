from string_dist_utils import multi_lcs

__author__ = 'anton'

import sys

def generate_ids(lines):
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
