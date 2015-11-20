############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

__author__ = 'anton'
import sys

def calculate_dist_table(s1, s2):
    n1 = len(s1)
    n2 = len(s2)
    t = []
    t.append(range(n2 + 1))
    for i in range(1, n1 + 1):
        t_line = [i]
        for j in range(1, n2 + 1):
            if s1[i - 1] == s2[j - 1]:
                t_line.append(t[i - 1][j - 1])
            else:
                t_line.append(1 + min(t[i - 1][j - 1], t_line[j - 1], t[i - 1][j]))
        t.append(t_line)
    return t

def calculate_lcs_table(s1, s2):
    n1 = len(s1)
    n2 = len(s2)
    t = []
    t.append([0] * (n2 + 1))
    for i in range(1, n1 + 1):
        t_line = [0]
        for j in range(1, n2 + 1):
            if s1[i - 1] == s2[j - 1]:
                t_line.append(t[i - 1][j - 1] + 1)
            else:
                t_line.append(max(t[i - 1][j - 1], t_line[j - 1], t[i - 1][j]))
            t.append(t_line)
    return t

def lcs(s1, s2):
    t = calculate_dist_table(s1, s2)
    i = len(s1)
    j = len(s2)
    res = ""
    while i > 0 and j > 0:
        if t[i][j] == t[i - 1][j -1] + 1:
            i -= 1
            j -= 1
        elif t[i][j] == t[i - 1][j] + 1:
            i -= 1
        elif t[i][j] == t[i][j - 1] + 1:
            j -= 1
        else:
            i -= 1
            j -= 1
            res = s1[i] + res
    return res

def dist(s1, s2):
    return calculate_dist_table(s1, s2)[len(s1)][len(s2)]


def multi_lcs(ids):
    all_lcs = ids[0]
    for s in ids:
        all_lcs = lcs(all_lcs, s)
    return all_lcs
