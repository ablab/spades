import os.path
import sys

from id_generation import generate_ids
from string_dist_utils import lcs, dist


__author__ = 'anton'


class Barcode:
    def __init__(self, libs, id = None):
        self.id = id
        self.libs = libs
        if id == None:
            self.set_lcs_id()

    def set_lcs_id(self):
        result = None
        for lib in self.libs:
            for reads in lib:
                if None != result:
                    result = lcs(result, reads)
                else:
                    result = reads
        self.id = result

    def add_ps(self, prefix, suffix):
        for lib in self.libs:
            for i in range(len(lib)):
                lib[i] = os.path.abspath(prefix + lib[i] + suffix)

    def __str__(self):
        return self.id + " " + " ".join([" ".join(lib) for lib in self.libs])

def calculate_best_dist(files):
    best_dist = []
    n = len(files)
    for i in range(n):
        best = 1000
        for j in range(n):
            if i != j:
                tmp = dist(files[i], files[j])
                if best == -1 or tmp < best:
                    best = tmp
        best_dist.append(best)
    return best_dist

def normalizeR(s):
    pos = s.find("_R")
    if pos == -1 or len(s) <= pos + 2 or (s[pos + 2] != '1' and s[pos + 2] != '2'):
        return s
    return s[:pos] + s[pos + 3:]

def normalizeLR(s):
    s1 = s
    s = normalizeR(s)
    s2 = s
    pos = s.find("_L")
    if pos != -1:
        rpos = pos + 2
        while rpos < len(s) and s[rpos].isdigit():
            rpos += 1
        s = s[:pos] + s[rpos:]
    return s

def generate_dataset(barcode):
    datasets = dict()
    for file in barcode:
        datasets[normalizeR(file)] = []
    for file in barcode:
        datasets[normalizeR(file)].append(file)
    return Barcode([value for key, value in datasets.iteritems()])

def find_barcodesLR(files):
    barcodes = list()
    for file in files:
        found = False
        for i in range(len(barcodes)):
            if normalizeLR(barcodes[i][0]) == normalizeLR(file):
                found = True
                barcodes[i].append(file)
        if not found:
            barcodes.append([file])
    tmp = 0
    for barcode in barcodes:
        if tmp == 0:
            tmp = len(barcode)
        else:
            if tmp != len(barcode):
                return None
    if tmp % 2 == 1:
        return None
    return [generate_dataset(barcode) for barcode in barcodes]


def distance_based_pairing(files):
    best_dist = calculate_best_dist(files)
    pairs = []
    paired = set()
    n = len(files)
    for i in range(n):
        if i not in paired:
            for j in range(n):
                if i != j and j not in paired:
                    tmp = dist(files[i], files[j])
                    if tmp == best_dist[i] and tmp == best_dist[j]:
                        paired.add(i)
                        paired.add(j)
                        pairs.append(Barcode([[files[i],files[j]]]))
    single = [files[i] for i in range(n) if i not in paired]
    if len(single) != 0:
        return None
    return pairs

def check_int_ids(ids):
    for id in ids:
        if not id.isdigit():
            return False
    return True

def generate_barcode_ids(barcodes):
    ids = generate_ids([barcode.id for barcode in barcodes])
    for bcid in range(len(barcodes)):
        barcodes[bcid].id = "BC_" + ids[bcid]
    if check_int_ids(ids):
        return sorted(barcodes, key=lambda barcode: int(barcode.id[3:]))
    return barcodes

def prepare_barcodes(result, prefix, suffix):
    result = generate_barcode_ids(result)
    for barcode in result:
        barcode.add_ps(prefix, suffix)
    return result

def common_prefix(s1, s2):
    n = 0
    while n < len(s1) and n < len(s2) and s1[n] == s2[n]:
        n += 1
    return n

def common_suffix(s1, s2):
    n = 0
    while n < len(s1) and n < len(s2) and s1[-n - 1] == s2[-n - 1]:
        n += 1
    return n

def cut_common(files):
    prefix_len = len(files[0])
    for f in files:
        prefix_len = min(prefix_len, common_prefix(f, files[0]))
    prefix = files[0][:prefix_len]
    files = [f[prefix_len:] for f in files]
    suffix_len = len(files[0])
    for f in files:
        suffix_len = min(suffix_len, common_suffix(f, files[0]))
    suffix = files[0][-suffix_len:]
    return [f[:-suffix_len] for f in files], prefix, suffix

#todo: write better code
def extract_barcodes(files):
    files, prefix, suffix = cut_common(files)
    result = find_barcodesLR(files)
    if None != result:
        return prepare_barcodes(result, prefix, suffix)
#    result = distance_based_pairing(files)
#    if None != result:
#        return prepare_barcodes(result, prefix, suffix)
    return None
