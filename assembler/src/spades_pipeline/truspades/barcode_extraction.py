############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os.path
import sys
import logging

from id_generation import generate_ids
from string_dist_utils import lcs, dist


__author__ = 'anton'


class Barcode:
    def __init__(self, id, libs):
        self.id = id
        self.libs = list(libs)
        if id == None:
            self.set_lcs_id()

    def add_ps(self, prefix, suffix):
        for lib in self.libs:
            for i in range(len(lib)):
                lib[i] = os.path.abspath(prefix + lib[i] + suffix)

    def __str__(self):
        return self.id + " " + " ".join([" ".join(lib) for lib in self.libs])

def RemoveLabel(s, code, code_range):
    for pos in range(len(s)):
        if s[pos:].startswith(code):
            for i in code_range:
                new_pos = pos + len(code)
                tmp = str(i)
                if new_pos + len(tmp) <= len(s) and s[new_pos:].startswith(tmp):
                    return s[:pos] + s[new_pos + len(tmp):]
    return s

def NormalizeR(s):
    return RemoveLabel(s, "R", [1,2])

def NormalizeLR(s):
    s = NormalizeR(s)
    return RemoveLabel(s, "L", range(1, 20))

def check_int_ids(ids):
    for id in ids:
        if not id[1].isdigit():
            return False
    return True

def generate_barcode_list(barcodes):
    ids = list(zip(barcodes, generate_ids(barcodes)))
    if check_int_ids(ids):
        ids = sorted(ids, key=lambda barcode: int(barcode[1]))
    return [(bid, "BC_" + short_id) for bid, short_id in ids]

def Normalize(file_path):
    return NormalizeLR(os.path.basename(file_path))

def GroupBy(norm, l):
    result = dict()
    for line in l:
        key = norm(line)
        if not key in result:
            result[key] = []
        result[key].append(line)
    return result

def CheckSameSize(iter, size = -1):
    for vl in iter:
        if size == -1:
            size = len(vl)
        if size != len(vl):
            return False
    return True

#todo: write better code
def ExtractBarcodes(dirs):
    files = []
    for dir in dirs:
        for file in [os.path.abspath(os.path.join(dir, file)) for file in os.listdir(dir) if os.path.isfile(os.path.join(dir, file))]:
            files.append(file)
    barcode_dict = GroupBy(Normalize, files)
    if not CheckSameSize(barcode_dict.values()):
        return None
    for bid in barcode_dict.keys():
        barcode_dict[bid] = GroupBy(NormalizeR, barcode_dict[bid]).values()
        if not CheckSameSize(barcode_dict[bid], 2):
            return None
    short_barcodes = generate_barcode_list(list(barcode_dict.keys()))
    return [Barcode(short, barcode_dict[bid]) for bid, short in short_barcodes]

def ReadDataset(file, log = logging.getLogger("ReadDataset")):
    log.info("Reading dataset from " + file + "\n")
    if os.path.exists(file) and os.path.isfile(file):
        result = []
        f = open(file, "r")
        lines = f.xreadlines()
        for line in lines:
            line = line.strip()
            if line == "":
                continue
            split = line.split()
            id = split[0]
            datasets = []
            for i in range(1, len(split), 2):
                datasets.append([split[i], split[i + 1]])
            result.append(Barcode(id, datasets))
        f.close()
        return result
    else:
        log.info("Error: Dataset file does not exist\n" + file + "\n")
        sys.exit(1)

def print_dataset(dataset, output_file, log):
    log.info("Printing dataset to " + output_file)
    open(output_file, "w").write("\n".join([str(line).strip() for line in dataset]) + "\n")

