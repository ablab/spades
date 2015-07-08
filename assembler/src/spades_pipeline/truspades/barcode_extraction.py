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

def generate_dataset(bid, barcode):
    datasets = dict()
    for file in barcode:
        datasets[normalizeR(file)] = []
    for file in barcode:
        datasets[normalizeR(file)].append(file)
    return Barcode([value for key, value in datasets.iteritems()], bid)

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

def find_barcodesLR(files):
    barcodes = dict()
    for file in files:
        barcodes[normalizeLR(file)] = []
    for file in files:
        barcodes[normalizeLR(file)].append(file)
    prefix, suffix = FindCommon([bid for bid, files in barcodes.iteritems()])
    tmp = 0
    for bid, barcode in barcodes.iteritems():
        if tmp == 0:
            tmp = len(barcode)
        else:
            if tmp != len(barcode):
                return None
    if tmp % 2 == 1:
        return None
    return [generate_dataset(bid[prefix: len(bid) - suffix], barcode) for bid, barcode in barcodes.iteritems()]


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

#todo: write better code
def extract_barcodes(files):
    result = find_barcodesLR(files)
    if None != result:
        return generate_barcode_ids(result)
#    result = distance_based_pairing(files)
#    if None != result:
#        return prepare_barcodes(result, prefix, suffix)
    return None

def ReadDataset(file):
    sys.stdout.write("Reading dataset from " + file + "\n")
    if os.path.exists(file) and os.path.isfile(file):
        result = []
        f = open(file, "r")
        lines = f.xreadlines()
        for line in lines:
            line = line.strip()
            if line == "":
                continue
            split = line.split(" ")
            id = split[0]
            datasets = []
            for i in range(1, len(split), 2):
                datasets.append([split[i], split[i + 1]])
            result.append(Barcode(datasets, id))
        f.close()
        return result
    else:
        sys.stderr.write("Error: Dataset file does not exist\n" + file + "\n")
        sys.exit(1)

def print_dataset(dataset, output_file):
    sys.stdout.write("Printing dataset to " + output_file + "\n")
    open(output_file, "w").write("\n".join([str(line).strip() for line in dataset]) + "\n")

