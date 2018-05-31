#!/usr/bin/env python3

############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os.path
from struct import Struct
import sys

from Bio.Seq import Seq

def read_int(file, size=None, signed=False):
    if size:
        bytes = file.read(size)
        if len(bytes) < size:
            raise EOFError
        return int.from_bytes(bytes, byteorder=sys.byteorder, signed=signed)

    byte = 0xFF
    result = 0
    shift = 0
    while byte & 0x80:
        bytes = file.read(1)
        if not len(bytes):
            raise EOFError
        byte = bytes[0]
        result = result | ((byte & 0x7F) << shift)
        shift += 7
    return result

def read_struct(file, struct):
    return struct.unpack(file.read(struct.size))

ST_SIZE = 8
ST_NUC = ST_SIZE * 4

def read_seq(file):
    def unpack_seq(bytes, length):
        nucs = ['A', 'C', 'G', 'T']
        for byte in bytes:
            tmp = byte
            for _ in range(min(4, length)):
                yield nucs[tmp & 3]
                tmp = tmp >> 2
            length -= 4
    length = read_int(file, 8)
    seq = "".join(unpack_seq(file.read((length + ST_NUC - 1) // ST_NUC * ST_SIZE), length))
    return Seq(seq)

#---- De Bruijn graph ----------------------------------------------------------
def show_grp(file, show_seq=False):
    _ = read_int(file) #max_id
    while True:
        try:
            start = read_int(file)
        except EOFError:
            break
        start_conj = read_int(file)
        count = read_int(file)
        for _ in range(count):
            edge = read_int(file)
            edge_conj = read_int(file)
            end = read_int(file)
            end_conj = read_int(file)
            seq = read_seq(file)

            if show_seq:
                print(">", edge, sep="")
                print(seq)
                if edge != edge_conj:
                    print(">", edge_conj, sep="")
                    print(seq.reverse_complement())
            else:
                print("Edge", edge, ":", start, "->", end, ", l =", len(seq), "~", edge_conj, ".")
                print("Edge", edge_conj, ":", end_conj, "->", start_conj, ", l =", len(seq), "~", edge, ".")

#---- Paired info --------------------------------------------------------------
def show_prd(file, clustered=False):
    size = read_int(file)
    point = Struct("fff" if clustered else "ff")
    for _ in range(size):
        e1 = read_int(file)
        inner_size = read_int(file)
        for _ in range(inner_size):
            e2 = read_int(file)
            hist_size = read_int(file)
            for _ in range(hist_size):
                p = read_struct(file, point)
                print(e1, e2, *p, ".")

#---- Edge sequences -----------------------------------------------------------
def show_sqn_old(file):
    while True:
        seq = read_seq(file)
        print(">", seq.id, "_length_", seq.length, sep="")
        print(seq.seq)

def show_sqn(file):
    show_grp(file, True)

#-------------------------------------------------------------------------------
showers = {".grp": show_grp, ".prd": show_prd, ".sqn": show_sqn}

basename, ext = os.path.splitext(sys.argv[1])
target = ext
if ext == ".sqn":
    target = ".grp"
with open(basename + target, "rb") as file:
    showers[ext](file)
