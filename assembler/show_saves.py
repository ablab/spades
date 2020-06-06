#!/usr/bin/env python3

############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os.path
from struct import Struct
import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

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

def show_kmidx(file):
    k = read_int(file)
    print("k: %d" % k)
    size = read_int(file, 8)
    print("Size: %d" % size)
    for _ in range(size):
        id, count, offset = read_int(file, 8), read_int(file, 4), read_int(file, 4)
        print(id, offset, count)


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
    _, _ = read_int(file), read_int(file) # max_vid, max_eid

    vertex_cnt = read_int(file)
    while True:
        try:
            start = read_int(file)
        except EOFError:
            break
        start_conj = read_int(file)
        while True:
            edge = read_int(file)
            if not edge: #null-term
                break
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
        while True:
            e2 = read_int(file)
            if not e2: #null-term
                break
            hist_size = read_int(file)
            for _ in range(hist_size):
                p = read_struct(file, point)
                print(e1, e2, *p, ".")

#---- Edge sequences -----------------------------------------------------------
def show_sqn(file):
    show_grp(file, True)

#---- Edge coverage ------------------------------------------------------------
def show_cvr(file):
    while True:
        edge = read_int(file)
        if not edge: # null-term
            break
        coverage = read_int(file)
        print(edge, coverage, ".")

#---- Long reads paths ---------------------------------------------------------
def show_mpr(file):
    size = read_int(file)
    for _ in range(size):
        count = read_int(file)
        print(count)
        for _ in range(count):
            weight = read_int(file)
            length = read_int(file)
            print(" Weight:", weight, "length:", length, end=' ')
            for _ in range(length):
                edge = read_int(file)
                print(edge, end=' ')
            print()
        print()

#-------------------------------------------------------------------------------
showers = {".grp" : show_grp,
           ".prd" : show_prd,
           ".sqn" : show_sqn,
           ".cvr" : show_cvr,
           ".flcvr" : show_cvr,
           ".mpr" : show_mpr,
           ".kmidx" : show_kmidx}

basename, ext = os.path.splitext(sys.argv[1])
target = ext
if ext in [".grp", ".sqn"]:
    target = ".grseq"
with open(basename + target, "rb") as file:
    showers[ext](file)
