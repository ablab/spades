#!/usr/bin/env python3

############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os.path
from struct import Struct
import sys

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

#---- De Bruijn graph ----------------------------------------------------------
def show_grp(file):
    vertex_count = read_int(file)
    edge_count = read_int(file)
    print(vertex_count, edge_count)
    vertex = Struct("qq")
    for _ in range(vertex_count):
        id = read_int(file)
        conj_id = read_int(file)
        #id, conj_id = read_struct(file, vertex)
        print("Vertex", id, "~", conj_id, ".")
    print()
    edge = Struct("qqqqq")
    for _ in range(edge_count):
        id = read_int(file)
        conj_id = read_int(file)
        start = read_int(file)
        end = read_int(file)
        length = read_int(file)
        #id, conj_id, start, end, length = read_struct(file, edge)
        print("Edge", id, ":", start, "->", end, ", l =", length, "~", conj_id, ".")

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
ST_SIZE = 8
ST_NUC = ST_SIZE * 4

def unpack_seq(bytes, length):
    nucs = ['A', 'C', 'G', 'T']
    for byte in bytes:
        tmp = byte
        for _ in range(min(4, length)):
            yield nucs[tmp & 3]
            tmp = tmp >> 2
        length -= 4

def show_sqn(file):
    while True:
        try:
            id = read_int(file)
        except EOFError:
            break
        length = read_int(file, 8)
        print(">", id, "_length_", length, sep="")
        print("".join(unpack_seq(file.read((length + ST_NUC - 1) // ST_NUC * ST_SIZE), length)))

#-------------------------------------------------------------------------------
showers = {".grp": show_grp, ".prd": show_prd, ".sqn": show_sqn}

filename = sys.argv[1]
with open(filename, "rb") as file:
    showers[os.path.splitext(filename)[1]](file)
