#!/usr/bin/env python
import os
import os.path
import sys
import shutil
import support

def write_k_file(k):
    kFile = os.path.join(os.getcwd(), "src/debruijn/k.hpp")
    fo = open(kFile, "w")
    fo.write("#pragma once\n\n")
    fo.write("namespace debruijn_graph {\n")
    fo.write("  const size_t K = " + str(k) + ";\n")
    fo.write("}\n")
    fo.close()

def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('Usage prebuild_spades.py <from> <to>')
        exit(1)

    prebuild_folder = './prebuild_spades'
    safe_mkdir(prebuild_folder)

    for i in range(int(sys.argv[1]), int(sys.argv[2]), 2):

        print('\n=== Buiding SPAdes with K ' + str(i) + ' ===\n')

        write_k_file(i)
        support.sys_call('make')

        new_folder = os.path.join(prebuild_folder, 'K') + str(i)
        safe_mkdir(new_folder)

        shutil.copy2('build/release/debruijn/spades', new_folder)
