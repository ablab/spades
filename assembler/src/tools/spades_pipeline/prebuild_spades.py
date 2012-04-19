#!/usr/bin/env python
import os
import os.path
import sys
import shutil

import support

def write_k_file(kFile, k):
    fo = open(kFile, "w")
    fo.write("#pragma once\n\n")
    fo.write("namespace debruijn_graph {\n")
    fo.write("  const size_t K = " + str(k) + ";\n")
    fo.write("}\n")
    fo.close()

def safe_mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('Usage prebuild_spades.py <from> <to>')
        exit(1)

    safe_mkdir('build/static')
    support.sys_call('cmake ../../src -Dstatic_build=\"\"', './build/static')

    prebuild_folder = 'prebuild_spades'
    safe_mkdir(prebuild_folder)

    #building SPAdes

    from_k = int(sys.argv[1])
    to_k   = int(sys.argv[2]) + 1

    for i in range(from_k, to_k, 2):

        print('\n=== Buiding SPAdes with K ' + str(i) + ' ===\n')

        write_k_file(os.path.join(os.getcwd(), "src/debruijn/k.hpp"), i)
        support.sys_call('make -C build/static/debruijn')

        new_folder = os.path.join(prebuild_folder, 'bin/K') + str(i)
        safe_mkdir(new_folder)

        shutil.copy2('build/static/debruijn/spades', new_folder)

    #building BayesHammer
    support.sys_call('make -C build/static/hammer hammer')
    new_folder = os.path.join(prebuild_folder, 'bayeshammer')
    safe_mkdir(new_folder)
    shutil.copy2('build/static/hammer/hammer', new_folder)




