#!/usr/bin/env python
import os
import os.path
import sys
import shutil

import support
import build

def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('Usage prebuild_spades.py <path to copy to> <from> <to>')
        exit(1)

    safe_mkdir('build')
    safe_mkdir('build/static')
    
    support.sys_call('cmake ../../src -Dstatic_build=\"\"', './build/static')

    prebuild_folder = sys.argv[1]
    safe_mkdir(prebuild_folder)

    from_k = int(sys.argv[2])
    to_k   = int(sys.argv[3]) + 2

    for i in range(from_k, to_k, 2):

        print('\n=== Buiding SPAdes with K ' + str(i) + ' ===\n')

        build.write_k_file(os.path.join(os.getcwd(), "src/debruijn/k.hpp"), i)
        support.sys_call('make -C build/static/debruijn')

        new_folder = os.path.join(prebuild_folder, 'K') + str(i)
        safe_mkdir(new_folder)

        shutil.copy2('build/static/debruijn/spades', new_folder)
