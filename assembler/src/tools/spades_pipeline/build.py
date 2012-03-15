#!/usr/bin/env python

import os
import support

def build_spades(dir):
    support.sys_call('cmake src', dir)
    support.sys_call('make', dir)

def build_k(spades_folder, k):

    import shutil

    build_folder_k = spades_folder + "build" + k + "/"
    if os.path.exists(build_folder_k) :
        shutil.rmtree(build_folder_k)

    shutil.copytree("/usr/share/spades/src", build_folder_k + "src")
    shutil.copytree("/usr/share/spades/ext", build_folder_k + "ext")
    kFile = build_folder_k + "src/debruijn/k.hpp"
    fo = open(kFile, "w")
    fo.write("#pragma once\n\n")
    fo.write("namespace debruijn_graph {\n")
    fo.write("	const size_t K = " + k + ";\n")
    fo.write("}\n")
    fo.close()

    print("\n== Compiling with K=" + k + " ==\n")
    build_spades(build_folder_k)


def build_spades_n_copy(cfg):

    import time

    print("\n== Compilation started. Don't start another instance of SPAdes before compilation ends ==\n")

    for K in cfg.iterative_K:
        precompiled_folder = os.getenv('HOME') + '/.spades/precompiled/'

        if not os.path.exists(precompiled_folder) :
            os.makedirs(precompiled_folder)

        k = str(K)

        compiledFlag = precompiled_folder + "compiled" + k
        lockFlag = precompiled_folder + "lock" + k

        if not os.path.exists(compiledFlag) :
            while os.path.exists(lockFlag) :
                try:
                    time.sleep(10)
                    os.unlink(lockFile)
                except:
                    print("Failed to delete " + lockFlag)
            if not os.path.exists(compiledFlag) :
                open(lockFlag, "w")
                try :
                    build_k(precompiled_folder, k)
                    compiledFile = open(compiledFlag, "w")
                    compiledFile.close()
                finally :
                    os.unlink(lockFlag)

    print("\n== Compilation finished successfully. Feel free to start another instance of SPAdes ==\n")