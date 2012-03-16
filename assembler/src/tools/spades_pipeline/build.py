#!/usr/bin/env python

import os
import support
import shutil
import time
from os import path


def build_spades(dir):
    support.sys_call('cmake src', dir)
    print(dir)
    support.sys_call('make debruijn', dir)

def build_k(spades_folder, str_k, spades_home):

    build_folder_k = path.join(spades_folder, "build" + str_k)
    if os.path.exists(build_folder_k) :
        shutil.rmtree(build_folder_k)

    shutil.copytree(path.join(spades_home, "src"  ), path.join(build_folder_k, "src"  ))
    shutil.copytree(path.join(spades_home, "ext"  ), path.join(build_folder_k, "ext"  ))
    shutil.copy    (path.join(spades_home, "gen_k"), path.join(build_folder_k, "gen_k"))

    os.system(path.join(build_folder_k, "gen_k") + " " + str_k)

    print("\n== Compiling with K=" + str_k + " ==\n")
    build_spades(build_folder_k)

def build_spades_n_copy(cfg, spades_home):

    precompiled_folder = path.join(os.getenv('HOME'), '.spades/precompiled/')

    if path.exists(precompiled_folder):
        if os.path.getmtime(precompiled_folder) < os.path.getmtime(__file__):
            shutil.rmtree(precompiled_folder)

    print("\n== Compilation started ==\n")

    for K in cfg.iterative_K:

        if not os.path.exists(precompiled_folder) :
            os.makedirs(precompiled_folder)

        str_k = str(K)

        compiledFlag = path.join(precompiled_folder, "compiled") + str_k
        lockFlag     = path.join(precompiled_folder, "lock")     + str_k

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
                    build_k(precompiled_folder, str_k, spades_home)
                    compiledFile = open(compiledFlag, "w")
                    compiledFile.close()
                finally :
                    os.unlink(lockFlag)

    print("\n== Compilation finished successfully ==\n")