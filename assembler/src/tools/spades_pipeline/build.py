#!/usr/bin/env python
import os
import support
import shutil
import fcntl


def build_spades(dir):
    if not os.path.exists(os.path.join(dir, "Makefile")):
        support.sys_call('cmake src', dir)

    support.sys_call('make spades', dir)

def syncFiles(src, dest):
    if os.path.isfile(src):
        if not os.path.exists(dest):
            shutil.copy2(src, dest)
        elif str(os.path.getmtime(src)) != str(os.path.getmtime(dest)):
            shutil.copy2(src, dest)
    else:
        if not os.path.exists(dest):
            os.makedirs(dest)
        for file in os.listdir(dest):
            if file != "k.hpp":
                destFile = os.path.join(dest, file)
                if os.path.isfile(destFile):
                    if not os.path.exists(os.path.join(src, file)):
                        os.unlink(destFile)

        if os.path.exists(src): # False if src is broken link
            for file in os.listdir(src):
                if file != "k.hpp":
                    syncFiles(os.path.join(src, file), os.path.join(dest, file))

def kFile_required(kFile, str_k):
    if not os.path.exists(kFile):
        return True
    input = open(kFile, "r")
    try:
        for line in input:
            if line.startswith("  const size_t K = " + str_k + ";"):
                return False
    finally:
        input.close()

    return True

def build_k(spades_folder, str_k, spades_home):

    build_folder_k = os.path.join(spades_folder, "build" + str_k)

    syncFiles(os.path.join(spades_home, "src"), os.path.join(build_folder_k, "src"))
    syncFiles(os.path.join(spades_home, "ext"), os.path.join(build_folder_k, "ext"))
    shutil.copy2(os.path.join(spades_home, "log4cxx.properties"), os.path.join(build_folder_k, "log4cxx.properties"))

    kFile = os.path.join(build_folder_k, "src/debruijn/k.hpp")
    if kFile_required(kFile, str_k):
        fo = open(kFile, "w")
        fo.write("#pragma once\n\n")
        fo.write("namespace debruijn_graph {\n")
        fo.write("  const size_t K = " + str_k + ";\n")
        fo.write("}\n")
        fo.close()

    print("\n== Compiling with K=" + str_k + " ==\n")
    build_spades(build_folder_k)

def build_spades_n_copy(cfg, spades_home):

    precompiled_folder = cfg.compilation_dir

    print("\n== Compilation started ==\n")

    for K in cfg.iterative_K:

        if not os.path.exists(precompiled_folder) :
            os.makedirs(precompiled_folder)

        str_k = str(K)

        lockFlag = os.path.join(precompiled_folder, "lock") + str_k

        fo = open(lockFlag, "w")
        fcntl.lockf(fo, fcntl.LOCK_EX)
        try :
            build_k(precompiled_folder, str_k, spades_home)
        finally :
            fcntl.lockf(fo, fcntl.LOCK_UN)

    print("\n== Compilation finished successfully ==\n")

def build_hammer(cfg, spades_home):

    precompiled_folder = cfg.compilation_dir

    print("\n== Compilation started ==\n")

    if not os.path.exists(precompiled_folder):
        os.makedirs(precompiled_folder)

    lockFlag = os.path.join(precompiled_folder, "lock_hammer")

    fo = open(lockFlag, "w")
    fcntl.lockf(fo, fcntl.LOCK_EX)
    try :
        build_folder = os.path.join(precompiled_folder, "build_hammer")

        syncFiles(os.path.join(spades_home, "src"), os.path.join(build_folder, "src"))
        syncFiles(os.path.join(spades_home, "ext"), os.path.join(build_folder, "ext"))
        
        if not os.path.exists(os.path.join(build_folder, "Makefile")):
            support.sys_call('cmake src', build_folder)

        support.sys_call('make hammer', build_folder)
    finally :
        fcntl.lockf(fo, fcntl.LOCK_UN)

    print("\n== Compilation finished successfully ==\n")
