#!/usr/bin/env python
import os
import support
import shutil
import fcntl
import spades_init

spades_version = spades_init.spades_version

def build_spades(dir):
    if not os.path.exists(os.path.join(dir, "Makefile")):
        support.sys_call('cmake -G "Unix Makefiles" src', dir)

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
            destFile = os.path.join(dest, file)
            if os.path.isfile(destFile):
                if not os.path.exists(os.path.join(src, file)):
                    os.unlink(destFile)

        if os.path.exists(src): # False if src is broken link
            for file in os.listdir(src):
                syncFiles(os.path.join(src, file), os.path.join(dest, file))


def build_runtime_k(spades_folder, spades_home):
    build_folder = os.path.join(spades_folder, "build")

    syncFiles(os.path.join(spades_home, "src"), os.path.join(build_folder, "src"))
    syncFiles(os.path.join(spades_home, "ext"), os.path.join(build_folder, "ext"))

    print("\n== Compiling ==\n")
    build_spades(build_folder)


def build_spades_runtime_k(cfg, spades_home):
    precompiled_folder = cfg.compilation_dir

    print("\n== Compilation started ==\n")

    if not os.path.exists(precompiled_folder):
        os.makedirs(precompiled_folder)

    binary_file = os.path.join(precompiled_folder, 'release' + spades_version, 'bin', 'spades')

    if os.path.isfile(binary_file):
        dest = os.path.join(precompiled_folder, 'build', 'debruijn')
        if not os.path.exists(dest):
            os.makedirs(dest)
        shutil.copy2(binary_file, dest)
        print("Downloaded SPAdes binary is used instead of compilation.\n")
        print("\n== Compilation finished successfully ==\n")
        return

    lockFlag = os.path.join(precompiled_folder, "lock")

    fo = open(lockFlag, "w")
    fcntl.lockf(fo, fcntl.LOCK_EX)
    try:
        build_runtime_k(precompiled_folder, spades_home)
    finally:
        fcntl.lockf(fo, fcntl.LOCK_UN)

    print("\n== Compilation finished successfully ==\n")


def build_hammer(cfg, spades_home):
    precompiled_folder = cfg.compilation_dir

    print("\n== Compilation started ==\n")

    if not os.path.exists(precompiled_folder):
        os.makedirs(precompiled_folder)

    binary_file = os.path.join(precompiled_folder, 'release' + spades_version, 'bayeshammer', 'hammer')
    if os.path.isfile(binary_file):
        dest = os.path.join(precompiled_folder, 'build_hammer', 'hammer')
        if not os.path.exists(dest):
            os.makedirs(dest)
        shutil.copy2(binary_file, dest)
        print("BayesHammer downloaded binary used instead of compilation.\n")
        return

    lockFlag = os.path.join(precompiled_folder, "lock_hammer")

    fo = open(lockFlag, "w")
    fcntl.lockf(fo, fcntl.LOCK_EX)
    try:
        build_folder = os.path.join(precompiled_folder, "build_hammer")

        syncFiles(os.path.join(spades_home, "src"), os.path.join(build_folder, "src"))
        syncFiles(os.path.join(spades_home, "ext"), os.path.join(build_folder, "ext"))

        if not os.path.exists(os.path.join(build_folder, "Makefile")):
            support.sys_call('cmake -G "Unix Makefiles" src', build_folder)

        support.sys_call('make hammer', build_folder)
    finally:
        fcntl.lockf(fo, fcntl.LOCK_UN)

    print("\n== Compilation finished successfully ==\n")
