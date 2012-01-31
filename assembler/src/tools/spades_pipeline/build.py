#!/usr/bin/env python

import os
import support

def build_spades(cfg, K):

    support.sys_call("./gen_k " + str(K))
    support.sys_call("make")


def build_spades_n_copy(cfg):

    import datetime
    import shutil

    print("\n== Compilation started. Don't start another instance of SPAdes before compilation ends ==\n")

    cfg_folder = "configs/debruijn"

    for K in cfg.iterative_K:
        build_folder_k = cfg.build_path + str(K) + r"/"

        os.makedirs(build_folder_k + "configs", 0755)
        shutil.copytree(cfg_folder, build_folder_k + cfg_folder)

        print("\n== Compiling with K=" + str(K) + " ==\n")
        build_spades(cfg, K)
        shutil.copy("build/release/debruijn/debruijn", build_folder_k)

    print("\n== Compilation finished successfully. Feel free to start another instance of SPAdes ==\n")
