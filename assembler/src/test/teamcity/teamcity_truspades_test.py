#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import os
import shutil
import sys


def run_test():
    shutil.rmtree('bin', True)
    shutil.rmtree('build', True)
    shutil.rmtree('build_spades', True)
    ecode = os.system('./prepare_cfg')
    if ecode != 0:
        print("Preparing configuration files finished abnormally with exit code " + str(ecode))
        sys.exit(2)
    ecode = os.system('./spades_compile.sh')
    if ecode != 0:
        print("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)
    cmd = "./truspades.py --test"
    ecode = os.system(cmd)
    output_dir = "spades_test_truspades"
    os.system("chmod -R 777 " + output_dir)
    if ecode != 0:
        print("SPAdes finished abnormally with exit code " + str(ecode))
        sys.exit(4)
    log = open(os.path.join(output_dir, "truspades.log")).readlines()
    if log[-1].strip() != "TruSPAdes test passed correctly":
        sys.exit(5)

if __name__ == "__main__":
    run_test()
