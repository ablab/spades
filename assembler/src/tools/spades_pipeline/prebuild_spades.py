#!/usr/bin/env python
import os
import shutil

import support


def safe_mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

if __name__ == "__main__":

    safe_mkdir('build/static')
    support.sys_call('cmake ../../src -Dstatic_build=\"\"', './build/static')

    prebuild_folder = 'prebuild_spades'

    safe_mkdir(prebuild_folder)

    #building SPAdes

    print('\n=== Buiding SPAdes with ===\n')

    support.sys_call('make -C build/static/debruijn')

    safe_mkdir(prebuild_folder)

    shutil.copy2('build/static/debruijn/spades', prebuild_folder)

    #building BayesHammer
    support.sys_call('make -C build/static/hammer hammer')
    shutil.copy2('build/static/hammer/hammer', prebuild_folder)




