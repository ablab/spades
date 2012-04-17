#!/usr/bin/env python

import os
import sys
import urllib2
import spades_init

spades_init.init()
spades_version = spades_init.spades_version

import support

for k in sys.argv[1:]:
    print("\n======= Binary download for " + k +  " started.\n")
    data = urllib2.urlopen('http://spades.bioinf.spbau.ru/release' + spades_version + '/bin/K' + k + '/spades')
    dir = os.path.join(os.getenv('HOME'), '.spades', 'release' + spades_version, 'bin', 'K' + k)
    if not os.path.exists(dir):
        os.makedirs(dir)
    file = os.path.join(dir, 'spades')
    support.save_data_to_file(data, file)
    print("\n======= Binary download for " + k +  " finished.\n")
