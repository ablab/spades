#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import urllib2
import spades_init

spades_init.init()
spades_version = spades_init.spades_version
spades_bin_dir = os.path.join(spades_init.spades_home, 'bin')

import support

print("\n======= Binaries download started.\n")
if not os.path.exists(spades_bin_dir):
    os.makedirs(spades_bin_dir)

print("\n======= BayesHammer download started.\n")
data = urllib2.urlopen('http://spades.bioinf.spbau.ru/release' + spades_version + '/hammer')
file = os.path.join(spades_bin_dir, 'hammer')
support.save_data_to_file(data, file)
print("\n======= BayesHammer download finished.\n")

print("\n======= SPAdes download started.\n")
data = urllib2.urlopen('http://spades.bioinf.spbau.ru/release' + spades_version + '/spades')
file = os.path.join(spades_bin_dir, 'spades')
support.save_data_to_file(data, file)
print("\n======= SPAdes download finished.\n")

print("\n======= BWA download started.\n")
data = urllib2.urlopen('http://spades.bioinf.spbau.ru/release' + spades_version + '/bwa-spades')
file = os.path.join(spades_bin_dir, 'bwa-spades')
support.save_data_to_file(data, file)
print("\n======= BWA download finished.\n")

print("\n======= Binaries download finished.\n")
