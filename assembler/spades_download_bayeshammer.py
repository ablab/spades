#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import urllib2
import spades_init

spades_init.init()
spades_version = spades_init.spades_version
spades_build_dir = spades_init.spades_build_dir

import support

print("\n======= BayesHammer download started.\n")
data = urllib2.urlopen('http://spades.bioinf.spbau.ru/release' + spades_version + '/bayeshammer/hammer')
dir = os.path.join(spades_build_dir, 'release' + spades_version, 'bayeshammer')
if not os.path.exists(dir):
    os.makedirs(dir)
file = os.path.join(dir, 'hammer')
support.save_data_to_file(data, file)
print("\n======= BayesHammer download finished.\n")

