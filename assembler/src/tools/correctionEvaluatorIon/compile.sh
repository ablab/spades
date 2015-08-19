#!/usr/bin/env bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

wget https://dl.dropbox.com/u/7916095/const/dmd_2.062_linux_64.tar.bz2
tar xjvf dmd_2.062_linux_64.tar.bz2
export PATH=`pwd`/dmd2/linux/bin64:$PATH
git clone https://github.com/biod/BioD
cd BioD
git checkout d109ffc38
cd ..
make -j5
wget https://dl.dropbox.com/u/7916095/const/sambamba_10_03_2013.bz2
bunzip2 sambamba_10_03_2013.bz2
mv sambamba_10_03_2013 sambamba
chmod +x sambamba
