#!/bin/bash

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
if [ ! -e input ];
then
    ln -s /tmp/data/input input
fi
