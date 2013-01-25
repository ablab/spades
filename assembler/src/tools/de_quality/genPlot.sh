#!/bin/sh

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

mkdir -p data/output/plot
java -cp build/ PlotFPR -s
