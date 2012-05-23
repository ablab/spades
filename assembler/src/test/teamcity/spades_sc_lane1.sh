#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -f spades_output/ECOLI_SC_LANE_1_BH_woHUMAN/latest
./spades.py src/test/teamcity/spades_config_sc_lane1.info
python src/test/teamcity/assess.py spades_output/ECOLI_SC_LANE_1_BH_woHUMAN/latest/quality_results/all.tsv 80000 3
popd
