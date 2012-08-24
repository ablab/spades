#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf spades_output/ECOLI_SC_LANE_1_BH_woHUMAN
./spades.py --config-file src/test/teamcity/spades_config_sc_lane1.info
python src/test/teamcity/assess.py spades_output/ECOLI_SC_LANE_1_BH_woHUMAN/quality_results/transposed_report.tsv 85000 6
popd
