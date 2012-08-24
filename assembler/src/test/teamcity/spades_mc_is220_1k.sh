#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf spades_output/ECOLI_IS220_QUAKE_1K
./spades.py --config-file src/test/teamcity/spades_config_mc_is220_1k.info
python src/test/teamcity/assess.py spades_output/ECOLI_IS220_QUAKE_1K/quality_results/transposed_report.tsv 1000 0
popd
