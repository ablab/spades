#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -f spades_output/ECOLI_IS220_QUAKE/latest
./spades.py src/test/teamcity/spades_config_mc_is220.info
python src/test/teamcity/assess.py spades_output/ECOLI_IS220_QUAKE/latest/quality_results/all.tsv 80000 3
popd
