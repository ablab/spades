#!/bin/bash
./generate_configs.sh ECOLI_IS220_BH_iter4.info ECOLI_SC_LANE_1_BH_woHUMAN.info SAUREUS_SC_LANE_7_BH.info ECOLI_JGI_0012_BH.info HMP_LANE_3_BH.info GEBA_1_BH.info
cd _generic
echo ";;;;;;;;;;;;;; only K = 55 ;;;;;;;;;;;;;;;;;" >> ECOLI_JGI_0012_BH.info
mcedit ECOLI_JGI_0012_BH.info
for f in *
do
    mv $f $f.dataset
done
