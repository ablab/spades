
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#!/bin/sh

PBSIMPATH=
REF=

$PBSIMPATH/src/pbsim $REF --model_qc $PBSIMPATH/data/model_qc_clr --seed 115249 --prefix sim_pacbio
