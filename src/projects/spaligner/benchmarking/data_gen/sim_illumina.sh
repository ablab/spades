
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#!/bin/sh

ARTSIMPATH=
REF=

$ARTSIMPATH/art_illumina -ss HS25 -sam -i $REF -p -l 150 -f 20 -m 200 -s 10 -rs 115249 -o paired_dat
