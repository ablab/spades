############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

./build/abruijn/abruijn \
---input ./data/input/E.Coli.K12.MG1655/EAS20_8/quaked/cropped/s_6.first400000_1.fastq.gz \
--input ./data/input/E.Coli.K12.MG1655/EAS20_8/quaked/s_6_1.cor.fastq.gz \
---input ./data/input/E.Coli.K12.MG1655/emul/MG1655-K12_emul1.fastq.gz \
--output ./data/abruijn/earmark \
--take 1 \
--mode 1 \
---cut 1000 \
