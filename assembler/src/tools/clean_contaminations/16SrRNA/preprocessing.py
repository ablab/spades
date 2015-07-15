#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Makes changes to arb-silva's fasta file:
# * RNA -> DNA (means U -> T)
# * remove spaces in sequences
# * change spaces to underscores in names

# Good for preprocessing before bowtie-build or bwa index.

# Usage example: silva_RNAfasta_to_DNAfasta.py < SSURef_108_tax_silva_trunc.fasta > SSURef_108_tax_silva_trunc_DNA.fasta

import sys

for line in sys.stdin:
        if line[0] != '>':
                line = line.translate(None, ' ')
                line = line.replace('U', 'T')
        else:
                line = line.replace(' ', '_')
        print line,
