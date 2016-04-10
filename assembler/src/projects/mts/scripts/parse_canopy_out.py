#!/usr/bin/env python

import os
import sys
from itertools import izip

print " ".join(sys.argv)
if len(sys.argv) < 3: 
    print "Usage: %s <canopy output file> <annotation output root>" % sys.argv[0]
    sys.exit(1) 

out_root=sys.argv[2]

samples_annotation=dict()

with open(sys.argv[1], "r") as input_file:
    for l in input_file:
        annotation_str = l.split()
        bin_id = annotation_str[0]
        sample_contig = annotation_str[1].split('-')
        sample = sample_contig[0]
        contig = sample_contig[1]

        if sample not in samples_annotation:
            samples_annotation[sample] = dict()

        annotation = samples_annotation[sample]

        if contig not in annotation:
            annotation[contig] = list()

        annotation[contig].append(bin_id)

    
for sample in samples_annotation:
    with open(out_root + '/' + sample + '.ann', 'w') as sample_out:
        annotation = samples_annotation[sample]
        
        for contig in annotation:
            print >>sample_out, contig + ' : ' + ' '.join(annotation[contig])
