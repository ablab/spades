#!/bin/bash

NANOSIMPATH=
REF=
REALREADS=

# Profiling stage, make sure to set the mode of read_analysis.py to -r-x or above
python2 $NANOSIMPATH/read_analysis.py -i $REALREADS -r $REF -o sim_nanopore

# Simulation stage, suppose the genome to be simulated is called test.fasta and make sure to provide the correct path to it
python2 $NANOSIMPATH/simulator.py circular --seed 115249 -n 40000 -r $REF  -c sim_nanopore
