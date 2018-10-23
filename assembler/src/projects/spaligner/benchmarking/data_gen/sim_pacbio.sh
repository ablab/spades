#!/bin/sh

PBSIMPATH=/home/tdvorkina/soft/pbsim/PBSIM-PacBio-Simulator
REF=ref/MG1655-K12.fasta

$PBSIMPATH/src/pbsim $REF --model_qc $PBSIMPATH/data/model_qc_clr --seed 115249 --prefix sim_pacbio
