#!/bin/sh
mkdir -p data/output/hists
java -cp build/ GetErrors -f distance_estimation
