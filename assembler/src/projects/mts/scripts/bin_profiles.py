#!/usr/bin/env python
from __future__ import (print_function)

import pandas
from pandas import DataFrame
import sys

profiles_in = pandas.read_table(sys.argv[1], index_col=0)
binning_out = pandas.read_csv(sys.argv[2], index_col=0, names=["bin"], dtype=str)
table = profiles_in.join(binning_out) #, on="contig")
profiles = table.groupby("bin").median()
profiles.to_csv(sys.stdout, sep="\t", header=False)
