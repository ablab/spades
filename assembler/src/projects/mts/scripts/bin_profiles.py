#!/usr/bin/env python
from __future__ import (print_function)

import pandas
from pandas import DataFrame
import sys

jgi_mode = False
has_var = False
for arg in sys.argv[3:]:
    if arg == "--jgi":
        jgi_mode = True
    elif arg == "-v":
        has_var = True
profiles_in = pandas.read_table(sys.argv[1], index_col=0, header=0 if jgi_mode else None)
profiles_in = profiles_in[profiles_in.columns[1 if jgi_mode else 0::2 if has_var else 1]]
binning_out = pandas.read_table(sys.argv[2], index_col=0, names=["bin"], dtype=str)
table = profiles_in.join(binning_out)
profiles = table.groupby("bin").median()
profiles.to_csv(sys.stdout, sep="\t", header=False)
