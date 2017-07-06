#!/usr/bin/env python
from __future__ import (print_function)

import re
import sys

from common import contig_length
import numpy
import pandas
from pandas import DataFrame

in_fn = sys.argv[1]
d = pandas.read_table(sys.argv[1], names=["name", "bin"], dtype=str)
d["sample"] = d.apply(lambda row: re.findall("\\w+\\d+", row["name"])[0], axis=1)
d["length"] = d.apply(lambda row: contig_length(row["name"]), axis=1)
del d["name"]
info = d.groupby(["bin", "sample"], as_index=False).sum()
info = info.groupby("bin", as_index=False)["length"].max()
info = info[info["length"] > 500000]
info.to_csv(sys.stdout, sep="\t", header=False, index=False)
