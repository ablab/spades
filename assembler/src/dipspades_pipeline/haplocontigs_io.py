#!usr/bin/env python

import os
import sys

def write_haplocontigs_in_file(filename, haplocontigs):
	if os.path.exists(filename):
		os.remove(filename);
	hapfile = open(filename, 'a')
	for hapcontig in haplocontigs:
		hapfile.write(hapcontig + "\n")
