#!/usr/bin/python

import sys
import os
import shutil
import re
import subprocess
logfile_out = open('../../tools/quality/results/all.txt', 'r')
cur_metric_id = 0  
cur_line = 0;      
for line in logfile_out:
	if cur_line == 1:
		if line.split('|')[1] < 1000:
			sys.exit(1);
		if line.split('|')[9] > 0:
			sys.exit(1);			
		print(line.split('|')[1])
		print(line.split('|')[9])
	cur_line += 1;
	
logfile_out.close()
sys.exit(0);

