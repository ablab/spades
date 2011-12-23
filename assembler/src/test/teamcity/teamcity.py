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
		n50 = int(line.split('|')[1]);
		mis = int(line.split('|')[11]);
		print('n50 = ' + str(n50) + ' missasembled contigs ' + str(mis));
		if n50 < 75000:
			print('n50 too small')
			sys.exit(1);
		if (mis > 1):
			print('too much miss')
			sys.exit(1);			
		
	cur_line += 1;
	
logfile_out.close()
sys.exit(0);

