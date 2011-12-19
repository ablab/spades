#!/usr/bin/python

import sys
import os
import string

def readline(f):
    while 1:
	line = f.readline()
	if not line:
	    return ""
	line = line.split(';')[0].strip()
	if line: return line

def miss(f):
    if (not f) or (f.lower() == "n/a"):
	return False
    return not os.path.exists(f)

cfg = sys.argv[1]
if not os.path.exists(cfg):
    print "no such file:", cfg
cfg = open(cfg, 'r')
while 1:
    ds = readline(cfg)
    if not ds: break
    s = readline(cfg)
    if s != "{":
	print "'{' expected, but found", s
	exit(1)
    p = {}
    while 1:
	s = readline(cfg)
	if s == "}":
	    break
	s = s.split();
	if (len(s) != 2):
	    print "Invalid property line:", s
	    exit(2)
	p[s[0]] = s[1]
    m = miss(p.get("reference_genome")) or miss(p.get("first")) or miss(p.get("second")) or miss(p.get("single_first")) or miss(p.get("single_second")) or miss(p.get("jumping_first")) or miss(p.get("jumping_second")) or miss(p.get("jumping_single_first")) or miss(p.get("jumping_single_second"))
    print ds, not m
