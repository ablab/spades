#!/usr/bin/python

import sys
import os
import string
import re

def readline(f):
    while 1:
	line = f.readline()
	if not line:
	    return ""
	line = line.split(';')[0].strip()
	if line: return line

def missFile(f):
    if (not f) or (f.lower() == "n/a"):
	return False
    return not os.path.exists(f)

files = ["first", "second"]
files += ["single_" + x for x in files]
files += ["jumping_" + x for x in files]
files += ["reference_genome"]

def miss(ds):
    return reduce(lambda x, y: x or y, [missFile(ds.get(f)) for f in files])

def check(ds):
    print ds["name"], not miss(ds)

def tar(ds):
    if miss(ds):
	print "#####", ds["name"], "is missing!", "#####"
	return
    s = [ds.get(f) for f in files]
    s = filter(lambda x: x, s)
    s = reduce(lambda x, y: x + " " + y, s)
    print "tarring", ds["name"], "..."
    os.system("tar -cf " + ds["name"] + ".tar " + s)

def process(cfg, func, filt):
    if not os.path.exists(cfg):
        print "no such file:", cfg
    cfg = open(cfg, 'r')
    while 1:
        ds = readline(cfg)
        if not ds: break
        if ds.startswith("#include"): continue
        s = readline(cfg)
        if s != "{":
	    print "'{' expected, but found", s
	    exit(1)
        ds = {"name": ds}
        while 1:
	    s = readline(cfg)
	    if s == "}":
	        break
	    s = s.split();
	    if (len(s) != 2):
	        print "Invalid property line:", s
	        exit(2)
	    ds[s[0]] = s[1]
	if filt(ds): func(ds)

if sys.argv[1] == "check":
    process(sys.argv[2], check, lambda ds: True);
if sys.argv[1] == "tar":
    regexp = ("^.*" + sys.argv[3] + ".*$") if 3 < len(sys.argv) else ""
    filt = lambda ds: re.match(regexp, ds["name"])
    process(sys.argv[2], tar, filt);
