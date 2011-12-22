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

files = ["first", "second"]
files += ["single_" + x for x in files]
files += ["jumping_" + x for x in files]
files += ["reference_genome"]

def check(ds):
    m = reduce(lambda x, y: x or y, [miss(ds.get(f)) for f in files])
    print ds["name"], not m

def tar(ds):
    s = [ds.get(f) for f in files]
    s = filter(lambda x: x, s)
    s = reduce(lambda x, y: x + " " + y, s)
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
    sets = sys.argv[3:]
    process(sys.argv[2], tar, lambda ds: ds["name"] in sets);
