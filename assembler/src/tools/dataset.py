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
    return not os.path.isfile(f)

def presentFile(f):
    if (not f) or (f.lower() == "n/a"):
	return False
    return os.path.isfile(f)

files = ["first", "second"]
files += ["single_" + x for x in files]
files += ["jumping_" + x for x in files]
files += ["reference_genome"]

def check(ds):
    print ds["name"], "is present"

def tar(ds):
    s = filter(presentFile, [ds.get(f) for f in files])
    s = reduce(lambda x, y: x + " " + y, s)
    a = ds["name"] + ".tar"
    if os.path.exists(a):
	print "ERROR:", a, "already exists"
	exit(2)
    print "tarring", ds["name"], "..."
    os.system("tar -cf " + a + " " + s)

def md5(ds):
    s = filter(presentFile, [ds.get(f) for f in files])
    print ds["name"]
    for f in s:
	os.system("md5sum " + f)

def process(cfg, func, filt):
    if not os.path.isfile(cfg):
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
	if filt(ds):
	    if reduce(lambda x, y: x or y, [missFile(ds.get(f)) for f in files]):
		print ds["name"], "is missing!!!!!!!!!!!!!!!!!!!!"
	    else:
		func(ds)

def printDS(p):
    print "{"
    for (a, b) in p:
	print "\t" + a + "\t" + b
    print "}"

def hammer(prefix):
    os.system('((ls -1 ' + prefix + '* 2> /dev/null | grep left.cor | grep -v single) || echo "") | head -n 1')

if sys.argv[1] == "check":
    process(sys.argv[2], check, lambda ds: True);
if sys.argv[1] == "tar":
    regexp = ""
    if 3 < len(sys.argv): regexp = "^.*" + sys.argv[3] + ".*$"
    filt = lambda ds: re.match(regexp, ds["name"])
    process(sys.argv[2], tar, filt);
if sys.argv[1] == "md5":
    regexp = ""
    if 3 < len(sys.argv): regexp = "^.*" + sys.argv[3] + ".*$"
    filt = lambda ds: re.match(regexp, ds["name"])
    process(sys.argv[2], md5, filt);
if sys.argv[1] == "hammer":
    hammer(sys.argv[2])
