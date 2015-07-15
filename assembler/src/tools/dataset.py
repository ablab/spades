#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import string
import re
import subprocess
import datetime

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

def read_files():
    read_files = ["first", "second"]
    read_files += ["single_" + x for x in read_files]
    read_files += ["jumping_" + x for x in read_files]
    return read_files

def ref_files():
    return ["reference_genome"]

def misc_props():
    return ["RL", "IS", "jump_is", "single_cell"]

def files():
    return read_files() + ref_files()

def props():
    return read_files() + misc_props() + ref_files()

def check(ds):
    print ds["name"], "is present"

def tar(ds):
    s = filter(presentFile, [ds.get(f) for f in files()])
    s = reduce(lambda x, y: x + " " + y, s)
    a = ds["name"] + ".tar"
    if os.path.exists(a):
	print "ATTENTION:", a, "already exists"
	return
    print "tarring", ds["name"], "..."
    os.system("tar -cf " + a + " " + s)

def md5(ds):
    s = filter(presentFile, [ds.get(f) for f in files()])
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
	    if reduce(lambda x, y: x or y, [missFile(ds.get(f)) for f in files()]):
		print ds["name"], "is missing!!!!!!!!!!!!!!!!!!!!"
	    else:
		func(ds)

def hammer_interactive(prefix):
    output_dir = "bh" + datetime.date.today().strftime('%Y%m%d')
    if os.path.exists(output_dir):
	print output_dir, "already exists!", "Please enter another directory name:"
	output_dir = raw_input().strip()
    output_dir = os.path.abspath(output_dir)
    #left_cor = subprocess.check_output('((ls -1 ' + prefix + '* 2> /dev/null | grep left.cor | grep -v single) || echo "")', shell=True).strip()
    try:
	ls = subprocess.check_output('ls -1 ' + prefix + '*', shell=True)
    except:
	#print "Not found:", prefix + '*'
	return
    ls = filter(os.path.isfile, ls.split('\n'))
    for i in range(len(ls)):
	f = ls[i]
	print i, ":", os.path.basename(f)
    def askFile(prop):
	print 'Which file is "' + prop + '"?', "(enter number from 0 to", (len(ls) - 1), "or press Enter if none)"
	a = raw_input().strip()
	if not a:
	    return ""
	return ls[int(a)]
    def askProp(prop):
	print "Enter", prop, "(or press Enter if you don't know yet)"
	a = raw_input().strip()
	if not a:
	    return "TODO"
	return a
    asked_props = {}
    for prop in props():
	if prop in read_files():
	    a = askFile(prop)
	else:
	    a = askProp(prop)
	if a:
	    asked_props[prop] = a
    hammer(asked_props, output_dir, True)

def hammer(given_props, output_dir, verbose):
    dataset_entry = []
    cmd = ["mkdir " + output_dir]
    for prop in props():
	if not (prop in given_props):
	    continue
	val = given_props[prop]
	if prop in read_files():
	    newfile = output_dir + "/" + os.path.basename(val)
	    cmd += ["cp " + val + " " + output_dir + "/"]
	    cmd += ["gzip -9 " + newfile]
	    newfile += ".gz"
	    val = re.sub(".*input/", "", newfile)
	dataset_entry += [(prop, val)]
    dataset_entry = ["{"] + map(lambda (a, b): "\t" + a + "\t" + b, dataset_entry) + ["}"]
    dataset_entry = reduce(lambda x, y: x + "\n" + y, dataset_entry)
    if verbose:
        print dataset_entry
    for c in cmd:
	if verbose:
	    print c
	os.system(c)
    return dataset_entry

if __name__ == "__main__":
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
        hammer_interactive(sys.argv[2])
