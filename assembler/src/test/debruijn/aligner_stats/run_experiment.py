import subprocess
import datetime
import os
import sys


if len(sys.argv) < 4:
    print "Usage: run_experiment.py <exe-file> <gap_algo: bf/dijkstra> <datasets list: E.coli_synth,C.elegans,Rumen,Smarr>"
    exit(-1)
run_file = sys.argv[1]
gap_mode = sys.argv[2]
datasets = sys.argv[3].split(",")
feature_name = run_file.split("/")[-3]
exe_name = run_file.split("/")[-1]

resultsdir = "/home/tdvorkina/results/"
datasets_info = "/Sid/tdvorkina/gralign/readme.txt"

run_dir = resultsdir + "/" + feature_name
run_prefix = run_dir + "/" + exe_name


def extract_params(name, fl):
    reads = ""
    saves = ""
    tp = ""
    K = -1
    with open(fl, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith(name + "_lreads"):
                reads = ln.strip().split("=")[-1]
            if ln.startswith(name + "_saves"):
                saves = ln.strip().split("=")[-1]
            if ln.startswith(name + "_K"):
                K = ln.strip().split("=")[-1]
            if ln.startswith(name + "_type"):
                tp = ln.strip().split("=")[-1]
    res = {"reads": reads, "saves": saves, "K": K, "type": tp}
    return res


def run(exe, params, prefix, gap_mode):
    cmd = exe
    args = str(params[d]["K"]) + " " + params[d]["saves"] + " " + params[d]["reads"] + " " + params[d]["type"] + " " + gap_mode + " " + prefix
    print "Run " + cmd + " " + args
    process = subprocess.Popen(cmd + " " + args + " > " + prefix + ".log", shell=True, stderr = subprocess.PIPE)#, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    process.wait()
    output, error = process.communicate()
    if process.returncode:
        print error
        return
    fout = open(prefix + ".log", "a+")
    fout.write(cmd + args)
    fout.close()

params = {}
for d in datasets:
    params[d] = extract_params(d, datasets_info)

if not os.path.exists(resultsdir + "/" + feature_name):
    os.makedirs(resultsdir + "/" + feature_name)
for d in datasets:
    print "Dataset: ", d
    run(run_file, params, run_prefix + "_" + d, gap_mode)

