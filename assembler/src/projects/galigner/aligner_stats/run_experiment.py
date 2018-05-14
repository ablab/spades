import subprocess
import datetime
import os
import sys


if len(sys.argv) < 4:
    print "Usage: run_experiment.py <exe-file> <yaml-template> <datasets list: E.coli_synth,C.elegans,Rumen,Smarr>"
    exit(-1)
run_file = sys.argv[1]
yaml_template = sys.argv[2]
datasets = sys.argv[3].split(",")
feature_name = run_file.split("/")[-3]
exe_name = run_file.split("/")[-1]

resultsdir = "/home/tdvorkina/results/"
datasets_info = "/Sid/tdvorkina/gralign/readme.txt"
yaml_name = yaml_template.split("/")[-1][:-5]

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
    res = {"path_to_sequences": reads, "path_to_graphfile": saves, "k": K, "data_type": tp}
    return res

def prepare_cfg(yaml_template, params, prefix, yaml_name):
    tmp_cfg = prefix + "_" + yaml_name + "_cfg.yaml"
    with open(tmp_cfg, "w") as fout:
        with open(yaml_template, "r") as fin:
            for ln in fin.readlines():
                print_ln = ln
                for p in params.keys():
                    if ln.startswith(p):
                        print_ln = p + ": " + params[p] + "\n"
                        break
                fout.write(print_ln)
    return tmp_cfg


def run(exe, params, prefix, yaml_name):
    cfg = prepare_cfg(yaml_template, params, prefix, yaml_name)
    cmd = exe
    print "Run " + cmd + " " + cfg  + " -o " + prefix + "_" + yaml_name + " > " + prefix + "_" + yaml_name + ".log"
    process = subprocess.Popen(cmd + " " + cfg + " -o " + prefix  + "_" + yaml_name + " > " + prefix + "_" + yaml_name + ".log", shell=True, stderr = subprocess.PIPE)#, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    process.wait()
    output, error = process.communicate()
    if process.returncode:
        print error
        return
    fout = open(prefix + "_" + yaml_name + ".log", "a+")
    fout.write(cmd + " " + cfg + " -o " + prefix  + "_" + yaml_name)
    fout.close()

params = {}
for d in datasets:
    params[d] = extract_params(d, datasets_info)

if not os.path.exists(resultsdir + "/" + feature_name):
    os.makedirs(resultsdir + "/" + feature_name)
for d in datasets:
    print "Dataset: ", d
    run(run_file, params[d], run_prefix + "_" + d, yaml_name)

