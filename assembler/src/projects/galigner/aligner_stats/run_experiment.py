import subprocess
import datetime
import os
import sys


if len(sys.argv) < 4:
    print "Usage: run_experiment.py <exe-file> <gap_algo: bf/dijkstra> <datasets list: E.coli_synth,C.elegans,Rumen,Smarr>"
    exit(-1)
run_file = sys.argv[1]
yaml_template = sys.argv[2]
datasets = sys.argv[3].split(",")
feature_name = run_file.split("/")[-3]
exe_name = run_file.split("/")[-1]

resultsdir = "/home/tdvorkina/results/"
datasets_info = "/Sid/tdvorkina/gralign/readme.txt"
#yaml_template = "/home/tdvorkina/tmp/algorithmic-biology/assembler/src/projects/galigner/galigner.yaml"

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

def prepare_cfg(yaml_template, params, prefix, gap_mode):
    tmp_cfg = prefix + "_" + gap_mode + "_cfg.yaml"
    with open(tmp_cfg, "w") as fout:
        with open(yaml_template, "r") as fin:
            for ln in fin.readlines():
                print_ln = ln
                for p in params.keys():
                    if ln.startswith(p):
                        print_ln = p + ": " + params[p] + "\n"
                        break
                if ln.startswith("run_dijkstra"):
                    if gap_mode == "bf":
                        print_ln = "run_dijkstra: false\n"
                    else:
                        print_ln = "run_dijkstra: true\n"
                fout.write(print_ln)
    return tmp_cfg


def run(exe, params, prefix, gap_mode):
    cfg = prepare_cfg(yaml_template, params, prefix, gap_mode)
    cmd = exe
    print "Run " + cmd + " " + cfg  + " -o " + prefix + "_" + gap_mode + " > " + prefix + "_" + gap_mode + ".log"
    process = subprocess.Popen(cmd + " " + cfg + " -o " + prefix  + "_" + gap_mode + " > " + prefix + "_" + gap_mode + ".log", shell=True, stderr = subprocess.PIPE)#, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    process.wait()
    output, error = process.communicate()
    if process.returncode:
        print error
        return
    fout = open(prefix + "_" + gap_mode + ".log", "a+")
    fout.write(cmd + " " + cfg + " -o " + prefix  + "_" + gap_mode)
    fout.close()

params = {}
for d in datasets:
    params[d] = extract_params(d, datasets_info)

if not os.path.exists(resultsdir + "/" + feature_name):
    os.makedirs(resultsdir + "/" + feature_name)
for d in datasets:
    print "Dataset: ", d
    run(run_file, params[d], run_prefix + "_" + d, gap_mode)

