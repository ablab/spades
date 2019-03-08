#!/usr/bin/python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing SPAdes
# provide a path to .yaml file with test description

import sys
import os
import shutil
import argparse
import subprocess
from traceback import print_exc
import filecmp

#Log class, use it, not print
class Log:

    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)
        sys.stdout.flush()

    def warn(self, s):
        msg = "WARNING: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('info', metavar='CONFIG_FILE', type=str,  help='a path to .yaml file with test description')
    parser.add_argument("--run_name", "-n", help="output dir custom suffix", type=str)
    parser.add_argument("--spades_path", "-p", help="custom directory to spades.py", type=str)
    parser.add_argument("--no_cfg_and_compilation", dest="cfg_compilation",
                        default=True,
                        help="don't copy configs or compile SPAdes even if .info file states otherwise",
                        action='store_false')
    parser.add_argument("--spades_cfg_dir", "-c", help="SPAdes config directory", type=str)
    args = parser.parse_args()
    return args


def load_info(info_filename):
    sys.path.append("./ext/src/python_libs/")
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

    info = pyyaml.load(open(info_filename, 'r'))
    return info


# Create output folder
def create_output_dir(args, dataset_info, test):
    #make dirs and remembering history
    output_dir = os.path.join(dataset_info["output_dir"], dataset_info["name"])
    if ("name" in test):
        output_dir = os.path.join(output_dir, test["name"])

    if "build_agent" in dataset_info:
        output_dir += "_" + dataset_info["build_agent"]

    # add custom suffix if present
    if args.run_name:
        output_dir += "_" + args.run_name

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    return output_dir


# Save meta information about this teamcity.py run
def save_run_info(args, output_dir):
    run_info = open(os.path.join(output_dir, "test_run.info"), "w")
    run_info.write(".info file: " + args.info + "\n")
    if args.run_name:
        run_info.write("run name: " + args.run_name + "\n")
    if args.spades_path:
        run_info.write("path to spades.py: " + str(args.spades_path) + "\n")
    if args.spades_cfg_dir:
        run_info.write("spades config direrctory: " + args.spades_cfg_dir + "\n")
    run_info.close()


# Compile SPAdes
def compile_spades(args, dataset_info, working_dir):
    if not args.cfg_compilation:
        log.log("Forced to use current SPAdes build, will not compile SPAdes")
    elif 'spades_compile' not in dataset_info or dataset_info["spades_compile"]:
        comp_params = ' '
        if 'compilation_params' in dataset_info:
            comp_params = " ".join(dataset_info["compilation_params"])

        bin_dir = 'build_spades'
        if not os.path.exists(bin_dir):
            os.makedirs(bin_dir)
        os.chdir(bin_dir)

        #Compilation
        err_code = os.system('cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=' + working_dir + ' ' + os.path.join(working_dir, 'src') + comp_params)
        err_code = err_code | os.system('make -j 16')
        err_code = err_code | os.system('make install')

        os.chdir(working_dir)

        if err_code != 0:
            # Compile from the beginning if failed
            shutil.rmtree('bin', True)
            shutil.rmtree('build_spades', True)
            return os.system('./spades_compile.sh ' + comp_params)
    return 0


#Create SPAdes command line
def make_spades_cmd(args, dataset_info, test, spades_params_list, spades_dir, output_dir):
    #Correct paths in params
    input_file_option_list = ['-1', '-2', "--12", "-s", "--merged",
                              "--pe", "--s", "--mp",
                              "--hqmp", "--nxmate",
                              "--sanger", "--pacbio", "--nanopore",
                              "--tslr", "--trusted-contigs", "--untrusted-contigs"]

    def contain_substr(str, substr_list):
        for substr in substr_list:
            if substr in str:
                if (substr == str) or not str[len(substr)].isalpha():
                    return True
        return False

    spades_params = []
    for option in (spades_params_list):
        if not isinstance(option, list):
            spades_params.append(option)
        else:
            if contain_substr(option[0], input_file_option_list):
                for i in range(0, len(option) - 1):
                    spades_params.append(option[i])

                parts = option[-1].split(':')
                if option[-1][0] != '/':
                    spades_params.append(os.path.join(dataset_info["dataset_path"], parts[-1]))
                else:
                    spades_params.append(parts[-1])

                if (len(parts) > 1):
                    spades_params[-1] = parts[0] + ":" + spades_params[-1]
            else:
                for elem in option:
                    spades_params.append(elem)

    if args.spades_cfg_dir:
        spades_params.append("--configs-dir")
        spades_params.append(args.spades_cfg_dir)

    spades_exec = "spades.py"

    cmd = os.path.join(spades_dir, spades_exec) + " " + " ".join(spades_params)
    if not("supress_output_dir" in test and test["supress_output_dir"]) and \
            not ("supress_output_dir" in dataset_info and dataset_info["supress_output_dir"]):
        cmd += " -o " + output_dir
    cmd += " > /dev/null"
    return cmd


def contain_prohib(prohib_subs, s):
    for subs in prohib_subs:
        if (subs in s):
            return True
    return False


def cmp_folder(output_dir, etalon_dir, ignore, ignore_substr):
    log.log("cmp folder " + output_dir + " with etalon")
    os.system("find " + output_dir + " -type f -exec sed -i '/_dir/d' {} \;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/_dir/d' {} \;")

    os.system("find " + output_dir + " -type f -exec sed -i '/agent/d' {} \;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/agent/d' {} \;")

    dircmp = filecmp.dircmp(etalon_dir, output_dir, ignore=ignore)
    dircmp.diff_files = [x for x in dircmp.diff_files if
                         not contain_prohib(ignore_substr + ["fasta", "fastq", "fastg", "gfa", "paths"], x)]
    dircmp.right_only = [x for x in dircmp.right_only if not contain_prohib(ignore_substr, x)]
    dircmp.left_only = [x for x in dircmp.left_only if not contain_prohib(ignore_substr, x)]

    if (dircmp.diff_files != []):
        log.err(str(dircmp.diff_files) + " differ from etalon")
        return 12
    elif (dircmp.right_only != []):
        log.err(str(dircmp.right_only) + " present in output but don't present in etalon")
        return 12
    elif (dircmp.left_only != []):
        log.err(str(dircmp.left_only) + " present in etalon but don't present in output")
        return 12
    else:
        return 0


def cmp_folder_rec(output_dir, etalon_dir, ignore, ignore_substr):
    folder_ignore = ['tmp', 'saves', '.bin_reads', 'path_extend']
    subdirs_etalon = [x for x in os.listdir(etalon_dir)
                      if os.path.isdir(os.path.join(etalon_dir, x)) and x not in folder_ignore]

    subdirs_output = [x for x in os.listdir(output_dir)
                        if os.path.isdir(os.path.join(output_dir, x)) and x not in folder_ignore]

    errcode = cmp_folder(output_dir, etalon_dir, [
                             "params.txt", "spades.log",
                             "tmp", "saves", ".bin_reads", "warnings.log",
                             "input_dataset.yaml", "dont_correct.yaml",
                             "final.lib_data", "corrected.yaml", "run_spades.sh", "run_spades.yaml"] + ignore,
               ignore_substr)

    if (errcode != 0):
        return errcode

    subdirs_etalon.sort()
    subdirs_output.sort()
    if (subdirs_etalon != subdirs_output):
        log.err("Etalon dirs(" + str(subdirs_etalon) +") != output dirs (" + str(subdirs_output) + ")")
        return 12

    for subdir in subdirs_etalon:
        errcode = cmp_folder_rec(output_dir + "/" + subdir, etalon_dir + "/" + subdir, ignore, ignore_substr)
        if (errcode != 0):
            return errcode
    return 0


def cmp_with_etalon(output_dir, etalon_dir, ignore=None, ignore_substr=None):
    if (ignore == None):
        ignore = []
    if (ignore_substr == None):
        ignore_substr = []

    return cmp_folder_rec(output_dir, etalon_dir, ignore, ignore_substr)


def run_spades(spades_dir, args, dataset_info, test, spades_params, output_dir, finish_with_error):
    if args.spades_path:
        spades_dir = args.spades_path
        log.log("Different spades.py path specified: " + spades_dir)

    spades_cmd = make_spades_cmd(args, dataset_info, test, spades_params, spades_dir, output_dir)
    log.log("Launching: " + spades_cmd)

    ecode = os.system(spades_cmd)
    if finish_with_error:
        if ecode == 0:
            log.err("SPAdes finished without error")
            return 4
    else:
        if ecode != 0:
            log.err("SPAdes finished abnormally with exit code " + str(ecode))
            return 4
    return 0


def etalon_saves(dataset_info, test, output_dir):
    if 'etalon_saves' in dataset_info:
        log.log("Comparing etalon saves now")
        etalon_folder = dataset_info["etalon_saves"]
        if ("name" in test):
            etalon_folder += test["name"]

        ecode = cmp_with_etalon(output_dir, os.path.join(etalon_folder))

        if ecode != 0:
            log.err("Comparing etalon saves did not pass, exit code " + str(ecode))
            return 12
    return 0


def handle_one_test(test, args, dataset_info, working_dir):
    if "name" in test:
        log.log("Start test: " + test["name"])

    output_dir = create_output_dir(args, dataset_info, test)
    save_run_info(args, output_dir)

    finish_with_error = (("finish_with_error" in dataset_info) and dataset_info["finish_with_error"])
    spades_params = []

    ecode = 0
    if "phases" in test:
        if "phases" not in dataset_info:
            dataset_info["phases"] = [] * len(test["phases"])
        for i in range(len(test["phases"])):
            spades_params = []
            if "spades_params" in test["phases"][i]:
                spades_params += test["phases"][i]["spades_params"]
            if "spades_params" in dataset_info["phases"][i]:
                spades_params += dataset_info["phases"][i]["spades_params"]
            if "spades_params" in dataset_info:
                spades_params += dataset_info["spades_params"]
            if "spades_params" in test:
                spades_params += test["spades_params"]

            if "name" in test["phases"][i]:
                phase_name = test["phases"][i]["name"]
            else:
                phase_name = dataset_info["phases"][i]["name"]
                test["phases"][i]["name"] = dataset_info["phases"][i]["name"]

            phase_outputdir = os.path.join(output_dir, phase_name)
            if i == 0:
                os.makedirs(phase_outputdir)
            else:
                prev_phase_name = test["phases"][i - 1]["name"]
                shutil.copytree(os.path.join(output_dir, prev_phase_name), phase_outputdir)

            cur_finish_with_error = finish_with_error
            if i < len(test["phases"]) - 1:
                cur_finish_with_error = False

            local_ecode = run_spades(working_dir, args, dataset_info, test["phases"][i], spades_params,
                       phase_outputdir, cur_finish_with_error)
            if local_ecode != 0:
                ecode = 4
                break
    else:
        if "spades_params" in dataset_info:
            spades_params += dataset_info["spades_params"]
        if "spades_params" in test:
            spades_params += test["spades_params"]

        phase_outputdir = os.path.join(output_dir, "out")
        os.makedirs(phase_outputdir)

        local_ecode = run_spades(working_dir, args, dataset_info, test, spades_params,
                                 phase_outputdir, finish_with_error)
        if local_ecode != 0:
            ecode = local_ecode

    if ecode == 0:
        ecode = etalon_saves(dataset_info, test, output_dir)

    if ecode != 0:
        log.log("TEST NOT PASS")
    else:
        log.log("TEST PASS")

    return ecode

### main ###
try:
    if len(sys.argv) == 1:
        command = 'python {} -h'.format(sys.argv[0])
        subprocess.call(command, shell=True)
        sys.exit(1)

    sys.stderr = sys.stdout
    exit_code = 0
    args = parse_args()
    dataset_info = load_info(args.info)
    working_dir = os.getcwd()

    # compile
    ecode = compile_spades(args, dataset_info, working_dir)
    if ecode != 0:
        log.err("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)

    cnt_pass = 0
    for test in dataset_info["tests"]:
        local_ecode = handle_one_test(test, args, dataset_info, working_dir)
        if local_ecode == 0:
            cnt_pass += 1
        else:
            exit_code = local_ecode
    log.log(str(cnt_pass) + "/" + str(len(dataset_info["tests"])) + " TESTS PASSED")
    sys.exit(exit_code)

except SystemExit:
    raise

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(239)
