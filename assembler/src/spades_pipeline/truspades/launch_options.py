############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

__author__ = 'anton'
import getopt
import os
import sys
import options_storage

class Options:
    def set_default_options(self):
        self.threads = 8
        self.dataset_file = None
        self.input_dirs = None
        self.print_dataset = False
        self.print_commands = False
        self.output_dir = None
        self.command_list = None
        self.spades_options = ""
        self.continue_launch = False
        self.index = ""
        self.reference = ""
        self.mode = "run_truspades"
        self.possible_modes = ["run_truspades", "generate_dataset", "construct_subreferences"]
        self.test = False
        self.clean = False

    def __init__(self, argv, bin, home, version):
        if len(argv) == 1:
            print_usage_and_exit(1, version)
        long_params = "test clean help-hidden construct-dataset reference= reference-index= do= continue " \
                      "threads= help version dataset= input-dir= additional-options=".split(" ")
        short_params = "o:t:hv"
        self.set_default_options()
        self.bin = bin
        self.home = home
        self.version = version
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], short_params, long_params)
            if len(tmp) != 0:
                print_usage_and_exit(1, self.version)
        except getopt.GetoptError:
            _, exc, _ = sys.exc_info()
            sys.stderr.write(str(exc) + "\n")
            print_usage_and_exit(1, self.version)
        for (key, value) in options_list:
            if key == "--version" or key == "-v":
                print_version_and_exit(self.version)
            if key == "--help" or key == "-h":
                print_usage_and_exit(1, self.version)
            elif key == "--test":
                dir = os.path.abspath("spades_test") + "_truspades"
                self.output_dir = dir
                self.input_dirs = [os.path.join(self.home, "test_dataset_truspades")]
                self.test = True
            elif key == "--do":
                self.mode = value
            elif key == "--construct-dataset":
                self.mode = "generate_dataset"
            elif key == "--dataset":
                self.dataset_file = value
            elif key == "--input-dir":
                if self.input_dirs is None:
                    self.input_dirs = []
                self.input_dirs.append(value)
            elif key == "--run-truspades":
                self.mode = "run_truspades"
            elif key == "--reference-index":
                self.index = value
            elif key == "--reference":
                self.reference = value
            elif key == "--continue":
                self.continue_launch = True
            elif key == "--additional-options":
                self.spades_options = value
            elif key == "-o":
                self.output_dir = value
            elif key == "--threads" or key == "-t":
                self.threads = int(value)
            elif key == "--clean":
                self.clean = True
            elif key == "--help-hidden":
                print_usage_and_exit(0, self.version, show_hidden=True)
        if not self.mode in self.possible_modes:
            sys.stderr.write("Error: --do parameter can only have one of the following values: " + ", ".join(self.possible_modes) + "\n")
            print_usage_and_exit(1, self.version)
        if None == self.output_dir or os.path.isfile(self.output_dir):
            sys.stderr.write("Error: Please provide output directory\n")
            print_usage_and_exit(1, self.version)
        if self.continue_launch:
            return
        cnt = len([option for option in [self.dataset_file, self.input_dirs, self.command_list] if option != None])
        if cnt != 1:
            sys.stderr.write("Error: exactly one of dataset-file and input-dir must be specified\n")
            print_usage_and_exit(1, self.version)
        if self.mode == "construct_subreferences":
            if self.index == "":
                sys.stderr.write("Error: Please provide reference index for BWA")
                print_usage_and_exit(1, self.version)
            if self.reference == "":
                sys.stderr.write("Error: Please provide reference for subreference construction")
                print_usage_and_exit(1, self.version)


def print_usage_and_exit(code, version, show_hidden=False):
    sys.stderr.write("SPAdes genome assembler v" + str(version) + " [truSPAdes mode]\n\n")
    sys.stderr.write("Usage: " + str(sys.argv[0]) + " [options] -o <output_dir>" + "\n")
    sys.stderr.write("" + "\n")
    sys.stderr.write("Basic options:" + "\n")
    sys.stderr.write("-h/--help\t\t\tprints this usage message" + "\n")
    sys.stderr.write("-v/--version\t\t\tprints version" + "\n")
    sys.stderr.write("--test\t\t\t\trun truSPAdes on toy dataset" + "\n")
    sys.stderr.write("-o\t\t<output_dir>\tdirectory to store all the resulting files (required)" + "\n")
    sys.stderr.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")
    sys.stderr.write("--continue\t\t\tcontinue interrupted launch" + "\n")
    sys.stderr.write("--construct-dataset\t\tparse dataset from input folder" + "\n")
    sys.stderr.write("" + "\n")
    sys.stderr.write("Input options:" + "\n")
    sys.stderr.write("--input-dir\t<directory>\tdirectory with input data. Note that the directory should contain only files with reads. This option can be used several times to provide several input directories." + "\n")
    sys.stderr.write("--dataset\t<file>\t\tfile with dataset description" + "\n")
    if show_hidden:
        pass
        #ToDo
    # sys.stderr.write("" + "\n")
    # sys.stderr.write("Output options:" + "\n")
    # sys.stderr.write("--print-dataset\tprints file with dataset generated after analysis of input directory contents" + "\n")
    # sys.stderr.write("--print-commands\tprints file with truspades commands that would assemble barcodes from dataset" + "\n")
    # sys.stderr.write("--run-truspades\truns truSPAdes on all barcodes" + "\n")
    sys.stderr.flush()
    sys.exit(code)


def print_version_and_exit(version):
    options_storage.version(version, mode="tru")
    sys.exit(0)
