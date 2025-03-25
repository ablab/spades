#!/usr/bin/python

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing SPAdes
# provide a path to .yaml file with test description

import workflow_base
from workflow_base import log

import sys
import os
import filecmp

def is_contain_one_of_substring(subs, s):
    for sub in subs:
        if (sub in s):
            return True
    return False


def create_ignore_list_of_files_in_dir(dir, allow_sub):
    ignore_list = []
    for file in os.listdir(dir):
        if os.path.isfile(os.path.join(dir, file)):
            if (not is_contain_one_of_substring(allow_sub, file)):
                ignore_list.append(file)

    return ignore_list


def cmp_folder(output_dir, etalon_dir, ignore, allowed_substring):
    log.log("cmp folder " + output_dir + " with etalon")

    #delete lines with tmp folders from logs (tmp folders have a different name each time)
    os.system("find " + output_dir + " -type f -exec sed -i '/_dir/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/_dir/d' {} \\;")

    os.system("find " + output_dir + " -type f -exec sed -i '/corrector_/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/corrector_/d' {} \\;")

    os.system("find " + output_dir + " -type f -exec sed -i '/tmp/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/tmp/d' {} \\;")

    os.system("find " + output_dir + " -type f -exec sed -i '/agent/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/agent/d' {} \\;")

    os.system("find " + output_dir + " -type f -exec sed -i '/hammer_/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/hammer_/d' {} \\;")

    os.system("find " + output_dir + " -type f -exec sed -i '/spades_/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/spades_/d' {} \\;")

    os.system("find " + output_dir + " -type f -exec sed -i '/version/d' {} \\;")
    os.system("find " + etalon_dir + " -type f -exec sed -i '/version/d' {} \\;")

    ignore_list = ignore
    ignore_list += create_ignore_list_of_files_in_dir(output_dir, ignore + allowed_substring)
    ignore_list += create_ignore_list_of_files_in_dir(etalon_dir, ignore + allowed_substring)

    dircmp = filecmp.dircmp(etalon_dir, output_dir, ignore=ignore_list)
    dircmp.diff_files = [x for x in dircmp.diff_files if is_contain_one_of_substring(allowed_substring, x)]
    dircmp.right_only = [x for x in dircmp.right_only if is_contain_one_of_substring(allowed_substring, x)]
    dircmp.left_only = [x for x in dircmp.left_only if is_contain_one_of_substring(allowed_substring, x)]

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


def cmp_folder_rec(output_dir, etalon_dir, ignore, allowed_substr, print_info=False):
    folder_ignore = ['tmp', 'saves', '.bin_reads', 'path_extend']
    subdirs_etalon = [x for x in os.listdir(etalon_dir)
                      if os.path.isdir(os.path.join(etalon_dir, x)) and x not in folder_ignore]

    subdirs_output = [x for x in os.listdir(output_dir)
                        if os.path.isdir(os.path.join(output_dir, x)) and x not in folder_ignore]

    errcode = cmp_folder(output_dir, etalon_dir, ["tmp", "saves", ".bin_reads", 'test_run.info'] + ignore,
                         allowed_substr)

    if (errcode != 0):
        return errcode

    subdirs_etalon.sort()
    subdirs_output.sort()
    if (subdirs_etalon != subdirs_output):
        log.err("Etalon dirs(" + str(subdirs_etalon) +") != output dirs (" + str(subdirs_output) + ")")
        return 12
    elif print_info:
        log.log("%s == %s" % (str(subdirs_etalon), str(subdirs_output)))

    for subdir in subdirs_etalon:
        errcode = cmp_folder_rec(output_dir + "/" + subdir, etalon_dir + "/" + subdir, ignore, allowed_substr, print_info)
        if (errcode != 0):
            return errcode
    return 0


def cmp_with_etalon(output_dir, etalon_dir, ignore=None, allowed_substr=None, print_info=False):
    if print_info:
        log.log("Compare %s and %s" % (output_dir, etalon_dir))
    if (ignore == None):
        ignore = []
    if (allowed_substr == None):
        allowed_substr = []

    return cmp_folder_rec(output_dir, etalon_dir, ignore, allowed_substr, print_info)


def etalon_saves(dataset_info, test, output_dir, log):
    if 'etalon_saves' in dataset_info:
        log.log("Comparing etalon saves now")
        etalon_folder = dataset_info["etalon_saves"]
        if ("name" in test):
            etalon_folder += test["name"]

        ecode = cmp_with_etalon(output_dir, os.path.join(etalon_folder),
                                allowed_substr=[".yaml", ".sh", "params.txt"])

        if ecode != 0:
            log.err("Comparing etalon saves did not pass, exit code " + str(ecode))
            return 12
    return 0

if __name__ == "__main__":
    print("Workflow main")
    sys.stdout.flush()
    workflow_base.main(etalon_saves)
