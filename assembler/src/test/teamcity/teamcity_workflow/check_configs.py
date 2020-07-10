#!/usr/bin/python

############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing the correct choosing of kmers in SPAdes
# provide a path to .yaml file with test description

import workflow_base
from workflow_base import log
import os


def check_one_out_folder(dataset_info, test, line_in_config, output_dir):
    for line_to_check in (test["line_in_config"] + dataset_info["line_in_config"] + line_in_config):
        if "not present" in line_to_check:
            log.log("Checking \"" + line_to_check["line"] + "\" not present in " + line_to_check["config"])
        else:
            log.log("Checking \"" + line_to_check["line"] + "\" present in " + line_to_check["config"])

        line_present = False
        line_with = ""
        with open(os.path.join(output_dir, line_to_check["config"])) as f:
            for line in f:
                if line_to_check["line"] in line:
                    line_present = True
                    line_with = line

        if line_present == ("not present" in line_to_check):
            log.err("Checking fail!")
            if line_present:
                log.err("Line with pattern: " + line_with)
            return 12
    return 0


def check_lines_in_cfg(dataset_info, test, output_dir):
    if "line_in_config" not in test:
        test["line_in_config"] = []
    if "line_in_config" not in dataset_info:
        dataset_info["line_in_config"] = []

    if "phases" in test:
        for i in range(len(test["phases"])):
            if "name" in test["phases"][i]:
                phase_name = test["phases"][i]["name"]
            else:
                phase_name = dataset_info["phases"][i]["name"]

            phase_outputdir = os.path.join(output_dir, phase_name)
            phases_line = []
            if "line_in_config" in test["phases"][i]:
                phases_line = test["phases"][i]["line_in_config"]

            err_code = check_one_out_folder(dataset_info, test, phases_line, phase_outputdir)
            if err_code != 0:
                return err_code
        return 0
    else:
        output_dir = os.path.join(output_dir, "out")
        return check_one_out_folder(dataset_info, test, [], output_dir)


workflow_base.main(check_test=check_lines_in_cfg)