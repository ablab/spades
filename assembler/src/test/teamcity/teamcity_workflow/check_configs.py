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


def check_lines_in_cfg(dataset_info, test, output_dir):
    if "line_in_config" not in test:
        test["line_in_config"] = []
    if "line_in_config" not in dataset_info:
        dataset_info["line_in_config"] = []

    output_dir = os.path.join(output_dir, "out")
    for line_to_check in (test["line_in_config"] + dataset_info["line_in_config"]):
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


workflow_base.main(check_test=check_lines_in_cfg)