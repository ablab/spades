#!/usr/bin/python

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2021-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing SPAdes
# provide a path to .yaml file with test description

import sys
import workflow_base
import teamcity_workflow


def cmp_phases_output(dataset_info, test, output_dir, log):
    print("Cmp output")
    sys.stdout.flush()
    log.log("Start compare outputs")
    outputs = workflow_base.get_outdirs(dataset_info, test, output_dir)
    log.log("Outputs to compare: " + str(outputs))
    for i in range(1, len(outputs)):
        ecode = teamcity_workflow.cmp_with_etalon(outputs[0], outputs[i],
                                                  allowed_substr=[".fasta", ".paths", ".gfa", ".fastg"], print_info=True)
        if ecode != 0:
            log.err("Comparing outputs did not pass, exit code " + str(ecode))
            return 12
    return 0


print("Stable test")
sys.stdout.flush()
if __name__ == "__main__":
    workflow_base.main(check_test=cmp_phases_output)