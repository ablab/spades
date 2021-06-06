#!/usr/bin/python

############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing SPAdes
# provide a path to .yaml file with test description

import workflow_base
import teamcity_workflow
from workflow_base import log

def cmp_phases_output(dataset_info, test, output_dir):
    outputs = workflow_base.get_outdirs(dataset_info, test, output_dir)
    for i in range(1, len(outputs)):
        ecode = teamcity_workflow.cmp_with_etalon(outputs[0], outputs[i],
                                                  allowed_substr=[".yaml", ".sh", "params.txt"])
        if ecode != 0:
            log.err("Comparing outputs did not pass, exit code " + str(ecode))
            return 12
    return 0


workflow_base.main(cmp_phases_output)