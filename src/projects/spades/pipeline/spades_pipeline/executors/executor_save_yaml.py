#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import executors
import commands_parser
import options_storage


class Executor(executors.ExecutorBase):
    def __init__(self, log):
        super(Executor, self).__init__(log)

    def execute(self, commands):
        super(Executor, self).execute(commands)
        commands_parser.write_commands_to_sh(commands, os.path.join(options_storage.args.output_dir, "run_spades.sh"))
        commands_parser.write_commands_to_yaml(commands,
                                               os.path.join(options_storage.args.output_dir,
                                               "run_spades.yaml"))
        return None

    def dump_commands(self, commands, outputfile):
        commands_parser.write_commands_to_sh(commands, outputfile)

    def join(self, job_name):
        assert (job_name is None)

    def kill(self, job_name):
        assert (job_name is None)
