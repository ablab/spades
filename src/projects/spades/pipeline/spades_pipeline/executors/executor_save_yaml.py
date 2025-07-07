#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
from .executor_base import ExecutorBase
from ..commands_parser import write_commands_to_sh, write_commands_to_yaml
from ..options_storage import OptionStorage

options_storage = OptionStorage()


class Executor(ExecutorBase):
    def __init__(self, log):
        super(Executor, self).__init__(log)

    def execute(self, commands):
        super(Executor, self).execute(commands)
        write_commands_to_sh(commands, os.path.join(options_storage.args.output_dir, "run_spades.sh"))
        write_commands_to_yaml(commands,
                               os.path.join(options_storage.args.output_dir,
                                            "run_spades.yaml"))
        return None

    def dump_commands(self, commands, outputfile):
        write_commands_to_sh(commands, outputfile)

    def join(self, job_name):
        assert (job_name is None)

    def kill(self, job_name):
        assert (job_name is None)
