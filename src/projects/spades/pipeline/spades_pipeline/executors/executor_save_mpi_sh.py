#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
from . import executor_base
from .. import commands_parser
from ..options_storage import OptionStorage

options_storage = OptionStorage()


class Executor(executor_base.ExecutorBase):
    def __init__(self, log):
        super(Executor, self).__init__(log)

    def execute(self, commands):
        super(Executor, self).execute(commands)
        commands_parser.write_commands_to_mpi_sh(commands, os.path.join(options_storage.args.output_dir, "run_spades.sh"))
        commands_parser.write_commands_to_yaml(commands,
                                               os.path.join(options_storage.args.output_dir,
                                                            "run_spades.yaml"))
        return None

    def dump_commands(self, commands, outputfile):
        commands_parser.write_commands_to_mpi_sh(commands, outputfile)

    def join(self, job_name):
        assert (job_name is None)

    def kill(self, job_name):
        assert (job_name is None)
