#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
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

    def dump_commands(self, commands, outputfile):
        commands_parser.write_commands_to_sh(commands, outputfile)

    def touch_file(self, command):
        pass
