#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
from .. import support
from . import executor_base
from .. import commands_parser
from ..options_storage import OptionStorage

options_storage = OptionStorage()

class Executor(executor_base.ExecutorBase):
    def __init__(self, log):
        super(Executor, self).__init__(log)

    def execute(self, commands):
        for num, command in enumerate(commands):
            stage_checkpoint_path = options_storage.get_stage_filename(num, command.short_name)
            if options_storage.args.continue_mode:
                if os.path.isfile(stage_checkpoint_path) and \
                        ("_start" not in command.short_name) and \
                        ("_finish" not in command.short_name):
                    self.log.info("===== Skipping %s (already processed)" % command.stage)
                    continue

            if "_finish" not in command.short_name:
                self.log.info("\n===== %s started. \n" % command.stage)

            if command.__str__() != "true":
                self.log.info("\n== Running: %s\n" % command.__str__())
                command.run(self.log)

            self.rm_files(command)
            self.check_output(command)

            if "_start" not in command.short_name:
                self.log.info("\n===== %s finished. \n" % command.stage)

            self.touch_file(command, num)
            if options_storage.args.stop_after == command.short_name or \
                    ("_finish" in command.short_name and
                             options_storage.args.stop_after == command.short_name.split('_')[0]):
                self.log.info("\n======= Skipping the rest of SPAdes "
                              "pipeline (--stop-after was set to '%s'). "
                              "You can continue later with --continue or "
                              "--restart-from options\n" % options_storage.args.stop_after)
                return None
        return None

    def rm_files(self, command):
        if options_storage.args.no_clear_after:
            return

        for fpath in command.del_after:
            fpath_abs = os.path.join(options_storage.args.output_dir, fpath)
            if os.path.isdir(fpath_abs):
                shutil.rmtree(fpath_abs)
            elif os.path.isfile(fpath_abs):
                os.remove(fpath_abs)

    def check_output(self, command):
        for fpath in command.output_files:
            if not os.path.isfile(fpath):
                support.error(command.stage + " finished abnormally: %s not found!" % fpath)

    def dump_commands(self, commands, outputfile):
        commands_parser.write_commands_to_sh(commands, outputfile)

    def touch_file(self, command, num):
        path = options_storage.get_stage_filename(num, command.short_name)
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        open(path, 'a').close()

    def join(self, job_name):
        assert (job_name is None)

    def kill(self, job_name):
        assert (job_name is None)
