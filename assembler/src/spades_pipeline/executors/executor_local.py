#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import support
import executors
import commands_parser
import options_storage


class Executor(executors.ExecutorBase):
    def __init__(self, log):
        super(Executor, self).__init__(log)

    def execute(self, commands):
        for num in range(len(commands)):
            command = commands[num]
            if options_storage.args.continue_mode:
                stage_file_name = "stage_%d_%s" % (num, command.short_name)
                stage_checkpoint_path = os.path.join(options_storage.args.output_dir, stage_file_name)
                if os.path.isfile(stage_checkpoint_path) and \
                        ("_start" not in command.short_name) and \
                        ("_finish" not in command.short_name):
                    self.log.info("===== Skipping %s (already processed)" % command.STAGE)
                    continue

            if "_finish" not in command.short_name:
                self.log.info("\n===== %s started. \n" % command.STAGE)

            if command.__str__() != "true":
                self.log.info("\n== Running: %s\n" % command.__str__())
                command.run(self.log)

            self.rm_files(command)
            self.check_output(command)

            if "_start" not in command.short_name:
                self.log.info("\n===== %s finished. \n" % command.STAGE)

            self.touch_file(command, num)
            if options_storage.args.stop_after == command.short_name or \
                    ("_finish" in command.short_name and
                             options_storage.args.stop_after == command.short_name.split('_')[0]):
                self.log.info("\n======= Skipping the rest of SPAdes "
                              "pipeline (--stop-after was set to '%s'). "
                              "You can continue later with --continue or "
                              "--restart-from options\n" % options_storage.args.stop_after)
                break

    def rm_files(self, command):
        for fpath in command.del_after:
            if os.path.isdir(fpath):
                shutil.rmtree(fpath)
            elif os.path.isfile(fpath):
                os.remove(fpath)

    def check_output(self, command):
        for fpath in command.output_files:
            if not os.path.isfile(fpath):
                support.error(command.STAGE + " finished abnormally: %s not found!" % fpath)

    def dump_commands(self, commands, outputfile):
        commands_parser.write_commands_to_sh(commands, outputfile)

    def touch_file(self, command, num):
        stage_file_name = "stage_%d_%s" % (num, command.short_name)
        stage_checkpoint_path = os.path.join(options_storage.args.output_dir, stage_file_name)
        open(stage_checkpoint_path, 'a').close()
