#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import support
import executors
import commands_parser
import options_storage
import executor_local


class Executor(executor_local.Executor):
    def execute(self, commands):
        for num, command in enumerate(commands):
            stage_checkpoint_path = options_storage.get_stage_filename(num, command.short_name)

            if options_storage.args.continue_mode:
                if os.path.isfile(stage_checkpoint_path) and \
                        ("_start" not in command.short_name) and \
                        ("_finish" not in command.short_name):
                    self.log.info("===== Skipping %s (already processed)" % command.STAGE)
                    continue

            if "_finish" not in command.short_name:
                self.log.info("\n===== %s started. \n" % command.STAGE)

            # `true' command does nothing, it corresponds to an arbitrary stage
            # used for cleanup, restart-from, and other such stuff. We skip its
            # actual running for the sake of log purity and beauty
            if command.__str__() != "true":
                if (command.mpi_support):
                    # cmd = "mpiexec -np 4 xterm -e gdb -ex run --args " + command.__str__()
                    valgrind = "valgrind" if options_storage.args.grid_valgrind else ""
                    cmd = "mpiexec --bind-to none -np {NODES} {VALGRIND} ".format(NODES=options_storage.args.grid_nnodes, VALGRIND=valgrind) + command.mpi_str()
                    self.log.info("\n== Running: %s\n" % cmd)
                    support.sys_call(cmd, self.log)
                else:
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
                return None
        return None

    def dump_commands(self, commands, outputfile):
        commands_parser.write_commands_to_mpi_sh(commands, outputfile)
