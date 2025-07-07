#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from abc import ABCMeta, abstractmethod

from .. import support
from .. import commands_parser
from ..options_storage import OptionStorage

options_storage = OptionStorage()


class ExecutorBase(object):
    __metaclass__ = ABCMeta

    def __init__(self, log):
        self.log = log

    @abstractmethod
    def execute(self, commands):
        pass

    @abstractmethod
    def dump_commands(self, commands, outputfile):
        pass

    @abstractmethod
    def join(self, job_name):
        pass

    @abstractmethod
    def kill(self, job_name):
        pass


class ExecutorCluster(ExecutorBase):
    grid_engine = None
    grid_engine_submit_command = None
    grid_engine_name_option = None
    grid_engine_output_option = None
    grid_engine_err_output_option = None
    grid_engine_thread_option = None
    grid_engine_dependency_option = None
    grid_engine_memory_option = None
    grid_engine_queue = None
    grid_engine_minimum_node_mem = None
    grid_engine_mpi_runtime = None
    grid_engine_mpi_runtime_args = None
    grid_engine_wait_command = None
    grid_engine_kill_command = None

    def join(self, job_name):
        support.sys_call(self.grid_engine_wait_command.format(JOB_NAME=job_name), logger_instance=self.log)

    def kill(self, job_name):
        # support.sys_call(self.grid_engine_kill_command.format(JOB_NAME=job_name), log=self.log)
        os.system(self.grid_engine_kill_command.format(JOB_NAME=job_name))

    def run_cluster_command(self, cmd, uuid):
        support.sys_call(cmd, logger_instance=self.log)
        return uuid

    def execute(self, commands):
        jobs = []
        def prev_id():
            if not jobs:
                return ""
            else:
                return jobs[-1]

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

            # `true' command does nothing, it corresponds to an arbitrary stage
            # used for cleanup, restart-from, and other such stuff We skip its
            # actual running for the sake of log purity and beauty
            if command.__str__() != "true":
                self.log.info("\n==Submitting: %s\n" % command.__str__())
                if commands[num].mpi_support:
                    cmd = self.get_MPI_command(command, prev_id())
                else:
                    cmd = self.get_not_MPI_command(command, prev_id())
                jid = self.run_cluster_command(cmd, command.job_uuid)
                if "_start" not in command.short_name:
                    self.log.info("\n===== %s submitted. Job ID: %s \n" % (command.stage, jid))
                jobs.append(jid)

            touch_command = commands_parser.Command(command.stage + "_touch",
                    "touch",
                                                    [stage_checkpoint_path],
                    "touch",
                                                    job_uuid=command.job_uuid + "_touch")

            touch_jid = self.run_cluster_command(self.get_not_MPI_command(touch_command, prev_id()), touch_command.job_uuid)
            jobs.append(touch_jid)

            # FIXME implement
            # self.rm_files(command)
            # self.check_output(command)

            if options_storage.args.stop_after == command.short_name or \
                    ("_finish" in command.short_name and
                             options_storage.args.stop_after == command.short_name.split('_')[0]):
                self.log.info("\n======= Skipping the rest of SPAdes "
                              "pipeline (--stop-after was set to '%s'). "
                              "You can continue later with --continue or "
                              "--restart-from options\n" % options_storage.args.stop_after)
                break
        return jobs

    def dump_commands(self, commands, outputfile):
        with open(outputfile, 'w') as fw:
            fw.write(self.get_MPI_sh_preambula() + "\n")

            prev_id = ""
            for i in range(len(commands)):
                if (commands[i].mpi_support):
                    cmd = self.get_MPI_sh_command(commands[i], prev_id)
                else:
                    cmd = self.get_not_MPI_sh_command(commands[i], prev_id)
                fw.write(cmd + "\n")
                prev_id = commands[i].job_uuid
        self.log.info("Commands were saved to " + outputfile)

    def get_MPI_sh_preambula(self):
        preambula = ""
        log_file = options_storage.args.output_dir + "/spades.log"
        preambula += "LOG_OUT=\"" + self.grid_engine_output_option.format(OUT=log_file) + "\"\n"
        preambula += "ERR_OUT=\"" + self.grid_engine_err_output_option.format(ERR=log_file) + "\"\n"
        memory_in_kb = int(options_storage.args.memory * 1024 * 1024)
        preambula += "QUEUE=\"" + self.grid_engine_queue.format(QUEUE=options_storage.args.grid_queue)  + "\"\n"
        preambula += "CLUSTER_ARGS=\"$QUEUE " + \
                     self.grid_engine_memory_option.format(MEMORY=memory_in_kb, TOTAL_MEMORY=memory_in_kb * options_storage.args.grid_nnodes) + " " + \
                     self.grid_engine_thread_option.format(NNODES=options_storage.args.grid_nnodes,
                                                           NCPUS=options_storage.args.threads,
                                                           NPROCESSORS=options_storage.args.grid_nnodes * options_storage.args.threads) + " " + \
                     self.grid_engine_minimum_node_mem.format(MEMORY=memory_in_kb) + "\"\n"
        preambula += "MPIRUN_ARGS=\"" + self.grid_engine_mpi_runtime_args.format(
            NNODES=options_storage.args.grid_nnodes,
            NCPUS=options_storage.args.threads) + "\"\n"
        preambula += "PYTHON=\"" + sys.executable + "\"\n"

        return preambula

    def get_MPI_sh_command(self, command, prev_job_name=""):
        cmd = self.grid_engine_submit_command + " "
        cmd += self.grid_engine_name_option.format(JOB_NAME=command.job_uuid) + " "
        cmd += "$LOG_OUT "
        cmd += "$ERR_OUT "
        if prev_job_name != "":
            cmd += self.grid_engine_dependency_option.format(WAIT_TAG=prev_job_name) + " "
        cmd += "$CLUSTER_ARGS "
        cmd += self.grid_engine_mpi_runtime + " $MPIRUN_ARGS "
        cmd1 = cmd
        cmd = "# === STAGE " + command.stage + "(MPI) === \n"
        cmd += "CMD=\"" + command.mpi_sh_str() + "\"\n\n"
        cmd += cmd1
        cmd += "$CMD\n\n"
        return cmd

    def get_not_MPI_sh_command(self, command, prev_job_name=""):
        cmd = "#=== STAGE " + command.stage + " (not MPI) ===\n"
        cmd += "CMD=\"" + command.sh_str() + "\"\n\n"

        cmd += self.grid_engine_submit_command + " "
        cmd += self.grid_engine_name_option.format(JOB_NAME=command.job_uuid) + " "
        cmd += "$LOG_OUT "
        cmd += "$ERR_OUT "
        if prev_job_name != "":
            cmd += self.grid_engine_dependency_option.format(WAIT_TAG=prev_job_name) + " "
        cmd += "$QUEUE "
        cmd += "$CMD\n\n"
        return cmd


    def get_MPI_command(self, command, prev_job_name=""):
        cmd = self.grid_engine_submit_command + " "
        cmd += self.grid_engine_name_option.format(JOB_NAME=command.job_uuid) + " "
        log_file = options_storage.args.output_dir + "/spades.log"
        cmd += self.grid_engine_output_option.format(OUT=log_file) + "  "
        cmd += self.grid_engine_err_output_option.format(ERR=log_file) + "  "
        if prev_job_name != "":
            cmd += self.grid_engine_dependency_option.format(WAIT_TAG=prev_job_name) + " "
        cmd += self.grid_engine_queue.format(QUEUE=options_storage.args.grid_queue) + " "
        memory_in_kb = int(options_storage.args.memory * 1024 * 1024)
        cmd += self.grid_engine_memory_option.format(MEMORY=memory_in_kb, TOTAL_MEMORY=memory_in_kb * options_storage.args.grid_nnodes) + " "
        cmd += self.grid_engine_thread_option.format(NNODES=options_storage.args.grid_nnodes,
                NCPUS=options_storage.args.threads,
                NPROCESSORS=options_storage.args.grid_nnodes * options_storage.args.threads) + " "
        cmd += self.grid_engine_minimum_node_mem.format(MEMORY=memory_in_kb) + " "

        cmd += self.grid_engine_mpi_runtime + " " + self.grid_engine_mpi_runtime_args.format(
            NNODES=options_storage.args.grid_nnodes,
            NCPUS=options_storage.args.threads) + " "

        if options_storage.args.grid_profile:
            name = command.stage + "_" + command.short_name + "_" + command.job_uuid
            profile = options_storage.args.output_dir + "/" + name + ".prof"
            profile_line = " -x CPUPROFILE={PROFILE} ompi_profile_helper.sh ".format(PROFILE=profile)
        else:
            profile_line = ""

        if options_storage.args.grid_valgrind:
            valgrind_line = " valgrind  --track-origins=yes "
        else:
            valgrind_line = ""

        if options_storage.args.grid_coredump:
            coredump_line = "ulimit -c unlimited && "
        else:
            coredump_line = ""
        cmd = coredump_line + " " + cmd
        cmd += profile_line + " "
        cmd += valgrind_line + " "
        cmd += command.mpi_str()
        return cmd


    def get_not_MPI_command(self, command, prev_job_name=""):
        cmd = self.grid_engine_submit_command + " "
        cmd += self.grid_engine_name_option.format(JOB_NAME=command.job_uuid) + " "
        log_file = options_storage.args.output_dir + "/spades.log"
        cmd += self.grid_engine_output_option.format(OUT=log_file) + "  "
        cmd += self.grid_engine_err_output_option.format(ERR=log_file) + "  "
        if prev_job_name != "":
            cmd += self.grid_engine_dependency_option.format(WAIT_TAG=prev_job_name) + " "

        cmd += self.grid_engine_queue.format(QUEUE=options_storage.args.grid_queue) + " "

        cmd += command.__str__()
        return cmd
