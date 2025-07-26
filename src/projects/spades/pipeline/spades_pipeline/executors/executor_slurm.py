#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import sys

from . import executor_base
from ..options_storage import OptionStorage
from .. import support

options_storage = OptionStorage()


class Executor(executor_base.ExecutorCluster):
    grid_engine = "SLURM"
    grid_engine_submit_command = "sbatch"
    grid_engine_slurm_args = ("--hint=compute_bound --mem-bind=verbose,none --cpus-per-task {NCPUS} --open-mode=append "
                              "--kill-on-invalid-dep=yes --mem {MEMORY_MB}M --time {TIME} {EXTRA}")
    grid_engine_output_option = "-o {OUT}"
    grid_engine_err_output_option = "-e {ERR}"
    grid_engine_job_name = "--job-name {JOB_NAME}"
    grid_engine_set_command = "--wrap \"{COMMAND}\""
    grid_engine_dependency_option = "--dependency afterok:{WAIT_TAG}"
    grid_engine_nodes = "--nodes={NNODES} --ntasks={NNODES}"
    grid_engine_kill_command = "scancel {JOB_NAME}"
    grid_engine_srun_args = "--cpus-per-task {NCPUS}"

    def join(self, job_name):
        log_file = options_storage.args.output_dir + "/spades.log"
        cmd = self.grid_engine_submit_command.format(COMMAND="true", JOB_NAME="wait", OUT=log_file, ERR=log_file, NCPUS=1)
        cmd += " " + self.grid_engine_dependency_option.format(WAIT_TAG=job_name)
        cmd += " " + self.grid_engine_credentials.format(QUEUE=options_storage.args.grid_queue)
        cmd += " --wait"
        support.sys_call(cmd, logger_instance=self.log)

    def get_MPI_sh_preambula(self):
        memory_mb = int(options_storage.args.memory * 1024)
        preambula = "SLURM_ARGS=\"" + self.grid_engine_slurm_args.format(NCPUS=options_storage.args.threads,
                                                                         MEMORY_MB=memory_mb,
                                                                         TIME=options_storage.args.grid_time,
                                                                         EXTRA=options_storage.args.grid_extra,
                                                                         QUEUE=options_storage.args.grid_queue) + "\"\n"
        preambula += "SRUN_ARGS=\"" + self.grid_engine_srun_args.format(NCPUS=options_storage.args.threads) + "\"\n"
        log_file = options_storage.args.output_dir + "/spades.log"
        preambula += "LOG_OUT=\"" + self.grid_engine_output_option.format(OUT=log_file) + "\"\n"
        preambula += "ERR_OUT=\"" + self.grid_engine_err_output_option.format(ERR=log_file) + "\"\n"
        preambula += "PYTHON=\"" + sys.executable + "\"\n"
        return preambula

    def get_sh_command(self, command, prev_id, mpi):
        cmd_str = "#=== STAGE " + command.stage + (" (MPI) ===\n" if mpi else " (not MPI) ===\n")
        cmd_str += "CMD=\"" + command.mpi_sh_str() + "\"\n"
        cmd_str += "SID1=$(" + self.grid_engine_submit_command + " $SLURM_ARGS " + \
                   self.grid_engine_job_name.format(JOB_NAME=command.job_uuid) + " $LOG_OUT $ERR_OUT "
        if mpi:
            cmd_str += self.grid_engine_set_command.format(COMMAND="srun $SRUN_ARGS $CMD")
        else:
            cmd_str += self.grid_engine_set_command.format(COMMAND="$CMD")

        if prev_id != "":
            cmd_str += " " + self.grid_engine_dependency_option.format(WAIT_TAG="$SID1")

        if mpi:
            cmd_str += " " + self.grid_engine_nodes.format(NNODES=options_storage.args.grid_nnodes)
        cmd_str += ")\n"
        cmd_str += "SID1=\"${SID1##* }\"\n"
        return cmd_str

    def get_command(self, command, prev_id, mpi):
        log_file = options_storage.args.output_dir + "/spades.log"
        if mpi:
            if options_storage.args.grid_profile:
                name = command.stage + "_" + command.short_name + "_" + command.job_uuid
                profile = options_storage.args.output_dir + "/" + name + ".prof"
                profile_line = "-x CPUPROFILE={PROFILE} ompi_profile_helper.sh".format(PROFILE=profile)
            else:
                profile_line = ""

            if options_storage.args.grid_valgrind:
                # valgrind_line = "valgrind  --track-origins=yes --suppressions=/global/software/sl-7.x86_64/modules/gcc/7.4.0/openmpi/4.0.1-gcc/share/openmpi/openmpi-valgrind.supp"
                valgrind_line = "-x valgrind  --track-origins=yes"
            else:
                valgrind_line = ""

            if options_storage.args.grid_coredump:
                coredump_line = "ulimit -c unlimited;"
            else:
                coredump_line = ""
            command_line = ("{COREDUMP} srun ".format(COREDUMP=coredump_line) +
                            self.grid_engine_srun_args.format(NCPUS=options_storage.args.threads) +
                            " {VALGRIND} {PROFILE}".format(PROFILE=profile_line, VALGRIND=valgrind_line) +
                            " " + command.mpi_str())
        else:
            command_line = command.mpi_str()
        cmd = self.grid_engine_submit_command + " "
        memory_mb = int(options_storage.args.memory * 1024)
        cmd += self.grid_engine_slurm_args.format(NCPUS=options_storage.args.threads,
                                                  MEMORY_MB=memory_mb,
                                                  TIME=options_storage.args.grid_time,
                                                  EXTRA=options_storage.args.grid_extra,
                                                  QUEUE=options_storage.args.grid_queue) + " "
        cmd += self.grid_engine_job_name.format(JOB_NAME=command.job_uuid) + " "
        cmd += self.grid_engine_err_output_option.format(ERR=log_file) + " "
        cmd += self.grid_engine_output_option.format(OUT=log_file) + " "
        cmd += self.grid_engine_set_command.format(COMMAND=command_line) + " "

        if prev_id != "":
            cmd += " " + self.grid_engine_dependency_option.format(WAIT_TAG=prev_id)

        if mpi:
            cmd += " " + self.grid_engine_nodes.format(NNODES=options_storage.args.grid_nnodes)
        return cmd

    def get_MPI_command(self, command, prev_id=""):
        return self.get_command(command, prev_id=prev_id, mpi=True)

    def get_not_MPI_command(self, command, prev_id=""):
        return self.get_command(command, prev_id=prev_id, mpi=False)

    def get_MPI_sh_command(self, command, prev_id=""):
        return self.get_sh_command(command, prev_id=prev_id, mpi=True)

    def get_not_MPI_sh_command(self, command, prev_id=""):
        return self.get_sh_command(command, prev_id=prev_id, mpi=False)

    def run_cluster_command(self, cmd, uuid):
        import re
        import os
        self.log.info("Submit cluster job: " + cmd)
        output = os.popen(cmd).read()
        jobid_search = re.search(r"^Submitted batch job (\d+)$", output)
        assert jobid_search
        return jobid_search.group(1)
