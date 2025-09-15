#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random
import string
import sys

from . import support
from .options_storage import OptionStorage
options_storage = OptionStorage()


class Command(object):
    def __init__(self, stage, path, args, short_name, config_dir="",
                 mpi_support=False, job_uuid="",
                 del_after=None, output_files=None):
        self.stage = stage
        self.path = path
        self.args = args
        self.short_name = short_name
        self.mpi_support = mpi_support
        self.job_uuid = self.generate_job_uuid()
        if job_uuid != "":
            self.job_uuid = job_uuid
        self.config_dir = config_dir
        self.del_after = del_after
        if self.del_after is None:
            self.del_after = []
        self.output_files = output_files
        if self.output_files is None:
            self.output_files = []

    def to_list(self):
        return [self.path.format(spades_core="spades-core")] + self.args

    def to_sh_list(self):
        if self.path == sys.executable:
            return ["$PYTHON"] + self.args
        return [self.path.format(spades_core="spades-core")] + self.args

    def to_mpi_list(self):
        return [self.path.format(spades_core="spades-hpc")] + self.args

    def to_mpi_sh_list(self):
        if self.path == sys.executable:
            return ["$PYTHON"] + self.args
        return [self.path.format(spades_core="spades-hpc")] + self.args

    def __str__(self):
        return ' '.join(self.to_list())

    def sh_str(self):
        return ' '.join(self.to_sh_list())

    def mpi_str(self):
        return ' '.join(self.to_mpi_list())

    def mpi_sh_str(self):
        return ' '.join(self.to_mpi_sh_list())

    def run(self, log):
        support.sys_call(self.to_list(), log)

    def to_dict(self, spades_core="spades-core"):
        return {"STAGE": self.stage,
                "path": self.path.format(spades_core=spades_core),
                "args": self.args,
                "short_name": self.short_name,
                "mpi_support": self.mpi_support,
                "job_uuid": self.job_uuid,
                "config_dir": self.config_dir,
                "output_files": self.output_files,
                "del_after": self.del_after}

    def generate_job_uuid(self):
        return (('hpcSPAdes_' if self.mpi_support else 'SPAdes_') + self.stage.replace(' ', '_') + "_" +
                ''.join([random.choice(string.ascii_uppercase + string.digits) for _ in range(32)]))


def write_commands_to_sh(commands, output_file):
    with open(output_file, 'w') as fw:
        fw.write("set -e\n")
        for command in commands:
            fw.write(str(command) + "\n")


def write_commands_to_mpi_sh(commands, output_file):
    with open(output_file, 'w') as fw:
        fw.write("set -e\n")
        for command in commands:
            fw.write(command.mpi_str() + "\n")


def write_commands_to_yaml(commands, output_file):
    import pyyaml3 as yaml

    data = [command.to_dict() for command in commands]

    with open(output_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def write_commands_to_mpi_yaml(commands, output_file):
    import pyyaml3 as yaml

    data = [command.to_dict(spades_core="spades-hpc") for command in commands]

    with open(output_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def read_commands_from_yaml(yaml_fpath):
    import pyyaml3 as yaml

    with open(yaml_fpath) as stream:
        data = yaml.load(stream)
    commands = []
    for kwargs in data:
        commands.append(Command(**kwargs))
    return commands
