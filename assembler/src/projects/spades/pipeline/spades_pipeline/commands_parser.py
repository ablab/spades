#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random
import string

import support
import sys


class Command(object):
    def __init__(self, STAGE, path, args, short_name, config_dir="",
                 del_after=None, output_files=None):
        self.STAGE = STAGE
        self.path = path
        self.args = args
        self.short_name = short_name
        self.config_dir = config_dir
        self.del_after = del_after
        if self.del_after is None:
            self.del_after = []
        self.output_files = output_files
        if self.output_files is None:
            self.output_files = []

    def to_list(self):
        return [self.path] + self.args

    def __str__(self):
        return ' '.join(self.to_list())

    def run(self, log):
        support.sys_call(self.to_list(), log)

    def to_dict(self):
        return {"STAGE": self.STAGE,
                "path": self.path,
                "args": self.args,
                "short_name": self.short_name,
                "config_dir": self.config_dir,
                "output_files": self.output_files,
                "del_after": self.del_after}


def write_commands_to_sh(commands, output_file):
    with open(output_file, 'w') as fw:
        fw.write("set -e\n")
        for command in commands:
            fw.write(command.__str__() + "\n")


def write_commands_to_yaml(commands, output_file):
    if sys.version.startswith("2."):
        import pyyaml2 as yaml
    elif sys.version.startswith("3."):
        import pyyaml3 as yaml

    data = [command.to_dict() for command in commands]

    with open(output_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def read_commands_from_yaml(yaml_fpath):
    if sys.version.startswith("2."):
        import pyyaml2 as yaml
    elif sys.version.startswith("3."):
        import pyyaml3 as yaml

    with open(yaml_fpath) as stream:
        data = yaml.load(stream)
    commands = []
    for kwargs in data:
        commands.append(Command(**kwargs))
    return commands
