#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import support

class cfg_placeholder:
    pass


def file_lines(filename):
    return open(filename).readlines()


def skip_info_comment(line):
    return line.split(';')[0].strip()


def skip_double_quotes(line):
    line = line.strip()
    if line.endswith('"'):
        line = line[:-1].strip()
        line = line.replace('"', '', 1)
    return line


def check_property(prop_line):
    if len(prop_line.split()) > 1: # property is set, i.e. has value
        if prop_line.split()[1] != "N/A":
            return True
    return False


def bool_to_str(b):
    if b:
        return "true"
    return "false"


def process_spaces(str):
    return support.process_spaces(str)


def vars_from_lines(lines):
    class var_metadata:
        def __init__(self, value, line_num, indent):
            self.value = value
            self.line_num = line_num
            self.indent = indent

    def valid_var_name(name):
        for sym in name:
            if not sym.isalpha() and sym != '_':
                return False

            return True

    def var_from_line(line, line_num):
        l = skip_double_quotes(skip_info_comment(line)).split()
        if len(l) == 0 or not valid_var_name(l[0]):
            return None, None

        def indent(s):
            return s[: len(s) - len(s.lstrip())]

        return l[0], var_metadata(l[1:], line_num, indent(line))

    vars = dict()

    for i in range(len(lines)):
        var, meta = var_from_line(lines[i], i)
        if var is not None:
            vars[var] = meta

    return vars


def substitute_params(filename, var_dict, log):
    lines = file_lines(filename)
    vars_in_file = vars_from_lines(lines)

    for var, value in var_dict.items():
        if var not in vars_in_file:
            support.error("Couldn't find " + var + " in " + filename, log)

        meta = vars_in_file[var]
        lines[meta.line_num] = meta.indent + str(var) + " " + str(value) + "\n"

    file = open(filename, "w")
    file.writelines(lines)
    file.close()


# configs with more priority should go first in parameters
def merge_configs(*cfgs):
    res = cfg_placeholder()

    for cfg in reversed(cfgs):
        res.__dict__.update(cfg.__dict__)

    return res


def load_config_from_vars(cfg_vars):
    cfg = cfg_placeholder()

    def load_value(value):
        if value == 'True' or value == 'true':
            return True
        elif value == 'False' or value == 'false':
            return False
        elif value.isdigit():
            return int(value)
        else:
            return value #string as-is

    def load_value_list(value_list):
        if len(value_list) > 1:
            return [load_value(one_value) for one_value in value_list]

        if len(value_list) == 1:
            return load_value(value_list[0])

        return None

    for var, meta in cfg_vars.items():
        cfg.__dict__[var] = load_value_list(meta.value)

    return cfg


def empty_config():
    return load_config_from_vars(dict())


def load_config_from_file(filename):
    return load_config_from_vars(vars_from_lines(file_lines(filename)))


def load_config_from_info_file(filename):
    lines = file_lines(filename)
    blocks = dict()

    cur_block_name = "common"
    blocks[cur_block_name] = []
    for i in range(1, len(lines) + 1):
        prev_line = skip_info_comment(lines[i - 1])
        cur_line = ""
        if i < len(lines):
            cur_line = lines[i]

        if cur_line.startswith('{'):
            cur_block_name = prev_line
            blocks[cur_block_name] = []
        elif cur_line.startswith('}'):
            if check_property(prev_line):
                blocks[cur_block_name].append(prev_line)
            cur_block_name = "common"
        elif not lines[i - 1].startswith('{') and not lines[i - 1].startswith('}') and not lines[i - 1].strip() == '':
            if check_property(prev_line):
                blocks[cur_block_name].append(prev_line)

    cfg = dict()
    for block_name in blocks.iterkeys():
        cfg[block_name] = load_config_from_vars(vars_from_lines(blocks[block_name]))

    return cfg
