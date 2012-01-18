#!/usr/bin/env python

import sys
import support

def config_file_name():
    if len(sys.argv) != 2:
        support.error("Usage: <script_name>.py <config_file_name>", "")

    return sys.argv[1]

def file_lines(filename):
    return open(filename).readlines()

def vars_from_lines(lines):

    class var_metadata:
        def __init__(self, value, line_num, indent):
            self.value    = value
            self.line_num = line_num
            self.indent   = indent

    def valid_var_name(name):
        for sym in name:
            if not sym.isalpha() and sym != '_':
                return False

            return True

    def var_from_line(line, line_num):
        l = line.split()
        if len(l) == 0 or not valid_var_name(l[0]):
            return None, None

        def indent(s):
            return s[ : len(s) - len(s.lstrip())]

        return l[0], var_metadata(l[1:], line_num, indent(line))

    vars = dict()

    for i in xrange(len(lines)):
        var, meta = var_from_line(lines[i], i)
        if var is not None:
            vars[var] = meta

    return vars

def substitute_params(filename, var_dict):

    lines = file_lines(filename)
    vars_in_file = vars_from_lines(lines)

    for var, value in var_dict.items():
        if var not in vars_in_file:
            error("Couldn't find " + var + " in " + filename)

        meta = vars_in_file[var]
        lines[meta.line_num] = meta.indent + str(var) + " " + str(value) + "\n"

    file = open(filename, "w")
    file.writelines(lines)

def load_config_from_vars(cfg_vars):

    class cfg_placeholder:
        pass

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

        return load_value(value_list[0])

    for var, meta in cfg_vars.items():
        cfg.__dict__[var] = load_value_list(meta.value)

    return cfg

def load_config_from_file(filename):
    return load_config_from_vars(vars_from_lines(file_lines(filename)))

