#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
from abc import ABCMeta, abstractmethod

import options_storage
import support
import commands_parser


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
    def touch_file(self, command):
        pass
