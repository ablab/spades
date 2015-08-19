#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web_service.settings")

    from django.core.management import execute_from_command_line

    execute_from_command_line(sys.argv)
