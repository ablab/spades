############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Helper for testing.
"""

import sys
import warnings
import os.path


def warnings_to_stdout():
    """ Redirect all warnings to stdout.
    """
    showwarning_orig = warnings.showwarning

    def showwarning(msg, cat, fname, lno, file=None, line=0):
        showwarning_orig(msg, cat, os.path.basename(fname), line, sys.stdout)

    warnings.showwarning = showwarning
    #warnings.simplefilter('always')
