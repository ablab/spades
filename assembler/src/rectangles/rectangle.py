############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from diagonal import Diagonal

class Rectangle(object):
    def __init__(self, e1, e2):
        self.e1, self.e2 = e1, e2
        self.diagonals = {} # (D, pathset) -> Diagonal

    def add_diagonal(self, d, D, pathset=None):
        if (D, pathset) not in self.diagonals:
            self.diagonals[D, pathset] = Diagonal(self, d, D, pathset)

    def get_closest_diagonal(self, D):
        min_diff = 1e100
        closest = None
        for diag in self.diagonals.itervalues(): # TODO: optimize, linear -> log (me sure?)
            if abs(diag.D - D) < min_diff:
                min_diff = abs(diag.D - D)
                closest = diag
        return closest


