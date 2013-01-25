############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#########
# ENUMS #
#########

class Filter:
    normal = 0
    hist = 1
    closest = 2
    spades = 3
    pathsets = 4


##########
# CONFIG #
##########

filter = Filter.spades

min_overlap = 10 # for missing rectangles
