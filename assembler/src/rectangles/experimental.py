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
