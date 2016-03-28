############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(graph_support CXX)

add_library(graph_support STATIC
            components/connected_component.cpp paths/bidirectional_path.cpp graph_support/scaff_supplementary.cpp)

