############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(graph_support CXX)

add_library(graph_support STATIC
            connected_component.cpp bidirectional_path.cpp genome_consistance_checker.cpp scaff_supplementary.cpp)

