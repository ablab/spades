############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(abruijn CXX)

add_executable(abruijn
   main.cpp
   hash.cpp 
   abruijn_graph.cpp 
   graph_builder.cpp 
   toy_examples.cpp
)

target_link_libraries(abruijn ${COMMON_LIBRARIES} input)
