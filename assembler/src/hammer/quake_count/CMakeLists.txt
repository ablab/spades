############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake_count CXX)

add_executable(quake_count quake_count.cpp)

add_executable(quake_count_17 quake_count_17.cpp)
add_executable(quake_count_19 quake_count_19.cpp)
add_executable(quake_count_21 quake_count_21.cpp)
add_executable(quake_count_25 quake_count_25.cpp)
add_executable(quake_count_29 quake_count_29.cpp)
add_executable(quake_count_33 quake_count_33.cpp)
add_executable(quake_count_37 quake_count_37.cpp)
add_executable(quake_count_45 quake_count_45.cpp)
add_executable(quake_count_55 quake_count_55.cpp)
add_executable(quake_count_65 quake_count_65.cpp)
add_executable(quake_count_75 quake_count_75.cpp)


target_link_libraries(quake_count ${COMMON_LIBRARIES})
target_link_libraries(quake_count_17 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_19 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_21 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_25 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_29 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_33 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_37 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_45 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_55 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_65 ${COMMON_LIBRARIES})
target_link_libraries(quake_count_75 ${COMMON_LIBRARIES})




