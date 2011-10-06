project(quake_count CXX)

add_executable(correct
   bithash.cpp correct.cpp edit.cpp Read.cpp gzstream.C)

target_link_libraries(correct ${COMMON_LIBRARIES})




