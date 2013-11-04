project(cgce CXX)

include_directories(${CMAKE_SOURCE_DIR}/debruijn
                    ${CMAKE_SOURCE_DIR}/online_vis
                    ${CMAKE_SOURCE_DIR}/cap)

add_executable(cgce
  cgce.cpp)

target_link_libraries(cgce input mph_index yaml-cpp ${COMMON_LIBRARIES})
