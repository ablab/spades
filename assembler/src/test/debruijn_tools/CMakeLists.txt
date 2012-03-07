project(debruijn_tools CXX)

include_directories("${CMAKE_SOURCE_DIR}/debruijn")

add_executable(debruijn_tools
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 tools.cpp)

target_link_libraries(debruijn_tools ${COMMON_LIBRARIES} input)
