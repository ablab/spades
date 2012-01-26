project(include_test CXX)

add_executable(include_test
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 test.cpp)

target_link_libraries(include_test ${COMMON_LIBRARIES} input)

