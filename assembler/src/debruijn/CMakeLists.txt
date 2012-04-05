project(spades CXX)

add_executable(spades
 main.cpp)

target_link_libraries(spades ${COMMON_LIBRARIES} input)

