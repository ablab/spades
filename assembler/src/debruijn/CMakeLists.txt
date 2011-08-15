project(debruijn CXX)

add_executable(debruijn
 main.cpp)

target_link_libraries(debruijn ${COMMON_LIBRARIES} iolib)

