project(llvm-support)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

file(GLOB sources "[a-zA-Z]*.cpp" "[a-zA-Z]*.c")

add_library(llvm-support STATIC
            ${sources})
