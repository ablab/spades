project(yaml-cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

set(YAML_CPP_VERSION_MAJOR "0")
set(YAML_CPP_VERSION_MINOR "5")
set(YAML_CPP_VERSION_PATCH "0")
set(YAML_CPP_VERSION "${YAML_CPP_VERSION_MAJOR}.${YAML_CPP_VERSION_MINOR}.${YAML_CPP_VERSION_PATCH}")

file(GLOB sources "[a-zA-Z]*.cpp")

add_library(yaml-cpp STATIC
            ${sources})