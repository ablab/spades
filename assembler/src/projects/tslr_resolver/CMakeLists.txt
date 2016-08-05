project(tslr_resolver CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(tslr_resolver main.cpp)

target_link_libraries(tslr_resolver spades_modules ${COMMON_LIBRARIES} ssw)